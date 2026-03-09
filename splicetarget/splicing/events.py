#!/usr/bin/env python3
"""
splicetarget v0.1.0

events.py — Aberrant splicing event detection from classified isoforms.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum

from splicetarget.data.reference import Gene, Transcript
from splicetarget.isoforms.classify import ClassifiedIsoform, IsoformCategory

logger = logging.getLogger(__name__)


class EventType(str, Enum):
    """Types of aberrant splicing events."""

    EXON_SKIPPING = "exon_skipping"
    CRYPTIC_EXON = "cryptic_exon"
    INTRON_RETENTION = "intron_retention"
    ALT_5SS = "alt_5ss"             # alternative 5' splice site (donor shift)
    ALT_3SS = "alt_3ss"             # alternative 3' splice site (acceptor shift)
    COMPLEX = "complex"             # multiple simultaneous events

    @property
    def aso_amenable(self) -> bool:
        """Whether this event type is typically targetable by ASO/SSO therapy."""
        return self in (
            EventType.EXON_SKIPPING,
            EventType.CRYPTIC_EXON,
            EventType.INTRON_RETENTION,
            EventType.ALT_5SS,
            EventType.ALT_3SS,
        )


@dataclass
class SplicingEvent:
    """A single aberrant splicing event detected in patient data."""

    event_id: str
    event_type: EventType
    chrom: str
    start: int
    end: int
    strand: str
    gene_name: str

    # Event-specific details
    skipped_exons: list[tuple[int, int]] = field(default_factory=list)
    cryptic_exon: tuple[int, int] | None = None
    retained_intron: tuple[int, int] | None = None
    shifted_donor: tuple[int, int] | None = None     # (reference_pos, patient_pos)
    shifted_acceptor: tuple[int, int] | None = None   # (reference_pos, patient_pos)

    # Supporting evidence
    supporting_isoform_ids: list[str] = field(default_factory=list)
    total_read_support: int = 0
    fraction_of_gene_expression: float = 0.0

    # Reference context
    affected_ref_transcript: str = ""
    affected_exon_numbers: list[int] = field(default_factory=list)
    reading_frame_impact: str = ""  # "in-frame", "frameshift", "unknown"

    @property
    def event_length(self) -> int:
        return self.end - self.start

    @property
    def is_aso_target(self) -> bool:
        """Is this event a good candidate for ASO-based correction?"""
        return (
            self.event_type.aso_amenable
            and self.total_read_support >= 2
            and self.fraction_of_gene_expression >= 0.05
        )

    def summary(self) -> str:
        """One-line human-readable summary of the event."""
        summaries = {
            EventType.EXON_SKIPPING: (
                f"Exon skipping: {len(self.skipped_exons)} exon(s) skipped "
                f"at {self.chrom}:{self.start}-{self.end}"
            ),
            EventType.CRYPTIC_EXON: (
                f"Cryptic exon inclusion at {self.chrom}:{self.cryptic_exon[0]}-{self.cryptic_exon[1]}"
                if self.cryptic_exon else "Cryptic exon inclusion"
            ),
            EventType.INTRON_RETENTION: (
                f"Intron retention at {self.chrom}:{self.retained_intron[0]}-{self.retained_intron[1]}"
                if self.retained_intron else "Intron retention"
            ),
            EventType.ALT_5SS: (
                f"Alternative 5'SS: donor shifted by "
                f"{abs(self.shifted_donor[1] - self.shifted_donor[0])}bp"
                if self.shifted_donor else "Alternative 5' splice site"
            ),
            EventType.ALT_3SS: (
                f"Alternative 3'SS: acceptor shifted by "
                f"{abs(self.shifted_acceptor[1] - self.shifted_acceptor[0])}bp"
                if self.shifted_acceptor else "Alternative 3' splice site"
            ),
            EventType.COMPLEX: f"Complex splicing event at {self.chrom}:{self.start}-{self.end}",
        }
        support = f" [{self.total_read_support} reads, {self.fraction_of_gene_expression:.1%} of gene]"
        return summaries.get(self.event_type, str(self.event_type)) + support


def detect_aberrant_events(
    classified_isoforms: list[ClassifiedIsoform],
    gene: Gene,
    total_gene_reads: int | None = None,
    min_read_support: int = 2,
    junction_tolerance: int = 10,
) -> list[SplicingEvent]:
    """
    Detect aberrant splicing events by comparing patient isoforms to reference.

    Parameters
    ----------
    classified_isoforms : list[ClassifiedIsoform]
        Isoforms with structural classification.
    gene : Gene
        Reference gene model.
    total_gene_reads : int | None
        Total reads for this gene (for fraction calculation).
    min_read_support : int
        Minimum reads to report an event.
    junction_tolerance : int
        Position tolerance for junction matching.

    Returns
    -------
    list[SplicingEvent]
        Detected events sorted by read support (descending).
    """
    if total_gene_reads is None:
        total_gene_reads = sum(c.isoform.read_count for c in classified_isoforms)
    if total_gene_reads == 0:
        total_gene_reads = 1

    events: list[SplicingEvent] = []
    event_idx = 0
    canonical_tx = gene.get_canonical_transcript()

    for ci in classified_isoforms:
        if ci.isoform.read_count < min_read_support:
            continue

        frac = ci.isoform.read_count / total_gene_reads

        # Route to appropriate detection function based on category
        if ci.category == IsoformCategory.NNC:
            detected = _detect_from_nnc(ci, gene, canonical_tx, junction_tolerance)
        elif ci.category == IsoformCategory.IR:
            detected = _detect_intron_retention(ci, gene, canonical_tx)
        elif ci.category == IsoformCategory.NIC:
            detected = _detect_from_nic(ci, gene, canonical_tx, junction_tolerance)
        else:
            continue  # FSM/ISM are not aberrant

        for event in detected:
            event_idx += 1
            event.event_id = f"{gene.gene_name}.event{event_idx:04d}"
            event.supporting_isoform_ids.append(ci.isoform.isoform_id)
            event.total_read_support = ci.isoform.read_count
            event.fraction_of_gene_expression = frac
            events.append(event)

    events.sort(key=lambda e: e.total_read_support, reverse=True)

    logger.info(
        "Detected %d aberrant splicing events in %s (%d exon-skip, %d cryptic, %d IR, %d alt-SS)",
        len(events),
        gene.gene_name,
        sum(1 for e in events if e.event_type == EventType.EXON_SKIPPING),
        sum(1 for e in events if e.event_type == EventType.CRYPTIC_EXON),
        sum(1 for e in events if e.event_type == EventType.INTRON_RETENTION),
        sum(1 for e in events if e.event_type in (EventType.ALT_5SS, EventType.ALT_3SS)),
    )

    return events


def _detect_from_nnc(
    ci: ClassifiedIsoform,
    gene: Gene,
    canonical_tx: Transcript | None,
    tolerance: int,
) -> list[SplicingEvent]:
    """Detect events from NNC (Novel Not in Catalog) isoforms."""
    events = []
    iso = ci.isoform

    if canonical_tx is None:
        return events

    ref_introns = sorted(canonical_tx.introns)
    iso_introns = list(iso.intron_chain)
    ref_exon_chain = sorted((e.start, e.end) for e in canonical_tx.exons)

    # Detect exon skipping
    skipped = _find_skipped_exons(iso_introns, ref_exon_chain, ref_introns, tolerance)
    if skipped:
        frame_impact = _assess_frame_impact(skipped)
        events.append(SplicingEvent(
            event_id="",
            event_type=EventType.EXON_SKIPPING,
            chrom=iso.chrom,
            start=skipped[0][0],
            end=skipped[-1][1],
            strand=iso.strand,
            gene_name=gene.gene_name,
            skipped_exons=skipped,
            affected_ref_transcript=canonical_tx.transcript_id,
            reading_frame_impact=frame_impact,
        ))

    # Detect cryptic exon inclusion
    cryptic = _find_cryptic_exons(iso.exon_blocks, ref_exon_chain, ref_introns, tolerance)
    for ce in cryptic:
        events.append(SplicingEvent(
            event_id="",
            event_type=EventType.CRYPTIC_EXON,
            chrom=iso.chrom,
            start=ce[0],
            end=ce[1],
            strand=iso.strand,
            gene_name=gene.gene_name,
            cryptic_exon=ce,
            affected_ref_transcript=canonical_tx.transcript_id,
            reading_frame_impact=_assess_frame_impact_single(ce),
        ))

    # Detect alternative splice sites
    alt_ss = _find_alt_splice_sites(iso_introns, ref_introns, tolerance)
    events.extend([
        SplicingEvent(
            event_id="",
            event_type=ss_type,
            chrom=iso.chrom,
            start=min(ref_pos, patient_pos),
            end=max(ref_pos, patient_pos),
            strand=iso.strand,
            gene_name=gene.gene_name,
            shifted_donor=(ref_pos, patient_pos) if ss_type == EventType.ALT_5SS else None,
            shifted_acceptor=(ref_pos, patient_pos) if ss_type == EventType.ALT_3SS else None,
            affected_ref_transcript=canonical_tx.transcript_id,
        )
        for ss_type, ref_pos, patient_pos in alt_ss
    ])

    return events


def _detect_intron_retention(
    ci: ClassifiedIsoform,
    gene: Gene,
    canonical_tx: Transcript | None,
) -> list[SplicingEvent]:
    """Detect intron retention events."""
    events = []
    for intron in ci.retained_introns:
        events.append(SplicingEvent(
            event_id="",
            event_type=EventType.INTRON_RETENTION,
            chrom=ci.isoform.chrom,
            start=intron[0],
            end=intron[1],
            strand=ci.isoform.strand,
            gene_name=gene.gene_name,
            retained_intron=intron,
            affected_ref_transcript=canonical_tx.transcript_id if canonical_tx else "",
            reading_frame_impact="likely_frameshift",
        ))
    return events


def _detect_from_nic(
    ci: ClassifiedIsoform,
    gene: Gene,
    canonical_tx: Transcript | None,
    tolerance: int,
) -> list[SplicingEvent]:
    """Detect events from NIC (Novel In Catalog) isoforms — novel junction combos."""
    events = []
    iso = ci.isoform

    if canonical_tx is None:
        return events

    ref_exon_chain = sorted((e.start, e.end) for e in canonical_tx.exons)
    ref_introns = sorted(canonical_tx.introns)
    iso_introns = list(iso.intron_chain)

    # NIC often represents exon skipping with known donors/acceptors
    skipped = _find_skipped_exons(iso_introns, ref_exon_chain, ref_introns, tolerance)
    if skipped:
        events.append(SplicingEvent(
            event_id="",
            event_type=EventType.EXON_SKIPPING,
            chrom=iso.chrom,
            start=skipped[0][0],
            end=skipped[-1][1],
            strand=iso.strand,
            gene_name=gene.gene_name,
            skipped_exons=skipped,
            affected_ref_transcript=canonical_tx.transcript_id,
            reading_frame_impact=_assess_frame_impact(skipped),
        ))

    return events


# ── Helpers ──────────────────────────────────────────────────────

def _find_skipped_exons(
    iso_introns: list[tuple[int, int]],
    ref_exons: list[tuple[int, int]],
    ref_introns: list[tuple[int, int]],
    tolerance: int,
) -> list[tuple[int, int]]:
    """
    Find reference exons that are skipped in patient isoform.

    An exon is "skipped" if the patient has an intron that spans across
    where a reference exon should be.
    """
    skipped = []
    for iso_donor, iso_acceptor in iso_introns:
        # Check if any reference exons fall entirely within this patient intron
        for exon_start, exon_end in ref_exons:
            if exon_start >= iso_donor + tolerance and exon_end <= iso_acceptor - tolerance:
                # Verify this exon is flanked by introns in the reference
                # (i.e., it's an internal exon, not the first/last)
                is_internal = any(
                    abs(ri_end - exon_start) <= tolerance for ri_start, ri_end in ref_introns
                ) and any(
                    abs(ri_start - exon_end) <= tolerance for ri_start, ri_end in ref_introns
                )
                if is_internal:
                    skipped.append((exon_start, exon_end))

    return skipped


def _find_cryptic_exons(
    iso_exons: list[tuple[int, int]],
    ref_exons: list[tuple[int, int]],
    ref_introns: list[tuple[int, int]],
    tolerance: int,
) -> list[tuple[int, int]]:
    """
    Find novel exons in patient that don't correspond to any reference exon.

    A cryptic exon is a patient exon block that falls within a reference intron
    and does not match any known exon.
    """
    cryptic = []
    for iso_start, iso_end in iso_exons:
        # Check if this exon matches any reference exon
        matches_ref = any(
            abs(iso_start - ref_s) <= tolerance and abs(iso_end - ref_e) <= tolerance
            for ref_s, ref_e in ref_exons
        )
        if matches_ref:
            continue

        # Check if it falls within a reference intron (= cryptic exon)
        for intron_start, intron_end in ref_introns:
            if iso_start >= intron_start + tolerance and iso_end <= intron_end - tolerance:
                cryptic.append((iso_start, iso_end))
                break

    return cryptic


def _find_alt_splice_sites(
    iso_introns: list[tuple[int, int]],
    ref_introns: list[tuple[int, int]],
    tolerance: int,
) -> list[tuple[EventType, int, int]]:
    """
    Find alternative 5' and 3' splice sites.

    Detects cases where one end of an intron matches reference but the other
    end is shifted beyond tolerance — indicating use of a nearby cryptic
    splice site.
    """
    alt_sites: list[tuple[EventType, int, int]] = []
    min_shift = tolerance + 1   # must exceed tolerance to be "alternative"
    max_shift = 500             # cap: beyond this it's likely a different event

    for iso_donor, iso_acceptor in iso_introns:
        for ref_donor, ref_acceptor in ref_introns:
            donor_diff = abs(iso_donor - ref_donor)
            acceptor_diff = abs(iso_acceptor - ref_acceptor)

            # Alt 5'SS: acceptor matches, donor shifted
            if acceptor_diff <= tolerance and min_shift < donor_diff < max_shift:
                alt_sites.append((EventType.ALT_5SS, ref_donor, iso_donor))

            # Alt 3'SS: donor matches, acceptor shifted
            if donor_diff <= tolerance and min_shift < acceptor_diff < max_shift:
                alt_sites.append((EventType.ALT_3SS, ref_acceptor, iso_acceptor))

    return alt_sites


def _assess_frame_impact(exons: list[tuple[int, int]]) -> str:
    """Assess reading frame impact of skipped exons."""
    total_length = sum(end - start for start, end in exons)
    if total_length % 3 == 0:
        return "in-frame"
    return "frameshift"


def _assess_frame_impact_single(exon: tuple[int, int]) -> str:
    """Assess reading frame impact of a single cryptic/retained element."""
    length = exon[1] - exon[0]
    if length % 3 == 0:
        return "in-frame"
    return "frameshift"

# splicetarget v0.1.0
# Any usage is subject to this software's license.
