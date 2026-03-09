#!/usr/bin/env python3
"""
splicetarget v0.1.0

classify.py — SQANTI3-style isoform structural classification.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from enum import Enum

from splicetarget.data.reference import Gene, Transcript
from splicetarget.isoforms.collapse import CollapsedIsoform

logger = logging.getLogger(__name__)


class IsoformCategory(str, Enum):
    """Structural classification categories (SQANTI3-style)."""

    FSM = "FSM"   # Full Splice Match
    ISM = "ISM"   # Incomplete Splice Match
    NIC = "NIC"   # Novel In Catalog
    NNC = "NNC"   # Novel Not in Catalog
    GI = "GI"     # Genic Intron
    GG = "GG"     # Genic Genomic
    AS = "AS"     # Antisense
    IR = "IR"     # Intron Retention
    FUSION = "FUSION"

    @property
    def is_known(self) -> bool:
        return self in (IsoformCategory.FSM, IsoformCategory.ISM)

    @property
    def is_novel(self) -> bool:
        return self in (IsoformCategory.NIC, IsoformCategory.NNC, IsoformCategory.IR)

    @property
    def is_aberrant(self) -> bool:
        """Categories likely arising from aberrant splicing (therapeutic targets)."""
        return self in (IsoformCategory.NNC, IsoformCategory.IR)


@dataclass
class ClassifiedIsoform:
    """An isoform with structural classification relative to reference."""

    isoform: CollapsedIsoform
    category: IsoformCategory
    best_ref_transcript: str = ""
    subcategory: str = ""
    junction_match_pct: float = 0.0
    novel_junctions: list[tuple[int, int]] = None
    retained_introns: list[tuple[int, int]] = None
    novel_donors: list[int] = None
    novel_acceptors: list[int] = None

    def __post_init__(self) -> None:
        if self.novel_junctions is None:
            self.novel_junctions = []
        if self.retained_introns is None:
            self.retained_introns = []
        if self.novel_donors is None:
            self.novel_donors = []
        if self.novel_acceptors is None:
            self.novel_acceptors = []

    @property
    def is_therapy_candidate(self) -> bool:
        """Whether this isoform is a candidate for splice-switching therapy."""
        return self.category.is_aberrant and self.isoform.read_count >= 2


def classify_isoforms(
    isoforms: list[CollapsedIsoform],
    gene: Gene,
    junction_tolerance: int = 10,
) -> list[ClassifiedIsoform]:
    """
    Classify collapsed isoforms against a reference gene model.

    Parameters
    ----------
    isoforms : list[CollapsedIsoform]
        Collapsed patient isoforms.
    gene : Gene
        Reference gene with known transcripts.
    junction_tolerance : int
        Maximum bp wobble for junction matching (default: 10).

    Returns
    -------
    list[ClassifiedIsoform]
        Classified isoforms, same order as input.
    """
    ref_transcripts = list(gene.transcripts.values())

    # Pre-compute reference junction sets
    all_ref_donors: set[int] = set()
    all_ref_acceptors: set[int] = set()
    for tx in ref_transcripts:
        for donor, acceptor in tx.splice_junctions:
            all_ref_donors.add(donor)
            all_ref_acceptors.add(acceptor)

    classified = []
    for iso in isoforms:
        result = _classify_single(iso, ref_transcripts, all_ref_donors, all_ref_acceptors, junction_tolerance)
        classified.append(result)

    # Summary
    category_counts = {}
    for c in classified:
        category_counts[c.category.value] = category_counts.get(c.category.value, 0) + 1
    logger.info("Classification summary: %s", category_counts)

    return classified


def _classify_single(
    iso: CollapsedIsoform,
    ref_transcripts: list[Transcript],
    all_ref_donors: set[int],
    all_ref_acceptors: set[int],
    tolerance: int,
) -> ClassifiedIsoform:
    """Classify a single isoform against reference transcripts."""
    iso_introns = set(iso.intron_chain)

    # Handle mono-exon isoforms
    if not iso_introns:
        return _classify_mono_exon(iso, ref_transcripts)

    # Check for intron retention FIRST — an isoform whose exon block spans a
    # reference intron is IR even if its remaining junctions subset-match (ISM).
    retained = _detect_intron_retention(iso, ref_transcripts, tolerance)

    # Check against each reference transcript
    best_match_tx = ""
    best_match_pct = 0.0
    best_category = IsoformCategory.GG

    for tx in ref_transcripts:
        ref_introns = tx.splice_junctions
        if not ref_introns:
            continue

        matched, total = _count_junction_matches(iso_introns, ref_introns, tolerance)
        match_pct = matched / max(len(iso_introns), 1)

        if match_pct > best_match_pct:
            best_match_pct = match_pct
            best_match_tx = tx.transcript_id

        # FSM: all isoform junctions match all reference junctions
        if matched == len(iso_introns) == len(ref_introns):
            return ClassifiedIsoform(
                isoform=iso,
                category=IsoformCategory.FSM,
                best_ref_transcript=tx.transcript_id,
                junction_match_pct=1.0,
            )

        # ISM: all isoform junctions are a subset of reference
        # BUT only if no introns are retained (otherwise it's IR)
        if matched == len(iso_introns) and len(iso_introns) < len(ref_introns) and not retained:
            return ClassifiedIsoform(
                isoform=iso,
                category=IsoformCategory.ISM,
                best_ref_transcript=tx.transcript_id,
                subcategory=f"{matched}/{len(ref_introns)} junctions",
                junction_match_pct=match_pct,
            )

    # Intron retention takes priority when an exon block spans a reference intron
    if retained:
        return ClassifiedIsoform(
            isoform=iso,
            category=IsoformCategory.IR,
            best_ref_transcript=best_match_tx,
            junction_match_pct=best_match_pct,
            retained_introns=retained,
        )

    # NIC vs NNC: check if novel junctions use known donors/acceptors
    novel_junctions = []
    novel_donors = []
    novel_acceptors = []

    for donor, acceptor in iso_introns:
        donor_known = _position_in_set(donor, all_ref_donors, tolerance)
        acceptor_known = _position_in_set(acceptor, all_ref_acceptors, tolerance)

        if not donor_known or not acceptor_known:
            novel_junctions.append((donor, acceptor))
            if not donor_known:
                novel_donors.append(donor)
            if not acceptor_known:
                novel_acceptors.append(acceptor)

    if novel_junctions:
        # NNC: at least one completely novel splice site
        if novel_donors or novel_acceptors:
            return ClassifiedIsoform(
                isoform=iso,
                category=IsoformCategory.NNC,
                best_ref_transcript=best_match_tx,
                junction_match_pct=best_match_pct,
                novel_junctions=novel_junctions,
                novel_donors=novel_donors,
                novel_acceptors=novel_acceptors,
            )
    else:
        # NIC: novel combination of known donors and acceptors
        return ClassifiedIsoform(
            isoform=iso,
            category=IsoformCategory.NIC,
            best_ref_transcript=best_match_tx,
            junction_match_pct=best_match_pct,
        )

    return ClassifiedIsoform(
        isoform=iso,
        category=IsoformCategory.GG,
        best_ref_transcript=best_match_tx,
        junction_match_pct=best_match_pct,
    )


def _classify_mono_exon(
    iso: CollapsedIsoform,
    ref_transcripts: list[Transcript],
) -> ClassifiedIsoform:
    """Classify a mono-exon isoform."""
    # Check if it falls within a reference intron (GI)
    for tx in ref_transcripts:
        for intron_start, intron_end in tx.introns:
            if iso.start >= intron_start and iso.end <= intron_end:
                return ClassifiedIsoform(
                    isoform=iso,
                    category=IsoformCategory.GI,
                    best_ref_transcript=tx.transcript_id,
                    subcategory="mono-exon within intron",
                )

    return ClassifiedIsoform(
        isoform=iso,
        category=IsoformCategory.GG,
        subcategory="mono-exon",
    )


def _count_junction_matches(
    query: set[tuple[int, int]],
    reference: set[tuple[int, int]],
    tolerance: int,
) -> tuple[int, int]:
    """Count how many query junctions match reference junctions within tolerance."""
    matched = 0
    for q_donor, q_acceptor in query:
        for r_donor, r_acceptor in reference:
            if abs(q_donor - r_donor) <= tolerance and abs(q_acceptor - r_acceptor) <= tolerance:
                matched += 1
                break
    return matched, len(query)


def _position_in_set(pos: int, ref_set: set[int], tolerance: int) -> bool:
    """Check if a position is within tolerance of any position in the reference set."""
    return any(abs(pos - ref_pos) <= tolerance for ref_pos in ref_set)


def _detect_intron_retention(
    iso: CollapsedIsoform,
    ref_transcripts: list[Transcript],
    tolerance: int,
) -> list[tuple[int, int]]:
    """
    Detect reference introns that are retained (not spliced out) in the isoform.

    An intron is "retained" if a patient exon block spans across where the
    reference has a donor→acceptor junction.
    """
    retained = []

    for tx in ref_transcripts:
        for intron_start, intron_end in tx.introns:
            for exon_start, exon_end in iso.exon_blocks:
                # Exon block spans the entire intron
                if (exon_start <= intron_start - tolerance
                        and exon_end >= intron_end + tolerance):
                    retained.append((intron_start, intron_end))
                    break

    return list(set(retained))

# splicetarget v0.1.0
# Any usage is subject to this software's license.
