#!/usr/bin/env python3
"""
splicetarget v0.1.0

aso_design.py — ASO/SSO candidate design for splice-switching therapy.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum

import pysam

from splicetarget.splicing.events import EventType, SplicingEvent

logger = logging.getLogger(__name__)


class ASOStrategy(str, Enum):
    """Therapeutic strategy for the ASO candidate."""

    EXON_INCLUSION = "exon_inclusion"    # force skipped exon back in
    EXON_EXCLUSION = "exon_exclusion"    # block cryptic exon inclusion
    SPLICE_SITE_BLOCK = "splice_site_block"  # mask aberrant splice site
    INTRON_SPLICING = "intron_splicing"   # enhance splicing of retained intron

    @property
    def description(self) -> str:
        descs = {
            ASOStrategy.EXON_INCLUSION: "Block exonic splicing silencer to restore exon inclusion",
            ASOStrategy.EXON_EXCLUSION: "Mask cryptic splice site to prevent aberrant exon inclusion",
            ASOStrategy.SPLICE_SITE_BLOCK: "Block activated cryptic splice site to restore canonical splicing",
            ASOStrategy.INTRON_SPLICING: "Mask intron retention boundary to enhance splicing efficiency",
        }
        return descs[self]


@dataclass
class ASOCandidate:
    """A candidate antisense oligonucleotide binding site."""

    candidate_id: str
    event_id: str
    strategy: ASOStrategy

    # Target region on pre-mRNA
    chrom: str
    target_start: int
    target_end: int
    strand: str
    target_sequence: str    # sense strand (pre-mRNA sequence)
    aso_sequence: str       # antisense strand (the ASO itself, 5'→3')

    # Scoring components (0–1 scale, higher = better)
    gc_score: float = 0.0
    self_comp_score: float = 0.0
    accessibility_score: float = 0.0
    splice_proximity_score: float = 0.0
    length_score: float = 0.0

    # Composite
    composite_score: float = 0.0

    # Metadata
    aso_length: int = 0
    gc_content: float = 0.0
    max_self_complementarity: int = 0
    distance_to_splice_site: int = 0
    target_region_type: str = ""  # "exonic", "intronic", "junction"

    off_target_hits: int = -1  # -1 = not yet assessed

    @property
    def is_high_confidence(self) -> bool:
        return self.composite_score >= 0.6 and self.off_target_hits <= 3


# ── Design parameters ────────────────────────────────────────────

@dataclass
class DesignParams:
    """Configurable parameters for ASO candidate generation."""

    min_length: int = 18
    max_length: int = 25
    optimal_length: int = 20
    step_size: int = 1           # sliding window step
    gc_min: float = 0.40
    gc_max: float = 0.60
    gc_optimal: float = 0.50
    max_homopolymer: int = 4     # max consecutive identical bases
    splice_site_flank: int = 300  # bp to scan around splice site
    exon_flank: int = 150         # bp into exon to scan for ESE/ESS elements
    top_n: int = 20              # max candidates per event


# ── Main design function ─────────────────────────────────────────

def design_aso_candidates(
    events: list[SplicingEvent],
    reference_fasta: str,
    params: DesignParams | None = None,
) -> list[ASOCandidate]:
    """
    Generate and score ASO/SSO candidates for each aberrant splicing event.

    Parameters
    ----------
    events : list[SplicingEvent]
        Detected aberrant splicing events.
    reference_fasta : str
        Path to reference genome FASTA (indexed).
    params : DesignParams | None
        Design parameters. Uses defaults if None.

    Returns
    -------
    list[ASOCandidate]
        Scored candidates across all events, sorted by composite score.
    """
    if params is None:
        params = DesignParams()

    ref = pysam.FastaFile(reference_fasta)
    all_candidates: list[ASOCandidate] = []

    for event in events:
        if not event.is_aso_target:
            continue

        strategy, regions = _get_target_regions(event, params)
        candidates = []

        for region_start, region_end, region_type in regions:
            # Extract genomic sequence for the target region
            try:
                seq = ref.fetch(event.chrom, region_start, region_end).upper()
            except (ValueError, KeyError):
                logger.warning("Could not fetch sequence for %s:%d-%d", event.chrom, region_start, region_end)
                continue

            # Sliding window to generate candidate binding sites
            window_candidates = _generate_window_candidates(
                seq=seq,
                chrom=event.chrom,
                region_start=region_start,
                strand=event.strand,
                event=event,
                strategy=strategy,
                region_type=region_type,
                params=params,
            )
            candidates.extend(window_candidates)

        # Score and rank
        for cand in candidates:
            _score_candidate(cand, params)

        # Keep top N per event
        candidates.sort(key=lambda c: c.composite_score, reverse=True)
        top_candidates = candidates[:params.top_n]

        # Assign IDs
        for idx, cand in enumerate(top_candidates, 1):
            cand.candidate_id = f"{event.event_id}.aso{idx:03d}"

        all_candidates.extend(top_candidates)

    ref.close()

    all_candidates.sort(key=lambda c: c.composite_score, reverse=True)
    logger.info("Designed %d ASO candidates across %d events", len(all_candidates), len(events))

    return all_candidates


def _get_target_regions(
    event: SplicingEvent,
    params: DesignParams,
) -> tuple[ASOStrategy, list[tuple[int, int, str]]]:
    """
    Determine ASO strategy and target regions based on event type.

    Returns strategy and list of (start, end, region_type) to scan.
    """
    regions: list[tuple[int, int, str]] = []

    if event.event_type == EventType.EXON_SKIPPING:
        # Strategy: block exonic splicing silencers (ESS) to restore inclusion
        # Target: within skipped exon + flanking intronic regions
        strategy = ASOStrategy.EXON_INCLUSION
        for exon_start, exon_end in event.skipped_exons:
            # Exonic region (ESE/ESS elements)
            regions.append((exon_start, exon_end, "exonic"))
            # Flanking intronic regions (splice site modulation)
            regions.append((exon_start - params.splice_site_flank, exon_start, "intronic_upstream"))
            regions.append((exon_end, exon_end + params.splice_site_flank, "intronic_downstream"))

    elif event.event_type == EventType.CRYPTIC_EXON:
        # Strategy: mask the cryptic splice sites to block inclusion
        strategy = ASOStrategy.EXON_EXCLUSION
        if event.cryptic_exon:
            ce_start, ce_end = event.cryptic_exon
            # The cryptic exon itself
            regions.append((ce_start, ce_end, "cryptic_exonic"))
            # Splice sites flanking the cryptic exon
            regions.append((ce_start - 50, ce_start + 10, "cryptic_3ss"))
            regions.append((ce_end - 10, ce_end + 50, "cryptic_5ss"))

    elif event.event_type == EventType.INTRON_RETENTION:
        # Strategy: enhance splicing by masking retention-promoting elements
        strategy = ASOStrategy.INTRON_SPLICING
        if event.retained_intron:
            ir_start, ir_end = event.retained_intron
            # 5' splice site region
            regions.append((ir_start - 10, ir_start + params.exon_flank, "ir_5ss_region"))
            # 3' splice site region
            regions.append((ir_end - params.exon_flank, ir_end + 10, "ir_3ss_region"))
            # Branch point region (~18-40nt upstream of 3'SS)
            regions.append((ir_end - 50, ir_end - 15, "ir_branch_point"))

    elif event.event_type in (EventType.ALT_5SS, EventType.ALT_3SS):
        # Strategy: block the aberrant splice site
        strategy = ASOStrategy.SPLICE_SITE_BLOCK
        if event.shifted_donor:
            patient_pos = event.shifted_donor[1]
            regions.append((patient_pos - 25, patient_pos + 25, "alt_donor"))
        if event.shifted_acceptor:
            patient_pos = event.shifted_acceptor[1]
            regions.append((patient_pos - 25, patient_pos + 25, "alt_acceptor"))

    else:
        strategy = ASOStrategy.SPLICE_SITE_BLOCK
        regions.append((event.start, event.end, "complex"))

    # Sanitize: ensure no negative coordinates
    regions = [(max(0, s), e, t) for s, e, t in regions]

    return strategy, regions


def _generate_window_candidates(
    seq: str,
    chrom: str,
    region_start: int,
    strand: str,
    event: SplicingEvent,
    strategy: ASOStrategy,
    region_type: str,
    params: DesignParams,
) -> list[ASOCandidate]:
    """Generate candidate ASO sequences using a sliding window."""
    candidates = []

    for length in range(params.min_length, params.max_length + 1):
        for i in range(0, len(seq) - length + 1, params.step_size):
            target_seq = seq[i:i + length]

            # Pre-filter: skip sequences with long homopolymers
            if _has_homopolymer(target_seq, params.max_homopolymer):
                continue

            # Pre-filter: GC content must be in reasonable range (relaxed for tiling)
            gc = _gc_content(target_seq)
            if gc < params.gc_min - 0.10 or gc > params.gc_max + 0.10:
                continue

            aso_seq = _reverse_complement(target_seq)

            # Compute distance to nearest splice site
            target_start = region_start + i
            target_end = target_start + length
            dist = _distance_to_splice_site(target_start, target_end, event)

            candidates.append(ASOCandidate(
                candidate_id="",
                event_id=event.event_id,
                strategy=strategy,
                chrom=chrom,
                target_start=target_start,
                target_end=target_end,
                strand=strand,
                target_sequence=target_seq,
                aso_sequence=aso_seq,
                aso_length=length,
                gc_content=gc,
                max_self_complementarity=_max_self_comp(aso_seq),
                distance_to_splice_site=dist,
                target_region_type=region_type,
            ))

    return candidates


def _score_candidate(cand: ASOCandidate, params: DesignParams) -> None:
    """
    Compute multi-criteria composite score for an ASO candidate.

    Scoring components (each 0–1, higher = better):
        gc_score              — GC content proximity to optimal (0.50)
        self_comp_score       — Low self-complementarity preferred
        splice_proximity_score — Closer to splice site = higher score
        length_score          — Proximity to optimal length (20nt)
        accessibility_score   — Placeholder for RNA folding (default: 0.5)
    """
    # GC content score: bell curve around optimal
    gc_dev = abs(cand.gc_content - params.gc_optimal)
    cand.gc_score = max(0.0, 1.0 - (gc_dev / 0.20))  # 0 at ±20% from optimal

    # Self-complementarity score: penalize high self-comp
    max_comp = cand.max_self_complementarity
    cand.self_comp_score = max(0.0, 1.0 - (max_comp / 8.0))  # 0 at 8bp self-comp

    # Splice proximity score: exponential decay with distance
    dist = cand.distance_to_splice_site
    cand.splice_proximity_score = 1.0 / (1.0 + dist / 50.0)

    # Length score: bell curve around optimal length
    len_dev = abs(cand.aso_length - params.optimal_length)
    cand.length_score = max(0.0, 1.0 - (len_dev / 5.0))

    # RNA accessibility placeholder (would use ViennaRNA for real scoring)
    cand.accessibility_score = 0.5

    # Weighted composite
    cand.composite_score = (
        0.25 * cand.gc_score
        + 0.20 * cand.self_comp_score
        + 0.20 * cand.splice_proximity_score
        + 0.15 * cand.length_score
        + 0.20 * cand.accessibility_score
    )


# ── Sequence utilities ───────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGT", "TGCA")

def _reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def _gc_content(seq: str) -> float:
    """Compute GC fraction of a sequence."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GCgc")
    return gc / len(seq)


def _has_homopolymer(seq: str, max_run: int) -> bool:
    """Check if sequence contains a homopolymer run longer than max_run."""
    count = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            count += 1
            if count > max_run:
                return True
        else:
            count = 1
    return False


def _max_self_comp(seq: str) -> int:
    """
    Compute maximum self-complementarity length.

    Finds the longest stretch where the ASO can fold back on itself
    (form a hairpin). Simple sliding window approach.
    """
    rc = _reverse_complement(seq)
    max_match = 0

    for offset in range(1, len(seq)):
        match_len = 0
        for i in range(len(seq) - offset):
            if i + offset < len(rc) and seq[i] == rc[i + offset]:
                match_len += 1
                max_match = max(max_match, match_len)
            else:
                match_len = 0

    return max_match


def _distance_to_splice_site(start: int, end: int, event: SplicingEvent) -> int:
    """Compute minimum distance from ASO target to the nearest splice site."""
    sites = []

    if event.skipped_exons:
        for es, ee in event.skipped_exons:
            sites.extend([es, ee])
    if event.cryptic_exon:
        sites.extend(event.cryptic_exon)
    if event.retained_intron:
        sites.extend(event.retained_intron)
    if event.shifted_donor:
        sites.append(event.shifted_donor[1])
    if event.shifted_acceptor:
        sites.append(event.shifted_acceptor[1])

    if not sites:
        return 9999

    mid = (start + end) // 2
    return min(abs(mid - s) for s in sites)

# splicetarget v0.1.0
# Any usage is subject to this software's license.
