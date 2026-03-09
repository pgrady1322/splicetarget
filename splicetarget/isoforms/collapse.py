#!/usr/bin/env python3
"""
splicetarget v0.1.0

collapse.py — Isoform collapsing — merge redundant long reads into unique isoforms.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field

from splicetarget.data.io import ReadRecord

logger = logging.getLogger(__name__)


@dataclass
class CollapsedIsoform:
    """A unique isoform defined by its splice junction chain."""

    isoform_id: str
    chrom: str
    strand: str
    exon_blocks: list[tuple[int, int]]
    supporting_reads: list[str] = field(default_factory=list)
    representative_read: str = ""

    @property
    def read_count(self) -> int:
        return len(self.supporting_reads)

    @property
    def start(self) -> int:
        return self.exon_blocks[0][0] if self.exon_blocks else 0

    @property
    def end(self) -> int:
        return self.exon_blocks[-1][1] if self.exon_blocks else 0

    @property
    def intron_chain(self) -> tuple[tuple[int, int], ...]:
        """The defining feature: ordered intron (donor, acceptor) pairs."""
        introns = []
        for i in range(len(self.exon_blocks) - 1):
            introns.append((self.exon_blocks[i][1], self.exon_blocks[i + 1][0]))
        return tuple(introns)

    @property
    def exon_count(self) -> int:
        return len(self.exon_blocks)

    @property
    def total_exon_length(self) -> int:
        return sum(end - start for start, end in self.exon_blocks)

    def to_bed12(self) -> str:
        """Export isoform as BED12 format line."""
        block_sizes = ",".join(str(e - s) for s, e in self.exon_blocks)
        block_starts = ",".join(str(s - self.start) for s, _ in self.exon_blocks)
        return (
            f"{self.chrom}\t{self.start}\t{self.end}\t{self.isoform_id}\t"
            f"{self.read_count}\t{self.strand}\t{self.start}\t{self.end}\t"
            f"0,0,0\t{self.exon_count}\t{block_sizes}\t{block_starts}"
        )


def collapse_isoforms(
    reads: list[ReadRecord],
    junction_tolerance: int = 10,
    end_tolerance: int = 50,
    min_reads: int = 2,
    gene_name: str = "GENE",
) -> list[CollapsedIsoform]:
    """
    Collapse long reads into unique isoforms based on splice junction chains.

    Two reads are assigned to the same isoform if:
    1. All internal splice junctions match within `junction_tolerance` bp
    2. TSS and TTS (first/last exon ends) are within `end_tolerance` bp

    Parameters
    ----------
    reads : list[ReadRecord]
        Mapped, spliced reads with populated exon_blocks.
    junction_tolerance : int
        Maximum bp difference at each splice junction to still merge (default: 10).
    end_tolerance : int
        Maximum bp difference at 5'/3' transcript ends (default: 50).
    min_reads : int
        Minimum supporting reads to report an isoform (default: 2).
    gene_name : str
        Gene name prefix for isoform IDs.

    Returns
    -------
    list[CollapsedIsoform]
        Collapsed isoforms sorted by read count (descending).
    """
    if not reads:
        return []

    # Filter to multi-exon reads only
    spliced_reads = [r for r in reads if len(r.exon_blocks) >= 2 and r.is_mapped]
    mono_exon_reads = [r for r in reads if len(r.exon_blocks) == 1 and r.is_mapped]

    logger.info(
        "Collapsing %d spliced reads + %d mono-exon reads",
        len(spliced_reads), len(mono_exon_reads),
    )

    # Group by intron chain (with tolerance)
    clusters: list[list[ReadRecord]] = []

    for read in spliced_reads:
        read_introns = _get_intron_chain(read)
        placed = False

        for cluster in clusters:
            ref_introns = _get_intron_chain(cluster[0])
            if _chains_match(read_introns, ref_introns, junction_tolerance, end_tolerance, read, cluster[0]):
                cluster.append(read)
                placed = True
                break

        if not placed:
            clusters.append([read])

    # Handle mono-exon reads as separate isoforms (group by overlap)
    mono_clusters = _cluster_mono_exon(mono_exon_reads, end_tolerance)
    clusters.extend(mono_clusters)

    # Build CollapsedIsoform objects
    isoforms = []
    iso_idx = 0
    for cluster in clusters:
        if len(cluster) < min_reads:
            continue

        iso_idx += 1
        # Use the read with median length as representative
        sorted_by_len = sorted(cluster, key=lambda r: r.length)
        representative = sorted_by_len[len(sorted_by_len) // 2]

        # Consensus exon boundaries: median of each exon start/end
        consensus_exons = _consensus_exon_blocks(cluster)

        isoforms.append(CollapsedIsoform(
            isoform_id=f"{gene_name}.iso{iso_idx:04d}",
            chrom=representative.chrom,
            strand=representative.strand,
            exon_blocks=consensus_exons,
            supporting_reads=[r.name for r in cluster],
            representative_read=representative.name,
        ))

    isoforms.sort(key=lambda iso: iso.read_count, reverse=True)
    logger.info("Collapsed into %d isoforms (min_reads=%d)", len(isoforms), min_reads)
    return isoforms


def _get_intron_chain(read: ReadRecord) -> list[tuple[int, int]]:
    """Extract intron chain from a read's exon blocks."""
    introns = []
    for i in range(len(read.exon_blocks) - 1):
        introns.append((read.exon_blocks[i][1], read.exon_blocks[i + 1][0]))
    return introns


def _chains_match(
    chain_a: list[tuple[int, int]],
    chain_b: list[tuple[int, int]],
    junction_tol: int,
    end_tol: int,
    read_a: ReadRecord,
    read_b: ReadRecord,
) -> bool:
    """Check if two intron chains match within tolerance."""
    if len(chain_a) != len(chain_b):
        return False

    # Check internal junctions
    for (a_start, a_end), (b_start, b_end) in zip(chain_a, chain_b):
        if abs(a_start - b_start) > junction_tol:
            return False
        if abs(a_end - b_end) > junction_tol:
            return False

    # Check transcript ends (5' and 3')
    if abs(read_a.start - read_b.start) > end_tol:
        return False
    if abs(read_a.end - read_b.end) > end_tol:
        return False

    return True


def _cluster_mono_exon(
    reads: list[ReadRecord],
    tolerance: int,
) -> list[list[ReadRecord]]:
    """Cluster mono-exon reads by coordinate overlap."""
    if not reads:
        return []

    sorted_reads = sorted(reads, key=lambda r: (r.chrom, r.start))
    clusters: list[list[ReadRecord]] = [[sorted_reads[0]]]

    for read in sorted_reads[1:]:
        last_cluster = clusters[-1]
        ref = last_cluster[0]

        if (read.chrom == ref.chrom
                and abs(read.start - ref.start) <= tolerance
                and abs(read.end - ref.end) <= tolerance):
            last_cluster.append(read)
        else:
            clusters.append([read])

    return clusters


def _consensus_exon_blocks(cluster: list[ReadRecord]) -> list[tuple[int, int]]:
    """
    Compute consensus exon boundaries from a cluster of reads.

    Uses median of each exon start/end across all reads in the cluster.
    """
    n_exons = len(cluster[0].exon_blocks)

    # Group positions by exon index
    starts: dict[int, list[int]] = defaultdict(list)
    ends: dict[int, list[int]] = defaultdict(list)

    for read in cluster:
        for i, (s, e) in enumerate(read.exon_blocks):
            if i < n_exons:
                starts[i].append(s)
                ends[i].append(e)

    consensus = []
    for i in range(n_exons):
        if i in starts and i in ends:
            med_start = sorted(starts[i])[len(starts[i]) // 2]
            med_end = sorted(ends[i])[len(ends[i]) // 2]
            consensus.append((med_start, med_end))

    return consensus

# splicetarget v0.1.0
# Any usage is subject to this software's license.
