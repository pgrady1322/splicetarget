#!/usr/bin/env python3
"""
splicetarget v0.1.0

io.py — I/O helpers for sequencing data (BAM, FASTQ, FASTA) and annotations.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import pysam

logger = logging.getLogger(__name__)


# ── Data classes ──────────────────────────────────────────────────

@dataclass
class ReadRecord:
    """A single long-read transcript record."""

    name: str
    sequence: str
    quality: str | None
    chrom: str | None = None
    start: int | None = None
    end: int | None = None
    strand: str | None = None
    cigar: str | None = None
    mapping_quality: int = 0
    is_mapped: bool = False
    exon_blocks: list[tuple[int, int]] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def intron_blocks(self) -> list[tuple[int, int]]:
        """Derive introns from consecutive exon blocks."""
        introns = []
        for i in range(len(self.exon_blocks) - 1):
            intron_start = self.exon_blocks[i][1]
            intron_end = self.exon_blocks[i + 1][0]
            introns.append((intron_start, intron_end))
        return introns


@dataclass
class SpliceJunction:
    """A splice junction extracted from a long-read alignment."""

    chrom: str
    start: int       # donor (5' end of intron)
    end: int         # acceptor (3' end of intron)
    strand: str
    read_count: int = 1
    canonical: bool = True
    donor_motif: str = ""     # e.g., "GT"
    acceptor_motif: str = ""  # e.g., "AG"

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def key(self) -> tuple[str, int, int, str]:
        return (self.chrom, self.start, self.end, self.strand)


# ── BAM readers ──────────────────────────────────────────────────

def iter_aligned_reads(
    bam_path: str | Path,
    region: str | None = None,
    min_mapq: int = 0,
    require_splice: bool = False,
) -> Iterator[ReadRecord]:
    """
    Iterate over aligned reads from a BAM file.

    Parameters
    ----------
    bam_path : str | Path
        Path to an indexed BAM file.
    region : str | None
        Genomic region string (e.g., "chr1:1000-2000"). None for whole file.
    min_mapq : int
        Minimum mapping quality filter (default: 0).
    require_splice : bool
        If True, only yield reads containing 'N' CIGAR ops (spliced reads).

    Yields
    ------
    ReadRecord
    """
    bam_path = Path(bam_path)
    _validate_bam(bam_path)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        fetch_kwargs = {"region": region} if region else {}

        for read in bam.fetch(**fetch_kwargs):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            exon_blocks = _extract_exon_blocks(read)

            if require_splice and len(exon_blocks) < 2:
                continue

            strand = "-" if read.is_reverse else "+"

            yield ReadRecord(
                name=read.query_name,
                sequence=read.query_sequence or "",
                quality=_phred_to_str(read.query_qualities),
                chrom=read.reference_name,
                start=read.reference_start,
                end=read.reference_end,
                strand=strand,
                cigar=read.cigarstring,
                mapping_quality=read.mapping_quality,
                is_mapped=True,
                exon_blocks=exon_blocks,
            )


def extract_splice_junctions(
    bam_path: str | Path,
    region: str | None = None,
    min_mapq: int = 20,
    reference_fasta: str | Path | None = None,
) -> list[SpliceJunction]:
    """
    Extract and deduplicate splice junctions from a BAM file.

    For each junction, checks donor/acceptor dinucleotide motifs against the
    reference to assess canonicality (GT-AG, GC-AG, AT-AC).

    Parameters
    ----------
    bam_path : str | Path
        Path to an indexed BAM file.
    region : str | None
        Genomic region filter.
    min_mapq : int
        Minimum MAPQ for reads to consider.
    reference_fasta : str | Path | None
        Reference FASTA for motif extraction. If None, motifs are left empty.

    Returns
    -------
    list[SpliceJunction]
        Deduplicated junctions sorted by (chrom, start).
    """
    junction_counts: dict[tuple[str, int, int, str], int] = {}

    for read in iter_aligned_reads(bam_path, region=region, min_mapq=min_mapq, require_splice=True):
        for intron_start, intron_end in read.intron_blocks:
            key = (read.chrom, intron_start, intron_end, read.strand)
            junction_counts[key] = junction_counts.get(key, 0) + 1

    # Build junction objects with optional motif annotation
    ref = None
    if reference_fasta is not None:
        ref = pysam.FastaFile(str(reference_fasta))

    junctions = []
    for (chrom, start, end, strand), count in junction_counts.items():
        donor_motif = ""
        acceptor_motif = ""
        canonical = True

        if ref is not None:
            try:
                donor_motif = ref.fetch(chrom, start, start + 2).upper()
                acceptor_motif = ref.fetch(chrom, end - 2, end).upper()
                canonical = _is_canonical_splice(donor_motif, acceptor_motif)
            except (ValueError, KeyError):
                logger.warning("Could not fetch motif for %s:%d-%d", chrom, start, end)

        junctions.append(SpliceJunction(
            chrom=chrom,
            start=start,
            end=end,
            strand=strand,
            read_count=count,
            canonical=canonical,
            donor_motif=donor_motif,
            acceptor_motif=acceptor_motif,
        ))

    if ref is not None:
        ref.close()

    junctions.sort(key=lambda j: (j.chrom, j.start))
    logger.info("Extracted %d unique splice junctions from %s", len(junctions), bam_path)
    return junctions


# ── FASTA / FASTQ readers ────────────────────────────────────────

def iter_fasta_records(fasta_path: str | Path) -> Iterator[ReadRecord]:
    """Iterate over records in a FASTA file."""
    with pysam.FastxFile(str(fasta_path)) as fh:
        for entry in fh:
            yield ReadRecord(
                name=entry.name,
                sequence=entry.sequence,
                quality=entry.quality,
            )


def iter_fastq_records(fastq_path: str | Path) -> Iterator[ReadRecord]:
    """Iterate over records in a FASTQ file."""
    return iter_fasta_records(fastq_path)  # pysam.FastxFile handles both


# ── Internal helpers ─────────────────────────────────────────────

def _validate_bam(bam_path: Path) -> None:
    """Ensure BAM file exists and has an index."""
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    index_paths = [
        bam_path.with_suffix(".bam.bai"),
        bam_path.parent / (bam_path.name + ".bai"),
    ]
    if not any(p.exists() for p in index_paths):
        logger.info("No BAM index found; creating index for %s", bam_path)
        pysam.index(str(bam_path))


def _extract_exon_blocks(read: pysam.AlignedSegment) -> list[tuple[int, int]]:
    """
    Extract exon coordinate blocks from a read's CIGAR string.

    Treats M/=/X as exonic, N as intron (splice), D as deletion within exon.
    """
    if read.cigartuples is None:
        return []

    blocks = []
    pos = read.reference_start
    block_start = pos

    # CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8
    CONSUME_REF = {0, 2, 3, 7, 8}  # M, D, N, =, X
    INTRON = 3  # N

    for op, length in read.cigartuples:
        if op == INTRON:
            # Close current exon block at intron boundary
            if block_start < pos:
                blocks.append((block_start, pos))
            pos += length
            block_start = pos
        elif op in CONSUME_REF:
            pos += length
        # I, S, H, P don't consume reference

    # Close final exon block
    if block_start < pos:
        blocks.append((block_start, pos))

    return blocks


def _phred_to_str(qualities: pysam.array_type | None) -> str | None:
    """Convert pysam quality array to a Phred+33 ASCII string."""
    if qualities is None:
        return None
    return "".join(chr(q + 33) for q in qualities)


def _is_canonical_splice(donor: str, acceptor: str) -> bool:
    """Check if splice motif is canonical (GT-AG, GC-AG, AT-AC)."""
    canonical_pairs = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
    return (donor, acceptor) in canonical_pairs

# splicetarget v0.1.0
# Any usage is subject to this software's license.
