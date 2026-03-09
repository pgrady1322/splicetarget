#!/usr/bin/env python3
"""
splicetarget v0.1.0

aligner.py — minimap2-based long-read splice-aware alignment.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

logger = logging.getLogger(__name__)


class ReadType(str, Enum):
    """Supported long-read sequencing platforms and data types."""

    ISOSEQ = "isoseq"          # PacBio IsoSeq (HiFi CCS reads)
    PACBIO_CDNA = "pb_cdna"    # PacBio cDNA (non-IsoSeq)
    ONT_DRNA = "ont_drna"      # ONT direct RNA
    ONT_CDNA = "ont_cdna"      # ONT cDNA (PCR or direct)

    @property
    def minimap2_preset(self) -> str:
        """Map read type to minimap2 splice-aware preset."""
        presets = {
            ReadType.ISOSEQ: "splice:hq",       # high-quality CCS
            ReadType.PACBIO_CDNA: "splice:hq",
            ReadType.ONT_DRNA: "splice",          # noisier ONT reads
            ReadType.ONT_CDNA: "splice",
        }
        return presets[self]

    @property
    def expected_error_rate(self) -> float:
        """Typical per-base error rate for this read type."""
        rates = {
            ReadType.ISOSEQ: 0.001,
            ReadType.PACBIO_CDNA: 0.005,
            ReadType.ONT_DRNA: 0.05,
            ReadType.ONT_CDNA: 0.03,
        }
        return rates[self]


@dataclass
class AlignmentStats:
    """Summary statistics from a minimap2 alignment run."""

    total_reads: int = 0
    mapped_reads: int = 0
    spliced_reads: int = 0          # reads with ≥1 intron (N cigar)
    mean_mapping_quality: float = 0.0
    median_read_length: int = 0

    @property
    def mapping_rate(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return self.mapped_reads / self.total_reads

    @property
    def splice_rate(self) -> float:
        """Fraction of mapped reads that are spliced."""
        if self.mapped_reads == 0:
            return 0.0
        return self.spliced_reads / self.mapped_reads


def align_long_reads(
    reads_path: str | Path,
    reference_path: str | Path,
    output_bam: str | Path,
    read_type: ReadType = ReadType.ISOSEQ,
    junction_bed: str | Path | None = None,
    threads: int = 4,
    extra_args: list[str] | None = None,
) -> Path:
    """
    Align long-read RNA sequences to a reference genome using minimap2.

    Produces a coordinate-sorted, indexed BAM file with splice-aware alignment.

    Parameters
    ----------
    reads_path : str | Path
        Input reads (FASTA, FASTQ, or unaligned BAM).
    reference_path : str | Path
        Reference genome FASTA (must be indexed with .fai).
    output_bam : str | Path
        Path for output sorted BAM file.
    read_type : ReadType
        Sequencing platform/data type (determines minimap2 preset).
    junction_bed : str | Path | None
        Optional BED file of known splice junctions to guide alignment.
        Strongly recommended for improving junction accuracy.
    threads : int
        Number of CPU threads for minimap2 and samtools.
    extra_args : list[str] | None
        Additional command-line arguments to pass to minimap2.

    Returns
    -------
    Path
        Path to the sorted, indexed output BAM file.
    """
    reads_path = Path(reads_path)
    reference_path = Path(reference_path)
    output_bam = Path(output_bam)

    _check_dependencies()
    _validate_inputs(reads_path, reference_path)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Build minimap2 command
    cmd = [
        "minimap2",
        "-ax", read_type.minimap2_preset,
        "-t", str(threads),
        "--secondary=no",           # suppress secondary alignments
        "-uf",                      # transcript strand: forward
        "--MD",                     # include MD tag for variant calling
    ]

    # Add known junction guidance if provided
    if junction_bed is not None:
        junction_bed = Path(junction_bed)
        if junction_bed.exists():
            cmd.extend(["--junc-bed", str(junction_bed)])
            logger.info("Using known junctions from %s", junction_bed)

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend([str(reference_path), str(reads_path)])

    # Pipe minimap2 → samtools sort → BAM
    unsorted_sam = output_bam.with_suffix(".unsorted.sam")

    logger.info(
        "Aligning %s reads (%s preset) to %s",
        read_type.value, read_type.minimap2_preset, reference_path.name,
    )
    logger.info("Command: %s", " ".join(cmd))

    # minimap2 alignment
    with open(unsorted_sam, "w") as sam_out:
        result = subprocess.run(
            cmd,
            stdout=sam_out,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    if result.returncode != 0:
        logger.error("minimap2 failed:\n%s", result.stderr)
        raise RuntimeError(f"minimap2 alignment failed: {result.stderr[:500]}")

    # Parse minimap2 stderr for stats
    _log_minimap2_stats(result.stderr)

    # Sort and index with samtools
    _sort_and_index(unsorted_sam, output_bam, threads)

    # Cleanup intermediate SAM
    unsorted_sam.unlink(missing_ok=True)

    logger.info("Alignment complete: %s", output_bam)
    return output_bam


def _sort_and_index(sam_path: Path, bam_path: Path, threads: int) -> None:
    """Sort SAM → BAM and create BAI index."""
    logger.info("Sorting alignment → %s", bam_path)

    sort_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-o", str(bam_path),
        str(sam_path),
    ]
    subprocess.run(sort_cmd, check=True, capture_output=True, text=True)

    index_cmd = ["samtools", "index", "-@", str(threads), str(bam_path)]
    subprocess.run(index_cmd, check=True, capture_output=True, text=True)

    logger.info("Indexed: %s.bai", bam_path)


def compute_alignment_stats(bam_path: str | Path) -> AlignmentStats:
    """
    Compute alignment summary statistics from a BAM file using samtools.

    Parameters
    ----------
    bam_path : str | Path
        Path to an indexed BAM file.

    Returns
    -------
    AlignmentStats
    """
    import pysam

    bam_path = Path(bam_path)
    stats = AlignmentStats()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        mapqs = []
        lengths = []

        for read in bam.fetch(until_eof=True):
            stats.total_reads += 1

            if read.is_unmapped:
                continue

            if read.is_secondary or read.is_supplementary:
                continue

            stats.mapped_reads += 1
            mapqs.append(read.mapping_quality)
            lengths.append(read.query_length or 0)

            # Check for spliced (N in CIGAR)
            if read.cigartuples and any(op == 3 for op, _ in read.cigartuples):
                stats.spliced_reads += 1

        if mapqs:
            stats.mean_mapping_quality = sum(mapqs) / len(mapqs)
        if lengths:
            sorted_lens = sorted(lengths)
            mid = len(sorted_lens) // 2
            stats.median_read_length = sorted_lens[mid]

    return stats


# ── Helpers ──────────────────────────────────────────────────────

def _check_dependencies() -> None:
    """Verify minimap2 and samtools are available on PATH."""
    for tool in ("minimap2", "samtools"):
        if shutil.which(tool) is None:
            raise EnvironmentError(
                f"'{tool}' not found on PATH. Install via: conda install -c bioconda {tool}"
            )


def _validate_inputs(reads_path: Path, reference_path: Path) -> None:
    """Check that input files exist."""
    if not reads_path.exists():
        raise FileNotFoundError(f"Reads file not found: {reads_path}")
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference file not found: {reference_path}")


def _log_minimap2_stats(stderr: str) -> None:
    """Parse and log key statistics from minimap2 stderr output."""
    for line in stderr.strip().split("\n"):
        if "mapped" in line.lower() or "peak" in line.lower():
            logger.info("[minimap2] %s", line.strip())

# splicetarget v0.1.0
# Any usage is subject to this software's license.
