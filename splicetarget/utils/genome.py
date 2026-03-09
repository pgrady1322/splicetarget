#!/usr/bin/env python3
"""
splicetarget v0.1.0

genome.py — Genome utility functions — coordinates, regions, sequence helpers.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass
class GenomicRegion:
    """A genomic interval."""

    chrom: str
    start: int
    end: int
    strand: str = "+"
    name: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2

    def overlaps(self, other: GenomicRegion) -> bool:
        """Check if two regions overlap on the same chromosome."""
        if self.chrom != other.chrom:
            return False
        return self.start < other.end and other.start < self.end

    def contains(self, other: GenomicRegion) -> bool:
        """Check if self fully contains other."""
        return (
            self.chrom == other.chrom
            and self.start <= other.start
            and self.end >= other.end
        )

    def distance_to(self, other: GenomicRegion) -> int:
        """Compute distance between non-overlapping regions (0 if overlap)."""
        if self.overlaps(other):
            return 0
        if self.end <= other.start:
            return other.start - self.end
        return self.start - other.end

    def expand(self, bp: int) -> GenomicRegion:
        """Return a new region expanded by bp in both directions."""
        return GenomicRegion(
            chrom=self.chrom,
            start=max(0, self.start - bp),
            end=self.end + bp,
            strand=self.strand,
            name=self.name,
        )

    def to_region_string(self) -> str:
        """Format as 'chr:start-end'."""
        return f"{self.chrom}:{self.start}-{self.end}"

    def to_bed(self) -> str:
        """Format as BED line."""
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t0\t{self.strand}"


def parse_region(region_str: str) -> GenomicRegion:
    """
    Parse a genomic region string into a GenomicRegion.

    Accepts formats:
        - chr1:1000-2000
        - chr1:1,000-2,000
        - chr1:1000-2000:+
    """
    region_str = region_str.replace(",", "")
    match = re.match(r"^(\w+):(\d+)-(\d+)(?::([+-]))?$", region_str)
    if not match:
        raise ValueError(f"Cannot parse region string: '{region_str}'")
    return GenomicRegion(
        chrom=match.group(1),
        start=int(match.group(2)),
        end=int(match.group(3)),
        strand=match.group(4) or "+",
    )


def reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def gc_content(seq: str) -> float:
    """Compute GC fraction."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in "GC")
    return gc / len(seq)


def format_bp(bp: int) -> str:
    """Human-readable base pair count: 1234567 → '1.23 Mb'."""
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.2f} Mb"
    elif bp >= 1_000:
        return f"{bp / 1_000:.1f} kb"
    return f"{bp} bp"

# splicetarget v0.1.0
# Any usage is subject to this software's license.
