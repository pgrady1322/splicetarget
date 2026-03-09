#!/usr/bin/env python3
"""
splicetarget v0.1.0

test_utils.py — Tests for genome utility functions.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import pytest

from splicetarget.utils.genome import (
    GenomicRegion,
    format_bp,
    gc_content,
    parse_region,
    reverse_complement,
)


class TestGenomicRegion:
    """Tests for GenomicRegion dataclass."""

    def test_length(self):
        r = GenomicRegion("chr1", 100, 200)
        assert r.length == 100

    def test_overlaps(self):
        a = GenomicRegion("chr1", 100, 200)
        b = GenomicRegion("chr1", 150, 250)
        c = GenomicRegion("chr1", 300, 400)
        d = GenomicRegion("chr2", 100, 200)

        assert a.overlaps(b) is True
        assert a.overlaps(c) is False
        assert a.overlaps(d) is False  # different chrom

    def test_contains(self):
        outer = GenomicRegion("chr1", 100, 500)
        inner = GenomicRegion("chr1", 200, 300)
        assert outer.contains(inner) is True
        assert inner.contains(outer) is False

    def test_distance_to(self):
        a = GenomicRegion("chr1", 100, 200)
        b = GenomicRegion("chr1", 300, 400)
        assert a.distance_to(b) == 100

        c = GenomicRegion("chr1", 150, 250)
        assert a.distance_to(c) == 0  # overlapping

    def test_expand(self):
        r = GenomicRegion("chr1", 100, 200)
        expanded = r.expand(50)
        assert expanded.start == 50
        assert expanded.end == 250

    def test_expand_no_negative(self):
        r = GenomicRegion("chr1", 10, 100)
        expanded = r.expand(50)
        assert expanded.start == 0  # clamped

    def test_to_region_string(self):
        r = GenomicRegion("chr1", 100, 200)
        assert r.to_region_string() == "chr1:100-200"


class TestParseRegion:
    """Tests for parse_region()."""

    def test_basic(self):
        r = parse_region("chr1:1000-2000")
        assert r.chrom == "chr1"
        assert r.start == 1000
        assert r.end == 2000

    def test_with_commas(self):
        r = parse_region("chr1:1,000-2,000")
        assert r.start == 1000

    def test_with_strand(self):
        r = parse_region("chrX:500-1000:-")
        assert r.strand == "-"

    def test_invalid(self):
        with pytest.raises(ValueError):
            parse_region("not_a_region")


class TestSequenceHelpers:
    """Tests for sequence utility functions."""

    def test_reverse_complement(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_gc_content(self):
        assert gc_content("GCGC") == pytest.approx(1.0)
        assert gc_content("ATAT") == pytest.approx(0.0)

    def test_format_bp(self):
        assert format_bp(500) == "500 bp"
        assert format_bp(1500) == "1.5 kb"
        assert format_bp(1500000) == "1.50 Mb"

# splicetarget v0.1.0
# Any usage is subject to this software's license.
