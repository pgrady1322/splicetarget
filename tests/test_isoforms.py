#!/usr/bin/env python3
"""
splicetarget v0.1.0

test_isoforms.py — Tests for isoform collapsing and classification.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import pytest

from splicetarget.isoforms.collapse import collapse_isoforms, CollapsedIsoform
from splicetarget.isoforms.classify import classify_isoforms, IsoformCategory
from splicetarget.isoforms.quantify import quantify_isoforms


class TestIsoformCollapsing:
    """Tests for collapse_isoforms()."""

    def test_collapse_identical_reads(self, mock_reads_normal):
        """Identical reads should collapse into a single isoform."""
        isoforms = collapse_isoforms(mock_reads_normal, min_reads=2, gene_name="TEST")
        assert len(isoforms) == 1
        assert isoforms[0].read_count == 10
        assert isoforms[0].exon_count == 5

    def test_collapse_different_structures(self, mock_reads_normal, mock_reads_exon_skip):
        """Normal + exon-skip reads should collapse into two isoforms."""
        all_reads = mock_reads_normal + mock_reads_exon_skip
        isoforms = collapse_isoforms(all_reads, min_reads=2, gene_name="TEST")
        assert len(isoforms) == 2

        # Sort by exon count to identify them
        by_exons = sorted(isoforms, key=lambda i: i.exon_count)
        assert by_exons[0].exon_count == 4  # exon skip
        assert by_exons[1].exon_count == 5  # normal

    def test_min_reads_filter(self, mock_reads_normal):
        """Isoforms below min_reads threshold should be filtered out."""
        # Only 1 read — should produce 0 isoforms with min_reads=2
        isoforms = collapse_isoforms(mock_reads_normal[:1], min_reads=2, gene_name="TEST")
        assert len(isoforms) == 0

    def test_collapse_empty_input(self):
        """Empty input should return empty list."""
        assert collapse_isoforms([], min_reads=1) == []

    def test_intron_chain_property(self, mock_reads_normal):
        """Collapsed isoform should have correct intron chain."""
        isoforms = collapse_isoforms(mock_reads_normal, min_reads=2, gene_name="TEST")
        iso = isoforms[0]
        assert len(iso.intron_chain) == 4  # 5 exons → 4 introns

    def test_bed12_output(self, mock_reads_normal):
        """BED12 export should produce valid format."""
        isoforms = collapse_isoforms(mock_reads_normal, min_reads=2, gene_name="TEST")
        bed_line = isoforms[0].to_bed12()
        fields = bed_line.split("\t")
        assert len(fields) == 12
        assert fields[0] == "chr1"


class TestIsoformClassification:
    """Tests for classify_isoforms()."""

    def test_classify_fsm(self, mock_reads_normal, mock_gene):
        """Reads matching reference should be classified as FSM."""
        isoforms = collapse_isoforms(mock_reads_normal, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        assert len(classified) == 1
        assert classified[0].category == IsoformCategory.FSM

    def test_classify_exon_skip_as_nic_or_nnc(self, mock_reads_exon_skip, mock_gene):
        """Exon-skipping reads should be NIC (known donors/acceptors, novel combo)."""
        isoforms = collapse_isoforms(mock_reads_exon_skip, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        assert len(classified) == 1
        # Exon skipping uses known splice sites → NIC
        assert classified[0].category in (IsoformCategory.NIC, IsoformCategory.NNC)

    def test_classify_intron_retention(self, mock_reads_intron_retention, mock_gene):
        """Intron retention reads should be classified as IR."""
        isoforms = collapse_isoforms(mock_reads_intron_retention, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        assert len(classified) == 1
        assert classified[0].category == IsoformCategory.IR

    def test_classify_cryptic_exon(self, mock_reads_cryptic_exon, mock_gene):
        """Cryptic exon inclusion should be classified as NNC."""
        isoforms = collapse_isoforms(mock_reads_cryptic_exon, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        assert len(classified) == 1
        assert classified[0].category == IsoformCategory.NNC

    def test_is_therapy_candidate(self, mock_reads_intron_retention, mock_gene):
        """IR and NNC isoforms should be flagged as therapy candidates."""
        isoforms = collapse_isoforms(mock_reads_intron_retention, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        assert classified[0].is_therapy_candidate


class TestIsoformQuantification:
    """Tests for quantify_isoforms()."""

    def test_quantify_single_isoform(self, mock_reads_normal, mock_gene):
        """Single isoform should get 100% of expression."""
        isoforms = collapse_isoforms(mock_reads_normal, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        expression = quantify_isoforms(classified, total_gene_reads=10)
        assert len(expression) == 1
        assert expression[0].fraction == pytest.approx(1.0)

    def test_quantify_mixed_isoforms(self, mock_reads_normal, mock_reads_exon_skip, mock_gene):
        """Mixed isoforms should have proportional expression."""
        all_reads = mock_reads_normal + mock_reads_exon_skip
        isoforms = collapse_isoforms(all_reads, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        expression = quantify_isoforms(classified, total_gene_reads=15)
        fractions = [e.fraction for e in expression]
        assert sum(fractions) == pytest.approx(1.0)

    def test_aberrant_flag(self, mock_reads_intron_retention, mock_gene):
        """IR isoforms should be flagged as aberrant in expression data."""
        isoforms = collapse_isoforms(mock_reads_intron_retention, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        expression = quantify_isoforms(classified)
        assert all(e.is_aberrant for e in expression)

# splicetarget v0.1.0
# Any usage is subject to this software's license.
