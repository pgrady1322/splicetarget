#!/usr/bin/env python3
"""
splicetarget v0.1.0

test_therapeutic.py — Tests for ASO candidate design and scoring.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import pytest

from splicetarget.therapeutic.aso_design import (
    ASOCandidate,
    ASOStrategy,
    DesignParams,
    _gc_content,
    _has_homopolymer,
    _max_self_comp,
    _reverse_complement,
)
from splicetarget.therapeutic.scoring import (
    ScoringWeights,
    candidates_to_dataframe,
    rescore_candidates,
)


class TestSequenceUtils:
    """Tests for low-level sequence utility functions."""

    def test_reverse_complement(self):
        assert _reverse_complement("ATCG") == "CGAT"
        assert _reverse_complement("AAAA") == "TTTT"
        assert _reverse_complement("GCTA") == "TAGC"

    def test_gc_content(self):
        assert _gc_content("GCGC") == pytest.approx(1.0)
        assert _gc_content("ATAT") == pytest.approx(0.0)
        assert _gc_content("ATGC") == pytest.approx(0.5)
        assert _gc_content("") == pytest.approx(0.0)

    def test_homopolymer_detection(self):
        assert _has_homopolymer("AAAAATCG", 4) is True
        assert _has_homopolymer("AATCG", 4) is False
        assert _has_homopolymer("GGGGGATCG", 4) is True
        assert _has_homopolymer("ATCGATCG", 4) is False

    def test_self_complementarity(self):
        # Palindromic sequence should have high self-comp
        comp = _max_self_comp("ATCGATCGAT")
        assert comp > 0
        # Random-ish sequence should have lower self-comp
        comp2 = _max_self_comp("ATGCTTACGA")
        assert isinstance(comp2, int)


class TestASODesignParams:
    """Tests for design parameter validation."""

    def test_default_params(self):
        params = DesignParams()
        assert params.min_length == 18
        assert params.max_length == 25
        assert params.gc_min < params.gc_max

    def test_scoring_weights_valid(self):
        weights = ScoringWeights()
        weights.validate()  # should not raise

    def test_scoring_weights_invalid(self):
        weights = ScoringWeights(gc_content=0.5, self_complementarity=0.5)
        with pytest.raises(ValueError, match="sum to 1.0"):
            weights.validate()


class TestASOCandidate:
    """Tests for ASOCandidate properties."""

    def test_high_confidence(self):
        cand = ASOCandidate(
            candidate_id="test.aso001",
            event_id="test.event001",
            strategy=ASOStrategy.EXON_INCLUSION,
            chrom="chr1",
            target_start=1000,
            target_end=1020,
            strand="+",
            target_sequence="ATCGATCGATCGATCGATCG",
            aso_sequence="CGATCGATCGATCGATCGAT",
            composite_score=0.75,
            off_target_hits=1,
        )
        assert cand.is_high_confidence is True

    def test_low_confidence_off_targets(self):
        cand = ASOCandidate(
            candidate_id="test.aso002",
            event_id="test.event001",
            strategy=ASOStrategy.EXON_INCLUSION,
            chrom="chr1",
            target_start=1000,
            target_end=1020,
            strand="+",
            target_sequence="ATCGATCGATCGATCGATCG",
            aso_sequence="CGATCGATCGATCGATCGAT",
            composite_score=0.75,
            off_target_hits=10,  # too many off-targets
        )
        assert cand.is_high_confidence is False


class TestScoring:
    """Tests for candidate scoring and ranking."""

    def _make_candidate(self, gc: float = 0.5, self_comp: int = 2, dist: int = 10) -> ASOCandidate:
        return ASOCandidate(
            candidate_id="test.aso001",
            event_id="test.event001",
            strategy=ASOStrategy.EXON_INCLUSION,
            chrom="chr1",
            target_start=1000,
            target_end=1020,
            strand="+",
            target_sequence="A" * 20,
            aso_sequence="T" * 20,
            aso_length=20,
            gc_content=gc,
            max_self_complementarity=self_comp,
            distance_to_splice_site=dist,
            gc_score=max(0, 1.0 - abs(gc - 0.5) / 0.2),
            self_comp_score=max(0, 1.0 - self_comp / 8.0),
            splice_proximity_score=1.0 / (1.0 + dist / 50.0),
            length_score=1.0,
            accessibility_score=0.5,
            off_target_hits=0,
        )

    def test_rescore_preserves_order(self):
        """Better candidates should rank higher after rescoring."""
        good = self._make_candidate(gc=0.5, self_comp=1, dist=5)
        bad = self._make_candidate(gc=0.8, self_comp=6, dist=200)
        bad.candidate_id = "test.aso002"

        rescored = rescore_candidates([bad, good])
        assert rescored[0].candidate_id == good.candidate_id

    def test_dataframe_output(self):
        """DataFrame should have expected columns."""
        cand = self._make_candidate()
        cand.composite_score = 0.75
        df = candidates_to_dataframe([cand])
        assert "rank" in df.columns
        assert "aso_sequence" in df.columns
        assert "composite_score" in df.columns
        assert len(df) == 1

# splicetarget v0.1.0
# Any usage is subject to this software's license.
