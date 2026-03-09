#!/usr/bin/env python3
"""
splicetarget v0.1.0

scoring.py — Multi-criteria ASO candidate scoring and ranking.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd

from splicetarget.therapeutic.aso_design import ASOCandidate

logger = logging.getLogger(__name__)


@dataclass
class ScoringWeights:
    """Adjustable weights for composite ASO scoring."""

    gc_content: float = 0.20
    self_complementarity: float = 0.15
    splice_proximity: float = 0.20
    length_optimality: float = 0.10
    rna_accessibility: float = 0.15
    off_target_penalty: float = 0.20

    def validate(self) -> None:
        total = (
            self.gc_content + self.self_complementarity + self.splice_proximity
            + self.length_optimality + self.rna_accessibility + self.off_target_penalty
        )
        if abs(total - 1.0) > 0.01:
            raise ValueError(f"Scoring weights must sum to 1.0, got {total:.2f}")


def rescore_candidates(
    candidates: list[ASOCandidate],
    weights: ScoringWeights | None = None,
) -> list[ASOCandidate]:
    """
    Re-score candidates with custom weights (e.g., after off-target assessment).

    Parameters
    ----------
    candidates : list[ASOCandidate]
        ASO candidates with individual component scores populated.
    weights : ScoringWeights | None
        Custom weights. If None, uses defaults.

    Returns
    -------
    list[ASOCandidate]
        Re-scored and re-sorted candidates.
    """
    if weights is None:
        weights = ScoringWeights()
    weights.validate()

    for cand in candidates:
        # Off-target penalty: 1.0 if no hits, decaying with more hits
        if cand.off_target_hits >= 0:
            ot_score = max(0.0, 1.0 - (cand.off_target_hits / 5.0))
        else:
            ot_score = 0.5  # unknown = neutral

        cand.composite_score = (
            weights.gc_content * cand.gc_score
            + weights.self_complementarity * cand.self_comp_score
            + weights.splice_proximity * cand.splice_proximity_score
            + weights.length_optimality * cand.length_score
            + weights.rna_accessibility * cand.accessibility_score
            + weights.off_target_penalty * ot_score
        )

    candidates.sort(key=lambda c: c.composite_score, reverse=True)
    return candidates


def candidates_to_dataframe(candidates: list[ASOCandidate]) -> pd.DataFrame:
    """
    Convert ASO candidates to a pandas DataFrame for tabular output.

    Returns a presentation-ready table with columns suitable for the
    therapeutic design team.
    """
    rows = []
    for rank, cand in enumerate(candidates, 1):
        rows.append({
            "rank": rank,
            "candidate_id": cand.candidate_id,
            "event_id": cand.event_id,
            "strategy": cand.strategy.value,
            "chrom": cand.chrom,
            "target_start": cand.target_start,
            "target_end": cand.target_end,
            "strand": cand.strand,
            "aso_sequence": cand.aso_sequence,
            "aso_length": cand.aso_length,
            "target_region": cand.target_region_type,
            "gc_content": round(cand.gc_content, 3),
            "self_comp_max": cand.max_self_complementarity,
            "splice_site_dist": cand.distance_to_splice_site,
            "off_target_hits": cand.off_target_hits,
            "gc_score": round(cand.gc_score, 3),
            "self_comp_score": round(cand.self_comp_score, 3),
            "proximity_score": round(cand.splice_proximity_score, 3),
            "length_score": round(cand.length_score, 3),
            "access_score": round(cand.accessibility_score, 3),
            "composite_score": round(cand.composite_score, 3),
            "high_confidence": cand.is_high_confidence,
        })

    df = pd.DataFrame(rows)
    logger.info(
        "Generated candidate table: %d total, %d high-confidence",
        len(df),
        df["high_confidence"].sum() if len(df) > 0 else 0,
    )
    return df

# splicetarget v0.1.0
# Any usage is subject to this software's license.
