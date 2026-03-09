#!/usr/bin/env python3
"""
splicetarget v0.1.0

quantify.py — Isoform expression quantification (counts, fractions, TPM).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from splicetarget.isoforms.classify import ClassifiedIsoform

logger = logging.getLogger(__name__)


@dataclass
class IsoformExpression:
    """Expression quantification for a single isoform."""

    isoform_id: str
    category: str
    read_count: int
    fraction: float           # proportion of total gene expression
    tpm_estimate: float = 0.0 # TPM-like normalized expression
    exon_count: int = 0
    total_exon_length: int = 0
    is_aberrant: bool = False

    @property
    def rpk(self) -> float:
        """Reads per kilobase of isoform length."""
        if self.total_exon_length == 0:
            return 0.0
        return self.read_count / (self.total_exon_length / 1000)


def quantify_isoforms(
    classified: list[ClassifiedIsoform],
    total_gene_reads: int | None = None,
) -> list[IsoformExpression]:
    """
    Compute expression quantification for classified isoforms.

    Calculates per-isoform read counts, gene-level fractions,
    and a TPM-like normalized expression metric.

    Parameters
    ----------
    classified : list[ClassifiedIsoform]
        Classified isoforms for a single gene.
    total_gene_reads : int | None
        Total reads mapping to this gene. If None, computed from isoform counts.

    Returns
    -------
    list[IsoformExpression]
        Expression data sorted by fraction (descending).
    """
    if not classified:
        return []

    if total_gene_reads is None:
        total_gene_reads = sum(c.isoform.read_count for c in classified)

    if total_gene_reads == 0:
        total_gene_reads = 1  # avoid division by zero

    # Compute RPK for TPM normalization
    expressions = []
    total_rpk = 0.0

    for c in classified:
        iso = c.isoform
        frac = iso.read_count / total_gene_reads

        expr = IsoformExpression(
            isoform_id=iso.isoform_id,
            category=c.category.value,
            read_count=iso.read_count,
            fraction=frac,
            exon_count=iso.exon_count,
            total_exon_length=iso.total_exon_length,
            is_aberrant=c.category.is_aberrant,
        )
        total_rpk += expr.rpk
        expressions.append(expr)

    # Normalize to TPM-like values
    if total_rpk > 0:
        for expr in expressions:
            expr.tpm_estimate = (expr.rpk / total_rpk) * 1e6

    expressions.sort(key=lambda e: e.fraction, reverse=True)

    # Log summary
    aberrant_frac = sum(e.fraction for e in expressions if e.is_aberrant)
    logger.info(
        "Quantified %d isoforms: %.1f%% gene expression from aberrant splice forms",
        len(expressions), aberrant_frac * 100,
    )

    return expressions

# splicetarget v0.1.0
# Any usage is subject to this software's license.
