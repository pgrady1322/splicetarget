#!/usr/bin/env python3
"""
splicetarget v0.1.0

sashimi.py — Sashimi-style splice visualization for long-read isoform analysis.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.collections import PatchCollection

from splicetarget.data.reference import Gene
from splicetarget.isoforms.classify import ClassifiedIsoform, IsoformCategory
from splicetarget.splicing.events import EventType, SplicingEvent
from splicetarget.therapeutic.aso_design import ASOCandidate

logger = logging.getLogger(__name__)

# ── Color scheme ─────────────────────────────────────────────────

COLORS = {
    "exon_ref": "#2C3E50",           # dark blue-gray
    "exon_patient_known": "#27AE60", # green
    "exon_patient_novel": "#E74C3C", # red
    "intron": "#BDC3C7",             # light gray
    "junction_known": "#3498DB",     # blue
    "junction_novel": "#E74C3C",     # red
    "exon_skip": "#F39C12",          # orange
    "cryptic_exon": "#9B59B6",       # purple
    "intron_retention": "#E67E22",   # dark orange
    "aso_target": "#1ABC9C",         # teal
    "background": "#FAFAFA",
}

CATEGORY_COLORS = {
    IsoformCategory.FSM: "#27AE60",
    IsoformCategory.ISM: "#2ECC71",
    IsoformCategory.NIC: "#F39C12",
    IsoformCategory.NNC: "#E74C3C",
    IsoformCategory.IR: "#E67E22",
    IsoformCategory.GG: "#95A5A6",
    IsoformCategory.GI: "#95A5A6",
}


def plot_sashimi(
    gene: Gene,
    classified_isoforms: list[ClassifiedIsoform],
    events: list[SplicingEvent] | None = None,
    aso_candidates: list[ASOCandidate] | None = None,
    output_path: str | Path | None = None,
    figsize: tuple[float, float] = (16, 10),
    max_isoforms: int = 15,
    title: str | None = None,
) -> plt.Figure:
    """
    Generate a sashimi-style splice plot.

    Parameters
    ----------
    gene : Gene
        Reference gene model.
    classified_isoforms : list[ClassifiedIsoform]
        Patient isoforms with classification.
    events : list[SplicingEvent] | None
        Detected aberrant splicing events to highlight.
    aso_candidates : list[ASOCandidate] | None
        Top ASO candidates to annotate on the plot.
    output_path : str | Path | None
        If provided, save figure to this path (PNG/PDF/SVG).
    figsize : tuple[float, float]
        Figure dimensions.
    max_isoforms : int
        Maximum isoforms to display (top by read count).
    title : str | None
        Custom plot title. Defaults to gene name.

    Returns
    -------
    matplotlib.figure.Figure
    """
    events = events or []
    aso_candidates = aso_candidates or []

    # Limit isoforms shown
    sorted_isoforms = sorted(classified_isoforms, key=lambda c: c.isoform.read_count, reverse=True)
    display_isoforms = sorted_isoforms[:max_isoforms]

    # Determine genomic range
    gene_start = gene.start
    gene_end = gene.end
    padding = int((gene_end - gene_start) * 0.05)
    view_start = gene_start - padding
    view_end = gene_end + padding

    n_tracks = 1 + len(display_isoforms)  # reference + patient isoforms
    if events:
        n_tracks += 1  # events track
    if aso_candidates:
        n_tracks += 1  # ASO track

    fig, axes = plt.subplots(
        n_tracks, 1,
        figsize=figsize,
        gridspec_kw={"height_ratios": _compute_height_ratios(n_tracks, events, aso_candidates)},
        sharex=True,
    )
    fig.patch.set_facecolor(COLORS["background"])

    if n_tracks == 1:
        axes = [axes]

    ax_idx = 0

    # ── Track 1: Reference gene model ────────────────────────────
    ax = axes[ax_idx]
    _draw_reference_track(ax, gene, view_start, view_end)
    ax_idx += 1

    # ── Tracks 2–N: Patient isoforms ─────────────────────────────
    for ci in display_isoforms:
        ax = axes[ax_idx]
        _draw_isoform_track(ax, ci, view_start, view_end)
        ax_idx += 1

    # ── Events highlight track ───────────────────────────────────
    if events:
        ax = axes[ax_idx]
        _draw_events_track(ax, events, view_start, view_end)
        ax_idx += 1

    # ── ASO candidates track ─────────────────────────────────────
    if aso_candidates:
        ax = axes[ax_idx]
        _draw_aso_track(ax, aso_candidates[:10], view_start, view_end)

    # ── Final formatting ─────────────────────────────────────────
    axes[-1].set_xlabel(f"Genomic position on {gene.chrom}", fontsize=11)
    axes[-1].set_xlim(view_start, view_end)

    if title is None:
        title = f"Splice Analysis: {gene.gene_name} ({gene.gene_id})"
    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info("Sashimi plot saved: %s", output_path)

    return fig


def _draw_reference_track(ax: plt.Axes, gene: Gene, view_start: int, view_end: int) -> None:
    """Draw the reference gene model track."""
    ax.set_facecolor(COLORS["background"])
    y_center = 0.5

    # Draw gene body (thin line)
    ax.plot([gene.start, gene.end], [y_center, y_center], color=COLORS["intron"], linewidth=1, zorder=1)

    # Draw exons from canonical transcript
    canonical = gene.get_canonical_transcript()
    if canonical:
        for exon in canonical.exons:
            rect = mpatches.FancyBboxPatch(
                (exon.start, y_center - 0.3),
                exon.end - exon.start,
                0.6,
                boxstyle="round,pad=0.02",
                facecolor=COLORS["exon_ref"],
                edgecolor="black",
                linewidth=0.5,
                zorder=2,
            )
            ax.add_patch(rect)

            # Exon number label
            mid = (exon.start + exon.end) / 2
            if exon.exon_number and (exon.end - exon.start) > (view_end - view_start) * 0.02:
                ax.text(
                    mid, y_center, str(exon.exon_number),
                    ha="center", va="center", fontsize=7, color="white", fontweight="bold", zorder=3,
                )

    ax.set_ylim(-0.2, 1.2)
    ax.set_ylabel("Reference", fontsize=9, rotation=0, labelpad=60, va="center")
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)

    # Strand arrow
    arrow_dir = "→" if gene.strand == "+" else "←"
    ax.text(
        gene.start, 1.0, f"{gene.gene_name} {arrow_dir}",
        fontsize=10, fontweight="bold", va="bottom",
    )


def _draw_isoform_track(
    ax: plt.Axes,
    ci: ClassifiedIsoform,
    view_start: int,
    view_end: int,
) -> None:
    """Draw a single patient isoform track."""
    iso = ci.isoform
    ax.set_facecolor(COLORS["background"])
    y_center = 0.5

    color = CATEGORY_COLORS.get(ci.category, "#95A5A6")

    # Draw introns (thin lines)
    for i in range(len(iso.exon_blocks) - 1):
        intron_start = iso.exon_blocks[i][1]
        intron_end = iso.exon_blocks[i + 1][0]
        # Arc for splice junction
        mid = (intron_start + intron_end) / 2
        arc_height = min(0.4, (intron_end - intron_start) / (view_end - view_start) * 3)
        arc_x = np.linspace(intron_start, intron_end, 50)
        arc_y = y_center + arc_height * np.sin(np.linspace(0, np.pi, 50))
        ax.plot(arc_x, arc_y, color=color, linewidth=1.5, alpha=0.7)

    # Draw exon blocks
    for exon_start, exon_end in iso.exon_blocks:
        rect = mpatches.FancyBboxPatch(
            (exon_start, y_center - 0.25),
            exon_end - exon_start,
            0.5,
            boxstyle="round,pad=0.01",
            facecolor=color,
            edgecolor="black",
            linewidth=0.5,
            alpha=0.85,
        )
        ax.add_patch(rect)

    # Label
    label = f"{iso.isoform_id} [{ci.category.value}] ({iso.read_count} reads)"
    ax.set_ylabel(label, fontsize=7, rotation=0, labelpad=120, va="center")
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)


def _draw_events_track(
    ax: plt.Axes,
    events: list[SplicingEvent],
    view_start: int,
    view_end: int,
) -> None:
    """Highlight aberrant splicing events."""
    ax.set_facecolor(COLORS["background"])

    event_colors = {
        EventType.EXON_SKIPPING: COLORS["exon_skip"],
        EventType.CRYPTIC_EXON: COLORS["cryptic_exon"],
        EventType.INTRON_RETENTION: COLORS["intron_retention"],
        EventType.ALT_5SS: "#3498DB",
        EventType.ALT_3SS: "#2980B9",
    }

    for i, event in enumerate(events[:8]):
        y = 0.1 + (i * 0.1) % 0.8
        color = event_colors.get(event.event_type, "#95A5A6")

        ax.axvspan(event.start, event.end, alpha=0.2, color=color, zorder=1)
        ax.plot([event.start, event.end], [y, y], color=color, linewidth=3, zorder=2)
        ax.text(
            (event.start + event.end) / 2, y + 0.05,
            f"{event.event_type.value} ({event.total_read_support}r)",
            ha="center", va="bottom", fontsize=6, color=color, fontweight="bold",
        )

    ax.set_ylim(0, 1)
    ax.set_ylabel("Events", fontsize=9, rotation=0, labelpad=60, va="center")
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)


def _draw_aso_track(
    ax: plt.Axes,
    candidates: list[ASOCandidate],
    view_start: int,
    view_end: int,
) -> None:
    """Draw ASO candidate binding sites."""
    ax.set_facecolor(COLORS["background"])

    for i, cand in enumerate(candidates):
        y = 0.1 + (i * 0.08) % 0.8
        alpha = 0.3 + 0.7 * cand.composite_score

        ax.barh(
            y, cand.target_end - cand.target_start,
            left=cand.target_start, height=0.06,
            color=COLORS["aso_target"], alpha=alpha,
            edgecolor="black", linewidth=0.3,
        )
        ax.text(
            cand.target_start, y + 0.04,
            f"#{i+1} ({cand.composite_score:.2f})",
            fontsize=5, va="bottom", color=COLORS["aso_target"],
        )

    ax.set_ylim(0, 1)
    ax.set_ylabel("ASO\ntargets", fontsize=9, rotation=0, labelpad=60, va="center")
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)


def _compute_height_ratios(
    n_tracks: int,
    events: list[SplicingEvent],
    aso_candidates: list[ASOCandidate],
) -> list[float]:
    """Compute subplot height ratios."""
    ratios = [2.0]  # reference track is taller
    n_isoform_tracks = n_tracks - 1 - (1 if events else 0) - (1 if aso_candidates else 0)
    ratios.extend([1.0] * n_isoform_tracks)
    if events:
        ratios.append(1.5)
    if aso_candidates:
        ratios.append(1.0)
    return ratios

# splicetarget v0.1.0
# Any usage is subject to this software's license.
