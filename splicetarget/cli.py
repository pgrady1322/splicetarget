#!/usr/bin/env python3
"""
splicetarget v0.1.0

cli.py — CLI entry point — orchestrates the full analysis pipeline.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import json
import logging
import sys
import time
from pathlib import Path

import click
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

console = Console()

# ── Logging setup ────────────────────────────────────────────────

def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        handlers=[RichHandler(console=console, rich_tracebacks=True, show_time=True)],
    )


# ── Main CLI group ───────────────────────────────────────────────

@click.group()
@click.version_option(package_name="splicetarget")
def main() -> None:
    """splicetarget — Long-read splice isoform analysis & ASO target nomination."""
    pass


# ── Full pipeline command ────────────────────────────────────────

@main.command()
@click.option("--reads", "-r", required=True, type=click.Path(exists=True),
              help="Input reads (BAM, FASTA, or FASTQ). Use aligned BAM to skip alignment.")
@click.option("--reference", "-ref", required=True, type=click.Path(exists=True),
              help="Reference genome FASTA (indexed with .fai).")
@click.option("--annotation", "-a", required=True, type=click.Path(exists=True),
              help="Gene annotation GTF/GFF3 (e.g., GENCODE v44).")
@click.option("--gene", "-g", required=True, type=str,
              help="Target gene symbol (e.g., DMD, SMN2, CFTR).")
@click.option("--outdir", "-o", default="results", type=click.Path(),
              help="Output directory [default: results/].")
@click.option("--read-type", type=click.Choice(["isoseq", "pb_cdna", "ont_drna", "ont_cdna"]),
              default="isoseq", help="Sequencing platform/data type [default: isoseq].")
@click.option("--transcriptome", type=click.Path(exists=True), default=None,
              help="Transcriptome FASTA for off-target assessment (optional).")
@click.option("--min-reads", type=int, default=2,
              help="Minimum reads to report an isoform [default: 2].")
@click.option("--junction-tolerance", type=int, default=10,
              help="Splice junction wobble tolerance in bp [default: 10].")
@click.option("--threads", "-t", type=int, default=4,
              help="CPU threads for alignment [default: 4].")
@click.option("--config", type=click.Path(exists=True), default=None,
              help="YAML config file (overrides CLI options).")
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
def run(
    reads: str,
    reference: str,
    annotation: str,
    gene: str,
    outdir: str,
    read_type: str,
    transcriptome: str | None,
    min_reads: int,
    junction_tolerance: int,
    threads: int,
    config: str | None,
    verbose: bool,
) -> None:
    """Run the full splicetarget analysis pipeline."""
    _setup_logging(verbose)
    logger = logging.getLogger("splicetarget")

    start_time = time.time()
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    reads_path = Path(reads)

    console.print(f"\n[bold cyan]splicetarget[/] — Splice Analysis Pipeline", highlight=False)
    console.print(f"  Gene:       [bold]{gene}[/]")
    console.print(f"  Reads:      {reads}")
    console.print(f"  Reference:  {reference}")
    console.print(f"  Annotation: {annotation}")
    console.print(f"  Output:     {outdir}\n")

    # ── Load config if provided ──────────────────────────────────
    if config:
        import yaml
        with open(config) as f:
            cfg = yaml.safe_load(f)
        logger.info("Loaded config: %s", config)
        min_reads = cfg.get("min_reads", min_reads)
        junction_tolerance = cfg.get("junction_tolerance", junction_tolerance)

    # ── Stage 1: Alignment (if needed) ───────────────────────────
    console.rule("[bold]Stage 1: Alignment")

    if reads_path.suffix in (".bam", ".cram"):
        logger.info("Input is aligned BAM — skipping alignment stage")
        aligned_bam = reads_path
    else:
        from splicetarget.alignment.aligner import ReadType, align_long_reads
        aligned_bam = outdir_path / f"{gene}_aligned.bam"
        align_long_reads(
            reads_path=reads_path,
            reference_path=reference,
            output_bam=aligned_bam,
            read_type=ReadType(read_type),
            threads=threads,
        )

    # Compute alignment stats
    from splicetarget.alignment.aligner import compute_alignment_stats
    stats = compute_alignment_stats(str(aligned_bam))
    console.print(f"  Mapped: {stats.mapped_reads}/{stats.total_reads} ({stats.mapping_rate:.1%})")
    console.print(f"  Spliced: {stats.spliced_reads} ({stats.splice_rate:.1%} of mapped)")

    # ── Stage 2: Reference loading ───────────────────────────────
    console.rule("[bold]Stage 2: Load Reference Transcriptome")

    from splicetarget.data.reference import ReferenceTranscriptome
    ref_db = ReferenceTranscriptome(annotation)
    gene_model = ref_db.get_gene(gene)

    if gene_model is None:
        console.print(f"[red bold]ERROR: Gene '{gene}' not found in annotation[/]")
        sys.exit(1)

    console.print(f"  Gene: {gene_model.gene_name} ({gene_model.gene_id})")
    console.print(f"  Region: {gene_model.chrom}:{gene_model.start}-{gene_model.end}")
    console.print(f"  Transcripts: {gene_model.transcript_count}")

    # ── Stage 3: Isoform analysis ────────────────────────────────
    console.rule("[bold]Stage 3: Isoform Collapsing & Classification")

    from splicetarget.data.io import iter_aligned_reads
    region = f"{gene_model.chrom}:{gene_model.start}-{gene_model.end}"
    reads_list = list(iter_aligned_reads(str(aligned_bam), region=region, require_splice=False))
    console.print(f"  Reads in gene region: {len(reads_list)}")

    from splicetarget.isoforms.collapse import collapse_isoforms
    isoforms = collapse_isoforms(
        reads_list,
        junction_tolerance=junction_tolerance,
        min_reads=min_reads,
        gene_name=gene,
    )
    console.print(f"  Collapsed isoforms: {len(isoforms)}")

    from splicetarget.isoforms.classify import classify_isoforms
    classified = classify_isoforms(isoforms, gene_model, junction_tolerance=junction_tolerance)

    # Classification summary table
    cat_table = Table(title="Isoform Classification Summary")
    cat_table.add_column("Category", style="bold")
    cat_table.add_column("Count", justify="right")
    cat_table.add_column("Total Reads", justify="right")
    cat_counts: dict[str, tuple[int, int]] = {}
    for ci in classified:
        key = ci.category.value
        prev_n, prev_r = cat_counts.get(key, (0, 0))
        cat_counts[key] = (prev_n + 1, prev_r + ci.isoform.read_count)
    for cat, (n, r) in sorted(cat_counts.items()):
        cat_table.add_row(cat, str(n), str(r))
    console.print(cat_table)

    # Quantification
    from splicetarget.isoforms.quantify import quantify_isoforms
    expression = quantify_isoforms(classified, total_gene_reads=len(reads_list))

    # ── Stage 4: Aberrant event detection ────────────────────────
    console.rule("[bold]Stage 4: Aberrant Splicing Event Detection")

    from splicetarget.splicing.events import detect_aberrant_events
    events = detect_aberrant_events(
        classified,
        gene_model,
        total_gene_reads=len(reads_list),
        min_read_support=min_reads,
        junction_tolerance=junction_tolerance,
    )

    if events:
        event_table = Table(title="Detected Aberrant Splicing Events")
        event_table.add_column("Event", style="bold")
        event_table.add_column("Type")
        event_table.add_column("Location")
        event_table.add_column("Reads", justify="right")
        event_table.add_column("Gene %", justify="right")
        event_table.add_column("Frame", justify="center")
        event_table.add_column("ASO?", justify="center")

        for ev in events:
            event_table.add_row(
                ev.event_id,
                ev.event_type.value,
                f"{ev.chrom}:{ev.start}-{ev.end}",
                str(ev.total_read_support),
                f"{ev.fraction_of_gene_expression:.1%}",
                ev.reading_frame_impact or "—",
                "✓" if ev.is_aso_target else "✗",
            )
        console.print(event_table)
    else:
        console.print("  [yellow]No aberrant splicing events detected above threshold[/]")

    # ── Stage 5: ASO candidate design ────────────────────────────
    console.rule("[bold]Stage 5: ASO/SSO Target Design")

    from splicetarget.therapeutic.aso_design import design_aso_candidates
    candidates = design_aso_candidates(events, reference)

    if candidates:
        console.print(f"  Generated {len(candidates)} ASO candidates")

        # Off-target assessment (if transcriptome provided)
        if transcriptome:
            from splicetarget.therapeutic.offtarget import assess_off_targets
            from splicetarget.therapeutic.scoring import rescore_candidates
            assess_off_targets(candidates, transcriptome)
            candidates = rescore_candidates(candidates)
            console.print(f"  Off-target assessment complete (transcriptome: {transcriptome})")

        # Top candidates table
        from splicetarget.therapeutic.scoring import candidates_to_dataframe
        df = candidates_to_dataframe(candidates)
        top_df = df.head(10)

        cand_table = Table(title="Top 10 ASO Candidates")
        for col in ["rank", "candidate_id", "strategy", "aso_sequence", "aso_length",
                     "gc_content", "off_target_hits", "composite_score", "high_confidence"]:
            cand_table.add_column(col.replace("_", " ").title(), justify="right" if col != "aso_sequence" else "left")

        for _, row in top_df.iterrows():
            cand_table.add_row(*[str(row[col]) for col in top_df.columns if col in
                                 ["rank", "candidate_id", "strategy", "aso_sequence", "aso_length",
                                  "gc_content", "off_target_hits", "composite_score", "high_confidence"]])
        console.print(cand_table)

        # Save full table
        csv_path = outdir_path / f"{gene}_aso_candidates.csv"
        df.to_csv(csv_path, index=False)
        console.print(f"  Full candidate table: {csv_path}")
    else:
        console.print("  [yellow]No ASO candidates generated (no suitable events)[/]")

    # ── Stage 6: Visualization & Report ──────────────────────────
    console.rule("[bold]Stage 6: Visualization & Report")

    from splicetarget.visualization.sashimi import plot_sashimi
    plot_path = outdir_path / f"{gene}_sashimi.png"
    plot_sashimi(
        gene=gene_model,
        classified_isoforms=classified,
        events=events,
        aso_candidates=candidates[:10] if candidates else None,
        output_path=plot_path,
    )
    console.print(f"  Sashimi plot: {plot_path}")

    # Save pipeline results as JSON
    results_json = {
        "gene": gene,
        "gene_id": gene_model.gene_id,
        "region": f"{gene_model.chrom}:{gene_model.start}-{gene_model.end}",
        "total_reads": len(reads_list),
        "isoforms": len(isoforms),
        "classification": cat_counts,
        "aberrant_events": len(events),
        "aso_candidates": len(candidates) if candidates else 0,
        "high_confidence_candidates": sum(1 for c in candidates if c.is_high_confidence) if candidates else 0,
        "runtime_seconds": round(time.time() - start_time, 1),
    }
    json_path = outdir_path / f"{gene}_results.json"
    with open(json_path, "w") as f:
        json.dump(results_json, f, indent=2, default=str)

    # ── Summary ──────────────────────────────────────────────────
    elapsed = time.time() - start_time
    console.rule("[bold green]Pipeline Complete")
    console.print(f"  Runtime: {elapsed:.1f}s")
    console.print(f"  Results: {outdir_path}/")
    console.print(f"  Key outputs:")
    console.print(f"    • {json_path.name}          — structured results")
    if candidates:
        console.print(f"    • {csv_path.name}  — ASO candidate table")
    console.print(f"    • {plot_path.name}        — sashimi visualization")
    console.print()


# ── Individual stage commands ────────────────────────────────────

@main.command()
@click.option("--reads", "-r", required=True, type=click.Path(exists=True))
@click.option("--reference", "-ref", required=True, type=click.Path(exists=True))
@click.option("--output", "-o", required=True, type=click.Path())
@click.option("--read-type", default="isoseq",
              type=click.Choice(["isoseq", "pb_cdna", "ont_drna", "ont_cdna"]))
@click.option("--threads", "-t", default=4)
@click.option("--verbose", "-v", is_flag=True)
def align(reads: str, reference: str, output: str, read_type: str, threads: int, verbose: bool) -> None:
    """Run splice-aware alignment only (minimap2)."""
    _setup_logging(verbose)
    from splicetarget.alignment.aligner import ReadType, align_long_reads
    align_long_reads(reads, reference, output, ReadType(read_type), threads=threads)
    console.print(f"[green]Alignment complete → {output}[/]")


@main.command()
@click.option("--bam", required=True, type=click.Path(exists=True))
@click.option("--verbose", "-v", is_flag=True)
def stats(bam: str, verbose: bool) -> None:
    """Compute alignment QC statistics."""
    _setup_logging(verbose)
    from splicetarget.alignment.aligner import compute_alignment_stats
    s = compute_alignment_stats(bam)
    console.print(f"Total reads:    {s.total_reads}")
    console.print(f"Mapped reads:   {s.mapped_reads} ({s.mapping_rate:.1%})")
    console.print(f"Spliced reads:  {s.spliced_reads} ({s.splice_rate:.1%} of mapped)")
    console.print(f"Mean MAPQ:      {s.mean_mapping_quality:.1f}")
    console.print(f"Median length:  {s.median_read_length}")


if __name__ == "__main__":
    main()

# splicetarget v0.1.0
# Any usage is subject to this software's license.
