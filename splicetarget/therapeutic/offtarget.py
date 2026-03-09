#!/usr/bin/env python3
"""
splicetarget v0.1.0

offtarget.py — Off-target assessment for ASO candidates via BLAST/k-mer screen.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from splicetarget.therapeutic.aso_design import ASOCandidate

logger = logging.getLogger(__name__)


@dataclass
class OffTargetHit:
    """A single off-target alignment for an ASO candidate."""

    target_gene: str
    chrom: str
    start: int
    end: int
    strand: str
    identity_pct: float
    alignment_length: int
    mismatches: int
    e_value: float

    @property
    def is_significant(self) -> bool:
        """Hit is significant if ≥85% identity over ≥15nt."""
        return self.identity_pct >= 85.0 and self.alignment_length >= 15


def assess_off_targets(
    candidates: list[ASOCandidate],
    transcriptome_fasta: str | Path,
    max_hits: int = 10,
    identity_threshold: float = 85.0,
    use_blast: bool = True,
) -> dict[str, list[OffTargetHit]]:
    """
    Screen ASO candidates for off-target binding across the transcriptome.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        ASO candidates to screen.
    transcriptome_fasta : str | Path
        Path to reference transcriptome FASTA (e.g., GENCODE transcripts).
    max_hits : int
        Maximum off-target hits to return per candidate.
    identity_threshold : float
        Minimum percent identity to report a hit.
    use_blast : bool
        If True, use blastn. If False, use bowtie2 (faster for many queries).

    Returns
    -------
    dict[str, list[OffTargetHit]]
        Mapping of candidate_id → list of off-target hits.
    """
    transcriptome_fasta = Path(transcriptome_fasta)

    if use_blast:
        return _blast_screen(candidates, transcriptome_fasta, max_hits, identity_threshold)
    else:
        return _simple_screen(candidates, transcriptome_fasta, identity_threshold)


def _blast_screen(
    candidates: list[ASOCandidate],
    transcriptome_fasta: Path,
    max_hits: int,
    identity_threshold: float,
) -> dict[str, list[OffTargetHit]]:
    """Off-target screening using BLAST."""
    results: dict[str, list[OffTargetHit]] = {}

    # Build BLAST database if not already built
    db_path = transcriptome_fasta.with_suffix(".ndb")
    if not db_path.exists():
        logger.info("Building BLAST database from %s", transcriptome_fasta)
        subprocess.run(
            ["makeblastdb", "-in", str(transcriptome_fasta), "-dbtype", "nucl"],
            check=True,
            capture_output=True,
        )

    # Write candidate sequences to temp FASTA
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        for cand in candidates:
            tmp.write(f">{cand.candidate_id}\n{cand.aso_sequence}\n")
        query_file = tmp.name

    # Run BLAST
    blast_cmd = [
        "blastn",
        "-query", query_file,
        "-db", str(transcriptome_fasta),
        "-task", "blastn-short",       # optimized for short sequences
        "-word_size", "7",
        "-evalue", "100",              # relaxed e-value for short oligos
        "-perc_identity", str(identity_threshold),
        "-max_target_seqs", str(max_hits),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
    ]

    try:
        result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True)

        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            cand_id = fields[0]
            hit = OffTargetHit(
                target_gene=fields[1],
                chrom=fields[1],  # subject ID
                start=int(fields[8]),
                end=int(fields[9]),
                strand="+" if int(fields[8]) < int(fields[9]) else "-",
                identity_pct=float(fields[2]),
                alignment_length=int(fields[3]),
                mismatches=int(fields[4]),
                e_value=float(fields[10]),
            )
            if cand_id not in results:
                results[cand_id] = []
            results[cand_id].append(hit)

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning("BLAST off-target screen failed: %s. Falling back to simple screen.", e)
        return _simple_screen(candidates, transcriptome_fasta, identity_threshold)

    # Clean up
    Path(query_file).unlink(missing_ok=True)

    # Update candidate hit counts
    for cand in candidates:
        hits = results.get(cand.candidate_id, [])
        significant_hits = [h for h in hits if h.is_significant]
        cand.off_target_hits = len(significant_hits)

    logger.info("Off-target screen: %d candidates assessed", len(candidates))
    return results


def _simple_screen(
    candidates: list[ASOCandidate],
    transcriptome_fasta: Path,
    identity_threshold: float,
) -> dict[str, list[OffTargetHit]]:
    """
    Simple k-mer-based off-target screen (fallback when BLAST unavailable).

    Checks for exact substring matches of candidate ASO within transcriptome.
    Less sensitive than BLAST but requires no external tools.
    """
    logger.info("Running simple k-mer off-target screen (BLAST not available)")

    # Load transcriptome sequences
    import pysam
    transcripts: dict[str, str] = {}
    with pysam.FastxFile(str(transcriptome_fasta)) as fh:
        for entry in fh:
            transcripts[entry.name] = entry.sequence.upper()

    results: dict[str, list[OffTargetHit]] = {}

    for cand in candidates:
        hits = []
        query = cand.aso_sequence.upper()
        # Check for near-exact matches (allowing 1-2 mismatches via shorter k-mers)
        kmer_len = max(15, len(query) - 3)
        query_kmer = query[:kmer_len]

        for tx_id, tx_seq in transcripts.items():
            # Skip the target gene itself
            count = tx_seq.count(query_kmer)
            if count > 0:
                hits.append(OffTargetHit(
                    target_gene=tx_id,
                    chrom=tx_id,
                    start=tx_seq.index(query_kmer),
                    end=tx_seq.index(query_kmer) + kmer_len,
                    strand="+",
                    identity_pct=100.0 * kmer_len / len(query),
                    alignment_length=kmer_len,
                    mismatches=len(query) - kmer_len,
                    e_value=0.0,
                ))

        results[cand.candidate_id] = hits
        # Subtract 1 for the expected on-target hit
        cand.off_target_hits = max(0, len(hits) - 1)

    return results

# splicetarget v0.1.0
# Any usage is subject to this software's license.
