#!/usr/bin/env python3
"""
splicetarget v0.1.0

reference.py — Reference transcriptome handling for GENCODE/Ensembl annotations.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import gffutils

logger = logging.getLogger(__name__)


# ── Data classes ──────────────────────────────────────────────────

@dataclass
class Exon:
    """A single exon within a transcript."""

    chrom: str
    start: int
    end: int
    strand: str
    exon_number: int | None = None

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Transcript:
    """A reference transcript with ordered exons."""

    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    biotype: str = ""
    exons: list[Exon] = field(default_factory=list)

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @property
    def introns(self) -> list[tuple[int, int]]:
        """Derive intron coordinates from sorted exon list."""
        sorted_exons = sorted(self.exons, key=lambda e: e.start)
        introns = []
        for i in range(len(sorted_exons) - 1):
            introns.append((sorted_exons[i].end, sorted_exons[i + 1].start))
        return introns

    @property
    def splice_junctions(self) -> set[tuple[int, int]]:
        """Set of (donor, acceptor) positions for this transcript."""
        return set(self.introns)

    @property
    def exon_chain(self) -> tuple[tuple[int, int], ...]:
        """Canonical exon chain as sorted tuple of (start, end) pairs."""
        return tuple(sorted((e.start, e.end) for e in self.exons))


@dataclass
class Gene:
    """A gene locus with associated transcripts."""

    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    biotype: str = ""
    transcripts: dict[str, Transcript] = field(default_factory=dict)

    @property
    def transcript_count(self) -> int:
        return len(self.transcripts)

    def get_canonical_transcript(self) -> Transcript | None:
        """Return the longest coding transcript as a proxy for canonical."""
        coding = [
            t for t in self.transcripts.values()
            if "protein_coding" in t.biotype
        ]
        if not coding:
            coding = list(self.transcripts.values())
        if not coding:
            return None
        return max(coding, key=lambda t: t.end - t.start)


# ── Reference database ───────────────────────────────────────────

class ReferenceTranscriptome:
    """
    Indexed reference transcriptome from GTF/GFF3 annotation.

    Uses gffutils to build an SQLite-backed database for fast gene/transcript
    lookups. The database is cached to disk for subsequent runs.

    Parameters
    ----------
    gtf_path : str | Path
        Path to GENCODE/Ensembl GTF or GFF3 file.
    db_path : str | Path | None
        Path for gffutils SQLite database. If None, derived from gtf_path.
    force_rebuild : bool
        If True, rebuild the database even if cached version exists.
    """

    def __init__(
        self,
        gtf_path: str | Path,
        db_path: str | Path | None = None,
        force_rebuild: bool = False,
    ) -> None:
        self.gtf_path = Path(gtf_path)
        self.db_path = Path(db_path) if db_path else self.gtf_path.with_suffix(".db")
        self.db = self._load_or_create_db(force_rebuild)
        self._gene_cache: dict[str, Gene] = {}

    def _load_or_create_db(self, force_rebuild: bool) -> gffutils.FeatureDB:
        """Load existing gffutils database or create from GTF."""
        if self.db_path.exists() and not force_rebuild:
            logger.info("Loading cached annotation database: %s", self.db_path)
            return gffutils.FeatureDB(str(self.db_path))

        logger.info("Building annotation database from %s (this may take a few minutes)...", self.gtf_path)
        db = gffutils.create_db(
            str(self.gtf_path),
            str(self.db_path),
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
            disable_infer_genes=False,
            disable_infer_transcripts=False,
        )
        logger.info("Annotation database built: %s", self.db_path)
        return db

    def get_gene(self, gene_name: str) -> Gene | None:
        """
        Get a Gene object by gene name (symbol) or gene ID.

        Parameters
        ----------
        gene_name : str
            Gene symbol (e.g., "DMD") or Ensembl ID (e.g., "ENSG00000198947").

        Returns
        -------
        Gene | None
        """
        if gene_name in self._gene_cache:
            return self._gene_cache[gene_name]

        # Try by gene_name attribute first, then by ID
        gene_features = list(self.db.all_features(featuretype="gene"))
        gene_feature = None

        for gf in gene_features:
            names = gf.attributes.get("gene_name", [])
            ids = gf.attributes.get("gene_id", [])
            if gene_name in names or gene_name in ids:
                gene_feature = gf
                break

        if gene_feature is None:
            logger.warning("Gene '%s' not found in annotation", gene_name)
            return None

        gene = self._build_gene(gene_feature)
        self._gene_cache[gene_name] = gene
        return gene

    def get_gene_region(self, gene_name: str) -> Gene | None:
        """
        Efficient gene lookup using region-based query.

        For large GTF files, this approach narrows the search space
        by first identifying the gene feature, then only loading
        child transcripts within that region.
        """
        # For large GTFs, prefer iterating by featuretype with a filter
        try:
            for gf in self.db.all_features(featuretype="gene"):
                attrs = gf.attributes
                if gene_name in attrs.get("gene_name", []) or gene_name in attrs.get("gene_id", []):
                    return self._build_gene(gf)
        except Exception:
            logger.error("Failed to look up gene '%s'", gene_name)
        return None

    def get_transcripts_for_gene(self, gene_name: str) -> list[Transcript]:
        """Return all reference transcripts for a gene."""
        gene = self.get_gene(gene_name)
        if gene is None:
            return []
        return list(gene.transcripts.values())

    def get_all_splice_junctions(self, gene_name: str) -> set[tuple[int, int]]:
        """Return the union of all reference splice junctions for a gene."""
        transcripts = self.get_transcripts_for_gene(gene_name)
        all_junctions: set[tuple[int, int]] = set()
        for tx in transcripts:
            all_junctions.update(tx.splice_junctions)
        return all_junctions

    def _build_gene(self, gene_feature: gffutils.Feature) -> Gene:
        """Build a Gene object with all child transcripts and exons."""
        gene_id = gene_feature.attributes.get("gene_id", [""])[0]
        gene_name = gene_feature.attributes.get("gene_name", [gene_id])[0]
        biotype = gene_feature.attributes.get("gene_biotype", gene_feature.attributes.get("gene_type", [""]))[0]

        gene = Gene(
            gene_id=gene_id,
            gene_name=gene_name,
            chrom=gene_feature.seqid,
            start=gene_feature.start,
            end=gene_feature.end,
            strand=gene_feature.strand,
            biotype=biotype,
        )

        # Iterate through child transcripts
        for tx_feature in self.db.children(gene_feature, featuretype="transcript"):
            tx_id = tx_feature.attributes.get("transcript_id", [""])[0]
            tx_biotype = tx_feature.attributes.get("transcript_biotype",
                         tx_feature.attributes.get("transcript_type", [""]))[0]

            transcript = Transcript(
                transcript_id=tx_id,
                gene_id=gene_id,
                gene_name=gene_name,
                chrom=tx_feature.seqid,
                start=tx_feature.start,
                end=tx_feature.end,
                strand=tx_feature.strand,
                biotype=tx_biotype,
            )

            # Add exons
            exon_num = 0
            for exon_feature in self.db.children(tx_feature, featuretype="exon"):
                exon_num += 1
                transcript.exons.append(Exon(
                    chrom=exon_feature.seqid,
                    start=exon_feature.start,
                    end=exon_feature.end,
                    strand=exon_feature.strand,
                    exon_number=exon_num,
                ))

            # Sort exons by position
            transcript.exons.sort(key=lambda e: e.start)
            gene.transcripts[tx_id] = transcript

        logger.info(
            "Loaded gene %s (%s): %d transcripts",
            gene_name, gene_id, gene.transcript_count,
        )
        return gene

# splicetarget v0.1.0
# Any usage is subject to this software's license.
