#!/usr/bin/env python3
"""
splicetarget v0.1.0

conftest.py — Shared pytest configuration and fixtures.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import pytest

from splicetarget.data.io import ReadRecord, SpliceJunction
from splicetarget.data.reference import Exon, Gene, Transcript
from splicetarget.isoforms.collapse import CollapsedIsoform
from splicetarget.isoforms.classify import ClassifiedIsoform, IsoformCategory
from splicetarget.splicing.events import EventType, SplicingEvent


# ── Test gene model (mini DMD-like gene) ─────────────────────────

@pytest.fixture
def mock_gene() -> Gene:
    """
    A simplified 5-exon gene model mimicking a DMD-like structure.

    Exon layout (chr1):
        Exon 1: 1000-1200
        Exon 2: 2000-2150  (target for exon skipping tests)
        Exon 3: 3000-3300
        Exon 4: 4000-4100
        Exon 5: 5000-5500
    """
    tx = Transcript(
        transcript_id="TX001",
        gene_id="GENE001",
        gene_name="CandidateGene",
        chrom="chr1",
        start=1000,
        end=5500,
        strand="+",
        biotype="protein_coding",
        exons=[
            Exon(chrom="chr1", start=1000, end=1200, strand="+", exon_number=1),
            Exon(chrom="chr1", start=2000, end=2150, strand="+", exon_number=2),
            Exon(chrom="chr1", start=3000, end=3300, strand="+", exon_number=3),
            Exon(chrom="chr1", start=4000, end=4100, strand="+", exon_number=4),
            Exon(chrom="chr1", start=5000, end=5500, strand="+", exon_number=5),
        ],
    )

    gene = Gene(
        gene_id="GENE001",
        gene_name="CandidateGene",
        chrom="chr1",
        start=1000,
        end=5500,
        strand="+",
        biotype="protein_coding",
        transcripts={"TX001": tx},
    )
    return gene


@pytest.fixture
def mock_reads_normal() -> list[ReadRecord]:
    """Reads matching the reference transcript (FSM)."""
    return [
        ReadRecord(
            name=f"read_normal_{i}",
            sequence="A" * 500,
            quality=None,
            chrom="chr1",
            start=1000,
            end=5500,
            strand="+",
            is_mapped=True,
            exon_blocks=[(1000, 1200), (2000, 2150), (3000, 3300), (4000, 4100), (5000, 5500)],
        )
        for i in range(10)
    ]


@pytest.fixture
def mock_reads_exon_skip() -> list[ReadRecord]:
    """Reads with exon 2 skipped (exon skipping event)."""
    return [
        ReadRecord(
            name=f"read_skip_{i}",
            sequence="A" * 400,
            quality=None,
            chrom="chr1",
            start=1000,
            end=5500,
            strand="+",
            is_mapped=True,
            # Missing exon 2 (2000-2150): jumps from exon 1 to exon 3
            exon_blocks=[(1000, 1200), (3000, 3300), (4000, 4100), (5000, 5500)],
        )
        for i in range(5)
    ]


@pytest.fixture
def mock_reads_cryptic_exon() -> list[ReadRecord]:
    """Reads with a cryptic exon included between exons 2 and 3."""
    return [
        ReadRecord(
            name=f"read_cryptic_{i}",
            sequence="A" * 550,
            quality=None,
            chrom="chr1",
            start=1000,
            end=5500,
            strand="+",
            is_mapped=True,
            # Cryptic exon at 2500-2600 (within intron 2-3)
            exon_blocks=[
                (1000, 1200), (2000, 2150), (2500, 2600),
                (3000, 3300), (4000, 4100), (5000, 5500),
            ],
        )
        for i in range(4)
    ]


@pytest.fixture
def mock_reads_intron_retention() -> list[ReadRecord]:
    """Reads retaining intron 2-3 (failure to splice)."""
    return [
        ReadRecord(
            name=f"read_ir_{i}",
            sequence="A" * 600,
            quality=None,
            chrom="chr1",
            start=1000,
            end=5500,
            strand="+",
            is_mapped=True,
            # Exon 2 extends through intron into exon 3 (retention)
            exon_blocks=[(1000, 1200), (2000, 3300), (4000, 4100), (5000, 5500)],
        )
        for i in range(3)
    ]


@pytest.fixture
def mock_splice_junction() -> SpliceJunction:
    """A canonical GT-AG splice junction."""
    return SpliceJunction(
        chrom="chr1",
        start=1200,
        end=2000,
        strand="+",
        read_count=10,
        canonical=True,
        donor_motif="GT",
        acceptor_motif="AG",
    )

# splicetarget v0.1.0
# Any usage is subject to this software's license.
