#!/usr/bin/env python3
"""
splicetarget v0.1.0

test_splicing.py — Tests for aberrant splicing event detection.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

from __future__ import annotations

import pytest

from splicetarget.isoforms.collapse import collapse_isoforms
from splicetarget.isoforms.classify import classify_isoforms, IsoformCategory
from splicetarget.splicing.events import detect_aberrant_events, EventType


class TestSplicingEventDetection:
    """Tests for detect_aberrant_events()."""

    def test_detect_exon_skipping(self, mock_reads_exon_skip, mock_gene):
        """Exon skipping reads should produce an EXON_SKIPPING event."""
        isoforms = collapse_isoforms(mock_reads_exon_skip, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        events = detect_aberrant_events(classified, mock_gene, min_read_support=2)

        skip_events = [e for e in events if e.event_type == EventType.EXON_SKIPPING]
        assert len(skip_events) >= 1
        assert len(skip_events[0].skipped_exons) >= 1

    def test_detect_intron_retention(self, mock_reads_intron_retention, mock_gene):
        """IR reads should produce an INTRON_RETENTION event."""
        isoforms = collapse_isoforms(mock_reads_intron_retention, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        events = detect_aberrant_events(classified, mock_gene, min_read_support=2)

        ir_events = [e for e in events if e.event_type == EventType.INTRON_RETENTION]
        assert len(ir_events) >= 1
        assert ir_events[0].retained_intron is not None

    def test_detect_cryptic_exon(self, mock_reads_cryptic_exon, mock_gene):
        """Cryptic exon reads should produce a CRYPTIC_EXON event."""
        isoforms = collapse_isoforms(mock_reads_cryptic_exon, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        events = detect_aberrant_events(classified, mock_gene, min_read_support=2)

        ce_events = [e for e in events if e.event_type == EventType.CRYPTIC_EXON]
        assert len(ce_events) >= 1
        assert ce_events[0].cryptic_exon is not None

    def test_no_events_for_normal_reads(self, mock_reads_normal, mock_gene):
        """Normal (FSM) reads should not produce any events."""
        isoforms = collapse_isoforms(mock_reads_normal, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        events = detect_aberrant_events(classified, mock_gene)
        assert len(events) == 0

    def test_event_read_support(self, mock_reads_exon_skip, mock_gene):
        """Events should carry correct read support counts."""
        isoforms = collapse_isoforms(mock_reads_exon_skip, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        events = detect_aberrant_events(classified, mock_gene, min_read_support=2)
        for event in events:
            assert event.total_read_support >= 2

    def test_event_aso_amenability(self, mock_reads_exon_skip, mock_gene):
        """Exon skipping events should be flagged as ASO-amenable."""
        isoforms = collapse_isoforms(mock_reads_exon_skip, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        events = detect_aberrant_events(
            classified, mock_gene,
            total_gene_reads=5,
            min_read_support=2,
        )
        skip_events = [e for e in events if e.event_type == EventType.EXON_SKIPPING]
        if skip_events:
            assert skip_events[0].is_aso_target

    def test_frame_impact_assessment(self, mock_reads_exon_skip, mock_gene):
        """Exon skipping events should have reading frame impact assessed."""
        isoforms = collapse_isoforms(mock_reads_exon_skip, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        events = detect_aberrant_events(classified, mock_gene, min_read_support=2)
        skip_events = [e for e in events if e.event_type == EventType.EXON_SKIPPING]
        if skip_events:
            assert skip_events[0].reading_frame_impact in ("in-frame", "frameshift")

    def test_event_summary(self, mock_reads_exon_skip, mock_gene):
        """Event summary string should be human-readable."""
        isoforms = collapse_isoforms(mock_reads_exon_skip, min_reads=2, gene_name="TEST")
        classified = classify_isoforms(isoforms, mock_gene)
        events = detect_aberrant_events(
            classified, mock_gene,
            total_gene_reads=5,
            min_read_support=2,
        )
        for event in events:
            summary = event.summary()
            assert isinstance(summary, str)
            assert len(summary) > 10

# splicetarget v0.1.0
# Any usage is subject to this software's license.
