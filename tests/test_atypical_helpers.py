"""Tests for shared atypical-warning interpretation helpers."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import atypical_helpers  # noqa: E402


class AtypicalHelpersTestCase(unittest.TestCase):
    """Cover shared atypical-warning classification rules."""

    def test_classify_returns_non_atypical_for_missing_values(self) -> None:
        """Treat absent or NA-like atypical warnings as non-atypical."""
        self.assertEqual(atypical_helpers.classify_atypical_warning(None), (False, False))
        self.assertEqual(atypical_helpers.classify_atypical_warning("NA"), (False, False))

    def test_classify_marks_non_exception_atypical_values(self) -> None:
        """Treat ordinary atypical-warning values as excluded atypical cases."""
        self.assertEqual(
            atypical_helpers.classify_atypical_warning("cell culture adapted"),
            (True, False),
        )

    def test_classify_marks_unverified_source_organism_as_exception(self) -> None:
        """Treat the locked NCBI phrase as an atypical inclusion exception."""
        self.assertEqual(
            atypical_helpers.classify_atypical_warning("unverified source organism"),
            (True, True),
        )

    def test_classify_matches_exception_case_insensitively(self) -> None:
        """Match the locked exception phrase regardless of capitalisation."""
        self.assertEqual(
            atypical_helpers.classify_atypical_warning("UnVeRiFiEd Source Organism"),
            (True, True),
        )

    def test_detect_reads_warning_value_from_metadata_row(self) -> None:
        """Read and classify the warning value via the metadata column helper."""
        metadata_row = {
            "Accession": "ACC1",
            "Atypical_Warnings": "unverified source organism",
        }

        flags = atypical_helpers.detect_atypical_flags(
            metadata_row,
            find_column_by_normalised_name=lambda header, target: (
                "Atypical_Warnings" if target == "Atypical_Warnings" else None
            ),
        )

        self.assertEqual(flags, (True, True))
