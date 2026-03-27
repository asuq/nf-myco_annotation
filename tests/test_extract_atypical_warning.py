"""Tests for metadata atypical-warning extraction."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import extract_atypical_warning  # noqa: E402


class ExtractAtypicalWarningTestCase(unittest.TestCase):
    """Cover per-accession atypical-warning lookup."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_extracts_raw_warning_from_tsv_metadata(self) -> None:
        """Return the raw atypical-warning value for one matching accession."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tAtypical_Warnings",
                        "ACC1\tunverified source organism",
                    ]
                )
                + "\n",
            )

            self.assertEqual(
                extract_atypical_warning.extract_atypical_warning(metadata, "ACC1"),
                "unverified source organism",
            )

    def test_returns_na_for_missing_column_or_accession(self) -> None:
        """Return `NA` when the warning column or accession is absent."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            metadata = self.write_text_file(
                tmpdir / "metadata.csv",
                "\n".join(
                    [
                        "Accession,Organism_Name",
                        "ACC1,One",
                    ]
                )
                + "\n",
            )

            self.assertEqual(
                extract_atypical_warning.extract_atypical_warning(metadata, "ACC2"),
                "NA",
            )
