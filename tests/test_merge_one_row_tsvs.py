"""Tests for deterministic merging of one-row TSV fragments."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import merge_one_row_tsvs  # noqa: E402


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header list and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return list(reader.fieldnames), list(reader)


class MergeOneRowTsvsTestCase(unittest.TestCase):
    """Cover deterministic merging and header validation."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write UTF-8 text to a temporary file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_main_merges_rows_once_with_accession_sort(self) -> None:
        """Write one canonical header and sort merged rows by accession."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            header = self.write_text_file(
                tmpdir / "header.tsv",
                "accession\tstatus\twarnings\n",
            )
            row_b = self.write_text_file(
                tmpdir / "row_b.tsv",
                "accession\tstatus\twarnings\nACC_B\tdone\t\n",
            )
            row_a = self.write_text_file(
                tmpdir / "row_a.tsv",
                "accession\tstatus\twarnings\nACC_A\tfailed\twarn_a\n",
            )
            output = tmpdir / "merged.tsv"

            exit_code = merge_one_row_tsvs.main(
                [
                    "--header",
                    str(header),
                    "--input",
                    str(row_b),
                    "--input",
                    str(row_a),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            header_row, rows = read_tsv(output)
            self.assertEqual(header_row, ["accession", "status", "warnings"])
            self.assertEqual(
                [row["accession"] for row in rows],
                ["ACC_A", "ACC_B"],
            )

    def test_main_rejects_header_mismatch(self) -> None:
        """Fail when a fragment header does not match the canonical header exactly."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            header = self.write_text_file(
                tmpdir / "header.tsv",
                "accession\tstatus\twarnings\n",
            )
            mismatched = self.write_text_file(
                tmpdir / "row.tsv",
                "accession\twarnings\tstatus\nACC_A\twarn_a\tfailed\n",
            )

            exit_code = merge_one_row_tsvs.main(
                [
                    "--header",
                    str(header),
                    "--input",
                    str(mismatched),
                    "--output",
                    str(tmpdir / "merged.tsv"),
                ]
            )

            self.assertEqual(exit_code, 1)


if __name__ == "__main__":
    unittest.main()
