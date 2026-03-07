"""Tests for BUSCO summary parsing."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import summarise_busco  # noqa: E402


def read_tsv_row(path: Path) -> dict[str, str]:
    """Read a single-row TSV file."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0]


class SummariseBuscoTestCase(unittest.TestCase):
    """Cover existing-string, numeric reconstruction, and failure cases."""

    def write_json(self, path: Path, payload: dict[str, object]) -> Path:
        """Write a JSON file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_main_preserves_existing_busco_string(self) -> None:
        """Reuse a precomputed one-line summary when it is already present."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = self.write_json(
                tmpdir / "short_summary.json",
                {
                    "results": {
                        "one_line_summary": (
                            "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200"
                        )
                    }
                },
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC1",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "bacillota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(
                row["BUSCO_bacillota_odb12"],
                "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
            )
            self.assertEqual(row["busco_status"], "done")
            self.assertEqual(row["warnings"], "")

    def test_main_builds_busco_string_from_numeric_percentages(self) -> None:
        """Construct a cluster_ani-compatible summary string from numeric JSON values."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = self.write_json(
                tmpdir / "short_summary.json",
                {
                    "results": {
                        "Complete percentage": 96.5,
                        "Single copy percentage": 95.0,
                        "Duplicated percentage": 1.5,
                        "Fragmented percentage": 2.0,
                        "Missing percentage": 1.5,
                        "n_markers": 220,
                    }
                },
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC2",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "mycoplasmatota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(
                row["BUSCO_mycoplasmatota_odb12"],
                "C:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:220",
            )
            self.assertEqual(row["busco_status"], "done")

    def test_main_marks_missing_json_as_failed(self) -> None:
        """Emit an NA row when the BUSCO summary file is unavailable."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = tmpdir / "missing.json"
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC3",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "bacillota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["BUSCO_bacillota_odb12"], "NA")
            self.assertEqual(row["busco_status"], "failed")
            self.assertEqual(row["warnings"], "busco_summary_failed")


if __name__ == "__main__":
    unittest.main()
