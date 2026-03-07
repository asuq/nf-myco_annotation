"""Tests for CRISPRCasFinder result summarisation."""

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

import summarise_ccfinder  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into a list of dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SummariseCCFinderTestCase(unittest.TestCase):
    """Cover retained arrays, evidence filtering, and failure handling."""

    def write_json(self, path: Path, payload: dict[str, object]) -> Path:
        """Write a JSON file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_main_builds_strain_contig_and_crispr_tables(self) -> None:
        """Ignore evidence-level 1 arrays and aggregate retained arrays correctly."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            result_json = self.write_json(
                tmpdir / "result.json",
                {
                    "Sequences": [
                        {
                            "Id": "contig1",
                            "Length": 1000,
                            "Crisprs": [
                                {
                                    "Name": "c1_a",
                                    "Evidence_Level": 4,
                                    "Start": 1,
                                    "End": 100,
                                    "Spacers": [{}, {}, {}],
                                },
                                {
                                    "Name": "c1_b",
                                    "Evidence_Level": 1,
                                    "Start": 200,
                                    "End": 260,
                                    "Spacers": [{}, {}],
                                },
                            ],
                        },
                        {
                            "Id": "contig2",
                            "Length": 500,
                            "Crisprs": [
                                {
                                    "Name": "c2_a",
                                    "Evidence_Level": 3,
                                    "Start": 50,
                                    "End": 99,
                                    "Spacers": 4,
                                }
                            ],
                        },
                    ]
                },
            )
            outdir = tmpdir / "out"

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC1",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_rows = read_tsv(outdir / "ccfinder_strains.tsv")
            contig_rows = read_tsv(outdir / "ccfinder_contigs.tsv")
            crispr_rows = read_tsv(outdir / "ccfinder_crisprs.tsv")

            self.assertEqual(len(strain_rows), 1)
            self.assertEqual(strain_rows[0]["CRISPRS"], "2")
            self.assertEqual(strain_rows[0]["SPACERS_SUM"], "7")
            self.assertEqual(strain_rows[0]["CRISPR_FRAC"], "0.1")
            self.assertEqual(strain_rows[0]["ccfinder_status"], "done")

            self.assertEqual(len(contig_rows), 2)
            self.assertEqual(contig_rows[0]["CRISPRS"], "1")
            self.assertEqual(contig_rows[1]["SPACERS_SUM"], "4")

            self.assertEqual(len(crispr_rows), 2)
            self.assertEqual(crispr_rows[0]["crispr_id"], "c1_a")
            self.assertEqual(crispr_rows[1]["crispr_length"], "50")

    def test_main_handles_samples_without_retained_crisprs(self) -> None:
        """Write zero-valued summaries when only evidence-level 1 arrays exist."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            result_json = self.write_json(
                tmpdir / "result.json",
                {
                    "Sequences": [
                        {
                            "Id": "contig1",
                            "Length": 1000,
                            "Crisprs": [
                                {
                                    "Evidence_Level": 1,
                                    "Start": 1,
                                    "End": 100,
                                    "Spacers": 2,
                                }
                            ],
                        }
                    ]
                },
            )
            outdir = tmpdir / "out"

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC2",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_row = read_tsv(outdir / "ccfinder_strains.tsv")[0]
            crispr_rows = read_tsv(outdir / "ccfinder_crisprs.tsv")
            self.assertEqual(strain_row["CRISPRS"], "0")
            self.assertEqual(strain_row["SPACERS_SUM"], "0")
            self.assertEqual(strain_row["CRISPR_FRAC"], "0")
            self.assertEqual(crispr_rows, [])

    def test_main_marks_invalid_json_as_failed(self) -> None:
        """Emit NA-like strain output and empty detail tables on failure."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            outdir = tmpdir / "out"
            result_json = tmpdir / "missing.json"

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC3",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_row = read_tsv(outdir / "ccfinder_strains.tsv")[0]
            self.assertEqual(strain_row["CRISPRS"], "NA")
            self.assertEqual(strain_row["ccfinder_status"], "failed")
            self.assertEqual(strain_row["warnings"], "ccfinder_summary_failed")
            self.assertEqual(read_tsv(outdir / "ccfinder_contigs.tsv"), [])
            self.assertEqual(read_tsv(outdir / "ccfinder_crisprs.tsv"), [])


if __name__ == "__main__":
    unittest.main()
