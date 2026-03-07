"""Tests for CheckM2 summarisation."""

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

import summarise_checkm2  # noqa: E402


def read_tsv_row(path: Path) -> dict[str, str]:
    """Read a single-row TSV file."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0]


class SummariseCheckM2TestCase(unittest.TestCase):
    """Cover successful and degraded CheckM2 merge behaviour."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def write_report(
        self,
        path: Path,
        *,
        completeness: str,
        contamination: str,
        coding_density: str = "0.9",
        average_gene_length: str = "900",
        total_coding_sequences: str = "800",
        genome_size: str = "1000000",
        gc_content: str = "0.3",
    ) -> Path:
        """Write a minimal one-row CheckM2 TSV."""
        content = "\n".join(
            [
                "\t".join(
                    [
                        "Name",
                        "Completeness",
                        "Contamination",
                        "Coding_Density",
                        "Average_Gene_Length",
                        "Total_Coding_Sequences",
                        "Genome_Size",
                        "GC_Content",
                    ]
                ),
                "\t".join(
                    [
                        "sample",
                        completeness,
                        contamination,
                        coding_density,
                        average_gene_length,
                        total_coding_sequences,
                        genome_size,
                        gc_content,
                    ]
                ),
            ]
        )
        return self.write_text_file(path, content + "\n")

    def test_main_assigns_gcode4_and_low_quality_false(self) -> None:
        """Prefer translation table 4 when completeness differs by more than 10."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="95", contamination="2")
            report11 = self.write_report(tmpdir / "g11.tsv", completeness="80", contamination="1")
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC1",
                    "--gcode4-report",
                    str(report4),
                    "--gcode11-report",
                    str(report11),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Gcode"], "4")
            self.assertEqual(row["Low_quality"], "false")
            self.assertEqual(row["checkm2_status"], "done")
            self.assertEqual(row["warnings"], "")

    def test_main_marks_ambiguous_gcode_as_na(self) -> None:
        """Emit gcode NA when completeness differences do not cross the threshold."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="88", contamination="4")
            report11 = self.write_report(tmpdir / "g11.tsv", completeness="82", contamination="2")
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC2",
                    "--gcode4-report",
                    str(report4),
                    "--gcode11-report",
                    str(report11),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Gcode"], "NA")
            self.assertEqual(row["Low_quality"], "NA")
            self.assertEqual(row["checkm2_status"], "done")
            self.assertEqual(row["warnings"], "gcode_na")

    def test_main_keeps_partial_metrics_when_one_report_is_missing(self) -> None:
        """Write a summary row even when one CheckM2 report cannot be parsed."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="92", contamination="3")
            missing_report = tmpdir / "g11.tsv"
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC3",
                    "--gcode4-report",
                    str(report4),
                    "--gcode11-report",
                    str(missing_report),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Completeness_gcode4"], "92")
            self.assertEqual(row["Completeness_gcode11"], "NA")
            self.assertEqual(row["Gcode"], "NA")
            self.assertEqual(row["Low_quality"], "NA")
            self.assertEqual(row["checkm2_status"], "failed")
            self.assertEqual(row["warnings"], "checkm2_gcode11_failed")

    def test_main_keeps_gcode_na_at_exact_threshold_difference(self) -> None:
        """Do not assign a gcode when completeness differs by exactly 10."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="90", contamination="2")
            report11 = self.write_report(tmpdir / "g11.tsv", completeness="80", contamination="1")
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC4",
                    "--gcode4-report",
                    str(report4),
                    "--gcode11-report",
                    str(report11),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Gcode"], "NA")
            self.assertEqual(row["warnings"], "gcode_na")

    def test_main_marks_low_quality_true_at_boundary_score(self) -> None:
        """Use the locked `<= 50` cutoff for low-quality calls."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="70", contamination="4")
            report11 = self.write_report(tmpdir / "g11.tsv", completeness="50", contamination="2")
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC5",
                    "--gcode4-report",
                    str(report4),
                    "--gcode11-report",
                    str(report11),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Gcode"], "4")
            self.assertEqual(row["Low_quality"], "true")

    def test_main_marks_inconsistent_shared_stats_as_failed(self) -> None:
        """Treat mismatched shared assembly statistics as a failed CheckM2 merge."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="95", contamination="1")
            report11 = self.write_report(
                tmpdir / "g11.tsv",
                completeness="80",
                contamination="1",
                genome_size="999999",
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC6",
                    "--gcode4-report",
                    str(report4),
                    "--gcode11-report",
                    str(report11),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Gcode"], "NA")
            self.assertEqual(row["checkm2_status"], "failed")
            self.assertEqual(row["warnings"], "inconsistent_shared_stats")


if __name__ == "__main__":
    unittest.main()
