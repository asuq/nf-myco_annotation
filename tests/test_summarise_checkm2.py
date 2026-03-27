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
        contig_n50: str = "50000",
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
                        "Contig_N50",
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
                        contig_n50,
                    ]
                ),
            ]
        )
        return self.write_text_file(path, content + "\n")

    def test_main_writes_stable_output_columns_and_paired_metrics(self) -> None:
        """Emit the locked column order and preserve both paired metric sets."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(
                tmpdir / "g4.tsv",
                completeness="95.5",
                contamination="2.5",
                coding_density="0.91",
                average_gene_length="901",
                total_coding_sequences="801",
            )
            report11 = self.write_report(
                tmpdir / "g11.tsv",
                completeness="80.25",
                contamination="1.5",
                coding_density="0.83",
                average_gene_length="851",
                total_coding_sequences="781",
            )
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
            with output.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames, list(summarise_checkm2.OUTPUT_COLUMNS))
                rows = list(reader)

            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["Completeness_gcode4"], "95.5")
            self.assertEqual(row["Contamination_gcode4"], "2.5")
            self.assertEqual(row["Coding_Density_gcode4"], "0.91")
            self.assertEqual(row["Average_Gene_Length_gcode4"], "901")
            self.assertEqual(row["Total_Coding_Sequences_gcode4"], "801")
            self.assertEqual(row["Completeness_gcode11"], "80.25")
            self.assertEqual(row["Contamination_gcode11"], "1.5")
            self.assertEqual(row["Coding_Density_gcode11"], "0.83")
            self.assertEqual(row["Average_Gene_Length_gcode11"], "851")
            self.assertEqual(row["Total_Coding_Sequences_gcode11"], "781")
            self.assertEqual(row["Gcode"], "4")
            self.assertEqual(row["Low_quality"], "false")
            self.assertEqual(row["checkm2_status"], "done")
            self.assertEqual(row["warnings"], "")

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

    def test_main_assigns_gcode11_when_table_11_is_clearly_better(self) -> None:
        """Prefer translation table 11 when its completeness wins by more than 10."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="78", contamination="1")
            report11 = self.write_report(tmpdir / "g11.tsv", completeness="92", contamination="3")
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC11",
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
            self.assertEqual(row["Gcode"], "11")
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

    def test_main_assigns_gcode11_for_valid_close_scores_with_delta_then_11(self) -> None:
        """Use the fallback rule to assign translation table 11 for close valid pairs."""
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
                    "--gcode-rule",
                    summarise_checkm2.DELTA_THEN_ELEVEN_RULE,
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Gcode"], "11")
            self.assertEqual(row["Low_quality"], "false")
            self.assertEqual(row["checkm2_status"], "done")
            self.assertEqual(row["warnings"], "")

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

    def test_main_marks_reports_without_complete_shared_stats_as_inconsistent(self) -> None:
        """Require the full shared-stat set before assigning gcode."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_text_file(
                tmpdir / "g4.tsv",
                "\n".join(
                    [
                        "Completeness\tContamination\tCoding_Density\tAverage_Gene_Length\tTotal_Coding_Sequences",
                        "95\t2\t0.9\t900\t800",
                    ]
                )
                + "\n",
            )
            report11 = self.write_text_file(
                tmpdir / "g11.tsv",
                "\n".join(
                    [
                        "Completeness\tContamination\tCoding_Density\tAverage_Gene_Length\tTotal_Coding_Sequences",
                        "80\t1\t0.8\t850\t780",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC_SHARED",
                    "--gcode4-report",
                    str(report4),
                    "--gcode11-report",
                    str(report11),
                    "--gcode-rule",
                    summarise_checkm2.DELTA_THEN_ELEVEN_RULE,
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Gcode"], "NA")
            self.assertEqual(row["Low_quality"], "NA")
            self.assertEqual(row["checkm2_status"], "failed")
            self.assertEqual(row["warnings"], "inconsistent_shared_stats")

    def test_main_keeps_gcode_na_at_exact_threshold_difference_under_strict_rule(self) -> None:
        """Keep the strict rule behaviour at an exact threshold difference."""
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

    def test_main_assigns_gcode11_at_exact_threshold_difference_with_fallback_rule(self) -> None:
        """Use translation table 11 when the fallback rule sees an exact threshold tie."""
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
                    "--gcode-rule",
                    summarise_checkm2.DELTA_THEN_ELEVEN_RULE,
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["Gcode"], "11")
            self.assertEqual(row["Low_quality"], "false")
            self.assertEqual(row["warnings"], "")

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

    def test_main_uses_chosen_gcode11_metrics_for_low_quality(self) -> None:
        """Compute low-quality status from the chosen table 11 metrics only."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="20", contamination="0")
            report11 = self.write_report(tmpdir / "g11.tsv", completeness="70", contamination="5")
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC7",
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
            self.assertEqual(row["Gcode"], "11")
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

    def test_main_marks_missing_shared_stats_as_failed(self) -> None:
        """Treat asymmetric shared-stat fields as an invalid paired report set."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            report4 = self.write_report(tmpdir / "g4.tsv", completeness="95", contamination="1")
            report11 = self.write_text_file(
                tmpdir / "g11.tsv",
                "\n".join(
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
                            ]
                        ),
                        "\t".join(
                            [
                                "sample",
                                "80",
                                "1",
                                "0.8",
                                "850",
                                "780",
                                "1000000",
                            ]
                        ),
                    ]
                )
                + "\n",
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_checkm2.main(
                [
                    "--accession",
                    "ACC8",
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
            self.assertEqual(row["checkm2_status"], "failed")
            self.assertEqual(row["warnings"], "inconsistent_shared_stats")


if __name__ == "__main__":
    unittest.main()
