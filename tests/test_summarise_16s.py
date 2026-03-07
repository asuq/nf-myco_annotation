"""Tests for Barrnap-derived 16S summarisation."""

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

import summarise_16s  # noqa: E402


def read_status_row(path: Path) -> dict[str, str]:
    """Read a single-row status TSV."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0]


class Summarise16STestCase(unittest.TestCase):
    """Test status calling and best-hit selection."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_main_prefers_lowest_score_then_longest_sequence(self) -> None:
        """Choose the best intact 16S hit with the locked tie-break order."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "\n".join(
                    [
                        "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA",
                        "contig1\tbarrnap\trRNA\t120\t260\t5.0\t+\t.\tName=16S_rRNA",
                        "contig1\tbarrnap\trRNA\t300\t360\t7.0\t+\t.\tName=23S_rRNA",
                    ]
                )
                + "\n",
            )
            fasta = self.write_text_file(
                tmpdir / "rrna.fa",
                "\n".join(
                    [
                        ">hit1 16S ribosomal RNA",
                        "A" * 100,
                        ">hit2 16S ribosomal RNA",
                        "C" * 141,
                        ">hit3 23S ribosomal RNA",
                        "G" * 61,
                    ]
                )
                + "\n",
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC1",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "Yes")
            self.assertEqual(status_row["best_16S_header"], "hit2 16S ribosomal RNA")
            self.assertEqual(status_row["best_16S_length"], "141")
            self.assertEqual(status_row["include_in_all_best_16S"], "true")
            best_fasta = (outdir / "best_16S.fna").read_text(encoding="utf-8")
            self.assertTrue(best_fasta.startswith(">hit2 16S ribosomal RNA\n"))

    def test_main_reports_partial_when_only_partial_hits_exist(self) -> None:
        """Emit `partial` when only partial 16S records are available."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t90\t3.0\t+\t.\tName=16S_rRNA;partial=true\n",
            )
            fasta = self.write_text_file(
                tmpdir / "rrna.fa",
                ">hit1 partial 16S ribosomal RNA\n" + "A" * 90 + "\n",
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC2",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                    "--is-atypical",
                    "true",
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "partial")
            self.assertEqual(status_row["include_in_all_best_16S"], "false")

    def test_main_returns_na_for_malformed_barrnap_outputs(self) -> None:
        """Gracefully degrade to NA when GFF and FASTA hit counts disagree."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA\n",
            )
            fasta = self.write_text_file(
                tmpdir / "rrna.fa",
                "\n".join(
                    [
                        ">hit1 16S ribosomal RNA",
                        "A" * 100,
                        ">hit2 23S ribosomal RNA",
                        "C" * 50,
                    ]
                )
                + "\n",
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC3",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "NA")
            self.assertEqual(status_row["warnings"], "invalid_barrnap_output")
            self.assertEqual(
                (outdir / "best_16S.fna").read_text(encoding="utf-8"),
                "",
            )

    def test_main_uses_gff_order_as_final_tie_break(self) -> None:
        """Prefer the first GFF hit when score and length remain tied."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "\n".join(
                    [
                        "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA",
                        "contig1\tbarrnap\trRNA\t200\t299\t5.0\t+\t.\tName=16S_rRNA",
                    ]
                )
                + "\n",
            )
            fasta = self.write_text_file(
                tmpdir / "rrna.fa",
                "\n".join(
                    [
                        ">hit1 16S ribosomal RNA",
                        "A" * 100,
                        ">hit2 16S ribosomal RNA",
                        "C" * 100,
                    ]
                )
                + "\n",
            )
            outdir = tmpdir / "out"

            summarise_16s.main(
                [
                    "--accession",
                    "ACC4",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["best_16S_header"], "hit1 16S ribosomal RNA")

    def test_main_reports_no_when_no_16s_hits_exist(self) -> None:
        """Emit `No` and no cohort inclusion when Barrnap found only non-16S rRNA."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=23S_rRNA\n",
            )
            fasta = self.write_text_file(
                tmpdir / "rrna.fa",
                ">hit1 23S ribosomal RNA\n" + "A" * 100 + "\n",
            )
            outdir = tmpdir / "out"

            summarise_16s.main(
                [
                    "--accession",
                    "ACC5",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "No")
            self.assertEqual(status_row["best_16S_header"], "NA")
            self.assertEqual(status_row["include_in_all_best_16S"], "false")


if __name__ == "__main__":
    unittest.main()
