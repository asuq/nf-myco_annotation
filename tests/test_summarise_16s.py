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

import concat_best_16s  # noqa: E402
import summarise_16s  # noqa: E402


def read_status_row(path: Path) -> dict[str, str]:
    """Read a single-row status TSV."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0]


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Summarise16STestCase(unittest.TestCase):
    """Test status calling and best-hit selection."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def barrnap_header(
        self,
        rrna_name: str,
        contig: str,
        start0: int,
        end0: int,
        strand: str,
    ) -> str:
        """Build one Barrnap-style FASTA header."""
        return f"{rrna_name}::{contig}:{start0}-{end0}({strand})"

    def write_barrnap_fasta(
        self,
        path: Path,
        records: list[tuple[str, str]],
    ) -> Path:
        """Write a Barrnap FASTA file from header and sequence pairs."""
        lines: list[str] = []
        for header, sequence in records:
            lines.extend((f">{header}", sequence))
        return self.write_text_file(path, "\n".join(lines) + "\n")

    def test_main_prefers_lowest_score_then_longest_sequence(self) -> None:
        """Choose the best intact 16S hit with the locked tie-break order."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            best_header = self.barrnap_header("16S_rRNA", "contig1", 119, 260, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "\n".join(
                    [
                        "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA",
                        "contig1\tbarrnap\trRNA\t120\t260\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA",
                        "contig1\tbarrnap\trRNA\t300\t360\t7.0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA",
                    ]
                )
                + "\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [
                    (self.barrnap_header("16S_rRNA", "contig1", 0, 100, "+"), "A" * 100),
                    (best_header, "C" * 141),
                    (self.barrnap_header("23S_rRNA", "contig1", 299, 360, "+"), "G" * 61),
                ],
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
            self.assertEqual(status_row["best_16S_header"], best_header)
            self.assertEqual(status_row["best_16S_length"], "141")
            self.assertEqual(status_row["include_in_all_best_16S"], "true")
            best_fasta = (outdir / "best_16S.fna").read_text(encoding="utf-8")
            self.assertTrue(best_fasta.startswith(f">{best_header}\n"))

    def test_main_reports_partial_when_only_partial_hits_exist(self) -> None:
        """Emit `partial` and keep the sample out of the intact cohort."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            partial_header = self.barrnap_header("16S_rRNA", "KB890753.1", 102, 1222, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                (
                    "KB890753.1\tbarrnap\trRNA\t103\t1222\t5.7e-275\t+\t.\t"
                    "Name=16S_rRNA;product=16S ribosomal RNA (partial);"
                    "note=aligned only 70 percent of the 16S ribosomal RNA\n"
                ),
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(partial_header, "A" * 1120)],
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
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "partial")
            self.assertEqual(status_row["best_16S_header"], partial_header)
            self.assertEqual(status_row["include_in_all_best_16S"], "false")

    def test_main_uses_atypical_warnings_to_exclude_from_cohort(self) -> None:
        """Derive atypical exclusion from the design-spec Atypical_Warnings field."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            intact_header = self.barrnap_header("16S_rRNA", "contig1", 0, 1500, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t1500\t2.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(intact_header, "A" * 1500)],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC2B",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                    "--atypical-warnings",
                    "cell culture adapted",
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "Yes")
            self.assertEqual(status_row["best_16S_header"], intact_header)
            self.assertEqual(status_row["include_in_all_best_16S"], "false")

    def test_main_allows_unverified_source_exception_in_intact_cohort(self) -> None:
        """Keep the locked unverified-source exception in the intact cohort."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            intact_header = self.barrnap_header("16S_rRNA", "contig1", 0, 1500, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t1500\t2.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(intact_header, "A" * 1500)],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC2C",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                    "--atypical-warnings",
                    "Unverified source organism",
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["include_in_all_best_16S"], "true")

    def test_main_uses_metadata_lookup_for_intact_exception(self) -> None:
        """Resolve the intact-cohort exception from CSV metadata when requested."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            intact_header = self.barrnap_header("16S_rRNA", "contig1", 0, 1500, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t1500\t2.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(intact_header, "A" * 1500)],
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.csv",
                "Accession,Atypical_Warnings\nACC2D,unverified source organism\n",
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC2D",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                    "--metadata",
                    str(metadata),
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["include_in_all_best_16S"], "true")

    def test_main_prefers_explicit_warnings_over_metadata_lookup(self) -> None:
        """Honour explicit atypical warnings before the metadata lookup."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            intact_header = self.barrnap_header("16S_rRNA", "contig1", 0, 1500, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t1500\t2.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(intact_header, "A" * 1500)],
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC2E\tunverified source organism\n",
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC2E",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                    "--metadata",
                    str(metadata),
                    "--atypical-warnings",
                    "cell culture adapted",
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["include_in_all_best_16S"], "false")

    def test_main_falls_back_when_metadata_does_not_resolve_sample(self) -> None:
        """Treat the sample as non-atypical when metadata lookup finds no row."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            intact_header = self.barrnap_header("16S_rRNA", "contig1", 0, 1500, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t1500\t2.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(intact_header, "A" * 1500)],
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nOTHER\tcell culture adapted\n",
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC2F",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                    "--metadata",
                    str(metadata),
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["include_in_all_best_16S"], "true")

    def test_main_returns_na_when_gff_has_16s_but_fasta_does_not(self) -> None:
        """Gracefully degrade to NA when a 16S GFF hit is missing from FASTA."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "\n".join(
                    [
                        "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA",
                        "contig1\tbarrnap\trRNA\t120\t220\t7.0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA",
                    ]
                )
                + "\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(self.barrnap_header("23S_rRNA", "contig1", 119, 220, "+"), "C" * 101)],
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
            first_header = self.barrnap_header("16S_rRNA", "contig1", 0, 100, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "\n".join(
                    [
                        "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA",
                        "contig1\tbarrnap\trRNA\t200\t299\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA",
                    ]
                )
                + "\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [
                    (first_header, "A" * 100),
                    (self.barrnap_header("16S_rRNA", "contig1", 199, 299, "+"), "C" * 100),
                ],
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
            self.assertEqual(status_row["best_16S_header"], first_header)

    def test_main_reports_no_when_no_16s_hits_exist(self) -> None:
        """Emit `No` and no cohort inclusion when Barrnap found only non-16S rRNA."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(self.barrnap_header("23S_rRNA", "contig1", 0, 100, "+"), "A" * 100)],
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

    def test_main_matches_mixed_rrna_order_by_key(self) -> None:
        """Match 16S records by Barrnap key even when mixed rRNA order differs."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            best_header = self.barrnap_header("16S_rRNA", "contig1", 319, 470, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "\n".join(
                    [
                        "contig1\tbarrnap\trRNA\t1\t100\t8.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA",
                        "contig1\tbarrnap\trRNA\t120\t220\t1.0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA",
                        "contig1\tbarrnap\trRNA\t240\t300\t1.0\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA",
                        "contig1\tbarrnap\trRNA\t320\t470\t4.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA",
                    ]
                )
                + "\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [
                    (self.barrnap_header("23S_rRNA", "contig1", 119, 220, "+"), "G" * 101),
                    (best_header, "C" * 151),
                    (self.barrnap_header("5S_rRNA", "contig1", 239, 300, "+"), "T" * 61),
                    (self.barrnap_header("16S_rRNA", "contig1", 0, 100, "+"), "A" * 100),
                ],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC6",
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
            self.assertEqual(status_row["best_16S_header"], best_header)
            best_fasta = (outdir / "best_16S.fna").read_text(encoding="utf-8")
            self.assertTrue(best_fasta.startswith(f">{best_header}\n"))

    def test_main_prefers_intact_hit_over_better_scoring_partial(self) -> None:
        """Prefer any intact 16S hit over a partial hit with a better score."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            intact_header = self.barrnap_header("16S_rRNA", "contig1", 1999, 3500, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "\n".join(
                    [
                        (
                            "contig1\tbarrnap\trRNA\t103\t1222\t1.0\t+\t.\t"
                            "Name=16S_rRNA;product=16S ribosomal RNA (partial);"
                            "note=aligned only 70 percent of the 16S ribosomal RNA"
                        ),
                        (
                            "contig1\tbarrnap\trRNA\t2000\t3500\t9.0\t+\t.\t"
                            "Name=16S_rRNA;product=16S ribosomal RNA"
                        ),
                    ]
                )
                + "\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [
                    (self.barrnap_header("16S_rRNA", "contig1", 102, 1222, "+"), "A" * 1120),
                    (intact_header, "C" * 1501),
                ],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC7",
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
            self.assertEqual(status_row["best_16S_header"], intact_header)
            self.assertEqual(status_row["include_in_all_best_16S"], "true")

    def test_main_accepts_negative_strand_barrnap_header(self) -> None:
        """Parse a Barrnap negative-strand 16S header and match it to GFF."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            minus_header = self.barrnap_header("16S_rRNA", "contig2", 499, 1650, "-")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                (
                    "contig2\tbarrnap\trRNA\t500\t1650\t2.0\t-\t.\t"
                    "Name=16S_rRNA;product=16S ribosomal RNA\n"
                ),
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(minus_header, "G" * 1151)],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC8",
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
            self.assertEqual(status_row["best_16S_header"], minus_header)

    def test_main_returns_na_when_fasta_has_16s_but_gff_does_not(self) -> None:
        """Gracefully degrade to NA when a 16S FASTA record is missing from GFF."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t120\t220\t7.0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(self.barrnap_header("16S_rRNA", "contig1", 0, 100, "+"), "A" * 100)],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC9",
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
            self.assertEqual((outdir / "best_16S.fna").read_text(encoding="utf-8"), "")

    def test_main_returns_na_on_duplicate_barrnap_fasta_key(self) -> None:
        """Gracefully degrade to NA when Barrnap FASTA contains duplicate keys."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            duplicate_header = self.barrnap_header("16S_rRNA", "contig1", 0, 100, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [
                    (duplicate_header, "A" * 100),
                    (duplicate_header, "C" * 100),
                ],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC10",
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
            self.assertEqual((outdir / "best_16S.fna").read_text(encoding="utf-8"), "")

    def test_concat_best_16s_excludes_partial_rows_from_intact_cohort(self) -> None:
        """Build the intact cohort from complete, explicitly included rows only."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)

            include_dir = tmpdir / "include"
            exclude_dir = tmpdir / "partial"
            include_dir.mkdir()
            exclude_dir.mkdir()

            include_status = self.write_text_file(
                include_dir / "16S_status.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings",
                        "ACC_OK\tYes\thit_ok\t100\ttrue\t",
                    ]
                )
                + "\n",
            )
            include_best = self.write_text_file(
                include_dir / "best_16S.fna",
                ">hit_ok\n" + "A" * 100 + "\n",
            )

            exclude_status = self.write_text_file(
                exclude_dir / "16S_status.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings",
                        "ACC_PART\tpartial\thit_partial\t90\tfalse\t",
                    ]
                )
                + "\n",
            )
            exclude_best = self.write_text_file(
                exclude_dir / "best_16S.fna",
                ">hit_partial\n" + "C" * 90 + "\n",
            )

            cohort_inputs = self.write_text_file(
                tmpdir / "cohort_inputs.tsv",
                "\n".join(
                    [
                        "accession\tstatus_tsv\tbest_16s_fasta",
                        f"ACC_OK\t{include_status}\t{include_best}",
                        f"ACC_PART\t{exclude_status}\t{exclude_best}",
                    ]
                )
                + "\n",
            )
            output_fasta = tmpdir / "all_best_16S.fna"
            output_manifest = tmpdir / "all_best_16S_manifest.tsv"

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--output-fasta",
                    str(output_fasta),
                    "--output-manifest",
                    str(output_manifest),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertEqual(
                output_fasta.read_text(encoding="utf-8"),
                ">hit_ok\n" + "A" * 100 + "\n",
            )
            manifest_rows = read_tsv_rows(output_manifest)
            self.assertEqual(len(manifest_rows), 1)
            self.assertEqual(manifest_rows[0]["accession"], "ACC_OK")

    def test_concat_best_16s_builds_partial_cohort_including_atypical_samples(self) -> None:
        """Build the partial cohort from all partial rows, including atypical ones."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)

            intact_dir = tmpdir / "intact"
            partial_dir = tmpdir / "partial"
            atypical_partial_dir = tmpdir / "atypical_partial"
            intact_dir.mkdir()
            partial_dir.mkdir()
            atypical_partial_dir.mkdir()

            intact_status = self.write_text_file(
                intact_dir / "16S_status.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings",
                        "ACC_OK\tYes\thit_ok\t100\ttrue\t",
                    ]
                )
                + "\n",
            )
            intact_best = self.write_text_file(
                intact_dir / "best_16S.fna",
                ">hit_ok\n" + "A" * 100 + "\n",
            )

            partial_status = self.write_text_file(
                partial_dir / "16S_status.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings",
                        "ACC_PART\tpartial\thit_partial\t90\tfalse\t",
                    ]
                )
                + "\n",
            )
            partial_best = self.write_text_file(
                partial_dir / "best_16S.fna",
                ">hit_partial\n" + "C" * 90 + "\n",
            )

            atypical_partial_status = self.write_text_file(
                atypical_partial_dir / "16S_status.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings",
                        "ACC_ATYP_PART\tpartial\thit_atyp_partial\t80\tfalse\tatypical",
                    ]
                )
                + "\n",
            )
            atypical_partial_best = self.write_text_file(
                atypical_partial_dir / "best_16S.fna",
                ">hit_atyp_partial\n" + "G" * 80 + "\n",
            )

            cohort_inputs = self.write_text_file(
                tmpdir / "cohort_inputs.tsv",
                "\n".join(
                    [
                        "accession\tstatus_tsv\tbest_16s_fasta",
                        f"ACC_OK\t{intact_status}\t{intact_best}",
                        f"ACC_PART\t{partial_status}\t{partial_best}",
                        (
                            "ACC_ATYP_PART\t"
                            f"{atypical_partial_status}\t{atypical_partial_best}"
                        ),
                    ]
                )
                + "\n",
            )
            output_fasta = tmpdir / "all_partial_16S.fna"
            output_manifest = tmpdir / "all_partial_16S_manifest.tsv"

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--cohort-kind",
                    "partial",
                    "--output-fasta",
                    str(output_fasta),
                    "--output-manifest",
                    str(output_manifest),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertEqual(
                output_fasta.read_text(encoding="utf-8"),
                ">hit_partial\n"
                + "C" * 90
                + "\n>hit_atyp_partial\n"
                + "G" * 80
                + "\n",
            )
            manifest_rows = read_tsv_rows(output_manifest)
            self.assertEqual(
                [row["accession"] for row in manifest_rows],
                ["ACC_PART", "ACC_ATYP_PART"],
            )

    def test_concat_best_16s_fails_when_included_fasta_is_empty(self) -> None:
        """Fail when a cohort-included sample has no FASTA content to append."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            status_path = self.write_text_file(
                tmpdir / "16S_status.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings",
                        "ACC_BAD\tYes\thit_bad\t100\ttrue\t",
                    ]
                )
                + "\n",
            )
            best_fasta = self.write_text_file(tmpdir / "best_16S.fna", "")
            cohort_inputs = self.write_text_file(
                tmpdir / "cohort_inputs.tsv",
                "\n".join(
                    [
                        "accession\tstatus_tsv\tbest_16s_fasta",
                        f"ACC_BAD\t{status_path}\t{best_fasta}",
                    ]
                )
                + "\n",
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 1)


if __name__ == "__main__":
    unittest.main()
