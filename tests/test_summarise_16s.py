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
        """Emit `partial` and keep the sample out of the intact cohort."""
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
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "partial")
            self.assertEqual(status_row["include_in_all_best_16S"], "false")

    def test_main_uses_atypical_warnings_to_exclude_from_cohort(self) -> None:
        """Derive atypical exclusion from the design-spec Atypical_Warnings field."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t1500\t2.0\t+\t.\tName=16S_rRNA\n",
            )
            fasta = self.write_text_file(
                tmpdir / "rrna.fa",
                ">hit1 16S ribosomal RNA\n" + "A" * 1500 + "\n",
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
