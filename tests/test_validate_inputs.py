"""Tests for the input-validation CLI."""

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

import validate_inputs  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ValidateInputsTestCase(unittest.TestCase):
    """Exercise success and failure cases for input validation."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_main_writes_expected_outputs_and_collision_safe_ids(self) -> None:
        """Validate manifest success, metadata fallback warnings, and collision suffixes."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_one = self.write_text_file(tmpdir / "genomes" / "a.fna", ">a\nACGT\n")
            fasta_two = self.write_text_file(tmpdir / "genomes" / "b.fna", ">b\nTGCA\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta,organism_name",
                        f"ACC-1,false,,{fasta_one},Known one",
                        f"ACC 1,true,Scaffold,{fasta_two},Candidate two",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "\n".join(
                    [
                        "accession,Tax_ID,Organism_Name",
                        "ACC-1,123,Known one",
                    ]
                )
                + "\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            validated_rows = read_tsv(outdir / "validated_samples.tsv")
            accession_map_rows = read_tsv(outdir / "accession_map.tsv")
            warning_rows = read_tsv(outdir / "validation_warnings.tsv")

            self.assertEqual(len(validated_rows), 2)
            self.assertEqual(validated_rows[0]["is_new"], "false")
            self.assertEqual(validated_rows[1]["is_new"], "true")
            self.assertEqual(validated_rows[0]["assembly_level"], "NA")
            self.assertEqual(validated_rows[1]["assembly_level"], "Scaffold")
            self.assertEqual(validated_rows[0]["genome_fasta"], str(fasta_one.resolve()))
            self.assertEqual(validated_rows[1]["genome_fasta"], str(fasta_two.resolve()))

            internal_ids = {row["accession"]: row["internal_id"] for row in accession_map_rows}
            self.assertEqual(len(set(internal_ids.values())), 2)
            self.assertTrue(internal_ids["ACC-1"].startswith("ACC_1_"))
            self.assertTrue(internal_ids["ACC 1"].startswith("ACC_1_"))

            warning_codes = {
                (row["accession"], row["warning_code"]) for row in warning_rows
            }
            self.assertIn(
                ("ACC 1", "missing_metadata_for_new_sample"),
                warning_codes,
            )
            self.assertIn(
                ("ACC-1", "internal_id_collision_resolved"),
                warning_codes,
            )
            self.assertIn(
                ("ACC 1", "internal_id_collision_resolved"),
                warning_codes,
            )

    def test_main_rejects_existing_samples_missing_metadata(self) -> None:
        """Fail when an existing sample has no matching metadata row."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nOTHER\t1\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse((outdir / "validated_samples.tsv").exists())

    def test_main_rejects_accessions_with_path_separators(self) -> None:
        """Fail when an accession is unsafe for published output folders."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"BAD/ACC,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nBAD/ACC,1\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse((outdir / "accession_map.tsv").exists())


if __name__ == "__main__":
    unittest.main()
