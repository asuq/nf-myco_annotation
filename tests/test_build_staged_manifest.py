"""Tests for the staged-genome manifest builder CLI."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "build_staged_manifest.py"


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class BuildStagedManifestTestCase(unittest.TestCase):
    """Cover staged-manifest building and validation."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write text to a temporary test file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def run_helper(self, *, input_path: Path, output_path: Path) -> subprocess.CompletedProcess[str]:
        """Run the staged-manifest builder CLI."""
        return subprocess.run(
            [sys.executable, str(SCRIPT), "--input", str(input_path), "--output", str(output_path)],
            check=False,
            capture_output=True,
            text=True,
            cwd=ROOT,
        )

    def test_writes_header_and_accession_sorted_rows(self) -> None:
        """Write one manifest header and sort rows by accession."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            input_path = self.write_text_file(
                tmpdir / "staged_manifest_rows.tsv",
                "\n".join(
                    [
                        "ACC_B\tid_b\tACC_B.fasta",
                        "ACC_A\tid_a\tACC_A.fasta",
                    ]
                )
                + "\n",
            )
            output_path = tmpdir / "staged_genomes.tsv"

            result = self.run_helper(input_path=input_path, output_path=output_path)

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                output_path.read_text(encoding="utf-8").splitlines()[0],
                "accession\tinternal_id\tstaged_filename",
            )
            self.assertEqual(
                [row["accession"] for row in read_tsv(output_path)],
                ["ACC_A", "ACC_B"],
            )

    def test_rejects_duplicate_accessions(self) -> None:
        """Fail when the raw input repeats an accession."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            input_path = self.write_text_file(
                tmpdir / "staged_manifest_rows.tsv",
                "\n".join(
                    [
                        "ACC_A\tid_a\tACC_A.fasta",
                        "ACC_A\tid_a2\tACC_A_2.fasta",
                    ]
                )
                + "\n",
            )
            output_path = tmpdir / "staged_genomes.tsv"

            result = self.run_helper(input_path=input_path, output_path=output_path)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("duplicate accession", result.stderr)

    def test_rejects_empty_required_fields(self) -> None:
        """Fail when one required field is empty."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            input_path = self.write_text_file(
                tmpdir / "staged_manifest_rows.tsv",
                "ACC_A\t\tACC_A.fasta\n",
            )
            output_path = tmpdir / "staged_genomes.tsv"

            result = self.run_helper(input_path=input_path, output_path=output_path)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("empty internal_id", result.stderr)

    def test_rejects_malformed_rows(self) -> None:
        """Fail when a row does not have exactly three TSV fields."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            input_path = self.write_text_file(
                tmpdir / "staged_manifest_rows.tsv",
                "ACC_A\tid_a\n",
            )
            output_path = tmpdir / "staged_genomes.tsv"

            result = self.run_helper(input_path=input_path, output_path=output_path)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("expected 3 tab-delimited fields", result.stderr)


if __name__ == "__main__":
    unittest.main()
