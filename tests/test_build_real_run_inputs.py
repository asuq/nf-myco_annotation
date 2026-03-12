"""Tests for the medium real-run input builder."""

from __future__ import annotations

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "build_real_run_inputs.py"
REQUIRED_METADATA_HEADER = [
    "Accession",
    "Tax_ID",
    "Organism_Name",
    "Assembly_Level",
    "N50",
    "Scaffolds",
    "Genome_Size",
    "Atypical_Warnings",
]
REQUIRED_SAMPLE_HEADER = [
    "accession",
    "is_new",
    "assembly_level",
    "genome_fasta",
]


def write_tsv(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    """Write a TSV file for one test fixture."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read one CSV file and return its header and rows."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        return list(reader.fieldnames), list(reader)


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read one TSV file and return its header and rows."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return list(reader.fieldnames), list(reader)


class BuildRealRunInputsTestCase(unittest.TestCase):
    """Exercise the medium real-run input generator."""

    maxDiff = None

    def run_builder(self, *args: str) -> subprocess.CompletedProcess[str]:
        """Run the builder CLI with captured text output."""
        return subprocess.run(
            ["python3", str(SCRIPT), *args],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def make_row(
        self,
        accession: str,
        phylum: str,
        *,
        role_tags: str,
        is_new: str = "false",
        include_metadata: str = "true",
    ) -> dict[str, str]:
        """Build one candidate TSV row."""
        return {
            "accession": accession,
            "genome_fasta": f"/data/{accession}.fna",
            "is_new": is_new,
            "assembly_level": "Complete Genome",
            "include_metadata": include_metadata,
            "tax_id": f"TAX_{accession}",
            "organism_name": f"Organism {accession}",
            "n50": "1000",
            "scaffolds": "1",
            "genome_size": "1000000",
            "atypical_warnings": "NA",
            "phylum": phylum,
            "role_tags": role_tags,
        }

    def make_valid_rows(self) -> list[dict[str, str]]:
        """Build one valid candidate manifest with extra edge-case coverage."""
        rows = [
            self.make_row(
                "MYCO_G4",
                "Mycoplasmatota",
                role_tags="gcode4_candidate;crispr_positive_candidate",
            ),
            self.make_row(
                "BACI_G11",
                "Bacillota",
                role_tags="gcode11_candidate;crispr_negative_candidate",
            ),
            self.make_row(
                "MYCO_METALESS",
                "Mycoplasmatota",
                role_tags="missing_metadata_case",
                is_new="true",
                include_metadata="false",
            ),
            self.make_row(
                "MYCO_ATYPICAL",
                "Mycoplasmatota",
                role_tags="atypical_excluded_candidate",
            ),
            self.make_row(
                "BACI_EXCEPTION",
                "Bacillota",
                role_tags="atypical_exception_candidate",
            ),
            self.make_row(
                "MYCO_ANI_1",
                "Mycoplasmatota",
                role_tags="ani_cluster_candidate",
            ),
            self.make_row(
                "BACI_ANI_1",
                "Bacillota",
                role_tags="ani_cluster_candidate",
            ),
            self.make_row(
                "BACI-PAIR",
                "Bacillota",
                role_tags="collision_candidate",
            ),
            self.make_row(
                "BACI PAIR",
                "Bacillota",
                role_tags="collision_candidate",
            ),
        ]
        for index in range(1, 10):
            rows.append(
                self.make_row(
                    f"MYCO_FILL_{index}",
                    "Mycoplasmatota",
                    role_tags="filler",
                )
            )
        for index in range(1, 10):
            rows.append(
                self.make_row(
                    f"BACI_FILL_{index}",
                    "Bacillota",
                    role_tags="filler",
                )
            )
        rows.append(
            self.make_row(
                "OTHER_FILL",
                "Pseudomonadota",
                role_tags="filler",
            )
        )
        return rows

    def test_builder_writes_valid_sample_and_metadata_outputs(self) -> None:
        """Write the expected sample-sheet and metadata headers and rows."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            candidate_path = root / "candidates.tsv"
            outdir = root / "generated"
            write_tsv(candidate_path, list(self.make_valid_rows()[0].keys()), self.make_valid_rows())

            result = self.run_builder(
                "--candidate-tsv",
                str(candidate_path),
                "--outdir",
                str(outdir),
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            sample_header, sample_rows = read_csv_rows(outdir / "sample_sheet.csv")
            metadata_header, metadata_rows = read_tsv_rows(outdir / "metadata.tsv")
            _selection_header, selection_rows = read_tsv_rows(outdir / "selection_report.tsv")
            _coverage_header, coverage_rows = read_tsv_rows(outdir / "coverage_report.tsv")

            self.assertEqual(sample_header, REQUIRED_SAMPLE_HEADER)
            self.assertEqual(metadata_header, REQUIRED_METADATA_HEADER)
            self.assertEqual(len(sample_rows), 20)
            self.assertTrue(sample_rows)
            self.assertTrue(selection_rows)
            self.assertTrue(coverage_rows)

    def test_builder_omits_metadata_rows_when_include_metadata_is_false(self) -> None:
        """Skip selected rows whose metadata should be omitted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            candidate_path = root / "candidates.tsv"
            outdir = root / "generated"
            write_tsv(candidate_path, list(self.make_valid_rows()[0].keys()), self.make_valid_rows())

            result = self.run_builder(
                "--candidate-tsv",
                str(candidate_path),
                "--outdir",
                str(outdir),
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            _metadata_header, metadata_rows = read_tsv_rows(outdir / "metadata.tsv")
            metadata_accessions = {row["Accession"] for row in metadata_rows}
            self.assertNotIn("MYCO_METALESS", metadata_accessions)

    def test_builder_rejects_invalid_is_new_and_include_metadata_combination(self) -> None:
        """Reject manifest rows that omit metadata for non-new samples."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            candidate_path = root / "candidates.tsv"
            rows = [
                self.make_row(
                    "ACC1",
                    "Mycoplasmatota",
                    role_tags="gcode4_candidate",
                    is_new="false",
                    include_metadata="false",
                ),
                self.make_row(
                    "ACC2",
                    "Bacillota",
                    role_tags="gcode11_candidate",
                ),
            ]
            write_tsv(candidate_path, list(rows[0].keys()), rows)

            result = self.run_builder("--candidate-tsv", str(candidate_path))

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("cannot omit metadata when is_new=false", result.stderr)

    def test_builder_requires_both_allowed_phyla(self) -> None:
        """Reject candidate manifests that do not contain both phyla."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            candidate_path = root / "candidates.tsv"
            rows = [
                self.make_row(f"MYCO_{index}", "Mycoplasmatota", role_tags="filler")
                for index in range(1, 21)
            ]
            write_tsv(candidate_path, list(rows[0].keys()), rows)

            result = self.run_builder("--candidate-tsv", str(candidate_path))

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("must include both allowed phyla", result.stderr)

    def test_builder_selection_is_deterministic(self) -> None:
        """Produce the same sample-sheet rows on repeated runs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            candidate_path = root / "candidates.tsv"
            rows = self.make_valid_rows()
            write_tsv(candidate_path, list(rows[0].keys()), rows)

            first_outdir = root / "first"
            second_outdir = root / "second"
            first = self.run_builder(
                "--candidate-tsv",
                str(candidate_path),
                "--outdir",
                str(first_outdir),
            )
            second = self.run_builder(
                "--candidate-tsv",
                str(candidate_path),
                "--outdir",
                str(second_outdir),
            )

            self.assertEqual(first.returncode, 0, msg=first.stderr)
            self.assertEqual(second.returncode, 0, msg=second.stderr)
            self.assertEqual(
                (first_outdir / "sample_sheet.csv").read_text(encoding="utf-8"),
                (second_outdir / "sample_sheet.csv").read_text(encoding="utf-8"),
            )

    def test_builder_selects_collision_pair_with_same_sanitised_base(self) -> None:
        """Select both collision candidates when a sanitised-base pair is available."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            candidate_path = root / "candidates.tsv"
            rows = self.make_valid_rows()
            write_tsv(candidate_path, list(rows[0].keys()), rows)

            result = self.run_builder(
                "--candidate-tsv",
                str(candidate_path),
                "--outdir",
                str(root / "generated"),
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            _header, sample_rows = read_csv_rows(root / "generated" / "sample_sheet.csv")
            accessions = {row["accession"] for row in sample_rows}
            self.assertIn("BACI-PAIR", accessions)
            self.assertIn("BACI PAIR", accessions)

    def test_builder_selects_at_least_two_ani_candidates_when_available(self) -> None:
        """Keep the medium cohort ANI coverage when the manifest provides it."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            candidate_path = root / "candidates.tsv"
            rows = self.make_valid_rows()
            write_tsv(candidate_path, list(rows[0].keys()), rows)

            result = self.run_builder(
                "--candidate-tsv",
                str(candidate_path),
                "--outdir",
                str(root / "generated"),
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            _header, coverage_rows = read_tsv_rows(root / "generated" / "coverage_report.tsv")
            ani_row = next(row for row in coverage_rows if row["coverage_key"] == "ani_cluster_candidate")
            self.assertEqual(ani_row["selected_count"], "2")
            self.assertEqual(ani_row["status"], "satisfied")


if __name__ == "__main__":
    unittest.main()
