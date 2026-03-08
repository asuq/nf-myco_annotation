"""Tests for the layered acceptance-test harness."""

from __future__ import annotations

import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import master_table_contract  # noqa: E402
import run_acceptance_tests  # noqa: E402


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read one TSV file into a header and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return reader.fieldnames, list(reader)


class RunAcceptanceTestsTestCase(unittest.TestCase):
    """Exercise cohort preparation and stable assertion helpers."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write text to a file and return the path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def write_gzip_file(self, path: Path, content: str) -> Path:
        """Write gzipped text and return the path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(content)
        return path

    def write_tsv_rows(
        self,
        path: Path,
        header: list[str],
        rows: list[dict[str, str]],
    ) -> Path:
        """Write TSV rows with a fixed header."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return path

    def write_source_catalog(self, path: Path, source_urls: dict[str, str]) -> Path:
        """Write one minimal source catalog for testing."""
        rows = [
            {
                "source_accession": "SRC_MYCO_A",
                "organism_name": "Myco A",
                "tax_id": "101",
                "assembly_level": "Complete Genome",
                "source_url": source_urls["SRC_MYCO_A"],
            },
            {
                "source_accession": "SRC_MYCO_B",
                "organism_name": "Myco B",
                "tax_id": "102",
                "assembly_level": "Complete Genome",
                "source_url": source_urls["SRC_MYCO_B"],
            },
            {
                "source_accession": "SRC_ECOLI",
                "organism_name": "E. coli",
                "tax_id": "103",
                "assembly_level": "Complete Genome",
                "source_url": source_urls["SRC_ECOLI"],
            },
            {
                "source_accession": "SRC_STREP",
                "organism_name": "Strep",
                "tax_id": "104",
                "assembly_level": "Complete Genome",
                "source_url": source_urls["SRC_STREP"],
            },
        ]
        return self.write_tsv_rows(
            path,
            ["source_accession", "organism_name", "tax_id", "assembly_level", "source_url"],
            rows,
        )

    def write_cohort_plan(self, path: Path) -> Path:
        """Write one role-complete cohort plan for testing."""
        rows = [
            {
                "accession": "SRC_MYCO_A",
                "source_accession": "SRC_MYCO_A",
                "is_new": "false",
                "include_metadata": "true",
                "atypical_warnings": "NA",
                "role_tags": "gcode4_candidate;crispr_negative_candidate",
            },
            {
                "accession": "SRC_MYCO_B",
                "source_accession": "SRC_MYCO_B",
                "is_new": "false",
                "include_metadata": "true",
                "atypical_warnings": "NA",
                "role_tags": "gcode4_candidate",
            },
            {
                "accession": "SRC_ECOLI",
                "source_accession": "SRC_ECOLI",
                "is_new": "false",
                "include_metadata": "true",
                "atypical_warnings": "NA",
                "role_tags": "gcode11_candidate",
            },
            {
                "accession": "SRC_STREP",
                "source_accession": "SRC_STREP",
                "is_new": "false",
                "include_metadata": "true",
                "atypical_warnings": "NA",
                "role_tags": "gcode11_candidate;crispr_positive_candidate",
            },
            {
                "accession": "SRC_MYCO_A_NEW",
                "source_accession": "SRC_MYCO_A",
                "is_new": "true",
                "include_metadata": "false",
                "atypical_warnings": "NA",
                "role_tags": "missing_metadata_case",
            },
            {
                "accession": "MYCO-ATYPICAL",
                "source_accession": "SRC_MYCO_A",
                "is_new": "false",
                "include_metadata": "true",
                "atypical_warnings": "cell culture adapted",
                "role_tags": "atypical_excluded_candidate;gcode4_candidate",
            },
            {
                "accession": "MYCO_EXCEPTION",
                "source_accession": "SRC_MYCO_A",
                "is_new": "false",
                "include_metadata": "true",
                "atypical_warnings": "unverified source organism",
                "role_tags": "atypical_exception_candidate;gcode4_candidate",
            },
            {
                "accession": "ECOLI-PAIR",
                "source_accession": "SRC_ECOLI",
                "is_new": "false",
                "include_metadata": "true",
                "atypical_warnings": "NA",
                "role_tags": "ani_cluster_candidate;collision_candidate;gcode11_candidate",
            },
            {
                "accession": "ECOLI PAIR",
                "source_accession": "SRC_ECOLI",
                "is_new": "false",
                "include_metadata": "true",
                "atypical_warnings": "NA",
                "role_tags": "ani_cluster_candidate;collision_candidate;gcode11_candidate",
            },
        ]
        return self.write_tsv_rows(
            path,
            [
                "accession",
                "source_accession",
                "is_new",
                "include_metadata",
                "atypical_warnings",
                "role_tags",
            ],
            rows,
        )

    def make_versions_table(
        self,
        path: Path,
        *,
        nextflow_version: str,
        checkm2_db_path: str,
    ) -> Path:
        """Write one versions table for comparison tests."""
        rows = [
            {
                "component": "nextflow",
                "kind": "runtime",
                "version": nextflow_version,
                "image_or_path": "NA",
                "notes": "workflow",
            },
            {
                "component": "pipeline",
                "kind": "pipeline",
                "version": "0.1.0",
                "image_or_path": "NA",
                "notes": "workflow manifest",
            },
            {
                "component": "checkm2_db",
                "kind": "database",
                "version": "v1",
                "image_or_path": checkm2_db_path,
                "notes": "CheckM2 database",
            },
            {
                "component": "python",
                "kind": "container",
                "version": "NA",
                "image_or_path": "python:3.12-slim",
                "notes": "params container reference",
            },
        ]
        return self.write_tsv_rows(
            path,
            ["component", "kind", "version", "image_or_path", "notes"],
            rows,
        )

    def test_prepare_cohort_builds_generated_inputs_and_checksums(self) -> None:
        """Prepare one hybrid cohort from local gzipped FASTA URLs."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source_urls = {
                "SRC_MYCO_A": self.write_gzip_file(
                    tmpdir / "downloads" / "myco_a.fna.gz",
                    ">a\nAAAA\n>b\nAAA\n",
                ).resolve().as_uri(),
                "SRC_MYCO_B": self.write_gzip_file(
                    tmpdir / "downloads" / "myco_b.fna.gz",
                    ">a\nAAAAAA\n",
                ).resolve().as_uri(),
                "SRC_ECOLI": self.write_gzip_file(
                    tmpdir / "downloads" / "ecoli.fna.gz",
                    ">a\nAAAAA\n>b\nAAAAA\n",
                ).resolve().as_uri(),
                "SRC_STREP": self.write_gzip_file(
                    tmpdir / "downloads" / "strep.fna.gz",
                    ">a\nAAAAAAAA\n",
                ).resolve().as_uri(),
            }
            source_catalog = self.write_source_catalog(tmpdir / "source_catalog.tsv", source_urls)
            cohort_plan = self.write_cohort_plan(tmpdir / "cohort_plan.tsv")

            prepared = run_acceptance_tests.prepare_cohort(
                work_root=tmpdir / "work",
                offline=False,
                source_catalog_path=source_catalog,
                cohort_plan_path=cohort_plan,
            )

            self.assertTrue(prepared.sample_csv.is_file())
            self.assertTrue(prepared.metadata_tsv.is_file())
            self.assertTrue(prepared.checksums_tsv.is_file())

            sample_header, sample_rows = run_acceptance_tests.read_csv(prepared.sample_csv)
            self.assertEqual(tuple(sample_header), run_acceptance_tests.SAMPLE_COLUMNS)
            self.assertEqual(len(sample_rows), 9)
            self.assertEqual(sample_rows[0]["genome_fasta"], str(prepared.source_stats["SRC_MYCO_A"].fasta_path))
            self.assertEqual(sample_rows[-1]["genome_fasta"], str(prepared.source_stats["SRC_ECOLI"].fasta_path))

            metadata_header, metadata_rows = read_tsv_rows(prepared.metadata_tsv)
            self.assertEqual(tuple(metadata_header), run_acceptance_tests.METADATA_COLUMNS)
            self.assertEqual(len(metadata_rows), 8)
            metadata_by_accession = {row["Accession"]: row for row in metadata_rows}
            self.assertNotIn("SRC_MYCO_A_NEW", metadata_by_accession)
            self.assertEqual(metadata_by_accession["MYCO-ATYPICAL"]["Atypical_Warnings"], "cell culture adapted")
            self.assertEqual(metadata_by_accession["MYCO_EXCEPTION"]["Atypical_Warnings"], "unverified source organism")

            self.assertEqual(prepared.source_stats["SRC_MYCO_A"].scaffolds, 2)
            self.assertEqual(prepared.source_stats["SRC_MYCO_A"].genome_size, 7)
            self.assertEqual(prepared.source_stats["SRC_MYCO_A"].n50, 4)

            checksum_header, checksum_rows = read_tsv_rows(prepared.checksums_tsv)
            self.assertEqual(tuple(checksum_header), run_acceptance_tests.DOWNLOAD_COLUMNS)
            self.assertEqual(len(checksum_rows), 4)

    def test_compare_versions_logically_ignores_runtime_only_drift(self) -> None:
        """Ignore runtime rows but detect stable provenance drift."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            local = self.make_versions_table(
                tmpdir / "local.tsv",
                nextflow_version="25.04.8",
                checkm2_db_path="/db/checkm2",
            )
            slurm = self.make_versions_table(
                tmpdir / "slurm.tsv",
                nextflow_version="25.04.9",
                checkm2_db_path="/db/checkm2",
            )

            run_acceptance_tests.compare_versions_logically(local, slurm)

            drifted = self.make_versions_table(
                tmpdir / "drifted.tsv",
                nextflow_version="25.04.9",
                checkm2_db_path="/db/other",
            )
            with self.assertRaisesRegex(
                run_acceptance_tests.AcceptanceTestError,
                "slurm_versions_logical_differs_from_local",
            ):
                run_acceptance_tests.compare_versions_logically(local, drifted)

    def test_assert_role_coverage_accepts_successful_positive_cohort(self) -> None:
        """Accept a synthetic positive cohort that proves each required role."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            cohort_plan_path = self.write_cohort_plan(tmpdir / "cohort_plan.tsv")
            source_catalog_path = self.write_source_catalog(
                tmpdir / "source_catalog.tsv",
                {
                    "SRC_MYCO_A": self.write_gzip_file(tmpdir / "a.fna.gz", ">a\nAAAA\n").resolve().as_uri(),
                    "SRC_MYCO_B": self.write_gzip_file(tmpdir / "b.fna.gz", ">a\nAAAA\n").resolve().as_uri(),
                    "SRC_ECOLI": self.write_gzip_file(tmpdir / "c.fna.gz", ">a\nAAAA\n").resolve().as_uri(),
                    "SRC_STREP": self.write_gzip_file(tmpdir / "d.fna.gz", ">a\nAAAA\n").resolve().as_uri(),
                },
            )
            catalog = run_acceptance_tests.load_source_catalog(source_catalog_path)
            plan = run_acceptance_tests.load_cohort_plan(cohort_plan_path, catalog)

            master_rows = {
                "SRC_MYCO_A": {"Gcode": "4", "CRISPRS": "0", "Cluster_ID": "NA", "Tax_ID": "101"},
                "SRC_MYCO_B": {"Gcode": "4", "CRISPRS": "0", "Cluster_ID": "NA", "Tax_ID": "102"},
                "SRC_ECOLI": {"Gcode": "11", "CRISPRS": "0", "Cluster_ID": "NA", "Tax_ID": "103"},
                "SRC_STREP": {"Gcode": "11", "CRISPRS": "2", "Cluster_ID": "NA", "Tax_ID": "104"},
                "SRC_MYCO_A_NEW": {"Gcode": "4", "CRISPRS": "0", "Cluster_ID": "NA", "Tax_ID": "NA"},
                "MYCO-ATYPICAL": {"Gcode": "4", "CRISPRS": "0", "Cluster_ID": "NA", "Tax_ID": "101"},
                "MYCO_EXCEPTION": {"Gcode": "4", "CRISPRS": "0", "Cluster_ID": "C000001", "Tax_ID": "101"},
                "ECOLI-PAIR": {"Gcode": "11", "CRISPRS": "0", "Cluster_ID": "C000002", "Tax_ID": "103"},
                "ECOLI PAIR": {"Gcode": "11", "CRISPRS": "0", "Cluster_ID": "C000002", "Tax_ID": "103"},
            }
            status_rows = {
                "SRC_MYCO_A": {"accession": "SRC_MYCO_A", "ani_included": "false", "ani_exclusion_reason": "low_quality", "warnings": "", "internal_id": "SRC_MYCO_A"},
                "SRC_MYCO_B": {"accession": "SRC_MYCO_B", "ani_included": "false", "ani_exclusion_reason": "low_quality", "warnings": "", "internal_id": "SRC_MYCO_B"},
                "SRC_ECOLI": {"accession": "SRC_ECOLI", "ani_included": "true", "ani_exclusion_reason": "", "warnings": "", "internal_id": "SRC_ECOLI"},
                "SRC_STREP": {"accession": "SRC_STREP", "ani_included": "true", "ani_exclusion_reason": "", "warnings": "", "internal_id": "SRC_STREP"},
                "SRC_MYCO_A_NEW": {"accession": "SRC_MYCO_A_NEW", "ani_included": "false", "ani_exclusion_reason": "partial_16s", "warnings": "missing_metadata_for_new_sample", "internal_id": "SRC_MYCO_A_NEW"},
                "MYCO-ATYPICAL": {"accession": "MYCO-ATYPICAL", "ani_included": "false", "ani_exclusion_reason": "atypical", "warnings": "", "internal_id": "MYCO_ATYPICAL"},
                "MYCO_EXCEPTION": {"accession": "MYCO_EXCEPTION", "ani_included": "true", "ani_exclusion_reason": "", "warnings": "", "internal_id": "MYCO_EXCEPTION"},
                "ECOLI-PAIR": {"accession": "ECOLI-PAIR", "ani_included": "true", "ani_exclusion_reason": "", "warnings": "internal_id_collision_resolved", "internal_id": "ECOLI_PAIR_aaa"},
                "ECOLI PAIR": {"accession": "ECOLI PAIR", "ani_included": "true", "ani_exclusion_reason": "", "warnings": "internal_id_collision_resolved", "internal_id": "ECOLI_PAIR_bbb"},
            }

            run_acceptance_tests.assert_role_coverage(
                plan=plan,
                master_rows=master_rows,
                status_rows=status_rows,
            )

    def test_assert_metadata_contract_uses_locked_append_columns(self) -> None:
        """Require the metadata prefix and append-column contract exactly."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            metadata_path = self.write_tsv_rows(
                tmpdir / "metadata.tsv",
                list(run_acceptance_tests.METADATA_COLUMNS),
                [
                    {
                        "Accession": "ACC1",
                        "Tax_ID": "101",
                        "Organism_Name": "Genome",
                        "Assembly_Level": "Complete Genome",
                        "N50": "100",
                        "Scaffolds": "1",
                        "Genome_Size": "100",
                        "Atypical_Warnings": "NA",
                    }
                ],
            )
            master_header = list(run_acceptance_tests.METADATA_COLUMNS) + master_table_contract.read_append_columns_asset()
            master_path = self.write_tsv_rows(
                tmpdir / "master_table.tsv",
                master_header,
                [
                    {column: "NA" for column in master_header} | {"Accession": "ACC1"}
                ],
            )

            run_acceptance_tests.assert_metadata_contract(master_path, metadata_path)


if __name__ == "__main__":
    unittest.main()
