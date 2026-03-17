"""Tests for the HPC matrix validation helper."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "validate_hpc_matrix.py"
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location("validate_hpc_matrix", SCRIPT)
assert SPEC is not None
assert SPEC.loader is not None
validate_hpc_matrix = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validate_hpc_matrix)


def write_text(path: Path, text: str) -> None:
    """Write one text file, creating parent directories first."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


class ValidateHpcMatrixTestCase(unittest.TestCase):
    """Exercise the HPC matrix validation helper."""

    def run_helper(self, *args: str) -> subprocess.CompletedProcess[str]:
        """Run the helper with text capture enabled."""
        return subprocess.run(
            ["python3", str(SCRIPT), *args],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_dbprep_validation_accepts_minimal_prepared_tree(self) -> None:
        """Validate a minimal successful dbprep case."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results_dir = root / "results"
            write_text(
                results_dir / "runtime_database_report.tsv",
                (
                    "component\tstatus\tsource\tdestination\tdetails\n"
                    f"taxdump\tprepared\tnative_download\t{root / 'db' / 'taxdump'}\tfiles=2\n"
                ),
            )
            write_text(
                results_dir / "nextflow_args.txt",
                f"--taxdump {root / 'db' / 'taxdump'}\n",
            )
            for name in ("trace.tsv", "report.html", "timeline.html", "dag.html"):
                write_text(results_dir / "pipeline_info" / name, "ok\n")

            taxdump_dir = root / "db" / "taxdump"
            write_text(taxdump_dir / "names.dmp", "names\n")
            write_text(taxdump_dir / "nodes.dmp", "nodes\n")
            write_text(taxdump_dir / ".nf_myco_ready.json", "{}\n")

            result = self.run_helper(
                "dbprep",
                "--results-dir",
                str(results_dir),
                "--taxdump",
                str(taxdump_dir),
                "--expected-component",
                "taxdump",
                "--expected-status",
                "prepared",
                "--expected-arg=--taxdump",
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)

    def test_medium_run_validation_accepts_allowed_phyla(self) -> None:
        """Validate a minimal successful medium real-data run."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            outdir = root / "out"
            write_text(
                outdir / "tables" / "master_table.tsv",
                (
                    "Accession\tphylum\n"
                    "ACC1\tMycoplasmatota\n"
                    "ACC2\tBacillota\n"
                ),
            )
            write_text(
                outdir / "tables" / "sample_status.tsv",
                "accession\tvalidation_status\nACC1\tdone\nACC2\tdone\n",
            )
            write_text(
                outdir / "tables" / "tool_and_db_versions.tsv",
                (
                    "component\tkind\tversion\timage_or_path\tnotes\n"
                    f"eggnog_mapper\ttool\t2.1.13\tNA\tNA\n"
                    f"taxdump\tdatabase\t20240914\t{(root / 'db' / 'taxdump').resolve()}\tNA\n"
                ),
            )
            for name in ("trace.tsv", "report.html", "timeline.html", "dag.html"):
                write_text(outdir / "pipeline_info" / name, "ok\n")

            result = self.run_helper(
                "medium-run",
                "--outdir",
                str(outdir),
                "--db-root",
                str(root / "db"),
                "--sample-count",
                "2",
                "--allowed-phylum",
                "Mycoplasmatota",
                "--allowed-phylum",
                "Bacillota",
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)

    def test_medium_phylum_helper_allows_missing_metadata_case_na(self) -> None:
        """Allow `phylum=NA` only for the deliberate missing-metadata case."""
        master_rows = {
            "MYCO_NEW_1": {"phylum": "NA"},
            "ACC1": {"phylum": "Mycoplasmatota"},
            "ACC2": {"phylum": "Bacillota"},
        }

        validate_hpc_matrix.assert_medium_phyla_allowed(
            master_rows,
            allowed_phyla={"Mycoplasmatota", "Bacillota"},
            missing_metadata_case_accessions={"MYCO_NEW_1"},
        )

    def test_medium_phylum_helper_rejects_unexpected_na(self) -> None:
        """Reject `phylum=NA` for accessions outside the missing-metadata role."""
        master_rows = {
            "ACC_BAD": {"phylum": "NA"},
            "ACC1": {"phylum": "Mycoplasmatota"},
            "ACC2": {"phylum": "Bacillota"},
        }

        with self.assertRaises(validate_hpc_matrix.ValidateHpcMatrixError) as context:
            validate_hpc_matrix.assert_medium_phyla_allowed(
                master_rows,
                allowed_phyla={"Mycoplasmatota", "Bacillota"},
                missing_metadata_case_accessions=set(),
            )

        self.assertIn("Observed unexpected phylum=NA rows in medium run: ACC_BAD", str(context.exception))


if __name__ == "__main__":
    unittest.main()
