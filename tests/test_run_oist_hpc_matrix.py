"""Tests for the OIST HPC matrix wrapper."""

from __future__ import annotations

import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "run_oist_hpc_matrix.sh"


class RunOistHpcMatrixScriptTestCase(unittest.TestCase):
    """Exercise the OIST HPC matrix wrapper."""

    def run_wrapper(self, *args: str) -> subprocess.CompletedProcess[str]:
        """Run the wrapper with text capture enabled."""
        return subprocess.run(
            ["bash", str(SCRIPT), *args],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_wrapper_help_lists_supported_modes(self) -> None:
        """Print wrapper help for the HPC matrix runner."""
        result = self.run_wrapper("--help")

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "<prepare|medium-prepare|db-full|db-reuse|db-matrix|p1|p2|all>",
            result.stdout,
        )
        self.assertIn("Run the OIST HPC validation campaign", result.stdout)
        self.assertIn("--hpc-root PATH", result.stdout)
        self.assertNotIn("--medium-candidates-tsv", result.stdout)

    def test_prepare_dry_run_prints_cohort_command(self) -> None:
        """Print the tracked-cohort preparation command."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "prepare",
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "bin/run_pipeline_test.sh prepare --work-root /tmp/hpc/acceptance",
            result.stdout,
        )

    def test_medium_prepare_dry_run_prints_medium_prepare_command(self) -> None:
        """Print the fixed medium-cohort preparation command."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "medium-prepare",
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "python3 bin/run_acceptance_tests.py prepare --work-root /tmp/hpc/medium",
            result.stdout,
        )
        self.assertIn(
            "--source-catalog /Users/asuq/Documents/Lab/Coding/nf-myco_annotation/assets/testdata/medium/source_catalog.tsv",
            result.stdout,
        )
        self.assertIn(
            "--cohort-plan /Users/asuq/Documents/Lab/Coding/nf-myco_annotation/assets/testdata/medium/cohort_plan.tsv",
            result.stdout,
        )

    def test_db_full_dry_run_uses_wrapper_and_validator(self) -> None:
        """Print the full DB-prep gate and its validation command."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "db-full",
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn("bin/run_pipeline_test.sh dbprep-slurm", result.stdout)
        self.assertIn("--dbprep-profile oist", result.stdout)
        self.assertIn("python3 bin/validate_hpc_matrix.py dbprep", result.stdout)
        self.assertNotIn("--padloc-db", result.stdout)

    def test_p1_dry_run_uses_profile_defaults(self) -> None:
        """Avoid explicit max-resource overrides in the real tracked run."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "p1",
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn("nextflow run . -profile oist", result.stdout)
        self.assertIn("--sample_csv /tmp/hpc/acceptance/generated/sample_sheet.csv", result.stdout)
        self.assertIn("--eggnog_db /tmp/hpc/db/Eggnog_db/Eggnog_Diamond_db", result.stdout)
        self.assertNotIn("--max_cpus", result.stdout)
        self.assertNotIn("--max_memory", result.stdout)
        self.assertNotIn("--max_time", result.stdout)
        self.assertNotIn("--padloc_db", result.stdout)

    def test_p2_auto_prepares_medium_inputs(self) -> None:
        """Generate the fixed medium cohort automatically when inputs are missing."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "p2",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn(
            "python3 bin/run_acceptance_tests.py prepare --work-root /tmp/hpc/medium",
            result.stdout,
        )
        self.assertIn(
            "--sample_csv /tmp/hpc/medium/generated/sample_sheet.csv",
            result.stdout,
        )
        self.assertIn(
            "--metadata /tmp/hpc/medium/generated/metadata.tsv",
            result.stdout,
        )

    def test_all_auto_prepares_medium_inputs(self) -> None:
        """Run the medium prepare step automatically during the full campaign."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "all",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn(
            "python3 bin/run_acceptance_tests.py prepare --work-root /tmp/hpc/medium",
            result.stdout,
        )

    def test_all_accepts_medium_inputs_after_mode(self) -> None:
        """Accept medium inputs even when they appear after the mode token."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "all",
            "--medium-sample-csv",
            "/tmp/medium_sample_sheet.csv",
            "--medium-metadata",
            "/tmp/medium_metadata.tsv",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn(
            "--sample_csv /tmp/medium_sample_sheet.csv",
            result.stdout,
        )
        self.assertIn(
            "--metadata /tmp/medium_metadata.tsv",
            result.stdout,
        )

    def test_wrapper_has_valid_bash_syntax(self) -> None:
        """Pass Bash syntax validation."""
        result = subprocess.run(
            ["bash", "-n", str(SCRIPT)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertEqual(result.stderr, "")


if __name__ == "__main__":
    unittest.main()
