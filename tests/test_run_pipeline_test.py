"""Tests for the manual pipeline test wrapper."""

from __future__ import annotations

import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "run_pipeline_test.sh"


class RunPipelineTestScriptTestCase(unittest.TestCase):
    """Exercise the Bash wrapper around the acceptance harness."""

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
        """Print wrapper help before a mode is chosen."""
        result = self.run_wrapper("--help")

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "bin/run_pipeline_test.sh [--dry-run] <prepare|unit|stub|local|slurm|all>",
            result.stdout,
        )
        self.assertIn("Manual wrapper around bin/run_acceptance_tests.py.", result.stdout)

    def test_wrapper_dry_run_stub_prints_delegated_command(self) -> None:
        """Print the delegated stub command without executing it."""
        result = self.run_wrapper("--dry-run", "stub")

        self.assertEqual(result.returncode, 0)
        self.assertEqual(
            result.stdout.strip(),
            "python3 bin/run_acceptance_tests.py stub",
        )
        self.assertEqual(result.stderr, "")

    def test_wrapper_dry_run_local_preserves_forwarded_arguments(self) -> None:
        """Forward real-run arguments unchanged during dry-run output."""
        result = self.run_wrapper(
            "--dry-run",
            "local",
            "--taxdump",
            "/tmp/taxdump",
            "--checkm2-db",
            "/tmp/checkm2",
            "--busco-download-dir",
            "/tmp/busco",
            "--eggnog-db",
            "/tmp/eggnog",
            "--padloc-db",
            "/tmp/padloc",
        )

        self.assertEqual(result.returncode, 0)
        self.assertEqual(
            result.stdout.strip(),
            "python3 bin/run_acceptance_tests.py local --taxdump /tmp/taxdump "
            "--checkm2-db /tmp/checkm2 --busco-download-dir /tmp/busco "
            "--eggnog-db /tmp/eggnog --padloc-db /tmp/padloc",
        )
        self.assertEqual(result.stderr, "")

    def test_wrapper_rejects_invalid_mode(self) -> None:
        """Fail with a clear mode-validation error."""
        result = self.run_wrapper("bogus")

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Error: unsupported mode bogus.", result.stderr)
        self.assertIn("Expected one of: prepare unit stub local slurm all", result.stderr)

    def test_wrapper_forwards_mode_help_to_acceptance_harness(self) -> None:
        """Delegate per-mode help text to the Python harness."""
        result = self.run_wrapper("local", "--help")

        self.assertEqual(result.returncode, 0)
        self.assertIn("Run the real-data local acceptance cohort.", result.stdout)
        self.assertIn("--checkm2-db CHECKM2_DB", result.stdout)

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
