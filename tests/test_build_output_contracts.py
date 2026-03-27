"""Tests for runtime output-contract generation."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import build_output_contracts  # noqa: E402


class BuildOutputContractsTestCase(unittest.TestCase):
    """Verify runtime BUSCO-aware contract generation."""

    def test_main_writes_dynamic_contracts_for_custom_lineages(self) -> None:
        """Write master-table and sample-status columns in configured lineage order."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            append_output = tmpdir / "master_table_append_columns.txt"
            status_output = tmpdir / "sample_status_columns.txt"

            exit_code = build_output_contracts.main(
                [
                    "--busco-lineage",
                    "custom_odb12",
                    "--busco-lineage",
                    "fallback_odb12",
                    "--append-columns-output",
                    str(append_output),
                    "--sample-status-columns-output",
                    str(status_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            append_columns = append_output.read_text(encoding="utf-8").splitlines()
            status_columns = status_output.read_text(encoding="utf-8").splitlines()
            self.assertIn("BUSCO_custom_odb12", append_columns)
            self.assertIn("BUSCO_fallback_odb12", append_columns)
            self.assertLess(
                append_columns.index("BUSCO_custom_odb12"),
                append_columns.index("BUSCO_fallback_odb12"),
            )
            self.assertIn("busco_custom_odb12_status", status_columns)
            self.assertIn("busco_fallback_odb12_status", status_columns)
            self.assertLess(
                status_columns.index("busco_custom_odb12_status"),
                status_columns.index("busco_fallback_odb12_status"),
            )

    def test_main_rejects_duplicate_lineages(self) -> None:
        """Fail when BUSCO lineage names repeat."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)

            exit_code = build_output_contracts.main(
                [
                    "--busco-lineage",
                    "custom_odb12",
                    "--busco-lineage",
                    "custom_odb12",
                    "--append-columns-output",
                    str(tmpdir / "append.txt"),
                    "--sample-status-columns-output",
                    str(tmpdir / "status.txt"),
                ]
            )

            self.assertEqual(exit_code, 1)


if __name__ == "__main__":
    unittest.main()
