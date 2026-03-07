"""Tests for the final tool and database provenance collector."""

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

import collect_versions  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class CollectVersionsTestCase(unittest.TestCase):
    """Cover version-file parsing and runtime/context reporting."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write test text to a file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_collects_versions_resources_and_containers(self) -> None:
        """Write a stable TSV from process versions plus config context."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            versions_a = self.write_text_file(
                tmpdir / "validate_versions.yml",
                '\n'.join(
                    [
                        '"VALIDATE_INPUTS":',
                        '  python: "3.12.1"',
                        '  script: "bin/validate_inputs.py"',
                    ]
                )
                + "\n",
            )
            versions_b = self.write_text_file(
                tmpdir / "checkm2_versions.yml",
                '\n'.join(
                    [
                        '"CHECKM2":',
                        '  checkm2: "1.0.2"',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(versions_a),
                    "--version-file",
                    str(versions_b),
                    "--nextflow-version",
                    "25.04.8",
                    "--pipeline-version",
                    "0.1.0",
                    "--git-commit",
                    "abc1234",
                    "--container-engine",
                    "apptainer",
                    "--use-biocontainers",
                    "true",
                    "--checkm2-db",
                    "/db/checkm2",
                    "--taxdump",
                    "/db/taxdump",
                    "--busco-lineage",
                    "bacillota_odb12",
                    "--busco-lineage",
                    "mycoplasmatota_odb12",
                    "--busco-download-dir",
                    "/db/busco",
                    "--eggnog-db",
                    "/db/eggnog",
                    "--container-ref",
                    "python=python:3.12-slim",
                    "--container-ref",
                    "checkm2=quay.io/biocontainers/checkm2:1.0.2",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            rows = read_tsv(output)
            row_map = {(row["category"], row["name"], row["source"]): row["value"] for row in rows}
            self.assertEqual(row_map[("runtime", "nextflow_version", "workflow")], "25.04.8")
            self.assertEqual(row_map[("resource", "busco_lineages", "params")], "bacillota_odb12;mycoplasmatota_odb12")
            self.assertEqual(row_map[("container", "python", "params")], "python:3.12-slim")
            self.assertEqual(row_map[("tool_version", "python", "VALIDATE_INPUTS")], "3.12.1")
            self.assertEqual(row_map[("script", "script", "VALIDATE_INPUTS")], "bin/validate_inputs.py")
            self.assertEqual(row_map[("tool_version", "checkm2", "CHECKM2")], "1.0.2")

    def test_missing_version_file_fails(self) -> None:
        """Fail cleanly when an input versions file is missing."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(tmpdir / "missing.yml"),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
