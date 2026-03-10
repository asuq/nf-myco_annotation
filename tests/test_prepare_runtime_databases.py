"""Tests for the runtime database preparation helper."""

from __future__ import annotations

import csv
import io
import json
import logging
import os
import sys
import tarfile
import tempfile
import unittest
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import prepare_runtime_databases  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read one TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PrepareRuntimeDatabasesTestCase(unittest.TestCase):
    """Cover in-place runtime database preparation and reuse."""

    def reset_logging(self) -> None:
        """Reset root logging handlers so stderr capture stays predictable."""
        root_logger = logging.getLogger()
        for handler in list(root_logger.handlers):
            root_logger.removeHandler(handler)

    def run_cli(self, args: list[str]) -> tuple[int, str, str]:
        """Run the helper CLI and capture stdout plus stderr."""
        self.reset_logging()
        stdout = io.StringIO()
        stderr = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            exit_code = prepare_runtime_databases.main(args)
        return exit_code, stdout.getvalue(), stderr.getvalue()

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write one UTF-8 test file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def read_marker(self, path: Path) -> dict[str, object]:
        """Read one ready marker JSON object."""
        return json.loads(
            (path / prepare_runtime_databases.MARKER_FILE_NAME).read_text(
                encoding="utf-8"
            )
        )

    def create_taxdump_dir(self, path: Path) -> Path:
        """Create one valid taxdump source directory."""
        self.write_text_file(path / "names.dmp", "1\t|\troot\t|\t\t|\tscientific name\t|\n")
        self.write_text_file(path / "nodes.dmp", "1\t|\t1\t|\tno rank\t|\n")
        return path

    def create_checkm2_dir(self, path: Path) -> Path:
        """Create one valid CheckM2 source directory."""
        self.write_text_file(path / "CheckM2_database.dmnd", "fake-dmnd\n")
        return path

    def create_busco_lineage_dir(self, path: Path) -> Path:
        """Create one valid BUSCO lineage source directory."""
        self.write_text_file(path / "dataset.cfg", "name=fake\n")
        return path

    def create_eggnog_dir(self, path: Path) -> Path:
        """Create one valid eggNOG source directory."""
        self.write_text_file(path / "eggnog.db", "sqlite-placeholder\n")
        self.write_text_file(path / "eggnog_proteins.dmnd", "diamond-placeholder\n")
        return path

    def create_padloc_dir(self, path: Path) -> Path:
        """Create one valid PADLOC source directory."""
        self.write_text_file(path / "hmm" / "padlocdb.hmm", "hmm-placeholder\n")
        return path

    def create_tar_archive(self, source: Path, archive_path: Path, root_name: str) -> Path:
        """Create one tar archive containing the source directory under one root."""
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive_path, "w:gz") as handle:
            handle.add(source, arcname=root_name)
        return archive_path

    def create_zip_archive(self, source: Path, archive_path: Path, root_name: str) -> Path:
        """Create one zip archive containing the source directory under one root."""
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive_path, "w") as handle:
            for child in source.rglob("*"):
                if child.is_dir():
                    continue
                relative = child.relative_to(source)
                handle.write(child, Path(root_name) / relative)
        return archive_path

    def test_main_prepares_all_directory_sources_and_writes_report(self) -> None:
        """Prepare all runtime database classes from directory sources."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            taxdump_source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            checkm2_source = self.create_checkm2_dir(tmpdir / "sources" / "checkm2")
            bacillota_source = self.create_busco_lineage_dir(
                tmpdir / "sources" / "busco" / "bacillota_odb12"
            )
            myco_source = self.create_busco_lineage_dir(
                tmpdir / "sources" / "busco" / "mycoplasmatota_odb12"
            )
            eggnog_source = self.create_eggnog_dir(tmpdir / "sources" / "eggnog")
            padloc_source = self.create_padloc_dir(tmpdir / "sources" / "padloc")
            report_path = tmpdir / "report.tsv"

            exit_code, stdout, _stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(taxdump_source),
                    "--taxdump-dest",
                    str(tmpdir / "prepared" / "taxdump"),
                    "--checkm2-source",
                    str(checkm2_source),
                    "--checkm2-dest",
                    str(tmpdir / "prepared" / "checkm2"),
                    "--busco-lineage-source",
                    f"bacillota_odb12={bacillota_source}",
                    "--busco-lineage-source",
                    f"mycoplasmatota_odb12={myco_source}",
                    "--busco-dest-root",
                    str(tmpdir / "prepared" / "busco"),
                    "--eggnog-source",
                    str(eggnog_source),
                    "--eggnog-dest",
                    str(tmpdir / "prepared" / "eggnog"),
                    "--padloc-source",
                    str(padloc_source),
                    "--padloc-dest",
                    str(tmpdir / "prepared" / "padloc"),
                    "--report",
                    str(report_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertIn("--taxdump", stdout)
            self.assertIn("--checkm2_db", stdout)
            self.assertIn("--busco_download_dir", stdout)
            self.assertIn("--eggnog_db", stdout)
            self.assertIn("--padloc_db", stdout)

            rows = read_tsv(report_path)
            self.assertEqual(len(rows), 7)
            row_map = {row["component"]: row for row in rows}
            self.assertEqual(row_map["taxdump"]["status"], "prepared")
            self.assertEqual(row_map["checkm2"]["status"], "prepared")
            self.assertEqual(row_map["busco:bacillota_odb12"]["status"], "prepared")
            self.assertEqual(row_map["busco:mycoplasmatota_odb12"]["status"], "prepared")
            self.assertEqual(
                row_map["busco_root"]["destination"],
                str((tmpdir / "prepared" / "busco").resolve()),
            )

            taxdump_dest = tmpdir / "prepared" / "taxdump"
            checkm2_dest = tmpdir / "prepared" / "checkm2"
            bacillota_dest = tmpdir / "prepared" / "busco" / "bacillota_odb12"
            myco_dest = tmpdir / "prepared" / "busco" / "mycoplasmatota_odb12"
            eggnog_dest = tmpdir / "prepared" / "eggnog"
            padloc_dest = tmpdir / "prepared" / "padloc"

            self.assertTrue((taxdump_dest / "names.dmp").is_file())
            self.assertTrue((checkm2_dest / "CheckM2_database.dmnd").is_file())
            self.assertTrue((bacillota_dest / "dataset.cfg").is_file())
            self.assertTrue((myco_dest / "dataset.cfg").is_file())
            self.assertTrue((eggnog_dest / "eggnog.db").is_file())
            self.assertTrue((padloc_dest / "hmm" / "padlocdb.hmm").is_file())

            self.assertEqual(self.read_marker(taxdump_dest)["component"], "taxdump")
            self.assertEqual(self.read_marker(checkm2_dest)["component"], "checkm2")
            self.assertEqual(
                self.read_marker(bacillota_dest)["component"],
                "busco:bacillota_odb12",
            )
            self.assertEqual(self.read_marker(eggnog_dest)["component"], "eggnog")
            self.assertEqual(self.read_marker(padloc_dest)["component"], "padloc")

    def test_main_accepts_checkm2_raw_dmnd_source(self) -> None:
        """Materialise one raw CheckM2 dmnd file into a prepared directory."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            dmnd_source = self.write_text_file(tmpdir / "sources" / "checkm2.dmnd", "fake\n")
            destination = tmpdir / "prepared" / "checkm2"

            exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--checkm2-source",
                    str(dmnd_source),
                    "--checkm2-dest",
                    str(destination),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue((destination / "checkm2.dmnd").is_file())
            self.assertEqual(self.read_marker(destination)["component"], "checkm2")

    def test_main_skips_valid_existing_destination_with_marker(self) -> None:
        """Reuse one ready destination rather than rebuilding it."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            destination = tmpdir / "prepared" / "taxdump"

            first_exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(source),
                    "--taxdump-dest",
                    str(destination),
                ]
            )
            second_exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(source),
                    "--taxdump-dest",
                    str(destination),
                    "--report",
                    str(tmpdir / "report.tsv"),
                ]
            )

            self.assertEqual(first_exit_code, 0)
            self.assertEqual(second_exit_code, 0)
            rows = read_tsv(tmpdir / "report.tsv")
            self.assertEqual(rows[0]["status"], "present")

    def test_main_reports_incomplete_destination_without_marker(self) -> None:
        """Fail when one destination is valid on disk but lacks the ready marker."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            destination = self.create_taxdump_dir(tmpdir / "prepared" / "taxdump")

            exit_code, _stdout, stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(source),
                    "--taxdump-dest",
                    str(destination),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertIn(prepare_runtime_databases.MARKER_FILE_NAME, stderr)

    def test_main_fails_on_invalid_existing_destination_without_force(self) -> None:
        """Fail cleanly when one destination exists but does not validate."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_padloc_dir(tmpdir / "sources" / "padloc")
            self.write_text_file(tmpdir / "prepared" / "padloc" / "broken.txt", "broken\n")

            exit_code, _stdout, stderr = self.run_cli(
                [
                    "--padloc-source",
                    str(source),
                    "--padloc-dest",
                    str(tmpdir / "prepared" / "padloc"),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertIn("invalid or partially prepared", stderr)

    def test_main_force_rebuilds_invalid_destination(self) -> None:
        """Replace one invalid destination when force is requested."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_padloc_dir(tmpdir / "sources" / "padloc")
            self.write_text_file(tmpdir / "prepared" / "padloc" / "broken.txt", "broken\n")

            exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--padloc-source",
                    str(source),
                    "--padloc-dest",
                    str(tmpdir / "prepared" / "padloc"),
                    "--force",
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue(
                (tmpdir / "prepared" / "padloc" / "hmm" / "padlocdb.hmm").is_file()
            )
            self.assertFalse((tmpdir / "prepared" / "padloc" / "broken.txt").exists())

    def test_main_honours_preparation_lock(self) -> None:
        """Fail fast when another preparation lock is already present."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            destination = tmpdir / "prepared" / "taxdump"
            lock_path = prepare_runtime_databases.lock_path(destination)
            self.write_text_file(lock_path, "{}\n")

            exit_code, _stdout, stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(source),
                    "--taxdump-dest",
                    str(destination),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertIn("Preparation lock already exists", stderr)

    def test_main_prepares_archive_sources_without_using_system_tmp(self) -> None:
        """Extract archive sources beside the destination parent rather than via system tmp."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)

            taxdump_archive = self.create_zip_archive(
                self.create_taxdump_dir(tmpdir / "sources" / "taxdump_dir"),
                tmpdir / "archives" / "taxdump.zip",
                "taxdump_payload",
            )
            busco_archive = self.create_tar_archive(
                self.create_busco_lineage_dir(tmpdir / "sources" / "busco_dir"),
                tmpdir / "archives" / "busco.tar.gz",
                "busco_payload",
            )
            eggnog_archive = self.create_tar_archive(
                self.create_eggnog_dir(tmpdir / "sources" / "eggnog_dir"),
                tmpdir / "archives" / "eggnog.tar.gz",
                "eggnog_payload",
            )
            padloc_archive = self.create_zip_archive(
                self.create_padloc_dir(tmpdir / "sources" / "padloc_dir"),
                tmpdir / "archives" / "padloc.zip",
                "padloc_payload",
            )

            with mock.patch.dict(os.environ, {"TMPDIR": str(tmpdir / "blocked-tmp")}):
                exit_code, _stdout, _stderr = self.run_cli(
                    [
                        "--taxdump-source",
                        str(taxdump_archive),
                        "--taxdump-dest",
                        str(tmpdir / "prepared" / "taxdump"),
                        "--busco-lineage-source",
                        f"bacillota_odb12={busco_archive}",
                        "--busco-dest-root",
                        str(tmpdir / "prepared" / "busco"),
                        "--eggnog-source",
                        str(eggnog_archive),
                        "--eggnog-dest",
                        str(tmpdir / "prepared" / "eggnog"),
                        "--padloc-source",
                        str(padloc_archive),
                        "--padloc-dest",
                        str(tmpdir / "prepared" / "padloc"),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertEqual(list((tmpdir / "prepared").glob(".prepare-*")), [])
            self.assertTrue((tmpdir / "prepared" / "taxdump" / "names.dmp").is_file())
            self.assertTrue(
                (tmpdir / "prepared" / "busco" / "bacillota_odb12" / "dataset.cfg").is_file()
            )
            self.assertTrue((tmpdir / "prepared" / "eggnog" / "eggnog.db").is_file())
            self.assertTrue(
                (tmpdir / "prepared" / "padloc" / "hmm" / "padlocdb.hmm").is_file()
            )


if __name__ == "__main__":
    unittest.main()
