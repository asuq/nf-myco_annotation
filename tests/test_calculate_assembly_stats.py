"""Tests for in-house assembly-stat calculation from staged FASTA files."""

from __future__ import annotations

import csv
import os
import stat
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "calculate_assembly_stats.sh"


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class CalculateAssemblyStatsTestCase(unittest.TestCase):
    """Cover the seqtk-based assembly-stats helper script."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write text to a temporary test file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def install_fake_seqtk(self, directory: Path) -> Path:
        """Install a fake seqtk executable that supports `seqtk comp`."""
        seqtk_path = directory / "seqtk"
        seqtk_path.write_text(
            textwrap.dedent(
                """\
                #!/usr/bin/env python3
                import sys
                from pathlib import Path

                def read_fasta(path: Path):
                    header = None
                    chunks = []
                    for raw_line in path.read_text(encoding="utf-8").splitlines():
                        line = raw_line.strip()
                        if not line:
                            continue
                        if line.startswith(">"):
                            if header is not None:
                                yield header, "".join(chunks)
                            header = line[1:]
                            chunks = []
                            continue
                        if header is None:
                            print("invalid FASTA", file=sys.stderr)
                            return
                        chunks.append(line)
                    if header is not None:
                        yield header, "".join(chunks)

                if len(sys.argv) != 3 or sys.argv[1] != "comp":
                    print("unsupported arguments", file=sys.stderr)
                    sys.exit(1)

                for name, sequence in read_fasta(Path(sys.argv[2])):
                    print(f"{name}\\t{len(sequence)}\\t0\\t0\\t0\\t0")
                """
            ),
            encoding="utf-8",
        )
        seqtk_path.chmod(seqtk_path.stat().st_mode | stat.S_IXUSR)
        return seqtk_path

    def run_helper(
        self,
        *,
        staged_manifest: Path,
        output: Path,
        path_prefix: Path,
    ) -> subprocess.CompletedProcess[str]:
        """Run the helper script with a controlled PATH."""
        environment = os.environ.copy()
        environment["PATH"] = f"{path_prefix}:{environment['PATH']}"
        return subprocess.run(
            ["bash", str(SCRIPT), "--staged-manifest", str(staged_manifest), "--output", str(output)],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )

    def test_calculates_single_contig_stats(self) -> None:
        """Emit expected stats for a single-contig FASTA."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            fasta = self.write_text_file(tmpdir / "ACC1.fasta", ">contig1\nAACCGGTT\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC1\tid_1\tACC1.fasta\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_tsv(output),
                [
                    {
                        "accession": "ACC1",
                        "n50": "8",
                        "scaffolds": "1",
                        "genome_size": "8",
                    }
                ],
            )

    def test_calculates_multi_contig_stats(self) -> None:
        """Emit the correct N50, scaffold count, and genome size for contigs."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(
                tmpdir / "ACC2.fasta",
                ">contig1\nAAAAAAAAAA\n>contig2\nAAAAAA\n>contig3\nAAA\n",
            )
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC2\tid_2\tACC2.fasta\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_tsv(output),
                [
                    {
                        "accession": "ACC2",
                        "n50": "10",
                        "scaffolds": "3",
                        "genome_size": "19",
                    }
                ],
            )

    def test_sorts_output_rows_by_accession(self) -> None:
        """Write assembly statistics in canonical accession order."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC_B.fasta", ">contig1\nAAAA\n")
            self.write_text_file(tmpdir / "ACC_A.fasta", ">contig1\nAAAAAA\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "\n".join(
                    [
                        "accession\tinternal_id\tstaged_filename",
                        "ACC_B\tid_b\tACC_B.fasta",
                        "ACC_A\tid_a\tACC_A.fasta",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                [row["accession"] for row in read_tsv(output)],
                ["ACC_A", "ACC_B"],
            )

    def test_fails_for_malformed_manifest_rows(self) -> None:
        """Fail with a direct row-shape error before duplicate checking."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC1.fasta", ">contig1\nAACCGGTT\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "\n".join(
                    [
                        "accession\tinternal_id\tstaged_filename",
                        "ACC1\tid_1\tACC1.fasta",
                        "ACC2",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Staged manifest row 3 is malformed", result.stderr)
            self.assertIn("expected at least 3 tab-delimited fields", result.stderr)

    def test_fails_for_duplicate_accessions(self) -> None:
        """Fail cleanly when the staged manifest repeats an accession."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC1.fasta", ">contig1\nAACCGGTT\n")
            self.write_text_file(tmpdir / "ACC1_copy.fasta", ">contig1\nTTGGCCAA\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "\n".join(
                    [
                        "accession\tinternal_id\tstaged_filename",
                        "ACC1\tid_1\tACC1.fasta",
                        "ACC1\tid_1_copy\tACC1_copy.fasta",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("duplicate accession values: ACC1", result.stderr)

    def test_fails_for_invalid_or_empty_fasta(self) -> None:
        """Fail cleanly when seqtk cannot produce contig lengths."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC3.fasta", "not-a-fasta\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC3\tid_3\tACC3.fasta\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Failed to calculate assembly statistics for ACC3", result.stderr)


if __name__ == "__main__":
    unittest.main()
