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
                    upper = sequence.upper()
                    print(
                        f"{name}\\t{len(sequence)}\\t"
                        f"{upper.count('A')}\\t{upper.count('C')}\\t"
                        f"{upper.count('G')}\\t{upper.count('T')}"
                    )
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
        jobs: int = 1,
        environment_overrides: dict[str, str] | None = None,
        cwd: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        """Run the helper script with a controlled PATH."""
        environment = os.environ.copy()
        environment["PATH"] = f"{path_prefix}:{environment['PATH']}"
        if environment_overrides:
            environment.update(environment_overrides)
        return subprocess.run(
            [
                "bash",
                str(SCRIPT),
                "--staged-manifest",
                str(staged_manifest),
                "--jobs",
                str(jobs),
                "--output",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
            cwd=cwd,
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
                        "gc_content": "50",
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
                        "gc_content": "0",
                    }
                ],
            )

    def test_accepts_crlf_manifest_rows(self) -> None:
        """Handle CRLF line endings in the staged manifest without losing file paths."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC4.fasta", ">contig1\nAACCGGTT\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\r\nACC4\tid_4\tACC4.fasta\r\n",
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
                        "accession": "ACC4",
                        "n50": "8",
                        "scaffolds": "1",
                        "genome_size": "8",
                        "gc_content": "50",
                    }
                ],
            )

    def test_excludes_ambiguous_bases_from_gc_denominator(self) -> None:
        """Calculate GC from canonical bases only."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC5.fasta", ">contig1\nACGTNNNN\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC5\tid_5\tACC5.fasta\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(read_tsv(output)[0]["gc_content"], "50")

    def test_uses_output_directory_when_tmpdir_is_unwritable(self) -> None:
        """Avoid implicit mktemp locations that can be read-only in containers."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC7.fasta", ">contig1\nAACCGGTT\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC7\tid_7\tACC7.fasta\n",
            )
            output_dir = tmpdir / "output"
            output_dir.mkdir()
            output = output_dir / "assembly_stats.tsv"
            blocked_tmp = tmpdir / "blocked_tmp"
            blocked_tmp.mkdir()
            blocked_tmp.chmod(stat.S_IREAD | stat.S_IEXEC)

            try:
                result = self.run_helper(
                    staged_manifest=manifest,
                    output=output,
                    path_prefix=tmpdir,
                    environment_overrides={"TMPDIR": str(blocked_tmp)},
                    cwd=blocked_tmp,
                )
            finally:
                blocked_tmp.chmod(stat.S_IREAD | stat.S_IWRITE | stat.S_IEXEC)

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(read_tsv(output)[0]["accession"], "ACC7")

    def test_parallel_jobs_preserve_manifest_order(self) -> None:
        """Emit deterministic row order while processing multiple genomes in parallel."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC2.fasta", ">contig1\nAAAAAAAAAA\n")
            self.write_text_file(tmpdir / "ACC1.fasta", ">contig1\nAAAA\n")
            self.write_text_file(tmpdir / "ACC3.fasta", ">contig1\nAAAAAAA\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                (
                    "accession\tinternal_id\tstaged_filename\n"
                    "ACC2\tid_2\tACC2.fasta\n"
                    "ACC1\tid_1\tACC1.fasta\n"
                    "ACC3\tid_3\tACC3.fasta\n"
                ),
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
                jobs=2,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                [row["accession"] for row in read_tsv(output)],
                ["ACC2", "ACC1", "ACC3"],
            )

    def test_parallel_jobs_report_failing_accession(self) -> None:
        """Fail the stage when one parallel worker cannot derive contig lengths."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC1.fasta", ">contig1\nAACCGGTT\n")
            self.write_text_file(tmpdir / "ACC2.fasta", "not-a-fasta\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                (
                    "accession\tinternal_id\tstaged_filename\n"
                    "ACC1\tid_1\tACC1.fasta\n"
                    "ACC2\tid_2\tACC2.fasta\n"
                ),
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
                jobs=2,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Failed to calculate assembly statistics for ACC2", result.stderr)

    def test_parallel_jobs_accept_paths_with_spaces(self) -> None:
        """Pass FASTA paths with spaces safely through the worker queue."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            fasta_dir = tmpdir / "with spaces"
            fasta_dir.mkdir(parents=True, exist_ok=True)
            self.write_text_file(fasta_dir / "ACC 5.fasta", ">contig1\nAACCGGTTAA\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC5\tid_5\twith spaces/ACC 5.fasta\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
                jobs=2,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(read_tsv(output)[0]["accession"], "ACC5")

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

    def test_fails_when_fasta_contains_no_canonical_bases(self) -> None:
        """Reject FASTA input that contains no A, C, G, or T bases."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            self.write_text_file(tmpdir / "ACC6.fasta", ">contig1\nNNNNNN\n")
            manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC6\tid_6\tACC6.fasta\n",
            )
            output = tmpdir / "assembly_stats.tsv"

            result = self.run_helper(
                staged_manifest=manifest,
                output=output,
                path_prefix=tmpdir,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Failed to calculate assembly statistics for ACC6", result.stderr)


if __name__ == "__main__":
    unittest.main()
