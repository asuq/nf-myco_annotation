"""Tests for rescuing ANI outputs from published pipeline results."""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Sequence
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import rescue_ani_from_results  # noqa: E402

SAMPLE_STATUS_COLUMNS = [
    line.strip()
    for line in (
        ROOT / "assets" / "tables" / "contracts" / "sample_status_columns.txt"
    ).read_text(encoding="utf-8").splitlines()
    if line.strip()
]
SCIPY_AVAILABLE = importlib.util.find_spec("scipy") is not None


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return reader.fieldnames, list(reader)


class RescueAniFromResultsTestCase(unittest.TestCase):
    """Cover the end-to-end ANI rescue CLI."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write one UTF-8 text file and return the path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def write_tsv_rows(
        self,
        path: Path,
        header: Sequence[str],
        rows: Sequence[dict[str, str]],
    ) -> Path:
        """Write one TSV file with a fixed header."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return path

    def make_initial_status_row(
        self,
        *,
        accession: str,
        internal_id: str,
        is_new: str = "false",
    ) -> dict[str, str]:
        """Build one seed status row matching the locked sample-status contract."""
        row: dict[str, str] = {}
        for column in SAMPLE_STATUS_COLUMNS:
            if column.endswith("_status") or column == "ani_included":
                row[column] = "na"
            elif column in {"gcode", "low_quality"}:
                row[column] = "NA"
            else:
                row[column] = ""
        row["accession"] = accession
        row["internal_id"] = internal_id
        row["is_new"] = is_new
        row["validation_status"] = "done"
        return row

    def barrnap_header(
        self,
        rrna_name: str,
        contig: str,
        start0: int,
        end0: int,
        strand: str,
    ) -> str:
        """Build one Barrnap FASTA header."""
        return f"{rrna_name}::{contig}:{start0}-{end0}({strand})"

    def write_barrnap_outputs(
        self,
        source_outdir: Path,
        *,
        accession: str,
        mode: str,
    ) -> None:
        """Write one sample's published Barrnap outputs."""
        barrnap_dir = source_outdir / "samples" / accession / "barrnap"
        if mode == "yes":
            header = self.barrnap_header("16S_rRNA", "contig1", 0, 1500, "+")
            gff_text = (
                "contig1\tbarrnap\trRNA\t1\t1500\t5.0\t+\t.\t"
                "Name=16S_rRNA;product=16S ribosomal RNA\n"
            )
            fasta_text = f">{header}\n" + ("A" * 1500) + "\n"
        elif mode == "partial":
            header = self.barrnap_header("16S_rRNA", "contig1", 100, 1200, "+")
            gff_text = (
                "contig1\tbarrnap\trRNA\t101\t1200\t5.0\t+\t.\t"
                "Name=16S_rRNA;product=16S ribosomal RNA (partial)\n"
            )
            fasta_text = f">{header}\n" + ("C" * 1100) + "\n"
        else:
            raise AssertionError(f"Unsupported Barrnap mode: {mode}")

        self.write_text_file(barrnap_dir / "rrna.gff", gff_text)
        self.write_text_file(barrnap_dir / "rrna.fa", fasta_text)

    def write_checkm2_report(
        self,
        source_outdir: Path,
        *,
        accession: str,
        translation_table: int,
        completeness: str,
        contamination: str,
    ) -> None:
        """Write one published CheckM2 quality report."""
        report_dir = (
            source_outdir / "samples" / accession / f"checkm2_gcode{translation_table}"
        )
        self.write_text_file(
            report_dir / "quality_report.tsv",
            "\n".join(
                [
                    "Name\tCompleteness\tContamination\tCoding_Density\tAverage_Gene_Length\tTotal_Coding_Sequences",
                    f"{accession}\t{completeness}\t{contamination}\t0.9\t900\t800",
                ]
            )
            + "\n",
        )

    def write_busco_summary(
        self,
        source_outdir: Path,
        *,
        accession: str,
        lineage: str,
        payload: str | None,
    ) -> None:
        """Write one published BUSCO JSON summary when provided."""
        if payload is None:
            return
        busco_dir = source_outdir / "samples" / accession / "busco" / lineage
        self.write_text_file(busco_dir / "short_summary.json", payload)

    def write_staged_fasta(
        self,
        source_outdir: Path,
        *,
        accession: str,
        internal_id: str,
        sequence: str = "ACGTACGTACGT",
    ) -> Path:
        """Write one published staged FASTA and return the path."""
        staged_path = source_outdir / "samples" / accession / "staged" / f"{internal_id}.fasta"
        return self.write_text_file(staged_path, f">{internal_id}\n{sequence}\n")

    def fake_run_subprocess(
        self,
        stats_by_accession: dict[str, tuple[str, str, str]],
        fastani_similarity: float = 97.5,
    ) -> mock.Mock:
        """Build a subprocess double for assembly-stats and FastANI commands."""

        def _run(
            command: Sequence[str | Path],
            *,
            cwd: Path,
            env: dict[str, str] | None = None,
        ) -> subprocess.CompletedProcess[str]:
            command_name = Path(str(command[0])).name
            if command_name == "calculate_assembly_stats.sh":
                _ = env
                manifest_path = Path(str(command[2]))
                output_path = Path(str(command[4]))
                _header, rows = read_tsv_rows(manifest_path)
                stats_rows = [
                    {
                        "accession": row["accession"],
                        "n50": stats_by_accession[row["accession"]][0],
                        "scaffolds": stats_by_accession[row["accession"]][1],
                        "genome_size": stats_by_accession[row["accession"]][2],
                    }
                    for row in rows
                ]
                self.write_tsv_rows(
                    output_path,
                    ("accession", "n50", "scaffolds", "genome_size"),
                    stats_rows,
                )
                return subprocess.CompletedProcess(
                    [str(part) for part in command],
                    0,
                    "",
                    "",
                )

            if str(command[0]) == "fake-fastani":
                path_list = cwd / "fastani_paths.txt"
                paths = [
                    line.strip()
                    for line in path_list.read_text(encoding="utf-8").splitlines()
                    if line.strip()
                ]
                matrix_lines = [str(len(paths))]
                if paths:
                    matrix_lines.append(paths[0])
                for index in range(1, len(paths)):
                    values = ["80.0000"] * (index - 1)
                    similarity = "100.0000" if index == 0 else f"{fastani_similarity:.4f}"
                    values.append(similarity)
                    matrix_lines.append(f"{paths[index]} {' '.join(values)}")
                self.write_text_file(cwd / "fastani.tsv", "")
                self.write_text_file(
                    cwd / "fastani.tsv.matrix",
                    "\n".join(matrix_lines) + ("\n" if matrix_lines else ""),
                )
                return subprocess.CompletedProcess(
                    [str(part) for part in command],
                    0,
                    "",
                    "",
                )

            raise AssertionError(f"Unexpected subprocess command: {command}")

        return mock.Mock(side_effect=_run)

    def test_parse_command_prefix_supports_container_wrappers(self) -> None:
        """Split a wrapped FastANI command into argv tokens."""
        self.assertEqual(
            rescue_ani_from_results.parse_command_prefix(
                'singularity exec /tmp/fastani.sif fastANI'
            ),
            ["singularity", "exec", "/tmp/fastani.sif", "fastANI"],
        )

    def test_build_tool_wrapper_env_creates_a_seqtk_shim(self) -> None:
        """Expose a wrapped seqtk command through a temporary PATH shim."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            env = rescue_ani_from_results.build_tool_wrapper_env(
                "singularity exec /tmp/seqtk.sif seqtk",
                executable_name="seqtk",
                wrapper_dir=tmpdir,
            )

            assert env is not None
            self.assertTrue(env["PATH"].startswith(f"{tmpdir}{os.pathsep}"))
            wrapper_text = (tmpdir / "seqtk").read_text(encoding="utf-8")
            self.assertIn("exec singularity exec /tmp/seqtk.sif seqtk \"$@\"", wrapper_text)

    def test_parse_args_accepts_one_busco_option_with_multiple_values(self) -> None:
        """Accept one BUSCO option followed by multiple lineage values."""
        args = rescue_ani_from_results.parse_args(
            [
                "--source-outdir",
                "/tmp/source",
                "--metadata",
                "/tmp/metadata.tsv",
                "--outdir",
                "/tmp/rescued",
                "--busco-lineage",
                "bacillota_odb12",
                "mycoplasmatota_odb12",
            ]
        )

        self.assertEqual(args.busco_lineage, ["bacillota_odb12", "mycoplasmatota_odb12"])

    def test_parse_args_rejects_repeated_busco_lineage_option(self) -> None:
        """Fail clearly when the BUSCO lineage option is repeated."""
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as excinfo:
                rescue_ani_from_results.parse_args(
                    [
                        "--source-outdir",
                        "/tmp/source",
                        "--metadata",
                        "/tmp/metadata.tsv",
                        "--outdir",
                        "/tmp/rescued",
                        "--busco-lineage",
                        "bacillota_odb12",
                        "--busco-lineage",
                        "mycoplasmatota_odb12",
                    ]
                )

        self.assertEqual(excinfo.exception.code, 2)
        self.assertIn("--busco-lineage may only be supplied once.", stderr.getvalue())

    @unittest.skipUnless(SCIPY_AVAILABLE, "SciPy is required for ANI rescue integration tests.")
    def test_main_recovers_ani_outputs_and_partial_final_tables(self) -> None:
        """Rescue ANI outputs from published per-sample Barrnap, CheckM2, and BUSCO artefacts."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source_outdir = tmpdir / "source_results"
            rescue_outdir = tmpdir / "rescued_results"
            metadata = tmpdir / "metadata.tsv"

            validated_rows = [
                {
                    "accession": "ACC1",
                    "is_new": "false",
                    "assembly_level": "NA",
                    "genome_fasta": str(tmpdir / "raw" / "ACC1.fasta"),
                    "internal_id": "ID1",
                },
                {
                    "accession": "ACC2",
                    "is_new": "false",
                    "assembly_level": "NA",
                    "genome_fasta": str(tmpdir / "raw" / "ACC2.fasta"),
                    "internal_id": "ID2",
                },
                {
                    "accession": "ACC3",
                    "is_new": "false",
                    "assembly_level": "NA",
                    "genome_fasta": str(tmpdir / "raw" / "ACC3.fasta"),
                    "internal_id": "ID3",
                },
                {
                    "accession": "ACC4",
                    "is_new": "false",
                    "assembly_level": "NA",
                    "genome_fasta": str(tmpdir / "raw" / "ACC4.fasta"),
                    "internal_id": "ID4",
                },
            ]
            initial_status_rows = [
                self.make_initial_status_row(accession="ACC1", internal_id="ID1"),
                self.make_initial_status_row(accession="ACC2", internal_id="ID2"),
                self.make_initial_status_row(accession="ACC3", internal_id="ID3"),
                self.make_initial_status_row(accession="ACC4", internal_id="ID4"),
            ]
            self.write_tsv_rows(
                source_outdir / "tables" / "validated_samples.tsv",
                ("accession", "is_new", "assembly_level", "genome_fasta", "internal_id"),
                validated_rows,
            )
            self.write_tsv_rows(
                source_outdir / "tables" / "sample_status.tsv",
                SAMPLE_STATUS_COLUMNS,
                initial_status_rows,
            )
            self.write_text_file(
                metadata,
                "\n".join(
                    [
                        "Accession\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings",
                        "ACC1\tGenome One\tChromosome\tNA\tNA\tNA\tNA",
                        "ACC2\tGenome Two\tScaffold\tNA\tNA\tNA\tNA",
                        "ACC3\tGenome Three\tScaffold\tNA\tNA\tNA\tNA",
                        "ACC4\tGenome Four\tScaffold\tNA\tNA\tNA\tUnverified source organism",
                    ]
                )
                + "\n",
            )

            for sample in validated_rows:
                self.write_staged_fasta(
                    source_outdir,
                    accession=sample["accession"],
                    internal_id=sample["internal_id"],
                )

            self.write_barrnap_outputs(source_outdir, accession="ACC1", mode="yes")
            self.write_barrnap_outputs(source_outdir, accession="ACC2", mode="partial")
            self.write_barrnap_outputs(source_outdir, accession="ACC3", mode="yes")
            self.write_barrnap_outputs(source_outdir, accession="ACC4", mode="yes")

            self.write_checkm2_report(
                source_outdir, accession="ACC1", translation_table=4, completeness="80", contamination="2"
            )
            self.write_checkm2_report(
                source_outdir, accession="ACC1", translation_table=11, completeness="95", contamination="1"
            )
            self.write_checkm2_report(
                source_outdir, accession="ACC2", translation_table=4, completeness="82", contamination="2"
            )
            self.write_checkm2_report(
                source_outdir, accession="ACC2", translation_table=11, completeness="94", contamination="1"
            )
            self.write_checkm2_report(
                source_outdir, accession="ACC3", translation_table=4, completeness="70", contamination="7"
            )
            self.write_checkm2_report(
                source_outdir, accession="ACC3", translation_table=11, completeness="90", contamination="8"
            )
            self.write_checkm2_report(
                source_outdir, accession="ACC4", translation_table=4, completeness="95", contamination="1"
            )
            self.write_checkm2_report(
                source_outdir, accession="ACC4", translation_table=11, completeness="80", contamination="2"
            )

            busco_payload_a = '{"results": {"C": 98.0, "S": 98.0, "D": 0.0, "F": 1.0, "M": 1.0, "n": 200}}'
            busco_payload_b = '{"results": {"C": 97.0, "S": 97.0, "D": 0.0, "F": 2.0, "M": 1.0, "n": 180}}'
            for accession in ("ACC1", "ACC2", "ACC3", "ACC4"):
                self.write_busco_summary(
                    source_outdir,
                    accession=accession,
                    lineage="bacillota_odb12",
                    payload=busco_payload_a,
                )
                self.write_busco_summary(
                    source_outdir,
                    accession=accession,
                    lineage="mycoplasmatota_odb12",
                    payload=busco_payload_b,
                )

            with mock.patch.object(
                rescue_ani_from_results,
                "run_subprocess",
                self.fake_run_subprocess(
                    {
                        "ACC1": ("50000", "2", "800000"),
                        "ACC2": ("45000", "3", "780000"),
                        "ACC3": ("42000", "4", "760000"),
                        "ACC4": ("47000", "2", "790000"),
                    }
                ),
            ):
                exit_code = rescue_ani_from_results.main(
                    [
                        "--source-outdir",
                        str(source_outdir),
                        "--metadata",
                        str(metadata),
                        "--outdir",
                        str(rescue_outdir),
                        "--busco-lineage",
                        "bacillota_odb12",
                        "mycoplasmatota_odb12",
                        "--fastani-binary",
                        "fake-fastani",
                    ]
                )

            self.assertEqual(exit_code, 0)

            _sixteen_s_header, sixteen_s_rows = read_tsv_rows(
                rescue_outdir / "tables" / "16s_statuses.tsv"
            )
            sixteen_s_by_accession = {row["accession"]: row for row in sixteen_s_rows}
            self.assertEqual(sixteen_s_by_accession["ACC1"]["16S"], "Yes")
            self.assertEqual(sixteen_s_by_accession["ACC2"]["16S"], "partial")
            self.assertEqual(sixteen_s_by_accession["ACC3"]["16S"], "Yes")

            _metadata_header, ani_metadata_rows = read_tsv_rows(
                rescue_outdir / "cohort" / "fastani" / "ani_metadata.tsv"
            )
            self.assertEqual(
                [row["accession"] for row in ani_metadata_rows],
                ["ACC1", "ACC4"],
            )

            _exclusions_header, exclusion_rows = read_tsv_rows(
                rescue_outdir / "cohort" / "fastani" / "ani_exclusions.tsv"
            )
            exclusions_by_accession = {row["accession"]: row for row in exclusion_rows}
            self.assertEqual(exclusions_by_accession["ACC1"]["ani_included"], "true")
            self.assertEqual(exclusions_by_accession["ACC4"]["ani_included"], "true")
            self.assertEqual(exclusions_by_accession["ACC2"]["ani_included"], "false")
            self.assertEqual(
                exclusions_by_accession["ACC2"]["ani_exclusion_reason"],
                "partial_16s",
            )
            self.assertEqual(exclusions_by_accession["ACC3"]["ani_included"], "false")
            self.assertEqual(
                exclusions_by_accession["ACC3"]["ani_exclusion_reason"],
                "low_quality",
            )

            _cluster_header, cluster_rows = read_tsv_rows(
                rescue_outdir / "cohort" / "ani_clusters" / "cluster.tsv"
            )
            self.assertEqual(
                [row["Accession"] for row in cluster_rows],
                ["ACC1", "ACC4"],
            )

            _summary_header, summary_rows = read_tsv_rows(
                rescue_outdir / "cohort" / "ani_clusters" / "ani_summary.tsv"
            )
            summary_by_accession = {row["Accession"]: row for row in summary_rows}
            self.assertEqual(summary_by_accession["ACC1"]["Cluster_ID"], "C000001")
            self.assertEqual(summary_by_accession["ACC4"]["Cluster_ID"], "C000001")

            _representatives_header, representative_rows = read_tsv_rows(
                rescue_outdir / "cohort" / "ani_clusters" / "ani_representatives.tsv"
            )
            self.assertEqual(len(representative_rows), 1)
            self.assertEqual(representative_rows[0]["Representative_Accession"], "ACC1")

            _master_header, master_rows = read_tsv_rows(
                rescue_outdir / "tables" / "master_table.tsv"
            )
            master_by_accession = {row["Accession"]: row for row in master_rows}
            self.assertEqual(master_by_accession["ACC1"]["16S"], "Yes")
            self.assertEqual(master_by_accession["ACC1"]["Cluster_ID"], "C000001")
            self.assertEqual(master_by_accession["ACC1"]["Codetta_Genetic_Code"], "NA")
            self.assertEqual(master_by_accession["ACC1"]["CRISPRS"], "NA")

            _status_header, status_rows = read_tsv_rows(
                rescue_outdir / "tables" / "sample_status.tsv"
            )
            status_by_accession = {row["accession"]: row for row in status_rows}
            self.assertEqual(status_by_accession["ACC1"]["barrnap_status"], "done")
            self.assertEqual(status_by_accession["ACC1"]["taxonomy_status"], "na")
            self.assertEqual(status_by_accession["ACC1"]["codetta_status"], "na")
            self.assertEqual(status_by_accession["ACC1"]["ani_included"], "true")
            self.assertEqual(status_by_accession["ACC2"]["ani_included"], "false")
            self.assertEqual(
                status_by_accession["ACC2"]["ani_exclusion_reason"],
                "partial_16s",
            )
            self.assertEqual(status_by_accession["ACC3"]["ani_included"], "false")
            self.assertEqual(
                status_by_accession["ACC3"]["ani_exclusion_reason"],
                "low_quality",
            )

            _provenance_header, provenance_rows = read_tsv_rows(
                rescue_outdir / "tables" / "rescue_provenance.tsv"
            )
            self.assertEqual(len(provenance_rows), 1)
            self.assertEqual(provenance_rows[0]["source_commit"], "595edef")
            self.assertEqual(provenance_rows[0]["ani_included_sample_count"], "2")
            self.assertEqual(provenance_rows[0]["ani_excluded_sample_count"], "2")

    def test_main_fails_when_no_staged_or_fallback_genome_exists(self) -> None:
        """Name the accession and expected staged path when FASTA resolution fails."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source_outdir = tmpdir / "source_results"
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings\n"
                "ACC1\tGenome One\tScaffold\tNA\tNA\tNA\tNA\n",
            )
            self.write_tsv_rows(
                source_outdir / "tables" / "validated_samples.tsv",
                ("accession", "is_new", "assembly_level", "genome_fasta", "internal_id"),
                [
                    {
                        "accession": "ACC1",
                        "is_new": "false",
                        "assembly_level": "NA",
                        "genome_fasta": str(tmpdir / "missing" / "ACC1.fasta"),
                        "internal_id": "ID1",
                    }
                ],
            )
            self.write_tsv_rows(
                source_outdir / "tables" / "sample_status.tsv",
                SAMPLE_STATUS_COLUMNS,
                [self.make_initial_status_row(accession="ACC1", internal_id="ID1")],
            )

            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                exit_code = rescue_ani_from_results.main(
                    [
                        "--source-outdir",
                        str(source_outdir),
                        "--metadata",
                        str(metadata),
                        "--outdir",
                        str(tmpdir / "rescued_results"),
                        "--busco-lineage",
                        "bacillota_odb12",
                        "mycoplasmatota_odb12",
                    ]
                )

            self.assertEqual(exit_code, 1)
            error_text = stderr.getvalue()
            self.assertIn("ACC1", error_text)
            self.assertIn("samples/ACC1/staged/ID1.fasta", error_text)

    @unittest.skipUnless(SCIPY_AVAILABLE, "SciPy is required for ANI rescue integration tests.")
    def test_main_uses_busco_lineage_order_for_primary_ani_gating(self) -> None:
        """Use the first supplied BUSCO lineage as the ANI primary column."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source_outdir = tmpdir / "source_results"
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings",
                        "ACC1\tGenome One\tScaffold\tNA\tNA\tNA\tNA",
                        "ACC2\tGenome Two\tScaffold\tNA\tNA\tNA\tNA",
                    ]
                )
                + "\n",
            )
            validated_rows = [
                {
                    "accession": "ACC1",
                    "is_new": "false",
                    "assembly_level": "NA",
                    "genome_fasta": str(tmpdir / "raw" / "ACC1.fasta"),
                    "internal_id": "ID1",
                },
                {
                    "accession": "ACC2",
                    "is_new": "false",
                    "assembly_level": "NA",
                    "genome_fasta": str(tmpdir / "raw" / "ACC2.fasta"),
                    "internal_id": "ID2",
                },
            ]
            self.write_tsv_rows(
                source_outdir / "tables" / "validated_samples.tsv",
                ("accession", "is_new", "assembly_level", "genome_fasta", "internal_id"),
                validated_rows,
            )
            self.write_tsv_rows(
                source_outdir / "tables" / "sample_status.tsv",
                SAMPLE_STATUS_COLUMNS,
                [
                    self.make_initial_status_row(accession="ACC1", internal_id="ID1"),
                    self.make_initial_status_row(accession="ACC2", internal_id="ID2"),
                ],
            )
            for sample in validated_rows:
                self.write_staged_fasta(
                    source_outdir,
                    accession=sample["accession"],
                    internal_id=sample["internal_id"],
                )
                self.write_barrnap_outputs(source_outdir, accession=sample["accession"], mode="yes")
                self.write_checkm2_report(
                    source_outdir,
                    accession=sample["accession"],
                    translation_table=4,
                    completeness="80",
                    contamination="2",
                )
                self.write_checkm2_report(
                    source_outdir,
                    accession=sample["accession"],
                    translation_table=11,
                    completeness="95",
                    contamination="1",
                )

            complete_payload = '{"results": {"C": 98.0, "S": 98.0, "D": 0.0, "F": 1.0, "M": 1.0, "n": 200}}'
            self.write_busco_summary(
                source_outdir,
                accession="ACC1",
                lineage="bacillota_odb12",
                payload=complete_payload,
            )
            self.write_busco_summary(
                source_outdir,
                accession="ACC2",
                lineage="bacillota_odb12",
                payload=None,
            )
            self.write_busco_summary(
                source_outdir,
                accession="ACC1",
                lineage="mycoplasmatota_odb12",
                payload=complete_payload,
            )
            self.write_busco_summary(
                source_outdir,
                accession="ACC2",
                lineage="mycoplasmatota_odb12",
                payload=complete_payload,
            )

            with mock.patch.object(
                rescue_ani_from_results,
                "run_subprocess",
                self.fake_run_subprocess(
                    {
                        "ACC1": ("50000", "2", "800000"),
                        "ACC2": ("49000", "2", "790000"),
                    }
                ),
            ):
                exit_code_one = rescue_ani_from_results.main(
                    [
                        "--source-outdir",
                        str(source_outdir),
                        "--metadata",
                        str(metadata),
                        "--outdir",
                        str(tmpdir / "rescued_order_one"),
                        "--busco-lineage",
                        "bacillota_odb12",
                        "mycoplasmatota_odb12",
                        "--fastani-binary",
                        "fake-fastani",
                    ]
                )
                exit_code_two = rescue_ani_from_results.main(
                    [
                        "--source-outdir",
                        str(source_outdir),
                        "--metadata",
                        str(metadata),
                        "--outdir",
                        str(tmpdir / "rescued_order_two"),
                        "--busco-lineage",
                        "mycoplasmatota_odb12",
                        "bacillota_odb12",
                        "--fastani-binary",
                        "fake-fastani",
                    ]
                )

            self.assertEqual(exit_code_one, 0)
            self.assertEqual(exit_code_two, 0)

            _order_one_header, order_one_rows = read_tsv_rows(
                tmpdir / "rescued_order_one" / "cohort" / "fastani" / "ani_metadata.tsv"
            )
            _order_two_header, order_two_rows = read_tsv_rows(
                tmpdir / "rescued_order_two" / "cohort" / "fastani" / "ani_metadata.tsv"
            )
            self.assertEqual([row["accession"] for row in order_one_rows], ["ACC1"])
            self.assertEqual(
                [row["accession"] for row in order_two_rows],
                ["ACC1", "ACC2"],
            )


if __name__ == "__main__":
    unittest.main()
