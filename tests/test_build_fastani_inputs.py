"""Tests for ANI eligibility and FastANI input preparation."""

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

import build_fastani_inputs  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class BuildFastAniInputsTestCase(unittest.TestCase):
    """Cover ANI eligibility gating and primary BUSCO selection."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write text to a test file and return the path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_builds_inputs_for_eligible_and_exception_samples(self) -> None:
        """Include valid samples and the locked atypical unverified-source exception."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            staged_a = self.write_text_file(tmpdir / "ACC1.fasta", ">a\nACGT\n")
            staged_b = self.write_text_file(tmpdir / "ACC2.fasta", ">b\nACGT\n")
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        f"ACC1\tfalse\tNA\t{staged_a}\tACC1",
                        f"ACC2\ttrue\tScaffold\t{staged_b}\tACC2",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings",
                        "ACC1\tGenome One\tChromosome\t50000\t2\t800000\tNA",
                        "ACC2\tGenome Two\tContig\t40000\t5\t700000\tUnverified source organism",
                    ]
                )
                + "\n",
            )
            staged_manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "\n".join(
                    [
                        "accession\tinternal_id\tstaged_filename",
                        "ACC1\tACC1\tACC1.fasta",
                        "ACC2\tACC2\tACC2.fasta",
                    ]
                )
                + "\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tGcode\tLow_quality",
                        "ACC1\tNA\t95\tNA\t1\t11\tfalse",
                        "ACC2\t90\tNA\t2\tNA\t4\tfalse",
                    ]
                )
                + "\n",
            )
            sixteen_s = self.write_text_file(
                tmpdir / "16s.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings",
                        "ACC1\tYes\th1\t1500\ttrue\t",
                        "ACC2\tYes\th2\t1490\tfalse\t",
                    ]
                )
                + "\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco_bacillota.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings",
                        "ACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t",
                        "ACC2\tbacillota_odb12\tC:96.0%[S:96.0%,D:0.0%],F:2.0%,M:2.0%,n:200\tdone\t",
                    ]
                )
                + "\n",
            )
            busco_mycoplasmatota = self.write_text_file(
                tmpdir / "busco_mycoplasmatota.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_mycoplasmatota_odb12\tbusco_status\twarnings",
                        "ACC1\tmycoplasmatota_odb12\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:180\tdone\t",
                    ]
                )
                + "\n",
            )
            outdir = tmpdir / "out"

            exit_code = build_fastani_inputs.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--staged-manifest",
                    str(staged_manifest),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(sixteen_s),
                    "--busco",
                    str(busco),
                    "--busco",
                    str(busco_mycoplasmatota),
                    "--primary-busco-column",
                    "BUSCO_bacillota_odb12",
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            metadata_rows = read_tsv(outdir / "ani_metadata.tsv")
            self.assertEqual([row["accession"] for row in metadata_rows], ["ACC1", "ACC2"])
            self.assertEqual(metadata_rows[0]["gcode"], "11")
            self.assertEqual(metadata_rows[0]["checkm2_completeness"], "95")
            self.assertEqual(metadata_rows[1]["assembly_level"], "Scaffold")
            self.assertTrue((outdir / "fastani_inputs" / "ACC1.fasta").exists())
            self.assertEqual(
                (outdir / "fastani_paths.txt").read_text(encoding="utf-8").splitlines(),
                ["fastani_inputs/ACC1.fasta", "fastani_inputs/ACC2.fasta"],
            )
            exclusion_rows = {
                row["accession"]: row for row in read_tsv(outdir / "ani_exclusions.tsv")
            }
            self.assertEqual(exclusion_rows["ACC1"]["ani_included"], "true")
            self.assertEqual(exclusion_rows["ACC2"]["ani_included"], "true")

    def test_excludes_low_quality_partial_and_non_exception_atypical_samples(self) -> None:
        """Record stable exclusion reasons for ineligible ANI samples."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            staged = self.write_text_file(tmpdir / "ACC3.fasta", ">a\nACGT\n")
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\n"
                f"ACC3\tfalse\tNA\t{staged}\tACC3\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings\n"
                "ACC3\tGenome Three\tScaffold\t50000\t10\t800000\tHost-associated\n",
            )
            staged_manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC3\tACC3\tACC3.fasta\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tGcode\tLow_quality\n"
                "ACC3\tNA\t94\tNA\t2\t11\ttrue\n",
            )
            sixteen_s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings\n"
                "ACC3\tpartial\th3\t900\tfalse\t\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\n"
                "ACC3\tbacillota_odb12\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:200\tdone\t\n",
            )
            outdir = tmpdir / "out"

            exit_code = build_fastani_inputs.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--staged-manifest",
                    str(staged_manifest),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(sixteen_s),
                    "--busco",
                    str(busco),
                    "--primary-busco-column",
                    "BUSCO_bacillota_odb12",
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertEqual(read_tsv(outdir / "ani_metadata.tsv"), [])
            exclusion_rows = read_tsv(outdir / "ani_exclusions.tsv")
            self.assertEqual(len(exclusion_rows), 1)
            self.assertEqual(exclusion_rows[0]["ani_included"], "false")
            self.assertEqual(
                exclusion_rows[0]["ani_exclusion_reason"],
                "low_quality;partial_16s;atypical",
            )
            self.assertTrue((outdir / "fastani_inputs").is_dir())
            self.assertEqual(list((outdir / "fastani_inputs").iterdir()), [])

    def test_uses_in_house_assembly_stats_when_metadata_metrics_are_missing(self) -> None:
        """Use computed stats for ANI eligibility and ANI metadata output."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            staged = self.write_text_file(tmpdir / "ACC4.fasta", ">a\nACGT\n")
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\n"
                f"ACC4\ttrue\tComplete Genome\t{staged}\tACC4\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings\n"
                "ACC4\tGenome Four\tNA\tNA\tNA\tNA\tNA\n",
            )
            staged_manifest = self.write_text_file(
                tmpdir / "staged_manifest.tsv",
                "accession\tinternal_id\tstaged_filename\nACC4\tACC4\tACC4.fasta\n",
            )
            assembly_stats = self.write_text_file(
                tmpdir / "assembly_stats.tsv",
                "accession\tn50\tscaffolds\tgenome_size\nACC4\t40000\t3\t120000\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tGcode\tLow_quality\n"
                "ACC4\tNA\t95\tNA\t1\t11\tfalse\n",
            )
            sixteen_s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings\n"
                "ACC4\tYes\th4\t1500\ttrue\t\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\n"
                "ACC4\tbacillota_odb12\tC:99.0%[S:99.0%,D:0.0%],F:0.0%,M:1.0%,n:200\tdone\t\n",
            )
            outdir = tmpdir / "out"

            exit_code = build_fastani_inputs.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--staged-manifest",
                    str(staged_manifest),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(sixteen_s),
                    "--busco",
                    str(busco),
                    "--primary-busco-column",
                    "BUSCO_bacillota_odb12",
                    "--assembly-stats",
                    str(assembly_stats),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            metadata_rows = read_tsv(outdir / "ani_metadata.tsv")
            self.assertEqual(len(metadata_rows), 1)
            self.assertEqual(metadata_rows[0]["n50"], "40000")
            self.assertEqual(metadata_rows[0]["scaffolds"], "3")
            self.assertEqual(metadata_rows[0]["genome_size"], "120000")
            exclusion_rows = read_tsv(outdir / "ani_exclusions.tsv")
            self.assertEqual(exclusion_rows[0]["ani_included"], "true")
            self.assertEqual(exclusion_rows[0]["ani_exclusion_reason"], "")


if __name__ == "__main__":
    unittest.main()
