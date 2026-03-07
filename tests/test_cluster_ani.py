"""Helper-level tests for ANI metadata parsing."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import cluster_ani  # noqa: E402


class ClusterAniHelperTestCase(unittest.TestCase):
    """Cover header normalization, BUSCO parsing, and ANI metadata gating."""

    def test_normalize_header_matches_nextflow_style(self) -> None:
        """Normalize mixed punctuation and case into lowercase underscore tokens."""
        self.assertEqual(cluster_ani.normalize_header(" Assembly Level "), "assembly_level")
        self.assertEqual(cluster_ani.normalize_header("BUSCO-bacillota_odb12"), "busco_bacillota_odb12")

    def test_try_parse_busco_extracts_complete_and_missing(self) -> None:
        """Read the compact BUSCO string format emitted by the BUSCO summariser."""
        parsed = cluster_ani.try_parse_busco(
            "C:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:220"
        )
        self.assertEqual(parsed, (96.5, 1.5))
        self.assertIsNone(cluster_ani.try_parse_busco("NA"))

    def test_build_genome_from_row_accepts_valid_metadata(self) -> None:
        """Produce a Genome record when all ANI-scoring metadata are present."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            genome_path = tmpdir / "genome.fna"
            genome_path.write_text(">a\nACGT\n", encoding="utf-8")
            genome, reasons = cluster_ani.build_genome_from_row(
                {
                    "path": str(genome_path),
                    "assembly_level": "Scaffold",
                    "gcode": "11",
                    "checkm2_completeness": "95.0",
                    "checkm2_contamination": "1.5",
                    "n50": "50000",
                    "scaffolds": "4",
                    "genome_size": "1000000",
                    "primary_busco": "C:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:220",
                    "organism_name": "Example",
                },
                acc="ACC1",
                busco_column="primary_busco",
            )

            self.assertEqual(reasons, [])
            assert genome is not None
            self.assertEqual(genome.Accession, "ACC1")
            self.assertEqual(genome.Gcode, 11)
            self.assertEqual(genome.Assembly_Level, "Scaffold")
            self.assertEqual(genome.BUSCO_C, 96.5)
            self.assertEqual(genome.BUSCO_M, 1.5)

    def test_build_genome_from_row_accumulates_exclusion_reasons(self) -> None:
        """Collect all missing or invalid ANI-scoring fields before exclusion."""
        genome, reasons = cluster_ani.build_genome_from_row(
            {
                "path": "missing.fna",
                "assembly_level": "unknown",
                "gcode": "NA",
                "checkm2_completeness": "NA",
                "checkm2_contamination": "bad",
                "n50": "1.5",
                "scaffolds": "",
                "genome_size": "NA",
                "primary_busco": "NA",
            },
            acc="ACC2",
            busco_column="primary_busco",
        )

        self.assertIsNone(genome)
        self.assertEqual(
            reasons,
            [
                "path_not_found",
                "invalid_assembly_level",
                "invalid_gcode",
                "missing_checkm2_completeness",
                "missing_checkm2_contamination",
                "missing_n50",
                "missing_scaffolds",
                "missing_genome_size",
                "missing_primary_busco",
            ],
        )

    def test_load_ani_metadata_supports_distinct_matrix_name_column(self) -> None:
        """Match PHYLIP matrix rows via a dedicated matrix-name column."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            genome_path = tmpdir / "genome.fna"
            genome_path.write_text(">a\nACGT\n", encoding="utf-8")
            metadata_tsv = tmpdir / "ani_metadata.tsv"
            metadata_tsv.write_text(
                "\n".join(
                    [
                        "accession\tmatrix_name\tpath\tassembly_level\tgcode\tcheckm2_completeness\tcheckm2_contamination\tn50\tscaffolds\tgenome_size\tBUSCO_bacillota_odb12",
                        f"ACC1\tfastani_inputs/ACC1.fasta\t{genome_path}\tScaffold\t11\t95\t1.5\t50000\t4\t1000000\tC:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:220",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            metadata, eligible_names = cluster_ani.load_ani_metadata(
                ani_metadata=metadata_tsv,
                matrix_names=["fastani_inputs/ACC1.fasta"],
                busco_column="BUSCO_bacillota_odb12",
                matrix_name_column="matrix_name",
            )

            self.assertEqual(eligible_names, ["fastani_inputs/ACC1.fasta"])
            self.assertIn("fastani_inputs/ACC1.fasta", metadata)
            self.assertEqual(metadata["fastani_inputs/ACC1.fasta"].Accession, "ACC1")


if __name__ == "__main__":
    unittest.main()
