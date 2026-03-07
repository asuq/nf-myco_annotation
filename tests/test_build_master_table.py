"""Tests for final master-table assembly."""

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

import build_master_table  # noqa: E402
import master_table_contract  # noqa: E402


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a CSV file into a header and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        return reader.fieldnames, list(reader)


class BuildMasterTableTestCase(unittest.TestCase):
    """Cover metadata preservation, supplemental fills, and duplicate detection."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_main_builds_master_table_with_preserved_metadata_and_na_fills(self) -> None:
        """Preserve metadata order and fill new-sample metadata from validated extras."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\torganism_name\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tKnown one\tid_1",
                        "ACC2\ttrue\tScaffold\t/path/two.fna\tCandidate two\tid_2",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings",
                        "ACC1\t123\tKnown one\tNA",
                    ]
                )
                + "\n",
            )
            taxonomy = self.write_text_file(
                tmpdir / "taxonomy.tsv",
                "\n".join(
                    [
                        "Tax_ID\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies",
                        "123\tBacteria\tFirmicutes\tMollicutes\tMycoplasmatales\tMycoplasmataceae\tMycoplasma\tsp1",
                    ]
                )
                + "\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tLow_quality\tcheckm2_status\twarnings",
                        "ACC1\t95\t80\t2\t1\t0.9\t0.8\t900\t850\t800\t780\t4\tfalse\tdone\t",
                        "ACC2\t70\t72\t5\t4\t0.7\t0.71\t700\t710\t600\t610\tNA\tNA\tdone\tgcode_na",
                    ]
                )
                + "\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\tinclude_in_all_best_16S\twarnings",
                        "ACC1\tYes\th1\t1500\ttrue\t",
                        "ACC2\tNo\tNA\tNA\tfalse\t",
                    ]
                )
                + "\n",
            )
            busco_one = self.write_text_file(
                tmpdir / "busco_b.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings",
                        "ACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t",
                        "ACC2\tbacillota_odb12\tNA\tfailed\tbusco_summary_failed",
                    ]
                )
                + "\n",
            )
            busco_two = self.write_text_file(
                tmpdir / "busco_m.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_mycoplasmatota_odb12\tbusco_status\twarnings",
                        "ACC1\tmycoplasmatota_odb12\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:180\tdone\t",
                        "ACC2\tmycoplasmatota_odb12\tNA\tfailed\tbusco_summary_failed",
                    ]
                )
                + "\n",
            )
            ccfinder = self.write_text_file(
                tmpdir / "ccfinder.tsv",
                "\n".join(
                    [
                        "accession\tCRISPRS\tSPACERS_SUM\tCRISPR_FRAC\tccfinder_status\twarnings",
                        "ACC1\t2\t7\t0.1\tdone\t",
                    ]
                )
                + "\n",
            )
            ani = self.write_text_file(
                tmpdir / "ani.tsv",
                "\n".join(
                    [
                        "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore",
                        "ACC1\tcluster_1\tyes\t100\t0.95",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "master_table.csv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(ROOT / "assets" / "master_table_append_columns.txt"),
                    "--taxonomy",
                    str(taxonomy),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco_one),
                    "--busco",
                    str(busco_two),
                    "--ccfinder-strains",
                    str(ccfinder),
                    "--ani",
                    str(ani),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            header, rows = read_csv_rows(output)
            expected_header = master_table_contract.build_master_table_columns(
                ["Accession", "Tax_ID", "Organism_Name", "Atypical_Warnings"]
            )
            self.assertEqual(header, expected_header)
            self.assertEqual(len(rows), 2)

            self.assertEqual(rows[0]["Accession"], "ACC1")
            self.assertEqual(rows[0]["Tax_ID"], "123")
            self.assertEqual(rows[0]["superkingdom"], "Bacteria")
            self.assertEqual(rows[0]["16S"], "Yes")
            self.assertEqual(rows[0]["CRISPRS"], "2")
            self.assertEqual(rows[0]["Cluster_ID"], "cluster_1")

            self.assertEqual(rows[1]["Accession"], "ACC2")
            self.assertEqual(rows[1]["Tax_ID"], "NA")
            self.assertEqual(rows[1]["Organism_Name"], "Candidate two")
            self.assertEqual(rows[1]["Atypical_Warnings"], "NA")
            self.assertEqual(rows[1]["superkingdom"], "NA")
            self.assertEqual(rows[1]["16S"], "No")
            self.assertEqual(rows[1]["CRISPRS"], "NA")
            self.assertEqual(rows[1]["Cluster_ID"], "NA")

    def test_main_rejects_duplicate_rows_in_derived_tables(self) -> None:
        """Fail when a derived table creates ambiguous duplicate accession joins."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nACC1\t123\n",
            )
            duplicate_checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tLow_quality\tcheckm2_status\twarnings",
                        "ACC1\t95\t80\t2\t1\t0.9\t0.8\t900\t850\t800\t780\t4\tfalse\tdone\t",
                        "ACC1\t95\t80\t2\t1\t0.9\t0.8\t900\t850\t800\t780\t4\tfalse\tdone\t",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "master_table.csv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(ROOT / "assets" / "master_table_append_columns.txt"),
                    "--checkm2",
                    str(duplicate_checkm2),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_main_rejects_duplicate_internal_ids_in_validated_samples(self) -> None:
        """Fail when the validated manifest maps two rows to the same internal ID."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tdup_id",
                        "ACC2\ttrue\tScaffold\t/path/two.fna\tdup_id",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nACC1\t123\n",
            )
            output = tmpdir / "master_table.csv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(ROOT / "assets" / "master_table_append_columns.txt"),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_main_rejects_metadata_with_both_accession_header_variants(self) -> None:
        """Fail when metadata contains both `accession` and `Accession` key columns."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "accession\tAccession\tTax_ID",
                        "ACC1\tACC1\t123",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "master_table.csv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(ROOT / "assets" / "master_table_append_columns.txt"),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_main_rejects_duplicate_busco_lineage_inputs(self) -> None:
        """Fail when the same BUSCO lineage column is supplied more than once."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nACC1\t123\n",
            )
            busco_one = self.write_text_file(
                tmpdir / "busco_1.tsv",
                "\n".join(
                    [
                        "accession\tBUSCO_bacillota_odb12",
                        "ACC1\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
                    ]
                )
                + "\n",
            )
            busco_two = self.write_text_file(
                tmpdir / "busco_2.tsv",
                "\n".join(
                    [
                        "accession\tBUSCO_bacillota_odb12",
                        "ACC1\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:200",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "master_table.csv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(ROOT / "assets" / "master_table_append_columns.txt"),
                    "--busco",
                    str(busco_one),
                    "--busco",
                    str(busco_two),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_main_rejects_append_contract_drift(self) -> None:
        """Fail when the provided append-column contract no longer matches the code contract."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nACC1\t123\n",
            )
            append_columns = self.write_text_file(
                tmpdir / "append_columns.txt",
                "\n".join(
                    master_table_contract.build_append_columns()[:-1]
                    + ["Unexpected_Column"]
                )
                + "\n",
            )
            output = tmpdir / "master_table.csv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(append_columns),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
