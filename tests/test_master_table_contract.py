"""Tests for the master-table output contract."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import master_table_contract  # noqa: E402


class MasterTableContractTestCase(unittest.TestCase):
    """Verify the locked derived-column order and asset consistency."""

    def test_default_contract_matches_repository_asset(self) -> None:
        """Keep the asset file aligned with the default code contract."""
        master_table_contract.validate_default_append_columns_asset()
        self.assertEqual(
            master_table_contract.read_append_columns_asset(),
            master_table_contract.build_append_columns(),
        )

    def test_build_master_table_columns_supports_custom_busco_order(self) -> None:
        """Insert BUSCO columns in caller-provided lineage order."""
        columns = master_table_contract.build_master_table_columns(
            metadata_columns=["Accession", "Tax_ID"],
            busco_lineages=["lineage_a", "lineage_b", "lineage_c"],
        )

        self.assertEqual(columns[:2], ["Accession", "Tax_ID"])
        self.assertIn("BUSCO_lineage_a", columns)
        self.assertIn("BUSCO_lineage_b", columns)
        self.assertIn("BUSCO_lineage_c", columns)
        self.assertLess(columns.index("BUSCO_lineage_a"), columns.index("BUSCO_lineage_b"))
        self.assertLess(columns.index("BUSCO_lineage_b"), columns.index("BUSCO_lineage_c"))
        self.assertEqual(columns[-4:], list(master_table_contract.ANI_COLUMNS))

    def test_codetta_columns_follow_gcode_immediately(self) -> None:
        """Keep the Codetta fields directly after the chosen gcode column."""
        columns = master_table_contract.build_append_columns()

        self.assertEqual(columns[columns.index("Gcode") + 1], "Codetta_Genetic_Code")
        self.assertEqual(
            columns[columns.index("Codetta_Genetic_Code") + 1],
            "Codetta_NCBI_Table_Candidates",
        )

    def test_normalise_busco_lineages_rejects_duplicates(self) -> None:
        """Fail when BUSCO lineage names repeat after normalisation."""
        with self.assertRaises(ValueError):
            master_table_contract.normalise_busco_lineages(
                ["bacillota_odb12", "bacillota_odb12"]
            )

    def test_build_master_table_columns_rejects_metadata_overlap(self) -> None:
        """Fail when metadata columns collide with derived output columns."""
        with self.assertRaises(ValueError):
            master_table_contract.build_master_table_columns(["Accession", "Gcode"])


if __name__ == "__main__":
    unittest.main()
