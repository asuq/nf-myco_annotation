"""Define the locked master-table column contract for v1."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence


DEFAULT_BUSCO_LINEAGES = ("bacillota_odb12", "mycoplasmatota_odb12")
TAXONOMY_COLUMNS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
CHECKM2_COLUMNS = (
    "Completeness_gcode4",
    "Completeness_gcode11",
    "Contamination_gcode4",
    "Contamination_gcode11",
    "Coding_Density_gcode4",
    "Coding_Density_gcode11",
    "Average_Gene_Length_gcode4",
    "Average_Gene_Length_gcode11",
    "Total_Coding_Sequences_gcode4",
    "Total_Coding_Sequences_gcode11",
)
GCODE_QC_COLUMNS = ("Gcode", "Low_quality", "16S")
CRISPR_COLUMNS = ("CRISPRS", "SPACERS_SUM", "CRISPR_FRAC")
ANI_COLUMNS = ("Cluster_ID", "Is_Representative", "ANI_to_Representative", "Score")
DEFAULT_APPEND_COLUMNS_ASSET = (
    Path(__file__).resolve().parents[1] / "assets" / "master_table_append_columns.txt"
)


def normalise_busco_lineages(
    busco_lineages: Sequence[str] | None = None,
) -> tuple[str, ...]:
    """Validate BUSCO lineage names and preserve their order."""
    raw_lineages = busco_lineages or DEFAULT_BUSCO_LINEAGES
    normalised = tuple(lineage.strip() for lineage in raw_lineages if lineage.strip())
    if not normalised:
        raise ValueError("At least one BUSCO lineage is required.")
    if len(set(normalised)) != len(normalised):
        raise ValueError("BUSCO lineages must be unique.")
    return normalised


def build_append_columns(busco_lineages: Sequence[str] | None = None) -> list[str]:
    """Return the appended derived columns in locked output order."""
    busco_columns = [
        f"BUSCO_{lineage}" for lineage in normalise_busco_lineages(busco_lineages)
    ]
    return [
        *TAXONOMY_COLUMNS,
        *CHECKM2_COLUMNS,
        *GCODE_QC_COLUMNS,
        *busco_columns,
        *CRISPR_COLUMNS,
        *ANI_COLUMNS,
    ]


def build_master_table_columns(
    metadata_columns: Sequence[str],
    busco_lineages: Sequence[str] | None = None,
) -> list[str]:
    """Return the full master-table header from metadata and derived contracts."""
    cleaned_metadata_columns = [column.strip() for column in metadata_columns]
    if not cleaned_metadata_columns:
        raise ValueError("Metadata columns cannot be empty.")
    if any(not column for column in cleaned_metadata_columns):
        raise ValueError("Metadata columns cannot contain empty names.")
    if len(set(cleaned_metadata_columns)) != len(cleaned_metadata_columns):
        raise ValueError("Metadata columns must be unique.")

    append_columns = build_append_columns(busco_lineages)
    overlap = sorted(set(cleaned_metadata_columns) & set(append_columns))
    if overlap:
        raise ValueError(
            "Metadata columns overlap with derived master-table columns: "
            + ", ".join(overlap)
        )
    return [*cleaned_metadata_columns, *append_columns]


def read_append_columns_asset(path: Path = DEFAULT_APPEND_COLUMNS_ASSET) -> list[str]:
    """Read the repository append-column asset file."""
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def validate_default_append_columns_asset(
    path: Path = DEFAULT_APPEND_COLUMNS_ASSET,
) -> None:
    """Raise when the default asset file drifts from the code-defined contract."""
    asset_columns = read_append_columns_asset(path)
    expected_columns = build_append_columns(DEFAULT_BUSCO_LINEAGES)
    if asset_columns != expected_columns:
        raise ValueError(
            f"Append-column asset {path} does not match the default v1 contract."
        )
