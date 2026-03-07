#!/usr/bin/env python3
"""Assemble the final master table from metadata and derived sample summaries."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import master_table_contract


LOGGER = logging.getLogger(__name__)
ROOT_DIR = Path(__file__).resolve().parents[1]

VALIDATED_SAMPLE_REQUIRED_COLUMNS = (
    "accession",
    "is_new",
    "assembly_level",
    "genome_fasta",
    "internal_id",
)
TAXONOMY_COLUMNS = tuple(master_table_contract.TAXONOMY_COLUMNS)
CHECKM2_COLUMNS = tuple(master_table_contract.CHECKM2_COLUMNS) + ("Gcode", "Low_quality")
SIXTEEN_S_COLUMNS = ("16S",)
CRISPR_COLUMNS = tuple(master_table_contract.CRISPR_COLUMNS)
ANI_COLUMNS = tuple(master_table_contract.ANI_COLUMNS)
DEFAULT_SAMPLE_STATUS_COLUMNS_ASSET = ROOT_DIR / "assets" / "sample_status_columns.txt"
MISSING_VALUE_TOKENS = {"", "na", "n/a", "null", "none"}


@dataclass(frozen=True)
class BuscoLineageSummary:
    """Store one BUSCO lineage summary row for one accession."""

    value: str
    status: str
    warnings: str


class MasterTableError(RuntimeError):
    """Raised when the final master table cannot be built unambiguously."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build the final one-row-per-sample master table."
    )
    parser.add_argument(
        "--validated-samples",
        required=True,
        type=Path,
        help="Path to the validated sample manifest TSV.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to the metadata table (CSV or TSV).",
    )
    parser.add_argument(
        "--append-columns",
        type=Path,
        default=master_table_contract.DEFAULT_APPEND_COLUMNS_ASSET,
        help="Path to the ordered derived-column contract file.",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        help="Optional taxonomy-expansion TSV keyed by Tax_ID.",
    )
    parser.add_argument(
        "--checkm2",
        type=Path,
        help="Optional per-sample CheckM2 summary TSV.",
    )
    parser.add_argument(
        "--16s-status",
        dest="sixteen_s_status",
        type=Path,
        help="Optional per-sample 16S status TSV.",
    )
    parser.add_argument(
        "--busco",
        action="append",
        default=[],
        type=Path,
        help="Optional per-lineage BUSCO summary TSV. May be supplied multiple times.",
    )
    parser.add_argument(
        "--ccfinder-strains",
        type=Path,
        help="Optional CRISPRCasFinder strain summary TSV.",
    )
    parser.add_argument(
        "--ani",
        type=Path,
        help="Optional ANI cluster summary TSV.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to the final master_table.tsv output.",
    )
    parser.add_argument(
        "--sample-status-output",
        type=Path,
        help="Optional path to the final sample_status.tsv output.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_key(key: str) -> str:
    """Normalise a column name for case-insensitive comparisons."""
    return "".join(character for character in key.lower() if character.isalnum())


def is_missing(value: str | None) -> bool:
    """Return True when a scalar should be treated as missing."""
    if value is None:
        return True
    return value.strip().lower() in MISSING_VALUE_TOKENS


def split_tokens(value: str | None) -> list[str]:
    """Split a semicolon-delimited field into stripped non-empty tokens."""
    if value is None:
        return []
    return [token.strip() for token in value.split(";") if token.strip()]


def join_tokens(tokens: Sequence[str]) -> str:
    """Join tokens into a stable semicolon-delimited string without duplicates."""
    return ";".join(dict.fromkeys(token for token in tokens if token))


def join_notes(notes: Sequence[str]) -> str:
    """Join note strings into a stable semicolon-delimited sentence list."""
    return "; ".join(dict.fromkeys(note for note in notes if note))


def sniff_delimiter(path: Path) -> str:
    """Guess whether a table is comma- or tab-separated."""
    sample = path.read_text(encoding="utf-8").splitlines()
    sample_text = "\n".join(sample[:10])
    if not sample_text.strip():
        raise MasterTableError(f"Input table is empty: {path}")
    try:
        dialect = csv.Sniffer().sniff(sample_text, delimiters=",\t")
    except csv.Error:
        return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    return dialect.delimiter


def read_table(path: Path, delimiter: str | None = None) -> tuple[list[str], list[dict[str, str]]]:
    """Read a delimited table into a header and row dictionaries."""
    if not path.is_file():
        raise MasterTableError(f"Missing input table: {path}")
    table_delimiter = delimiter or sniff_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=table_delimiter)
        rows = [row for row in reader if any(cell.strip() for cell in row)]

    if not rows:
        raise MasterTableError(f"Input table is empty: {path}")

    header = [column.strip() for column in rows[0]]
    if any(not column for column in header):
        raise MasterTableError(f"Input table has an empty header name: {path}")
    if len(set(header)) != len(header):
        raise MasterTableError(f"Input table contains duplicate headers: {path}")

    records: list[dict[str, str]] = []
    for line_number, values in enumerate(rows[1:], start=2):
        if len(values) != len(header):
            raise MasterTableError(
                f"Input table {path} has {len(values)} columns on line {line_number}, "
                f"expected {len(header)}."
            )
        records.append(
            {column: value.strip() for column, value in zip(header, values, strict=True)}
        )
    return header, records


def detect_key_column(header: Sequence[str], aliases: Sequence[str]) -> str:
    """Detect a join key column from a list of aliases."""
    matching_columns = [
        column
        for column in header
        if normalise_key(column) in {normalise_key(alias) for alias in aliases}
    ]
    if len(matching_columns) > 1:
        raise MasterTableError(
            f"Input table contains multiple possible key columns: {matching_columns}"
        )
    if not matching_columns:
        raise MasterTableError(
            f"Input table is missing a required key column. Expected one of: {aliases}"
        )
    return matching_columns[0]


def build_index(
    rows: Sequence[dict[str, str]],
    key_column: str,
    table_name: str,
) -> dict[str, dict[str, str]]:
    """Build a unique keyed row index and fail on duplicates."""
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get(key_column, "").strip()
        if not key:
            raise MasterTableError(f"{table_name} contains an empty key value.")
        if key in indexed:
            raise MasterTableError(f"{table_name} contains duplicate key {key!r}.")
        indexed[key] = row
    return indexed


def ensure_index_keys_are_valid(
    index: dict[str, dict[str, str]],
    validated_accessions: set[str],
    table_name: str,
) -> None:
    """Raise when an accession-keyed table contains rows outside the validated set."""
    unexpected = sorted(set(index) - validated_accessions)
    if unexpected:
        raise MasterTableError(
            f"{table_name} contains accession values not present in validated_samples.tsv: "
            + ", ".join(unexpected)
        )


def require_columns(
    header: Sequence[str],
    columns: Sequence[str],
    table_name: str,
) -> None:
    """Ensure a table contains a required set of columns."""
    missing = [column for column in columns if column not in header]
    if missing:
        raise MasterTableError(
            f"{table_name} is missing required columns: {', '.join(missing)}"
        )


def load_validated_samples(path: Path) -> list[dict[str, str]]:
    """Load and validate the canonical validated sample table."""
    header, rows = read_table(path, delimiter="\t")
    require_columns(header, VALIDATED_SAMPLE_REQUIRED_COLUMNS, "validated_samples.tsv")
    samples = build_index(rows, "accession", "validated_samples.tsv")
    internal_ids = [row["internal_id"] for row in samples.values()]
    if len(set(internal_ids)) != len(internal_ids):
        raise MasterTableError("validated_samples.tsv contains duplicate internal_id values.")
    return list(samples.values())


def load_sample_status_columns(
    path: Path = DEFAULT_SAMPLE_STATUS_COLUMNS_ASSET,
) -> list[str]:
    """Load the maintained sample-status column order asset."""
    if not path.is_file():
        raise MasterTableError(f"Missing sample-status column asset: {path}")
    columns = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not columns:
        raise MasterTableError(f"Sample-status column asset is empty: {path}")
    if len(set(columns)) != len(columns):
        raise MasterTableError(f"Sample-status column asset contains duplicate columns: {path}")
    return columns


def load_metadata(path: Path) -> tuple[list[str], str, dict[str, dict[str, str]]]:
    """Load metadata and preserve the original header order."""
    header, rows = read_table(path)
    key_column = detect_key_column(header, ("Accession", "accession"))
    index = build_index(rows, key_column, "metadata table")
    return header, key_column, index


def build_supplemental_metadata_map(sample_row: dict[str, str]) -> dict[str, str]:
    """Map validated sample extras onto normalised metadata-style names."""
    supplemental: dict[str, str] = {}
    for column, value in sample_row.items():
        if column in VALIDATED_SAMPLE_REQUIRED_COLUMNS:
            continue
        if not value or value == "NA":
            continue
        supplemental[normalise_key(column)] = value
    return supplemental


def build_metadata_row(
    sample_row: dict[str, str],
    metadata_header: Sequence[str],
    metadata_key_column: str,
    metadata_index: dict[str, dict[str, str]],
) -> dict[str, str]:
    """Build the preserved metadata block for one sample."""
    accession = sample_row["accession"]
    existing_row = metadata_index.get(accession)
    if existing_row is not None:
        return {column: existing_row.get(column, "NA") or "NA" for column in metadata_header}

    metadata_row = {column: "NA" for column in metadata_header}
    metadata_row[metadata_key_column] = accession
    if sample_row.get("is_new") == "true":
        supplemental = build_supplemental_metadata_map(sample_row)
        for column in metadata_header:
            normalised = normalise_key(column)
            if normalised in supplemental:
                metadata_row[column] = supplemental[normalised]
    return metadata_row


def load_taxonomy_index(path: Path | None) -> dict[str, dict[str, str]]:
    """Load taxonomy rows keyed by Tax_ID."""
    if path is None:
        return {}
    header, rows = read_table(path, delimiter="\t")
    key_column = detect_key_column(header, ("Tax_ID", "tax_id"))
    require_columns(header, TAXONOMY_COLUMNS, "taxonomy table")
    return build_index(rows, key_column, "taxonomy table")


def load_optional_accession_index(
    path: Path | None,
    required_columns: Sequence[str],
    table_name: str,
    validated_accessions: set[str],
) -> dict[str, dict[str, str]]:
    """Load an optional accession-keyed TSV with required columns."""
    if path is None:
        return {}
    header, rows = read_table(path, delimiter="\t")
    key_column = detect_key_column(header, ("accession", "Accession"))
    require_columns(header, required_columns, table_name)
    index = build_index(rows, key_column, table_name)
    ensure_index_keys_are_valid(index, validated_accessions, table_name)
    return index


def load_busco_index(
    paths: Sequence[Path],
    validated_accessions: set[str],
) -> tuple[dict[str, dict[str, BuscoLineageSummary]], set[str]]:
    """Load per-lineage BUSCO summaries into one accession-keyed mapping."""
    busco_index: dict[str, dict[str, BuscoLineageSummary]] = {}
    seen_busco_columns: set[str] = set()

    for path in paths:
        header, rows = read_table(path, delimiter="\t")
        key_column = detect_key_column(header, ("accession", "Accession"))
        busco_columns = [column for column in header if column.startswith("BUSCO_")]
        if len(busco_columns) != 1:
            raise MasterTableError(
                f"BUSCO table must contain exactly one BUSCO_<lineage> column: {path}"
            )
        busco_column = busco_columns[0]
        if busco_column in seen_busco_columns:
            raise MasterTableError(f"Duplicate BUSCO lineage column supplied: {busco_column}")
        seen_busco_columns.add(busco_column)

        table_index = build_index(rows, key_column, f"BUSCO table {path.name}")
        ensure_index_keys_are_valid(
            table_index,
            validated_accessions,
            f"BUSCO table {path.name}",
        )
        for accession, row in table_index.items():
            busco_index.setdefault(accession, {})
            busco_value = row.get(busco_column, "") or "NA"
            busco_status = row.get("busco_status", "")
            busco_warnings = row.get("warnings", "")
            if not busco_status:
                busco_status = "failed" if busco_value == "NA" else "done"
            if not busco_warnings and busco_status == "failed":
                busco_warnings = "busco_summary_failed"
            busco_index[accession][busco_column] = BuscoLineageSummary(
                value=busco_value,
                status=busco_status,
                warnings=busco_warnings,
            )

    return busco_index, seen_busco_columns


def load_append_columns(path: Path) -> list[str]:
    """Load and validate the append-column contract file."""
    columns = master_table_contract.read_append_columns_asset(path)
    if not columns:
        raise MasterTableError(f"Append-column contract is empty: {path}")
    if len(set(columns)) != len(columns):
        raise MasterTableError(f"Append-column contract contains duplicate columns: {path}")
    return columns


def derive_taxonomy_values(
    metadata_row: dict[str, str],
    taxonomy_index: dict[str, dict[str, str]],
) -> dict[str, str]:
    """Resolve taxonomy columns by Tax_ID or return NA-filled values."""
    tax_id_column = next(
        (column for column in metadata_row if normalise_key(column) == "taxid"),
        None,
    )
    blank_taxonomy = {column: "NA" for column in TAXONOMY_COLUMNS}
    if tax_id_column is None:
        return blank_taxonomy
    tax_id = metadata_row.get(tax_id_column, "NA")
    if not tax_id or tax_id == "NA":
        return blank_taxonomy
    taxonomy_row = taxonomy_index.get(tax_id)
    if taxonomy_row is None:
        return blank_taxonomy
    return {column: taxonomy_row.get(column, "NA") or "NA" for column in TAXONOMY_COLUMNS}


def find_column_by_normalised_name(
    columns: Sequence[str],
    target: str,
) -> str | None:
    """Find a column by case- and punctuation-insensitive name."""
    normalised_target = normalise_key(target)
    for column in columns:
        if normalise_key(column) == normalised_target:
            return column
    return None


def derive_taxonomy_status(
    metadata_row: dict[str, str],
    taxonomy_index: dict[str, dict[str, str]],
    taxonomy_requested: bool,
) -> tuple[str, list[str]]:
    """Return taxonomy status and any warning tokens for one sample."""
    if not taxonomy_requested:
        return "na", []
    tax_id_column = find_column_by_normalised_name(tuple(metadata_row), "tax_id")
    if tax_id_column is None:
        return "na", []
    tax_id = metadata_row.get(tax_id_column, "NA")
    if is_missing(tax_id):
        return "na", []
    if tax_id in taxonomy_index:
        return "done", []
    return "failed", ["taxonomy_missing"]


def has_non_na_values(row: dict[str, str], columns: Sequence[str]) -> bool:
    """Return True when any requested column contains a non-NA value."""
    return any(not is_missing(row.get(column)) for column in columns)


def derive_checkm2_statuses(
    accession: str,
    checkm2_index: dict[str, dict[str, str]],
    checkm2_requested: bool,
) -> tuple[str, str, str, list[str]]:
    """Return CheckM2/gcode status columns and warning tokens."""
    if not checkm2_requested:
        return "na", "na", "na", []
    row = checkm2_index.get(accession)
    if row is None:
        return "failed", "failed", "failed", ["missing_checkm2_summary"]

    warnings = split_tokens(row.get("warnings", ""))
    gcode4_columns = [column for column in CHECKM2_COLUMNS if column.endswith("_gcode4")]
    gcode11_columns = [column for column in CHECKM2_COLUMNS if column.endswith("_gcode11")]
    has_gcode4 = has_non_na_values(row, gcode4_columns)
    has_gcode11 = has_non_na_values(row, gcode11_columns)

    gcode4_status = "failed" if "checkm2_gcode4_failed" in warnings or not has_gcode4 else "done"
    gcode11_status = (
        "failed" if "checkm2_gcode11_failed" in warnings or not has_gcode11 else "done"
    )

    gcode_value = row.get("Gcode", "NA") or "NA"
    gcode_status = "done" if gcode_value in {"4", "11"} else "failed"
    return gcode4_status, gcode11_status, gcode_status, warnings


def derive_barrnap_status(
    accession: str,
    sixteen_s_index: dict[str, dict[str, str]],
    sixteen_s_requested: bool,
) -> tuple[str, list[str]]:
    """Return Barrnap status and warning tokens for one sample."""
    if not sixteen_s_requested:
        return "na", []
    row = sixteen_s_index.get(accession)
    if row is None:
        return "failed", ["missing_16s_summary"]
    warnings = split_tokens(row.get("warnings", ""))
    status_value = row.get("16S", "NA") or "NA"
    if status_value == "NA":
        return "failed", warnings or ["invalid_barrnap_output"]
    return "done", warnings


def derive_busco_statuses(
    accession: str,
    busco_index: dict[str, dict[str, BuscoLineageSummary]],
    provided_busco_columns: set[str],
    sample_status_columns: Sequence[str],
) -> tuple[dict[str, str], list[str]]:
    """Return BUSCO lineage status columns and warning tokens for one sample."""
    status_values = {
        column: "na"
        for column in sample_status_columns
        if column.startswith("busco_") and column.endswith("_status")
    }
    warnings: list[str] = []
    sample_busco = busco_index.get(accession, {})

    for status_column in status_values:
        lineage = status_column.removeprefix("busco_").removesuffix("_status")
        busco_column = f"BUSCO_{lineage}"
        if busco_column not in provided_busco_columns:
            continue
        summary = sample_busco.get(busco_column)
        if summary is None:
            status_values[status_column] = "failed"
            warnings.append(f"missing_{busco_column.lower()}_summary")
            continue
        status_values[status_column] = summary.status if summary.status else "done"
        warnings.extend(split_tokens(summary.warnings))
    return status_values, warnings


def derive_ccfinder_status(
    accession: str,
    ccfinder_index: dict[str, dict[str, str]],
    ccfinder_requested: bool,
    gcode_value: str,
) -> tuple[str, list[str]]:
    """Return CRISPRCasFinder status and warning tokens for one sample."""
    row = ccfinder_index.get(accession)
    if row is not None:
        warnings = split_tokens(row.get("warnings", ""))
        status = row.get("ccfinder_status", "") or "done"
        return status, warnings
    if gcode_value == "NA":
        return ("skipped", []) if ccfinder_requested else ("na", [])
    if not ccfinder_requested:
        return "na", []
    return "failed", ["missing_ccfinder_summary"]


def detect_atypical_flags(metadata_row: dict[str, str]) -> tuple[bool, bool]:
    """Return the atypical and exception flags from the preserved metadata block."""
    atypical_column = find_column_by_normalised_name(tuple(metadata_row), "Atypical_Warnings")
    if atypical_column is None:
        return False, False
    atypical_value = metadata_row.get(atypical_column, "")
    if is_missing(atypical_value):
        return False, False
    lowered = atypical_value.casefold()
    return True, "unverified source organism" in lowered


def derive_ani_decision(
    accession: str,
    master_row: dict[str, str],
    metadata_row: dict[str, str],
    ani_index: dict[str, dict[str, str]],
    ani_requested: bool,
    primary_busco_column: str | None,
) -> tuple[str, str]:
    """Return ANI inclusion status and exclusion reasons for one sample."""
    if not ani_requested:
        return "na", ""

    exclusion_reasons: list[str] = []
    gcode_value = master_row.get("Gcode", "NA") or "NA"
    low_quality_value = master_row.get("Low_quality", "NA") or "NA"
    sixteen_s_value = master_row.get("16S", "NA") or "NA"
    is_atypical, is_exception = detect_atypical_flags(metadata_row)

    if gcode_value == "NA":
        exclusion_reasons.append("gcode_na")
    if low_quality_value == "true":
        exclusion_reasons.append("low_quality")
    elif low_quality_value == "NA" and gcode_value in {"4", "11"}:
        exclusion_reasons.append("missing_low_quality")

    if sixteen_s_value == "partial":
        exclusion_reasons.append("partial_16s")
    elif sixteen_s_value == "No":
        exclusion_reasons.append("no_16s")
    elif sixteen_s_value == "NA":
        exclusion_reasons.append("16s_na")

    if is_atypical and not is_exception:
        exclusion_reasons.append("atypical")

    if primary_busco_column is not None and is_missing(master_row.get(primary_busco_column)):
        exclusion_reasons.append("missing_primary_busco")

    if accession in ani_index:
        if exclusion_reasons:
            raise MasterTableError(
                f"ANI summary contains ineligible accession {accession!r}: "
                + ";".join(exclusion_reasons)
            )
        return "true", ""

    if exclusion_reasons:
        return "false", join_tokens(exclusion_reasons)

    raise MasterTableError(
        f"ANI summary is missing eligible accession {accession!r}."
    )


def build_sample_status_row(
    sample_row: dict[str, str],
    metadata_present: bool,
    metadata_row: dict[str, str],
    master_row: dict[str, str],
    taxonomy_index: dict[str, dict[str, str]],
    taxonomy_requested: bool,
    checkm2_index: dict[str, dict[str, str]],
    checkm2_requested: bool,
    sixteen_s_index: dict[str, dict[str, str]],
    sixteen_s_requested: bool,
    busco_index: dict[str, dict[str, BuscoLineageSummary]],
    provided_busco_columns: set[str],
    ccfinder_index: dict[str, dict[str, str]],
    ccfinder_requested: bool,
    ani_index: dict[str, dict[str, str]],
    ani_requested: bool,
    sample_status_columns: Sequence[str],
) -> dict[str, str]:
    """Build one authoritative sample-status row from the keyed final joins."""
    row: dict[str, str] = {}
    for column in sample_status_columns:
        if column.endswith("_status") or column == "ani_included":
            row[column] = "na"
        elif column in {"gcode", "low_quality"}:
            row[column] = "NA"
        else:
            row[column] = ""

    accession = sample_row["accession"]
    row["accession"] = accession
    row["internal_id"] = sample_row["internal_id"]
    row["is_new"] = sample_row["is_new"]
    row["validation_status"] = "done"
    row["gcode"] = master_row.get("Gcode", "NA") or "NA"
    row["low_quality"] = master_row.get("Low_quality", "NA") or "NA"

    warnings: list[str] = []
    notes: list[str] = []

    taxonomy_status, taxonomy_warnings = derive_taxonomy_status(
        metadata_row=metadata_row,
        taxonomy_index=taxonomy_index,
        taxonomy_requested=taxonomy_requested,
    )
    row["taxonomy_status"] = taxonomy_status
    warnings.extend(taxonomy_warnings)

    barrnap_status, barrnap_warnings = derive_barrnap_status(
        accession=accession,
        sixteen_s_index=sixteen_s_index,
        sixteen_s_requested=sixteen_s_requested,
    )
    row["barrnap_status"] = barrnap_status
    warnings.extend(barrnap_warnings)

    (
        row["checkm2_gcode4_status"],
        row["checkm2_gcode11_status"],
        row["gcode_status"],
        checkm2_warnings,
    ) = derive_checkm2_statuses(
        accession=accession,
        checkm2_index=checkm2_index,
        checkm2_requested=checkm2_requested,
    )
    warnings.extend(checkm2_warnings)

    busco_status_values, busco_warnings = derive_busco_statuses(
        accession=accession,
        busco_index=busco_index,
        provided_busco_columns=provided_busco_columns,
        sample_status_columns=sample_status_columns,
    )
    row.update(busco_status_values)
    warnings.extend(busco_warnings)

    ccfinder_status, ccfinder_warnings = derive_ccfinder_status(
        accession=accession,
        ccfinder_index=ccfinder_index,
        ccfinder_requested=ccfinder_requested,
        gcode_value=row["gcode"],
    )
    row["ccfinder_status"] = ccfinder_status
    warnings.extend(ccfinder_warnings)

    if row["gcode"] == "NA":
        row["prokka_status"] = "skipped" if checkm2_requested else "na"
        row["padloc_status"] = "skipped" if checkm2_requested else "na"
        row["eggnog_status"] = "skipped" if checkm2_requested else "na"
    else:
        row["prokka_status"] = "na"
        row["padloc_status"] = "na"
        row["eggnog_status"] = "na"

    primary_busco_column = next(
        (column for column in master_row if column.startswith("BUSCO_")),
        None,
    )
    row["ani_included"], row["ani_exclusion_reason"] = derive_ani_decision(
        accession=accession,
        master_row=master_row,
        metadata_row=metadata_row,
        ani_index=ani_index,
        ani_requested=ani_requested,
        primary_busco_column=primary_busco_column,
    )

    if sample_row.get("is_new") == "true" and not metadata_present:
        warnings.append("missing_metadata_for_new_sample")
        notes.append("New sample metadata missing; unavailable fields filled with NA.")

    row["warnings"] = join_tokens(warnings)
    row["notes"] = join_notes(notes)
    return row


def build_master_row(
    sample_row: dict[str, str],
    metadata_header: Sequence[str],
    metadata_key_column: str,
    metadata_index: dict[str, dict[str, str]],
    taxonomy_index: dict[str, dict[str, str]],
    checkm2_index: dict[str, dict[str, str]],
    sixteen_s_index: dict[str, dict[str, str]],
    busco_index: dict[str, dict[str, BuscoLineageSummary]],
    ccfinder_index: dict[str, dict[str, str]],
    ani_index: dict[str, dict[str, str]],
    append_columns: Sequence[str],
) -> dict[str, str]:
    """Build one final master-table row for a validated sample."""
    metadata_row = build_metadata_row(
        sample_row,
        metadata_header,
        metadata_key_column,
        metadata_index,
    )
    row = {column: metadata_row[column] for column in metadata_header}

    derived_values = {column: "NA" for column in append_columns}
    derived_values.update(derive_taxonomy_values(metadata_row, taxonomy_index))

    accession = sample_row["accession"]
    for column in CHECKM2_COLUMNS:
        if accession in checkm2_index and column in checkm2_index[accession]:
            derived_values[column] = checkm2_index[accession][column] or "NA"

    if accession in sixteen_s_index and "16S" in sixteen_s_index[accession]:
        derived_values["16S"] = sixteen_s_index[accession]["16S"] or "NA"

    if accession in busco_index:
        for column, summary in busco_index[accession].items():
            if column in derived_values:
                derived_values[column] = summary.value or "NA"

    for column in CRISPR_COLUMNS:
        if accession in ccfinder_index and column in ccfinder_index[accession]:
            derived_values[column] = ccfinder_index[accession][column] or "NA"

    for column in ANI_COLUMNS:
        if accession in ani_index and column in ani_index[accession]:
            derived_values[column] = ani_index[accession][column] or "NA"

    row.update(derived_values)
    return row


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write a TSV with a fixed header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_build(args: argparse.Namespace) -> None:
    """Build the final master table from validated and derived inputs."""
    validated_samples = load_validated_samples(args.validated_samples)
    validated_accessions = {row["accession"] for row in validated_samples}
    metadata_header, metadata_key_column, metadata_index = load_metadata(args.metadata)
    append_columns = load_append_columns(args.append_columns)
    sample_status_columns = load_sample_status_columns()

    expected_header = master_table_contract.build_master_table_columns(metadata_header)
    if append_columns != expected_header[len(metadata_header) :]:
        raise MasterTableError(
            "Append-column contract does not match the default master-table column contract."
        )

    taxonomy_index = load_taxonomy_index(args.taxonomy)
    checkm2_index = load_optional_accession_index(
        args.checkm2,
        CHECKM2_COLUMNS,
        "CheckM2 summary table",
        validated_accessions,
    )
    sixteen_s_index = load_optional_accession_index(
        args.sixteen_s_status,
        SIXTEEN_S_COLUMNS,
        "16S status table",
        validated_accessions,
    )
    ccfinder_index = load_optional_accession_index(
        args.ccfinder_strains,
        CRISPR_COLUMNS,
        "CRISPR strain summary table",
        validated_accessions,
    )
    ani_index = load_optional_accession_index(
        args.ani,
        ANI_COLUMNS,
        "ANI summary table",
        validated_accessions,
    )
    busco_index, provided_busco_columns = load_busco_index(args.busco, validated_accessions)

    rows: list[dict[str, str]] = []
    sample_status_rows: list[dict[str, str]] = []
    for sample_row in validated_samples:
        metadata_present = sample_row["accession"] in metadata_index
        metadata_row = build_metadata_row(
            sample_row=sample_row,
            metadata_header=metadata_header,
            metadata_key_column=metadata_key_column,
            metadata_index=metadata_index,
        )
        master_row = build_master_row(
            sample_row=sample_row,
            metadata_header=metadata_header,
            metadata_key_column=metadata_key_column,
            metadata_index=metadata_index,
            taxonomy_index=taxonomy_index,
            checkm2_index=checkm2_index,
            sixteen_s_index=sixteen_s_index,
            busco_index=busco_index,
            ccfinder_index=ccfinder_index,
            ani_index=ani_index,
            append_columns=append_columns,
        )
        rows.append(master_row)
        sample_status_rows.append(
            build_sample_status_row(
                sample_row=sample_row,
                metadata_present=metadata_present,
                metadata_row=metadata_row,
                master_row=master_row,
                taxonomy_index=taxonomy_index,
                taxonomy_requested=args.taxonomy is not None,
                checkm2_index=checkm2_index,
                checkm2_requested=args.checkm2 is not None,
                sixteen_s_index=sixteen_s_index,
                sixteen_s_requested=args.sixteen_s_status is not None,
                busco_index=busco_index,
                provided_busco_columns=provided_busco_columns,
                ccfinder_index=ccfinder_index,
                ccfinder_requested=args.ccfinder_strains is not None,
                ani_index=ani_index,
                ani_requested=args.ani is not None,
                sample_status_columns=sample_status_columns,
            )
        )

    output_header = [*metadata_header, *append_columns]
    if len(rows) != len({row[metadata_key_column] for row in rows}):
        raise MasterTableError("Final master table contains duplicate accession rows.")

    write_tsv(args.output, output_header, rows)
    sample_status_output = args.sample_status_output or args.output.with_name("sample_status.tsv")
    write_tsv(sample_status_output, sample_status_columns, sample_status_rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the master-table builder CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_build(args)
    except MasterTableError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote master table to %s.", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
