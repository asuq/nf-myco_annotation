#!/usr/bin/env python3
"""Assemble the final master table from metadata and derived sample summaries."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Sequence

import master_table_contract


LOGGER = logging.getLogger(__name__)

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
        help="Path to the final master table CSV.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_key(key: str) -> str:
    """Normalise a column name for case-insensitive comparisons."""
    return "".join(character for character in key.lower() if character.isalnum())


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
) -> dict[str, dict[str, str]]:
    """Load an optional accession-keyed TSV with required columns."""
    if path is None:
        return {}
    header, rows = read_table(path, delimiter="\t")
    key_column = detect_key_column(header, ("accession", "Accession"))
    require_columns(header, required_columns, table_name)
    return build_index(rows, key_column, table_name)


def load_busco_index(paths: Sequence[Path]) -> dict[str, dict[str, str]]:
    """Load per-lineage BUSCO summaries into one accession-keyed mapping."""
    busco_index: dict[str, dict[str, str]] = {}
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
        for accession, row in table_index.items():
            busco_index.setdefault(accession, {})
            busco_index[accession][busco_column] = row[busco_column] or "NA"

    return busco_index


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


def build_master_row(
    sample_row: dict[str, str],
    metadata_header: Sequence[str],
    metadata_key_column: str,
    metadata_index: dict[str, dict[str, str]],
    taxonomy_index: dict[str, dict[str, str]],
    checkm2_index: dict[str, dict[str, str]],
    sixteen_s_index: dict[str, dict[str, str]],
    busco_index: dict[str, dict[str, str]],
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
        for column, value in busco_index[accession].items():
            if column in derived_values:
                derived_values[column] = value or "NA"

    for column in CRISPR_COLUMNS:
        if accession in ccfinder_index and column in ccfinder_index[accession]:
            derived_values[column] = ccfinder_index[accession][column] or "NA"

    for column in ANI_COLUMNS:
        if accession in ani_index and column in ani_index[accession]:
            derived_values[column] = ani_index[accession][column] or "NA"

    row.update(derived_values)
    return row


def write_csv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write the final master table as CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_build(args: argparse.Namespace) -> None:
    """Build the final master table from validated and derived inputs."""
    validated_samples = load_validated_samples(args.validated_samples)
    metadata_header, metadata_key_column, metadata_index = load_metadata(args.metadata)
    append_columns = load_append_columns(args.append_columns)

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
    )
    sixteen_s_index = load_optional_accession_index(
        args.sixteen_s_status,
        SIXTEEN_S_COLUMNS,
        "16S status table",
    )
    ccfinder_index = load_optional_accession_index(
        args.ccfinder_strains,
        CRISPR_COLUMNS,
        "CRISPR strain summary table",
    )
    ani_index = load_optional_accession_index(
        args.ani,
        ANI_COLUMNS,
        "ANI summary table",
    )
    busco_index = load_busco_index(args.busco)

    rows = [
        build_master_row(
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
        for sample_row in validated_samples
    ]

    output_header = [*metadata_header, *append_columns]
    if len(rows) != len({row[metadata_key_column] for row in rows}):
        raise MasterTableError("Final master table contains duplicate accession rows.")

    write_csv(args.output, output_header, rows)


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
