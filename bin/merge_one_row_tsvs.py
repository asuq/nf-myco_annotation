#!/usr/bin/env python3
"""Merge one-row TSV fragments deterministically against a canonical header."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Sequence


LOGGER = logging.getLogger(__name__)


class MergeOneRowTsvsError(RuntimeError):
    """Raised when one-row TSV fragments cannot be merged safely."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge one-row TSV fragments using one canonical header."
    )
    parser.add_argument(
        "--header",
        required=True,
        type=Path,
        help="Path to the canonical header TSV.",
    )
    parser.add_argument(
        "--input",
        action="append",
        dest="inputs",
        required=True,
        type=Path,
        help="Path to one one-row TSV fragment. May be supplied multiple times.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to the merged TSV output.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def read_canonical_header(path: Path) -> list[str]:
    """Read one canonical TSV header row."""
    if not path.is_file():
        raise MergeOneRowTsvsError(f"Missing canonical header file: {path}")
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) != 1:
        raise MergeOneRowTsvsError(
            f"Canonical header file must contain exactly one non-empty line: {path}"
        )
    return lines[0].split("\t")


def read_one_row_tsv(path: Path, canonical_header: Sequence[str]) -> dict[str, str]:
    """Read and validate one one-row TSV fragment."""
    if not path.is_file():
        raise MergeOneRowTsvsError(f"Missing TSV fragment: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise MergeOneRowTsvsError(f"TSV fragment is empty: {path}")
        if list(reader.fieldnames) != list(canonical_header):
            raise MergeOneRowTsvsError(
                f"TSV fragment header does not match the canonical header: {path}"
            )
        rows = list(reader)
    if len(rows) != 1:
        raise MergeOneRowTsvsError(
            f"TSV fragment must contain exactly one data row: {path}"
        )
    accession = rows[0].get("accession", "").strip()
    if not accession:
        raise MergeOneRowTsvsError(f"TSV fragment is missing accession data: {path}")
    return {column: rows[0].get(column, "") for column in canonical_header}


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write the merged TSV in canonical order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def merge_one_row_tsvs(
    header_path: Path,
    input_paths: Sequence[Path],
    output_path: Path,
) -> None:
    """Merge one-row TSV fragments sorted by accession."""
    canonical_header = read_canonical_header(header_path)
    merged_rows: list[dict[str, str]] = []
    seen_accessions: set[str] = set()

    for input_path in input_paths:
        row = read_one_row_tsv(input_path, canonical_header)
        accession = row["accession"].strip()
        if accession in seen_accessions:
            raise MergeOneRowTsvsError(
                f"Duplicate accession {accession!r} in TSV fragments."
            )
        seen_accessions.add(accession)
        merged_rows.append(row)

    merged_rows.sort(key=lambda row: row["accession"])
    write_tsv(output_path, canonical_header, merged_rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the deterministic TSV merger CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        merge_one_row_tsvs(
            header_path=args.header,
            input_paths=args.inputs,
            output_path=args.output,
        )
    except MergeOneRowTsvsError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote merged TSV to %s.", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
