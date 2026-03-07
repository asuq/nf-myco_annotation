#!/usr/bin/env python3
"""Concatenate per-sample best 16S FASTA files into a cohort FASTA."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Sequence


LOGGER = logging.getLogger(__name__)

REQUIRED_INPUT_COLUMNS = ("accession", "status_tsv", "best_16s_fasta")
REQUIRED_STATUS_COLUMNS = (
    "accession",
    "16S",
    "best_16S_header",
    "best_16S_length",
    "include_in_all_best_16S",
    "warnings",
)


class Cohort16SError(RuntimeError):
    """Raised when the cohort 16S inputs are malformed."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Concatenate cohort-eligible best 16S FASTA files."
    )
    parser.add_argument(
        "--inputs",
        required=True,
        type=Path,
        help="TSV listing accession, status_tsv, and best_16s_fasta paths.",
    )
    parser.add_argument(
        "--output-fasta",
        required=True,
        type=Path,
        help="Path to the concatenated all_best_16S.fna file.",
    )
    parser.add_argument(
        "--output-manifest",
        required=True,
        type=Path,
        help="Path to the cohort manifest TSV.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header list and row dictionaries."""
    if not path.is_file():
        raise Cohort16SError(f"Missing cohort 16S input file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        raise Cohort16SError(f"Empty TSV file: {path}")

    header = [value.strip() for value in rows[0]]
    data_rows: list[dict[str, str]] = []
    for line_number, values in enumerate(rows[1:], start=2):
        if len(values) != len(header):
            raise Cohort16SError(
                f"TSV file {path} has {len(values)} columns on line {line_number}, "
                f"expected {len(header)}."
            )
        data_rows.append(
            {column: value.strip() for column, value in zip(header, values, strict=True)}
        )
    return header, data_rows


def require_columns(
    header: Sequence[str],
    required_columns: Sequence[str],
    path: Path,
) -> None:
    """Raise when a TSV header is missing required columns."""
    missing = [column for column in required_columns if column not in header]
    if missing:
        raise Cohort16SError(
            f"TSV file {path} is missing required columns: {', '.join(missing)}"
        )


def load_cohort_rows(path: Path) -> list[dict[str, str]]:
    """Load and validate the cohort-input manifest rows."""
    header, rows = read_tsv(path)
    require_columns(header, REQUIRED_INPUT_COLUMNS, path)
    return rows


def load_status_row(path: Path) -> dict[str, str]:
    """Load and validate a single per-sample 16S status row."""
    header, rows = read_tsv(path)
    require_columns(header, REQUIRED_STATUS_COLUMNS, path)
    if len(rows) != 1:
        raise Cohort16SError(f"Status TSV must contain exactly one data row: {path}")
    return rows[0]


def append_fasta_content(output_handle, fasta_path: Path) -> None:
    """Append a non-empty FASTA file to the cohort output handle."""
    if not fasta_path.is_file():
        raise Cohort16SError(f"Missing best 16S FASTA file: {fasta_path}")
    content = fasta_path.read_text(encoding="utf-8")
    if not content.strip():
        raise Cohort16SError(
            f"Included cohort FASTA file is empty but marked for inclusion: {fasta_path}"
        )
    output_handle.write(content if content.endswith("\n") else content + "\n")


def write_manifest(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write the cohort manifest TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(REQUIRED_STATUS_COLUMNS),
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in REQUIRED_STATUS_COLUMNS})


def build_cohort_fasta(
    inputs_path: Path,
    output_fasta: Path,
    output_manifest: Path,
) -> None:
    """Build the cohort FASTA and manifest from per-sample 16S outputs."""
    manifest_rows: list[dict[str, str]] = []
    cohort_rows = load_cohort_rows(inputs_path)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as fasta_handle:
        for cohort_row in cohort_rows:
            status_path = Path(cohort_row["status_tsv"]).expanduser()
            best_fasta_path = Path(cohort_row["best_16s_fasta"]).expanduser()
            status_row = load_status_row(status_path)
            if status_row["accession"] != cohort_row["accession"]:
                raise Cohort16SError(
                    "Cohort input accession does not match the status TSV accession "
                    f"for {status_path}."
                )
            if status_row["include_in_all_best_16S"] != "true":
                continue
            append_fasta_content(fasta_handle, best_fasta_path)
            manifest_rows.append(status_row)

    write_manifest(output_manifest, manifest_rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the cohort 16S concatenation CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        build_cohort_fasta(
            inputs_path=args.inputs,
            output_fasta=args.output_fasta,
            output_manifest=args.output_manifest,
        )
    except Cohort16SError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote cohort 16S FASTA to %s.", args.output_fasta)
    return 0


if __name__ == "__main__":
    sys.exit(main())
