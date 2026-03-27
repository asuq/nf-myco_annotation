#!/usr/bin/env python3
"""Build a deterministic staged-genome manifest for ANI helpers."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


LOGGER = logging.getLogger(__name__)
MANIFEST_HEADER = ("accession", "internal_id", "staged_filename")


class BuildStagedManifestError(RuntimeError):
    """Raised when staged-manifest rows cannot be validated safely."""


@dataclass(frozen=True)
class StagedManifestRecord:
    """One validated staged-genome manifest row."""

    accession: str
    internal_id: str
    staged_filename: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a deterministic staged-genome manifest TSV."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to a headerless TSV of accession, internal_id, and staged_filename.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to write staged_genomes.tsv.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def parse_record(row: list[str], row_number: int) -> StagedManifestRecord:
    """Validate one raw TSV row and convert it to a record."""
    if len(row) != len(MANIFEST_HEADER):
        raise BuildStagedManifestError(
            "Staged manifest input row "
            f"{row_number} is malformed: expected {len(MANIFEST_HEADER)} tab-delimited "
            f"fields, found {len(row)}."
        )

    accession, internal_id, staged_filename = (field.strip() for field in row)
    if not accession:
        raise BuildStagedManifestError(
            f"Staged manifest input row {row_number} has an empty accession."
        )
    if not internal_id:
        raise BuildStagedManifestError(
            f"Staged manifest input row {row_number} has an empty internal_id."
        )
    if not staged_filename:
        raise BuildStagedManifestError(
            f"Staged manifest input row {row_number} has an empty staged_filename."
        )

    return StagedManifestRecord(
        accession=accession,
        internal_id=internal_id,
        staged_filename=staged_filename,
    )


def read_records(path: Path) -> list[StagedManifestRecord]:
    """Read, validate, and sort raw staged-manifest records."""
    if not path.is_file():
        raise BuildStagedManifestError(f"Missing staged manifest input rows: {path}")

    records: list[StagedManifestRecord] = []
    seen_accessions: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=1):
            if not row or all(not field.strip() for field in row):
                continue
            record = parse_record(row, row_number)
            if record.accession in seen_accessions:
                raise BuildStagedManifestError(
                    "Staged manifest input contains duplicate accession "
                    f"{record.accession!r}."
                )
            seen_accessions.add(record.accession)
            records.append(record)

    return sorted(records, key=lambda record: record.accession)


def write_manifest(path: Path, records: Sequence[StagedManifestRecord]) -> None:
    """Write the final staged-manifest TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(MANIFEST_HEADER)
        for record in records:
            writer.writerow(
                (record.accession, record.internal_id, record.staged_filename)
            )


def build_staged_manifest(input_path: Path, output_path: Path) -> None:
    """Build one deterministic staged-genome manifest TSV."""
    records = read_records(input_path)
    write_manifest(output_path, records)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the staged-manifest builder CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        build_staged_manifest(args.input, args.output)
    except BuildStagedManifestError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote staged manifest to %s.", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
