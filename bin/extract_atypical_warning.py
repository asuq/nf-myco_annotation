#!/usr/bin/env python3
"""Extract one sample's raw `Atypical_Warnings` value from a metadata table."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Sequence

from atypical_helpers import ATYPICAL_WARNINGS_COLUMN, is_missing_atypical_warning
from build_master_table import detect_key_column, find_column_by_normalised_name, read_table


LOGGER = logging.getLogger(__name__)
MISSING_VALUE = "NA"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract one accession's raw atypical-warning value from metadata."
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to the metadata table (CSV or TSV).",
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Original accession to look up in the metadata table.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def extract_atypical_warning(metadata_path: Path, accession: str) -> str:
    """Return the raw atypical-warning value for one accession or `NA`."""
    header, rows = read_table(metadata_path)
    accession_column = detect_key_column(header, ("Accession", "accession"))
    atypical_column = find_column_by_normalised_name(header, ATYPICAL_WARNINGS_COLUMN)
    if atypical_column is None:
        return MISSING_VALUE

    for row in rows:
        if row.get(accession_column, "").strip() != accession:
            continue
        value = row.get(atypical_column, "").strip()
        return MISSING_VALUE if is_missing_atypical_warning(value) else value
    return MISSING_VALUE


def main(argv: Sequence[str] | None = None) -> int:
    """Print one raw atypical-warning value for shell consumption."""
    args = parse_args(argv)
    configure_logging()
    try:
        print(extract_atypical_warning(args.metadata, args.accession))
    except Exception as error:
        LOGGER.error(str(error))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
