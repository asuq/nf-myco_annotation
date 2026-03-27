#!/usr/bin/env python3
"""Write runtime BUSCO-aware reporting contracts for the workflow."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Sequence

import master_table_contract


LOGGER = logging.getLogger(__name__)


class OutputContractError(RuntimeError):
    """Raised when runtime output contracts cannot be built safely."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build runtime master-table and sample-status column contracts."
    )
    parser.add_argument(
        "--busco-lineage",
        action="append",
        required=True,
        help="Configured BUSCO lineage name. May be supplied multiple times.",
    )
    parser.add_argument(
        "--append-columns-output",
        required=True,
        type=Path,
        help="Output path for master_table_append_columns.txt.",
    )
    parser.add_argument(
        "--sample-status-columns-output",
        required=True,
        type=Path,
        help="Output path for sample_status_columns.txt.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def write_columns(path: Path, columns: Sequence[str]) -> None:
    """Write one column name per line."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(columns) + "\n", encoding="utf-8")


def run_build(args: argparse.Namespace) -> None:
    """Build and write the requested output contracts."""
    try:
        busco_lineages = master_table_contract.normalise_busco_lineages(args.busco_lineage)
    except ValueError as error:
        raise OutputContractError(str(error)) from error

    append_columns = master_table_contract.build_append_columns(busco_lineages)
    sample_status_columns = master_table_contract.build_sample_status_columns(busco_lineages)
    write_columns(args.append_columns_output, append_columns)
    write_columns(args.sample_status_columns_output, sample_status_columns)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the output-contract builder CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_build(args)
    except OutputContractError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info(
        "Wrote runtime output contracts to %s and %s.",
        args.append_columns_output,
        args.sample_status_columns_output,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
