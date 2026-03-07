#!/usr/bin/env python3
"""Placeholder CLI for accession sanitisation and collision handling."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Placeholder CLI for building original-to-internal ID mappings."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to a TSV with accession values to sanitise.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to the accession mapping TSV to write.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure placeholder logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> int:
    """Run the placeholder CLI."""
    _ = parse_args()
    configure_logging()
    logging.error("sanitise_accessions.py is a placeholder and is not implemented yet.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
