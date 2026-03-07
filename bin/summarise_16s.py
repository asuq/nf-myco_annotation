#!/usr/bin/env python3
"""Placeholder CLI for deriving 16S status tables and FASTA outputs."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Placeholder CLI for Barrnap-derived 16S summarisation."
    )
    parser.add_argument(
        "--rrna-gff",
        required=True,
        type=Path,
        help="Path to the Barrnap GFF file.",
    )
    parser.add_argument(
        "--rrna-fasta",
        required=True,
        type=Path,
        help="Path to the Barrnap FASTA file.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for 16S summary outputs.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure placeholder logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> int:
    """Run the placeholder CLI."""
    _ = parse_args()
    configure_logging()
    logging.error("summarise_16s.py is a placeholder and is not implemented yet.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
