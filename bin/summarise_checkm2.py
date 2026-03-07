#!/usr/bin/env python3
"""Placeholder CLI for CheckM2 summarisation, gcode assignment, and QC calls."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Placeholder CLI for merging CheckM2 runs and assigning gcode."
    )
    parser.add_argument(
        "--gcode4-report",
        required=True,
        type=Path,
        help="Path to the CheckM2 quality report for translation table 4.",
    )
    parser.add_argument(
        "--gcode11-report",
        required=True,
        type=Path,
        help="Path to the CheckM2 quality report for translation table 11.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to the combined per-sample QC TSV.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure placeholder logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> int:
    """Run the placeholder CLI."""
    _ = parse_args()
    configure_logging()
    logging.error("summarise_checkm2.py is a placeholder and is not implemented yet.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
