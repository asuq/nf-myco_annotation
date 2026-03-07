#!/usr/bin/env python3
"""Placeholder CLI for validating sample and metadata inputs."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Placeholder CLI for validating the pipeline input tables."
    )
    parser.add_argument(
        "--sample-csv",
        required=True,
        type=Path,
        help="Path to the sample manifest CSV.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to the metadata table.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for validated manifest outputs.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure placeholder logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> int:
    """Run the placeholder CLI."""
    _ = parse_args()
    configure_logging()
    logging.error("validate_inputs.py is a placeholder and is not implemented yet.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
