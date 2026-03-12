#!/usr/bin/env python3
"""Patch the PADLOC launcher to use a task-local bootstrap data directory."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


BOOTSTRAP_DATA_MKDIR = 'mkdir -p "${SRC_DIR}/../data"'
BOOTSTRAP_DATA_MKDIR_REPLACEMENT = 'mkdir -p "${PADLOC_BOOTSTRAP_DATA}"'
DATA_INITIALISER = 'DATA=$(normpath "${SRC_DIR}/../data")'
DATA_INITIALISER_REPLACEMENT = 'DATA=$(normpath "${PADLOC_BOOTSTRAP_DATA}")'


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Patch the PADLOC launcher bootstrap data path.",
    )
    parser.add_argument(
        "launcher",
        type=Path,
        help="Path to the copied PADLOC launcher script.",
    )
    return parser


def patch_launcher(launcher_path: Path) -> None:
    """Rewrite the launcher bootstrap path initialisation in place."""
    launcher_text = launcher_path.read_text(encoding="utf-8")

    if BOOTSTRAP_DATA_MKDIR not in launcher_text:
        raise RuntimeError(
            f"Missing PADLOC bootstrap mkdir line in {launcher_path}."
        )
    if DATA_INITIALISER not in launcher_text:
        raise RuntimeError(
            f"Missing PADLOC data initialiser line in {launcher_path}."
        )

    launcher_text = launcher_text.replace(
        BOOTSTRAP_DATA_MKDIR,
        BOOTSTRAP_DATA_MKDIR_REPLACEMENT,
    )
    launcher_text = launcher_text.replace(
        DATA_INITIALISER,
        DATA_INITIALISER_REPLACEMENT,
    )
    launcher_path.write_text(launcher_text, encoding="utf-8")


def main() -> int:
    """Run the launcher patch command."""
    args = build_parser().parse_args()
    patch_launcher(args.launcher)
    return 0


if __name__ == "__main__":
    sys.exit(main())
