#!/usr/bin/env python3
"""Collect tool, container, and resource provenance into a final TSV report."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Sequence


LOGGER = logging.getLogger(__name__)
OUTPUT_COLUMNS = ("category", "name", "value", "source")
VERSION_CATEGORY_KEYS = {
    "script": "script",
    "transform": "helper",
}


class CollectVersionsError(RuntimeError):
    """Raised when version collection inputs are invalid."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Collect tool, container, and database provenance into a TSV."
    )
    parser.add_argument(
        "--version-file",
        action="append",
        default=[],
        type=Path,
        help="Path to one process versions.yml file. May be supplied multiple times.",
    )
    parser.add_argument(
        "--container-ref",
        action="append",
        default=[],
        help="Container reference as name=value.",
    )
    parser.add_argument(
        "--busco-lineage",
        action="append",
        default=[],
        help="Configured BUSCO lineage name. May be supplied multiple times.",
    )
    parser.add_argument(
        "--nextflow-version",
        default="NA",
        help="Nextflow version string.",
    )
    parser.add_argument(
        "--pipeline-version",
        default="NA",
        help="Pipeline manifest version string.",
    )
    parser.add_argument(
        "--git-commit",
        default="NA",
        help="Pipeline git commit hash or NA.",
    )
    parser.add_argument(
        "--container-engine",
        default="NA",
        help="Container engine name or NA.",
    )
    parser.add_argument(
        "--use-biocontainers",
        default="NA",
        help="Whether BioContainers are preferred for the run.",
    )
    parser.add_argument(
        "--checkm2-db",
        default="NA",
        help="CheckM2 database path.",
    )
    parser.add_argument(
        "--checkm2-db-label",
        default="NA",
        help="Optional CheckM2 database version or label.",
    )
    parser.add_argument(
        "--taxdump",
        default="NA",
        help="Pinned taxdump path.",
    )
    parser.add_argument(
        "--taxdump-label",
        default="NA",
        help="Optional pinned taxdump label.",
    )
    parser.add_argument(
        "--busco-download-dir",
        default="NA",
        help="BUSCO dataset root path.",
    )
    parser.add_argument(
        "--eggnog-db",
        default="NA",
        help="eggNOG database path.",
    )
    parser.add_argument(
        "--eggnog-db-label",
        default="NA",
        help="Optional eggNOG database label.",
    )
    parser.add_argument(
        "--padloc-db-label",
        default="NA",
        help="Optional PADLOC database label.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to tool_and_db_versions.tsv.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_value(value: str | None) -> str:
    """Normalise blank values to NA."""
    if value is None:
        return "NA"
    stripped = value.strip()
    return stripped if stripped else "NA"


def parse_name_value(token: str, argument_name: str) -> tuple[str, str]:
    """Parse one name=value token."""
    if "=" not in token:
        raise CollectVersionsError(
            f"{argument_name} entry must use name=value syntax: {token!r}"
        )
    name, value = token.split("=", 1)
    name = name.strip()
    if not name:
        raise CollectVersionsError(
            f"{argument_name} entry is missing the name before '=': {token!r}"
        )
    return name, normalise_value(value)


def parse_versions_file(path: Path) -> list[dict[str, str]]:
    """Parse one simple Nextflow versions.yml file into row dictionaries."""
    if not path.is_file():
        raise CollectVersionsError(f"Missing versions file: {path}")

    rows: list[dict[str, str]] = []
    current_source = "unknown"
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if not line.startswith((" ", "\t")):
                if not stripped.endswith(":"):
                    raise CollectVersionsError(
                        f"Malformed versions.yml header line in {path}: {stripped!r}"
                    )
                current_source = stripped[:-1].strip().strip("'\"")
                continue
            if ":" not in stripped:
                raise CollectVersionsError(
                    f"Malformed versions.yml entry in {path}: {stripped!r}"
                )
            name, value = stripped.split(":", 1)
            clean_name = name.strip()
            clean_value = value.strip().strip("'\"")
            category = VERSION_CATEGORY_KEYS.get(clean_name, "tool_version")
            rows.append(
                {
                    "category": category,
                    "name": clean_name,
                    "value": normalise_value(clean_value),
                    "source": current_source,
                }
            )
    return rows


def build_runtime_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    """Build runtime and resource report rows from CLI arguments."""
    rows = [
        {"category": "runtime", "name": "nextflow_version", "value": normalise_value(args.nextflow_version), "source": "workflow"},
        {"category": "runtime", "name": "pipeline_version", "value": normalise_value(args.pipeline_version), "source": "workflow"},
        {"category": "runtime", "name": "git_commit", "value": normalise_value(args.git_commit), "source": "workflow"},
        {"category": "runtime", "name": "container_engine", "value": normalise_value(args.container_engine), "source": "workflow"},
        {"category": "config", "name": "use_biocontainers", "value": normalise_value(args.use_biocontainers), "source": "params"},
        {"category": "resource", "name": "checkm2_db", "value": normalise_value(args.checkm2_db), "source": "params"},
        {"category": "resource", "name": "checkm2_db_label", "value": normalise_value(args.checkm2_db_label), "source": "params"},
        {"category": "resource", "name": "taxdump", "value": normalise_value(args.taxdump), "source": "params"},
        {"category": "resource", "name": "taxdump_label", "value": normalise_value(args.taxdump_label), "source": "params"},
        {"category": "resource", "name": "busco_download_dir", "value": normalise_value(args.busco_download_dir), "source": "params"},
        {"category": "resource", "name": "eggnog_db", "value": normalise_value(args.eggnog_db), "source": "params"},
        {"category": "resource", "name": "eggnog_db_label", "value": normalise_value(args.eggnog_db_label), "source": "params"},
        {"category": "resource", "name": "padloc_db_label", "value": normalise_value(args.padloc_db_label), "source": "params"},
        {
            "category": "resource",
            "name": "busco_lineages",
            "value": ";".join(args.busco_lineage) if args.busco_lineage else "NA",
            "source": "params",
        },
    ]
    for lineage in args.busco_lineage:
        rows.append(
            {
                "category": "resource",
                "name": f"busco_dataset_{lineage}",
                "value": (
                    str(Path(args.busco_download_dir) / lineage)
                    if normalise_value(args.busco_download_dir) != "NA"
                    else "NA"
                ),
                "source": "params",
            }
        )
    for token in args.container_ref:
        name, value = parse_name_value(token, "--container-ref")
        rows.append(
            {
                "category": "container",
                "name": name,
                "value": value,
                "source": "params",
            }
        )
    return rows


def deduplicate_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    """Deduplicate rows while keeping a stable, category-aware order."""
    seen: set[tuple[str, str, str, str]] = set()
    ordered_rows: list[dict[str, str]] = []
    category_order = {
        "runtime": 0,
        "config": 1,
        "resource": 2,
        "container": 3,
        "tool_version": 4,
        "script": 5,
        "helper": 6,
    }
    for row in sorted(
        rows,
        key=lambda row: (
            category_order.get(row["category"], 99),
            row["name"],
            row["source"],
            row["value"],
        ),
    ):
        key = (row["category"], row["name"], row["value"], row["source"])
        if key in seen:
            continue
        seen.add(key)
        ordered_rows.append(row)
    return ordered_rows


def write_tsv(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write the final versions TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_COLUMNS), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_collect_versions(args: argparse.Namespace) -> None:
    """Build the final tool and database versions TSV."""
    rows = build_runtime_rows(args)
    for version_file in args.version_file:
        rows.extend(parse_versions_file(version_file))
    write_tsv(args.output, deduplicate_rows(rows))


def main(argv: Sequence[str] | None = None) -> int:
    """Run the versions collector CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_collect_versions(args)
    except CollectVersionsError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote tool and database provenance to %s.", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
