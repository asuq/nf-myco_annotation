#!/usr/bin/env python3
"""Validate OIST HPC database-prep and real-run outputs."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Sequence

import run_acceptance_tests
from prepare_runtime_databases import (
    MARKER_FILE_NAME,
    PrepareRuntimeDatabasesError,
    build_busco_validator,
    validate_checkm2,
    validate_eggnog,
    validate_taxdump,
)


ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_COHORT_PLAN = ROOT_DIR / "assets" / "testdata" / "acceptance" / "cohort_plan.tsv"
DEFAULT_SOURCE_CATALOG = ROOT_DIR / "assets" / "testdata" / "acceptance" / "source_catalog.tsv"
PIPELINE_INFO_OUTPUTS = ("trace.tsv", "report.html", "timeline.html", "dag.html")


class ValidateHpcMatrixError(RuntimeError):
    """Raised when one HPC matrix validation fails."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate OIST HPC database-prep and real-run outputs."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    dbprep_parser = subparsers.add_parser(
        "dbprep",
        help="Validate one database-prep case.",
    )
    dbprep_parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Database-prep results directory.",
    )
    dbprep_parser.add_argument("--taxdump", type=Path, default=None)
    dbprep_parser.add_argument("--checkm2-db", type=Path, default=None)
    dbprep_parser.add_argument("--busco-db", type=Path, default=None)
    dbprep_parser.add_argument("--eggnog-db", type=Path, default=None)
    dbprep_parser.add_argument(
        "--expected-component",
        action="append",
        default=[],
        help="Expected report component. May be supplied multiple times.",
    )
    dbprep_parser.add_argument(
        "--expected-status",
        choices=("prepared", "present"),
        default=None,
        help="Expected status for all report rows.",
    )
    dbprep_parser.add_argument(
        "--expected-arg",
        action="append",
        default=[],
        help="Expected argument flag in nextflow_args.txt. May be supplied multiple times.",
    )

    tracked_parser = subparsers.add_parser(
        "tracked-run",
        help="Validate the tracked 9-sample real run.",
    )
    tracked_parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Published output directory.",
    )
    tracked_parser.add_argument(
        "--metadata-tsv",
        type=Path,
        required=True,
        help="Metadata TSV used for the run.",
    )
    tracked_parser.add_argument(
        "--db-root",
        type=Path,
        required=True,
        help="Root directory for the prepared runtime databases.",
    )
    tracked_parser.add_argument(
        "--cohort-plan",
        type=Path,
        default=DEFAULT_COHORT_PLAN,
        help="Tracked cohort plan TSV.",
    )
    tracked_parser.add_argument(
        "--source-catalog",
        type=Path,
        default=DEFAULT_SOURCE_CATALOG,
        help="Tracked source catalog TSV.",
    )

    medium_parser = subparsers.add_parser(
        "medium-run",
        help="Validate the medium real-data run.",
    )
    medium_parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Published output directory.",
    )
    medium_parser.add_argument(
        "--db-root",
        type=Path,
        required=True,
        help="Root directory for the prepared runtime databases.",
    )
    medium_parser.add_argument(
        "--sample-count",
        type=int,
        default=20,
        help="Expected number of rows in the master table.",
    )
    medium_parser.add_argument(
        "--allowed-phylum",
        action="append",
        default=["Mycoplasmatota", "Bacillota"],
        help="Allowed phylum. May be supplied multiple times.",
    )

    return parser.parse_args(argv)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read one TSV file into a header and row dictionaries."""
    if not path.is_file():
        raise ValidateHpcMatrixError(f"Missing TSV file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValidateHpcMatrixError(f"Empty TSV file: {path}")
        return list(reader.fieldnames), list(reader)


def read_rows_indexed_by(path: Path, key_column: str) -> dict[str, dict[str, str]]:
    """Read one TSV file and index rows by one unique key column."""
    _header, rows = read_tsv(path)
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get(key_column, "").strip()
        if not key:
            raise ValidateHpcMatrixError(f"Empty key {key_column} in {path}")
        if key in indexed:
            raise ValidateHpcMatrixError(f"Duplicate key {key!r} in {path}")
        indexed[key] = row
    return indexed


def require_pipeline_info(outdir: Path) -> None:
    """Require the standard pipeline-info artefacts for one run."""
    pipeline_info_dir = outdir / "pipeline_info"
    missing = [
        name for name in PIPELINE_INFO_OUTPUTS if not (pipeline_info_dir / name).is_file()
    ]
    if missing:
        raise ValidateHpcMatrixError(
            "Missing pipeline_info artefacts: " + ", ".join(sorted(missing))
        )


def assert_expected_components(
    rows: Sequence[dict[str, str]],
    expected_components: Sequence[str],
) -> None:
    """Assert that the report contains exactly the requested components."""
    if not expected_components:
        return
    seen = {row["component"] for row in rows}
    expected = set(expected_components)
    if seen != expected:
        raise ValidateHpcMatrixError(
            "Unexpected dbprep report components: "
            f"expected {sorted(expected)!r}, observed {sorted(seen)!r}"
        )


def assert_expected_status(
    rows: Sequence[dict[str, str]],
    expected_status: str | None,
) -> None:
    """Assert that all report rows share one expected status."""
    if expected_status is None:
        return
    bad_rows = [row["component"] for row in rows if row["status"] != expected_status]
    if bad_rows:
        raise ValidateHpcMatrixError(
            "Unexpected dbprep statuses for: " + ", ".join(sorted(bad_rows))
        )


def assert_expected_args(args_path: Path, expected_args: Sequence[str]) -> None:
    """Assert that the nextflow-args file contains exactly the expected DB flags."""
    if not expected_args:
        return
    contents = args_path.read_text(encoding="utf-8")
    db_flags = {"--taxdump", "--checkm2_db", "--busco_db", "--eggnog_db"}
    for flag in expected_args:
        if flag not in contents:
            raise ValidateHpcMatrixError(f"Missing expected flag {flag} in {args_path}")
    unexpected = sorted(
        flag for flag in db_flags if flag in contents and flag not in set(expected_args)
    )
    if unexpected:
        raise ValidateHpcMatrixError(
            "Unexpected database flags in nextflow_args.txt: " + ", ".join(unexpected)
        )


def validate_dbprep_component(
    component: str,
    destination: Path,
) -> None:
    """Validate one prepared database destination and its ready marker."""
    destination = destination.resolve()
    if not (destination / MARKER_FILE_NAME).is_file():
        raise ValidateHpcMatrixError(f"Missing ready marker for {component}: {destination}")
    if component == "taxdump":
        validate_taxdump(destination)
        return
    if component == "checkm2":
        validate_checkm2(destination)
        return
    if component == "busco_root":
        validator = build_busco_validator(
            run_acceptance_tests.DEFAULT_DBPREP_BUSCO_LINEAGES
        )
        validator(destination)
        return
    if component == "eggnog":
        validate_eggnog(destination)
        return
    raise ValidateHpcMatrixError(f"Unsupported dbprep component: {component}")


def validate_dbprep(args: argparse.Namespace) -> None:
    """Validate one database-prep case."""
    outputs = run_acceptance_tests.require_dbprep_outputs(args.results_dir.resolve())
    _header, rows = read_tsv(outputs["runtime_database_report.tsv"])
    assert_expected_components(rows, args.expected_component)
    assert_expected_status(rows, args.expected_status)
    assert_expected_args(outputs["nextflow_args.txt"], args.expected_arg)

    destinations = {
        "taxdump": args.taxdump,
        "checkm2": args.checkm2_db,
        "busco_root": args.busco_db,
        "eggnog": args.eggnog_db,
    }
    for row in rows:
        destination = destinations.get(row["component"])
        if destination is None:
            raise ValidateHpcMatrixError(
                f"Missing destination argument for component {row['component']}"
            )
        validate_dbprep_component(row["component"], destination)


def assert_no_failed_statuses(sample_status_path: Path) -> None:
    """Assert that the final sample-status table contains no failures."""
    _header, rows = read_tsv(sample_status_path)
    failed_rows = [
        row["accession"]
        for row in rows
        if any(value.strip() == "failed" for value in row.values())
    ]
    if failed_rows:
        raise ValidateHpcMatrixError(
            "Found failed sample-status rows: " + ", ".join(sorted(failed_rows))
        )


def assert_versions_table_clean(versions_path: Path, db_root: Path) -> None:
    """Assert that the final versions table is complete and points at the HPC DB root."""
    _header, rows = read_tsv(versions_path)
    if any(row["version"] == "unknown" for row in rows):
        raise ValidateHpcMatrixError("tool_and_db_versions.tsv contains unknown rows")

    eggnog_rows = [row for row in rows if row["component"] == "eggnog_mapper"]
    if len(eggnog_rows) != 1 or eggnog_rows[0]["version"] != "2.1.13":
        raise ValidateHpcMatrixError("eggnog_mapper version is not 2.1.13")

    database_rows = [row for row in rows if row["kind"] == "database"]
    expected_prefix = str(db_root.resolve())
    bad_rows = [
        row["component"]
        for row in database_rows
        if not row["image_or_path"].startswith(expected_prefix)
    ]
    if bad_rows:
        raise ValidateHpcMatrixError(
            "Database paths do not point at the HPC DB root: " + ", ".join(sorted(bad_rows))
        )


def load_tracked_plan(
    source_catalog_path: Path,
    cohort_plan_path: Path,
) -> tuple[run_acceptance_tests.CohortRecord, ...]:
    """Load the tracked cohort plan with the current harness logic."""
    catalog = run_acceptance_tests.load_source_catalog(source_catalog_path)
    return run_acceptance_tests.load_cohort_plan(cohort_plan_path, catalog)


def validate_tracked_run(args: argparse.Namespace) -> None:
    """Validate the tracked 9-sample run outputs."""
    require_pipeline_info(args.outdir.resolve())
    plan = load_tracked_plan(args.source_catalog, args.cohort_plan)
    run_acceptance_tests.assert_run_outputs(
        args.outdir.resolve(),
        plan,
        args.metadata_tsv.resolve(),
    )
    sample_status_path = args.outdir.resolve() / "tables" / "sample_status.tsv"
    versions_path = args.outdir.resolve() / "tables" / "tool_and_db_versions.tsv"
    assert_no_failed_statuses(sample_status_path)
    assert_versions_table_clean(versions_path, args.db_root)


def validate_medium_run(args: argparse.Namespace) -> None:
    """Validate the medium real-data run outputs."""
    require_pipeline_info(args.outdir.resolve())
    tables = run_acceptance_tests.require_output_tables(args.outdir.resolve())
    sample_status_path = tables["sample_status.tsv"]
    versions_path = tables["tool_and_db_versions.tsv"]
    master_table_path = tables["master_table.tsv"]

    assert_no_failed_statuses(sample_status_path)
    assert_versions_table_clean(versions_path, args.db_root)

    master_rows = read_rows_indexed_by(master_table_path, "Accession")
    if len(master_rows) != args.sample_count:
        raise ValidateHpcMatrixError(
            f"Expected {args.sample_count} samples, observed {len(master_rows)}"
        )

    allowed = {value.strip() for value in args.allowed_phylum if value.strip()}
    observed = {row["phylum"] for row in master_rows.values()}
    bad_phyla = sorted(phylum for phylum in observed if phylum not in allowed)
    if bad_phyla:
        raise ValidateHpcMatrixError(
            "Observed disallowed phyla in medium run: " + ", ".join(bad_phyla)
        )
    missing_allowed = sorted(phylum for phylum in allowed if phylum not in observed)
    if missing_allowed:
        raise ValidateHpcMatrixError(
            "Medium run is missing allowed phyla: " + ", ".join(missing_allowed)
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run one HPC matrix validation command."""
    args = parse_args(argv)
    try:
        if args.command == "dbprep":
            validate_dbprep(args)
        elif args.command == "tracked-run":
            validate_tracked_run(args)
        else:
            validate_medium_run(args)
    except (
        ValidateHpcMatrixError,
        run_acceptance_tests.AcceptanceTestError,
        PrepareRuntimeDatabasesError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
