#!/usr/bin/env python3
"""Prepare pinned runtime databases in operator-managed destinations."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import shlex
import shutil
import sys
import tarfile
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Callable, Iterator, Sequence


LOGGER = logging.getLogger(__name__)
MARKER_FILE_NAME = ".nf_myco_ready.json"
LOCK_FILE_SUFFIX = ".nf_myco_prepare.lock"
DEFAULT_BUSCO_LINEAGES = ("bacillota_odb12", "mycoplasmatota_odb12")
REPORT_COLUMNS = ("component", "status", "source", "destination", "details")
LINK_MODES = ("copy", "symlink", "hardlink")


class PrepareRuntimeDatabasesError(RuntimeError):
    """Raised when one runtime database cannot be prepared safely."""


@dataclass(frozen=True)
class ValidationResult:
    """Describe one validated runtime database layout."""

    required_paths: tuple[str, ...]
    details: dict[str, str]


@dataclass(frozen=True)
class PreparationRecord:
    """Summarise one prepared or reused database destination."""

    component: str
    status: str
    source: Path | None
    destination: Path
    details: str


Validator = Callable[[Path], ValidationResult]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Prepare pinned runtime databases in operator-managed destinations "
            "without relying on system tmp."
        )
    )
    parser.add_argument("--taxdump-source", type=Path, help="Pinned taxdump source path.")
    parser.add_argument("--taxdump-dest", type=Path, help="Final taxdump destination.")
    parser.add_argument("--checkm2-source", type=Path, help="CheckM2 database source path.")
    parser.add_argument("--checkm2-dest", type=Path, help="Final CheckM2 destination.")
    parser.add_argument(
        "--busco-lineage-source",
        action="append",
        default=[],
        help="BUSCO lineage mapping as LINEAGE=PATH. May be supplied multiple times.",
    )
    parser.add_argument(
        "--busco-dest-root",
        type=Path,
        help="Final BUSCO lineage root directory.",
    )
    parser.add_argument("--eggnog-source", type=Path, help="eggNOG database source path.")
    parser.add_argument("--eggnog-dest", type=Path, help="Final eggNOG destination.")
    parser.add_argument("--padloc-source", type=Path, help="PADLOC database source path.")
    parser.add_argument("--padloc-dest", type=Path, help="Final PADLOC destination.")
    parser.add_argument(
        "--scratch-root",
        type=Path,
        default=None,
        help=(
            "Optional scratch root for archive extraction. Defaults to the "
            "destination parent directory."
        ),
    )
    parser.add_argument(
        "--link-mode",
        choices=LINK_MODES,
        default="copy",
        help="How directory sources are materialised into final destinations.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace invalid or incomplete destinations.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Optional TSV report path describing the prepared destinations.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_path(path: Path) -> Path:
    """Expand one user path without requiring it to exist yet."""
    return path.expanduser().resolve(strict=False)


def timestamp_token() -> str:
    """Return one UTC timestamp token suitable for file names."""
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def marker_path(destination: Path) -> Path:
    """Return the ready-marker path for one prepared destination."""
    return destination / MARKER_FILE_NAME


def lock_path(destination: Path) -> Path:
    """Return the sidecar lock path for one destination."""
    return destination.parent / f".{destination.name}{LOCK_FILE_SUFFIX}"


def report_row(record: PreparationRecord) -> dict[str, str]:
    """Build one TSV report row."""
    return {
        "component": record.component,
        "status": record.status,
        "source": str(record.source) if record.source else "NA",
        "destination": str(record.destination),
        "details": record.details,
    }


def parse_name_path(token: str, argument_name: str) -> tuple[str, Path]:
    """Parse one NAME=PATH token."""
    if "=" not in token:
        raise PrepareRuntimeDatabasesError(
            f"{argument_name} entries must use LINEAGE=PATH syntax: {token!r}"
        )
    name, raw_path = token.split("=", 1)
    clean_name = name.strip()
    if not clean_name:
        raise PrepareRuntimeDatabasesError(
            f"{argument_name} entry is missing the lineage before '=': {token!r}"
        )
    clean_path = raw_path.strip()
    if not clean_path:
        raise PrepareRuntimeDatabasesError(
            f"{argument_name} entry is missing the source path after '=': {token!r}"
        )
    return clean_name, normalise_path(Path(clean_path))


def ensure_existing_path(path: Path, label: str) -> Path:
    """Require one existing source path."""
    source = normalise_path(path)
    if not source.exists():
        raise PrepareRuntimeDatabasesError(f"Missing {label}: {source}")
    return source


def ensure_parent_directory(path: Path) -> None:
    """Create the parent directory for one destination or report."""
    path.parent.mkdir(parents=True, exist_ok=True)


def remove_path(path: Path) -> None:
    """Remove one filesystem path regardless of type."""
    if not path.exists() and not path.is_symlink():
        return
    if path.is_symlink() or path.is_file():
        path.unlink()
        return
    shutil.rmtree(path)


def destination_is_empty(path: Path) -> bool:
    """Return whether one existing destination directory is empty."""
    return not any(path.iterdir())


def load_marker(destination: Path, component: str) -> dict[str, object] | None:
    """Load and validate one ready marker."""
    path = marker_path(destination)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise PrepareRuntimeDatabasesError(
            f"Marker file is not valid JSON for {component}: {path}"
        ) from error
    if not isinstance(data, dict):
        raise PrepareRuntimeDatabasesError(
            f"Marker file must contain a JSON object for {component}: {path}"
        )
    if data.get("component") != component:
        raise PrepareRuntimeDatabasesError(
            f"Marker component mismatch for {component}: {path}"
        )
    return data


def marker_exists(destination: Path) -> bool:
    """Return whether one destination already has a ready marker file."""
    return marker_path(destination).is_file()


def write_marker(
    *,
    component: str,
    source: Path,
    destination: Path,
    validation: ValidationResult,
) -> None:
    """Write one ready marker after validation has succeeded."""
    payload = {
        "component": component,
        "source": str(source),
        "destination": str(destination),
        "prepared_at": datetime.now(UTC).isoformat(),
        "required_paths": list(validation.required_paths),
        "details": validation.details,
    }
    marker_path(destination).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


@contextmanager
def preparation_lock(destination: Path, *, force: bool) -> Iterator[Path]:
    """Guard one destination against concurrent preparation."""
    sidecar = lock_path(destination)
    ensure_parent_directory(sidecar)
    if sidecar.exists():
        if not force:
            raise PrepareRuntimeDatabasesError(
                f"Preparation lock already exists for {destination}: {sidecar}"
            )
        LOGGER.warning("Removing existing preparation lock for %s", destination)
        remove_path(sidecar)
    payload = {
        "destination": str(destination),
        "pid": os.getpid(),
        "created_at": datetime.now(UTC).isoformat(),
    }
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    try:
        yield sidecar
    finally:
        sidecar.unlink(missing_ok=True)


def validate_taxdump(path: Path) -> ValidationResult:
    """Validate one pinned taxdump directory."""
    required = ("names.dmp", "nodes.dmp")
    missing = [name for name in required if not (path / name).is_file()]
    if missing:
        raise PrepareRuntimeDatabasesError(
            f"Taxdump directory is missing required files in {path}: {', '.join(missing)}"
        )
    return ValidationResult(
        required_paths=required,
        details={"files": str(len(required))},
    )


def validate_checkm2(path: Path) -> ValidationResult:
    """Validate one CheckM2 database directory."""
    candidates = sorted(file.name for file in path.glob("*.dmnd") if file.is_file())
    if len(candidates) != 1:
        raise PrepareRuntimeDatabasesError(
            f"CheckM2 destination must contain exactly one top-level .dmnd file: {path}"
        )
    return ValidationResult(
        required_paths=(candidates[0],),
        details={"dmnd": candidates[0]},
    )


def validate_eggnog(path: Path) -> ValidationResult:
    """Validate one eggNOG data directory."""
    required = ("eggnog.db", "eggnog_proteins.dmnd")
    missing = [name for name in required if not (path / name).is_file()]
    if missing:
        raise PrepareRuntimeDatabasesError(
            f"eggNOG destination is missing required files in {path}: {', '.join(missing)}"
        )
    return ValidationResult(
        required_paths=required,
        details={"files": str(len(required))},
    )


def validate_padloc(path: Path) -> ValidationResult:
    """Validate one PADLOC data directory."""
    required = ("hmm/padlocdb.hmm",)
    if not (path / required[0]).is_file():
        raise PrepareRuntimeDatabasesError(
            f"PADLOC destination must contain hmm/padlocdb.hmm: {path}"
        )
    return ValidationResult(
        required_paths=required,
        details={"hmm": "hmm/padlocdb.hmm"},
    )


def build_busco_lineage_validator(lineage: str) -> Validator:
    """Build a BUSCO lineage validator for one lineage directory."""

    def validate_busco_lineage(path: Path) -> ValidationResult:
        """Validate one BUSCO lineage directory."""
        required = ("dataset.cfg",)
        if not (path / required[0]).is_file():
            raise PrepareRuntimeDatabasesError(
                f"BUSCO lineage {lineage} is missing dataset.cfg in {path}"
            )
        return ValidationResult(
            required_paths=required,
            details={"lineage": lineage},
        )

    return validate_busco_lineage


def build_busco_validator(lineages: Sequence[str]) -> Validator:
    """Build a BUSCO root validator for one lineage set."""

    def validate_busco_root(path: Path) -> ValidationResult:
        """Validate one BUSCO lineage root directory."""
        required = tuple(f"{lineage}/dataset.cfg" for lineage in lineages)
        missing = [entry for entry in required if not (path / entry).is_file()]
        if missing:
            raise PrepareRuntimeDatabasesError(
                f"BUSCO destination is missing required lineage files in {path}: "
                + ", ".join(missing)
            )
        return ValidationResult(
            required_paths=required,
            details={"lineages": ";".join(lineages)},
        )

    return validate_busco_root


def resolve_directory_source(path: Path, validator: Validator, component: str) -> Path:
    """Resolve one directory source to the single valid content root."""
    candidates = [path]
    candidates.extend(sorted(child for child in path.iterdir() if child.is_dir()))
    matches: list[Path] = []
    for candidate in candidates:
        try:
            validator(candidate)
        except PrepareRuntimeDatabasesError:
            continue
        matches.append(candidate)
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise PrepareRuntimeDatabasesError(
            f"Could not find a valid {component} root inside {path}"
        )
    raise PrepareRuntimeDatabasesError(
        f"Found multiple valid {component} roots inside {path}; choose a more specific source path."
    )


def copy_directory_tree(source: Path, destination: Path) -> None:
    """Copy one directory tree into an existing destination directory."""
    for child in source.iterdir():
        target = destination / child.name
        if child.is_dir():
            shutil.copytree(child, target)
        else:
            shutil.copy2(child, target)


def hardlink_directory_tree(source: Path, destination: Path) -> None:
    """Materialise one directory tree with hardlinked files."""
    for child in source.iterdir():
        target = destination / child.name
        if child.is_dir():
            target.mkdir(parents=True, exist_ok=False)
            hardlink_directory_tree(child, target)
            continue
        os.link(child, target)


def symlink_directory_tree(source: Path, destination: Path) -> None:
    """Materialise one directory tree with top-level symlinks."""
    for child in source.iterdir():
        target = destination / child.name
        os.symlink(child, target, target_is_directory=child.is_dir())


def materialise_directory(source: Path, destination: Path, link_mode: str) -> None:
    """Materialise one directory source into the final destination."""
    destination.mkdir(parents=True, exist_ok=True)
    if not destination_is_empty(destination):
        raise PrepareRuntimeDatabasesError(
            f"Destination must be absent or empty before preparation: {destination}"
        )
    if link_mode == "copy":
        copy_directory_tree(source, destination)
        return
    if link_mode == "hardlink":
        hardlink_directory_tree(source, destination)
        return
    symlink_directory_tree(source, destination)


def materialise_checkm2_file(source: Path, destination: Path, link_mode: str) -> None:
    """Materialise one raw CheckM2 .dmnd source into the final destination."""
    destination.mkdir(parents=True, exist_ok=True)
    if not destination_is_empty(destination):
        raise PrepareRuntimeDatabasesError(
            f"Destination must be absent or empty before preparation: {destination}"
        )
    target = destination / source.name
    if link_mode == "copy":
        shutil.copy2(source, target)
        return
    if link_mode == "hardlink":
        os.link(source, target)
        return
    os.symlink(source, target)


def ensure_same_filesystem(scratch_root: Path, destination_parent: Path, component: str) -> None:
    """Require archive scratch and destination to share one filesystem."""
    if scratch_root.stat().st_dev != destination_parent.stat().st_dev:
        raise PrepareRuntimeDatabasesError(
            f"{component} archive scratch root must be on the same filesystem as "
            f"the destination parent for rename-based finalisation: "
            f"{scratch_root} -> {destination_parent}"
        )


def extract_zip(archive_path: Path, destination: Path) -> None:
    """Extract one zip archive safely."""
    destination_root = destination.resolve()
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            if not target.is_relative_to(destination_root):
                raise PrepareRuntimeDatabasesError(
                    f"Zip archive contains an unsafe path: {member.filename!r}"
                )
        archive.extractall(destination)


def extract_tar(archive_path: Path, destination: Path) -> None:
    """Extract one tar archive safely."""
    with tarfile.open(archive_path) as archive:
        archive.extractall(destination, filter="data")


def extract_archive(archive_path: Path, scratch_dir: Path) -> None:
    """Extract one supported archive into scratch space."""
    if zipfile.is_zipfile(archive_path):
        extract_zip(archive_path, scratch_dir)
        return
    if tarfile.is_tarfile(archive_path):
        extract_tar(archive_path, scratch_dir)
        return
    raise PrepareRuntimeDatabasesError(f"Unsupported archive format: {archive_path}")


def prepare_from_archive(
    *,
    component: str,
    source: Path,
    destination: Path,
    validator: Validator,
    scratch_root: Path | None,
) -> ValidationResult:
    """Prepare one database from an archive beside the final destination."""
    destination_parent = destination.parent
    destination_parent.mkdir(parents=True, exist_ok=True)
    scratch_base = normalise_path(scratch_root) if scratch_root else destination_parent
    scratch_base.mkdir(parents=True, exist_ok=True)
    ensure_same_filesystem(scratch_base, destination_parent, component)

    scratch_dir = scratch_base / f".prepare-{destination.name}-{timestamp_token()}"
    if scratch_dir.exists():
        raise PrepareRuntimeDatabasesError(
            f"Archive scratch directory already exists: {scratch_dir}"
        )

    try:
        scratch_dir.mkdir(parents=True, exist_ok=False)
        extract_archive(source, scratch_dir)
        resolved_root = resolve_directory_source(scratch_dir, validator, component)
        validation = validator(resolved_root)
        if destination.exists():
            if destination.is_dir() and destination_is_empty(destination):
                destination.rmdir()
            else:
                raise PrepareRuntimeDatabasesError(
                    f"Archive destination must be absent or empty before finalisation: {destination}"
                )
        resolved_root.rename(destination)
        if scratch_dir.exists():
            shutil.rmtree(scratch_dir)
        return validation
    except Exception:
        shutil.rmtree(scratch_dir, ignore_errors=True)
        raise


def describe_validation(validation: ValidationResult) -> str:
    """Render one validation summary into a short report string."""
    paths = ",".join(validation.required_paths)
    details = ";".join(f"{key}={value}" for key, value in sorted(validation.details.items()))
    return f"required_paths={paths};{details}"


def assess_destination(
    *,
    component: str,
    destination: Path,
    validator: Validator,
) -> tuple[str, ValidationResult | None]:
    """Classify one destination before preparation."""
    if not destination.exists():
        return "missing", None
    if not destination.is_dir():
        raise PrepareRuntimeDatabasesError(
            f"Destination must be a directory for {component}: {destination}"
        )
    try:
        validation = validator(destination)
    except PrepareRuntimeDatabasesError as error:
        if marker_exists(destination):
            raise PrepareRuntimeDatabasesError(
                f"Destination marker exists but validation failed for {component}: {error}"
            ) from error
        raise PrepareRuntimeDatabasesError(
            f"Destination is invalid or partially prepared for {component}: {error}"
        ) from error
    marker = load_marker(destination, component)
    if marker is None:
        raise PrepareRuntimeDatabasesError(
            f"Destination is valid but incomplete for {component}; "
            f"missing {MARKER_FILE_NAME}: {destination}"
        )
    return "ready", validation


def prepare_component(
    *,
    component: str,
    source: Path,
    destination: Path,
    validator: Validator,
    link_mode: str,
    scratch_root: Path | None,
    force: bool,
    allow_file_source: bool = False,
) -> PreparationRecord:
    """Prepare one non-BUSCO database destination."""
    source = ensure_existing_path(source, f"{component} source")
    destination = normalise_path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    try:
        _state, existing_validation = assess_destination(
            component=component,
            destination=destination,
            validator=validator,
        )
        if existing_validation is not None:
            LOGGER.info("Reusing ready %s destination: %s", component, destination)
            return PreparationRecord(
                component=component,
                status="present",
                source=source,
                destination=destination,
                details=describe_validation(existing_validation),
            )
    except PrepareRuntimeDatabasesError:
        if not destination.exists() or not force:
            raise
        LOGGER.warning("Replacing existing %s destination at %s", component, destination)
        remove_path(destination)

    with preparation_lock(destination, force=force):
        if source.is_dir():
            resolved_source = resolve_directory_source(source, validator, component)
            materialise_directory(resolved_source, destination, link_mode)
            validation = validator(destination)
        elif allow_file_source and source.is_file() and source.suffix == ".dmnd":
            materialise_checkm2_file(source, destination, link_mode)
            validation = validator(destination)
        else:
            validation = prepare_from_archive(
                component=component,
                source=source,
                destination=destination,
                validator=validator,
                scratch_root=scratch_root,
            )
        write_marker(
            component=component,
            source=source,
            destination=destination,
            validation=validation,
        )
    return PreparationRecord(
        component=component,
        status="prepared",
        source=source,
        destination=destination,
        details=describe_validation(validation),
    )


def determine_busco_lineages(mapping: dict[str, Path]) -> tuple[str, ...]:
    """Determine the BUSCO lineages that should exist after preparation."""
    if mapping:
        return tuple(sorted(mapping))
    return DEFAULT_BUSCO_LINEAGES


def prepare_busco_lineages(
    *,
    source_mapping: dict[str, Path],
    destination_root: Path,
    link_mode: str,
    scratch_root: Path | None,
    force: bool,
) -> list[PreparationRecord]:
    """Prepare the requested BUSCO lineage directories beneath one root."""
    destination_root = normalise_path(destination_root)
    destination_root.mkdir(parents=True, exist_ok=True)
    lineages = determine_busco_lineages(source_mapping)
    records: list[PreparationRecord] = []
    for lineage in lineages:
        source = source_mapping.get(lineage)
        destination = destination_root / lineage
        validator = build_busco_lineage_validator(lineage)
        component = f"busco:{lineage}"
        try:
            _state, existing_validation = assess_destination(
                component=component,
                destination=destination,
                validator=validator,
            )
            if existing_validation is not None:
                LOGGER.info("Reusing ready BUSCO lineage %s at %s", lineage, destination)
                records.append(
                    PreparationRecord(
                        component=component,
                        status="present",
                        source=source,
                        destination=destination,
                        details=describe_validation(existing_validation),
                    )
                )
                continue
        except PrepareRuntimeDatabasesError:
            if destination.exists() and not force:
                raise
            if destination.exists():
                LOGGER.warning("Replacing existing BUSCO lineage at %s", destination)
                remove_path(destination)

        if source is None:
            raise PrepareRuntimeDatabasesError(
                f"BUSCO lineage {lineage} is missing at {destination} and no source was supplied."
            )

        records.append(
            prepare_component(
                component=component,
                source=source,
                destination=destination,
                validator=validator,
                link_mode=link_mode,
                scratch_root=scratch_root,
                force=force,
            )
        )

    root_validation = build_busco_validator(lineages)(destination_root)
    records.append(
        PreparationRecord(
            component="busco_root",
            status="prepared",
            source=None,
            destination=destination_root,
            details=describe_validation(root_validation),
        )
    )
    return records


def validate_argument_pair(
    source: Path | None,
    destination: Path | None,
    label: str,
) -> None:
    """Require one source and destination to be supplied together."""
    if (source is None) == (destination is None):
        return
    raise PrepareRuntimeDatabasesError(
        f"{label} requires both source and destination arguments."
    )


def build_component_records(args: argparse.Namespace) -> list[PreparationRecord]:
    """Prepare the requested database destinations and return report records."""
    validate_argument_pair(args.taxdump_source, args.taxdump_dest, "Taxdump")
    validate_argument_pair(args.checkm2_source, args.checkm2_dest, "CheckM2")
    validate_argument_pair(args.eggnog_source, args.eggnog_dest, "eggNOG")
    validate_argument_pair(args.padloc_source, args.padloc_dest, "PADLOC")
    if args.busco_dest_root is None and args.busco_lineage_source:
        raise PrepareRuntimeDatabasesError(
            "BUSCO lineage sources require --busco-dest-root."
        )

    busco_sources = {
        lineage: ensure_existing_path(path, f"BUSCO lineage source for {lineage}")
        for lineage, path in (
            parse_name_path(token, "--busco-lineage-source")
            for token in args.busco_lineage_source
        )
    }

    records: list[PreparationRecord] = []
    if args.taxdump_source and args.taxdump_dest:
        records.append(
            prepare_component(
                component="taxdump",
                source=args.taxdump_source,
                destination=args.taxdump_dest,
                validator=validate_taxdump,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
            )
        )
    if args.checkm2_source and args.checkm2_dest:
        records.append(
            prepare_component(
                component="checkm2",
                source=args.checkm2_source,
                destination=args.checkm2_dest,
                validator=validate_checkm2,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
                allow_file_source=True,
            )
        )
    if args.busco_dest_root:
        records.extend(
            prepare_busco_lineages(
                source_mapping=busco_sources,
                destination_root=args.busco_dest_root,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
            )
        )
    if args.eggnog_source and args.eggnog_dest:
        records.append(
            prepare_component(
                component="eggnog",
                source=args.eggnog_source,
                destination=args.eggnog_dest,
                validator=validate_eggnog,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
            )
        )
    if args.padloc_source and args.padloc_dest:
        records.append(
            prepare_component(
                component="padloc",
                source=args.padloc_source,
                destination=args.padloc_dest,
                validator=validate_padloc,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
            )
        )
    if not records:
        raise PrepareRuntimeDatabasesError(
            "No database targets were requested. Supply at least one source/destination pair."
        )
    return records


def write_report(path: Path, records: Sequence[PreparationRecord]) -> None:
    """Write one TSV report describing the prepared destinations."""
    report_path = normalise_path(path)
    ensure_parent_directory(report_path)
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow(report_row(record))


def build_nextflow_arguments(records: Sequence[PreparationRecord]) -> list[tuple[str, str]]:
    """Build one copy-pastable Nextflow argument block."""
    mapping: dict[str, str] = {}
    for record in records:
        if record.component == "taxdump":
            mapping["--taxdump"] = str(record.destination)
        elif record.component == "checkm2":
            mapping["--checkm2_db"] = str(record.destination)
        elif record.component == "busco_root":
            mapping["--busco_download_dir"] = str(record.destination)
        elif record.component == "eggnog":
            mapping["--eggnog_db"] = str(record.destination)
        elif record.component == "padloc":
            mapping["--padloc_db"] = str(record.destination)
    return sorted(mapping.items())


def print_nextflow_arguments(records: Sequence[PreparationRecord]) -> None:
    """Print one ready-to-paste Nextflow command fragment."""
    arguments = build_nextflow_arguments(records)
    if not arguments:
        return
    print("Nextflow arguments:")
    print("nextflow run . \\")
    for index, (flag, value) in enumerate(arguments):
        suffix = " \\" if index < len(arguments) - 1 else ""
        print(f"  {flag} {shlex.quote(value)}{suffix}")


def main(argv: Sequence[str] | None = None) -> int:
    """Prepare requested runtime databases and print reusable paths."""
    args = parse_args(argv)
    configure_logging()
    try:
        records = build_component_records(args)
        if args.report:
            write_report(args.report, records)
        print_nextflow_arguments(records)
    except PrepareRuntimeDatabasesError as error:
        LOGGER.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
