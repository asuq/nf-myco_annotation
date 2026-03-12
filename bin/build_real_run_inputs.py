#!/usr/bin/env python3
"""Build medium real-run inputs for the OIST HPC campaign."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import run_acceptance_tests
import validate_inputs


DEFAULT_ALLOWED_PHYLA = ("Mycoplasmatota", "Bacillota")
DEFAULT_SAMPLE_COUNT = 20
REQUIRED_COLUMNS = (
    "accession",
    "genome_fasta",
    "is_new",
    "assembly_level",
    "include_metadata",
    "tax_id",
    "organism_name",
    "n50",
    "scaffolds",
    "genome_size",
    "atypical_warnings",
    "phylum",
    "role_tags",
)
EDGE_RULES = (
    ("gcode4_candidate", "tag", 1, False),
    ("gcode11_candidate", "tag", 1, False),
    ("crispr_positive_candidate", "tag", 1, False),
    ("crispr_negative_candidate", "tag", 1, False),
    ("ani_cluster_candidate", "tag", 2, False),
    ("atypical_excluded_candidate", "tag", 1, True),
    ("atypical_exception_candidate", "tag", 1, True),
    ("missing_metadata_case", "tag", 1, True),
    ("collision_candidate", "collision_pair", 2, True),
)


class BuildRealRunInputsError(RuntimeError):
    """Raised when the medium real-run inputs cannot be generated."""


@dataclass(frozen=True)
class CandidateRecord:
    """Store one candidate-manifest row."""

    accession: str
    genome_fasta: str
    is_new: str
    assembly_level: str
    include_metadata: bool
    tax_id: str
    organism_name: str
    n50: str
    scaffolds: str
    genome_size: str
    atypical_warnings: str
    phylum: str
    role_tags: tuple[str, ...]
    order: int


@dataclass(frozen=True)
class CoverageRow:
    """Store one medium-cohort coverage report row."""

    coverage_key: str
    required_count: int
    available_count: int
    selected_count: int
    status: str
    details: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Build sample_sheet.csv and metadata.tsv for the medium "
            "Mycoplasmatota/Bacillota real-data run."
        )
    )
    parser.add_argument(
        "--candidate-tsv",
        required=True,
        type=Path,
        help="Candidate manifest TSV with phylum and role_tags columns.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=default_outdir(),
        help="Output directory for generated inputs and reports.",
    )
    parser.add_argument(
        "--sample-count",
        type=int,
        default=DEFAULT_SAMPLE_COUNT,
        help="Target number of selected samples. Default: 20.",
    )
    return parser.parse_args(argv)


def default_outdir() -> Path:
    """Return the default generated-input output directory."""
    hpc_root = os.environ.get("HPC_ROOT", "").strip()
    if hpc_root:
        return Path(hpc_root) / "medium_inputs" / "generated"
    return Path("medium_inputs") / "generated"


def read_candidate_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read the candidate TSV into raw row dictionaries."""
    header, rows = run_acceptance_tests.read_tsv(path)
    missing = [column for column in REQUIRED_COLUMNS if column not in header]
    if missing:
        raise BuildRealRunInputsError(
            "Candidate TSV is missing required columns: " + ", ".join(missing)
        )
    return header, rows


def parse_candidate_rows(path: Path) -> list[CandidateRecord]:
    """Parse and validate candidate-manifest rows."""
    _header, rows = read_candidate_rows(path)
    records: list[CandidateRecord] = []
    seen_accessions: set[str] = set()

    for order, row in enumerate(rows):
        accession = row["accession"].strip()
        if not accession:
            raise BuildRealRunInputsError("Candidate TSV contains an empty accession.")
        if accession in seen_accessions:
            raise BuildRealRunInputsError(
                f"Candidate TSV contains duplicate accession {accession!r}."
            )
        seen_accessions.add(accession)
        validate_inputs.ensure_path_safe_accession(accession)

        is_new = (
            "true"
            if run_acceptance_tests.normalise_boolean(
                row["is_new"], field_name="is_new"
            )
            else "false"
        )
        include_metadata = run_acceptance_tests.normalise_boolean(
            row["include_metadata"],
            field_name="include_metadata",
        )
        if is_new == "false" and not include_metadata:
            raise BuildRealRunInputsError(
                f"Candidate {accession!r} cannot omit metadata when is_new=false."
            )

        role_tags = run_acceptance_tests.parse_role_tags(row["role_tags"])
        phylum = row["phylum"].strip()
        if not phylum:
            raise BuildRealRunInputsError(
                f"Candidate {accession!r} is missing a phylum value."
            )

        records.append(
            CandidateRecord(
                accession=accession,
                genome_fasta=row["genome_fasta"].strip(),
                is_new=is_new,
                assembly_level=run_acceptance_tests.normalise_optional(
                    row["assembly_level"]
                ),
                include_metadata=include_metadata,
                tax_id=run_acceptance_tests.normalise_optional(row["tax_id"]),
                organism_name=run_acceptance_tests.normalise_optional(
                    row["organism_name"]
                ),
                n50=run_acceptance_tests.normalise_optional(row["n50"]),
                scaffolds=run_acceptance_tests.normalise_optional(row["scaffolds"]),
                genome_size=run_acceptance_tests.normalise_optional(
                    row["genome_size"]
                ),
                atypical_warnings=run_acceptance_tests.normalise_optional(
                    row["atypical_warnings"]
                ),
                phylum=phylum,
                role_tags=role_tags,
                order=order,
            )
        )
    return records


def restrict_allowed_phyla(records: Sequence[CandidateRecord]) -> list[CandidateRecord]:
    """Restrict candidates to the supported phyla for the medium run."""
    allowed = [record for record in records if record.phylum in DEFAULT_ALLOWED_PHYLA]
    counts = Counter(record.phylum for record in allowed)
    if not allowed:
        raise BuildRealRunInputsError(
            "Candidate TSV does not contain any Mycoplasmatota or Bacillota genomes."
        )
    missing = [phylum for phylum in DEFAULT_ALLOWED_PHYLA if counts[phylum] == 0]
    if missing:
        raise BuildRealRunInputsError(
            "Candidate TSV must include both allowed phyla: " + ", ".join(missing)
        )
    return allowed


def build_target_counts(
    records: Sequence[CandidateRecord],
    sample_count: int,
) -> dict[str, int]:
    """Build preferred phylum counts for fill-stage balancing."""
    counts = Counter(record.phylum for record in records)
    if len(records) < sample_count:
        raise BuildRealRunInputsError(
            f"Need at least {sample_count} allowed candidates, observed {len(records)}."
        )

    target_counts = {phylum: min(sample_count // 2, counts[phylum]) for phylum in DEFAULT_ALLOWED_PHYLA}
    remaining = sample_count - sum(target_counts.values())
    while remaining > 0:
        candidates = [
            phylum
            for phylum in DEFAULT_ALLOWED_PHYLA
            if target_counts[phylum] < counts[phylum]
        ]
        if not candidates:
            break
        candidates.sort(
            key=lambda phylum: (
                target_counts[phylum],
                -(counts[phylum] - target_counts[phylum]),
                phylum,
            )
        )
        chosen = candidates[0]
        target_counts[chosen] += 1
        remaining -= 1

    if sum(target_counts.values()) != sample_count:
        raise BuildRealRunInputsError(
            "Unable to assign preferred phylum counts for the requested sample_count."
        )
    return target_counts


def selection_priority(
    record: CandidateRecord,
    current_counts: Counter[str],
    target_counts: dict[str, int],
) -> tuple[int, int, int]:
    """Return a stable selection priority for one candidate."""
    current = current_counts[record.phylum]
    target = target_counts[record.phylum]
    if current < target:
        return (0, current - target, record.order)
    return (1, current - target, record.order)


def select_ordered_candidates(
    candidates: Sequence[CandidateRecord],
    count: int,
    selected_accessions: set[str],
    current_counts: Counter[str],
    target_counts: dict[str, int],
) -> list[CandidateRecord]:
    """Select up to one requested number of unselected candidates."""
    remaining = [record for record in candidates if record.accession not in selected_accessions]
    remaining.sort(key=lambda record: selection_priority(record, current_counts, target_counts))
    return remaining[:count]


def select_collision_candidates(
    candidates: Sequence[CandidateRecord],
    selected_accessions: set[str],
    current_counts: Counter[str],
    target_counts: dict[str, int],
) -> tuple[list[CandidateRecord], int, str]:
    """Select one collision-candidate pair if a sanitised-base pair exists."""
    grouped: dict[str, list[CandidateRecord]] = defaultdict(list)
    for record in candidates:
        if record.accession in selected_accessions:
            continue
        if "collision_candidate" not in record.role_tags:
            continue
        grouped[validate_inputs.sanitise_accession(record.accession)].append(record)

    qualifying = {
        base: sorted(group, key=lambda record: record.order)
        for base, group in grouped.items()
        if len(group) >= 2
    }
    if not qualifying:
        return [], 0, "No collision_candidate pair with the same sanitised base is available."

    ordered_groups = sorted(
        qualifying.items(),
        key=lambda item: (
            min(record.order for record in item[1]),
            item[0],
        ),
    )
    base, group = ordered_groups[0]
    chosen = select_ordered_candidates(
        group,
        2,
        selected_accessions,
        current_counts,
        target_counts,
    )
    return chosen, len(group), f"Selected sanitised collision base {base!r}."


def add_selected_records(
    chosen: Sequence[CandidateRecord],
    *,
    selected: list[CandidateRecord],
    selected_accessions: set[str],
    selection_reasons: dict[str, str],
    current_counts: Counter[str],
    reason: str,
) -> None:
    """Add one ordered list of selected records to the final cohort."""
    for record in chosen:
        if record.accession in selected_accessions:
            continue
        selected.append(record)
        selected_accessions.add(record.accession)
        selection_reasons[record.accession] = reason
        current_counts[record.phylum] += 1


def build_coverage_row(
    key: str,
    required_count: int,
    available_count: int,
    selected_count: int,
    details: str,
) -> CoverageRow:
    """Build one coverage-report row."""
    if selected_count >= required_count:
        status = "satisfied"
    elif available_count == 0:
        status = "unavailable"
    else:
        status = "underfilled"
    return CoverageRow(
        coverage_key=key,
        required_count=required_count,
        available_count=available_count,
        selected_count=selected_count,
        status=status,
        details=details,
    )


def select_medium_cohort(
    records: Sequence[CandidateRecord],
    sample_count: int,
) -> tuple[list[CandidateRecord], dict[str, str], list[CoverageRow]]:
    """Select the medium cohort and coverage reports deterministically."""
    allowed = restrict_allowed_phyla(records)
    target_counts = build_target_counts(allowed, sample_count)
    selected: list[CandidateRecord] = []
    selected_accessions: set[str] = set()
    selection_reasons: dict[str, str] = {}
    current_counts: Counter[str] = Counter()
    coverage_rows: list[CoverageRow] = []

    for key, rule_type, required_count, optional in EDGE_RULES:
        if rule_type == "tag":
            tagged = [record for record in allowed if key in record.role_tags]
            chosen = select_ordered_candidates(
                tagged,
                required_count,
                selected_accessions,
                current_counts,
                target_counts,
            )
            add_selected_records(
                chosen,
                selected=selected,
                selected_accessions=selected_accessions,
                selection_reasons=selection_reasons,
                current_counts=current_counts,
                reason=key,
            )
            selected_count = sum(1 for record in selected if key in record.role_tags)
            details = "Optional edge-case category." if optional else "Core edge-case category."
            coverage_rows.append(
                build_coverage_row(
                    key,
                    required_count,
                    len(tagged),
                    selected_count,
                    details,
                )
            )
            continue

        chosen, available_count, details = select_collision_candidates(
            allowed,
            selected_accessions,
            current_counts,
            target_counts,
        )
        add_selected_records(
            chosen,
            selected=selected,
            selected_accessions=selected_accessions,
            selection_reasons=selection_reasons,
            current_counts=current_counts,
            reason=key,
        )
        selected_count = 0
        if available_count:
            selected_collision = [
                record
                for record in selected
                if "collision_candidate" in record.role_tags
            ]
            grouped = defaultdict(list)
            for record in selected_collision:
                grouped[validate_inputs.sanitise_accession(record.accession)].append(record)
            selected_count = max((len(group) for group in grouped.values()), default=0)
        coverage_rows.append(
            build_coverage_row(
                key,
                required_count,
                available_count,
                selected_count,
                details,
            )
        )

    remaining = [record for record in allowed if record.accession not in selected_accessions]
    while len(selected) < sample_count:
        if not remaining:
            raise BuildRealRunInputsError(
                "Not enough allowed candidates remain to fill the medium cohort."
            )
        remaining.sort(key=lambda record: selection_priority(record, current_counts, target_counts))
        chosen = remaining.pop(0)
        add_selected_records(
            [chosen],
            selected=selected,
            selected_accessions=selected_accessions,
            selection_reasons=selection_reasons,
            current_counts=current_counts,
            reason="balanced_fill",
        )
        remaining = [record for record in remaining if record.accession not in selected_accessions]

    for phylum in DEFAULT_ALLOWED_PHYLA:
        coverage_rows.append(
            build_coverage_row(
                f"phylum:{phylum}",
                1,
                sum(1 for record in allowed if record.phylum == phylum),
                sum(1 for record in selected if record.phylum == phylum),
                f"Preferred target count {target_counts[phylum]} for {phylum}.",
            )
        )

    return selected, selection_reasons, coverage_rows


def build_sample_rows(selected: Sequence[CandidateRecord]) -> list[dict[str, str]]:
    """Build sample_sheet.csv rows for the selected cohort."""
    return [
        {
            "accession": record.accession,
            "is_new": record.is_new,
            "assembly_level": record.assembly_level,
            "genome_fasta": record.genome_fasta,
        }
        for record in selected
    ]


def build_metadata_rows(selected: Sequence[CandidateRecord]) -> list[dict[str, str]]:
    """Build metadata.tsv rows for the selected cohort."""
    rows: list[dict[str, str]] = []
    for record in selected:
        if not record.include_metadata:
            continue
        rows.append(
            {
                "Accession": record.accession,
                "Tax_ID": record.tax_id,
                "Organism_Name": record.organism_name,
                "Assembly_Level": record.assembly_level,
                "N50": record.n50,
                "Scaffolds": record.scaffolds,
                "Genome_Size": record.genome_size,
                "Atypical_Warnings": record.atypical_warnings,
            }
        )
    return rows


def build_selection_rows(
    selected: Sequence[CandidateRecord],
    selection_reasons: dict[str, str],
) -> list[dict[str, str]]:
    """Build selection_report.tsv rows for the selected cohort."""
    rows: list[dict[str, str]] = []
    for index, record in enumerate(selected, start=1):
        rows.append(
            {
                "selection_order": str(index),
                "accession": record.accession,
                "phylum": record.phylum,
                "selection_reason": selection_reasons[record.accession],
                "role_tags": ";".join(record.role_tags),
                "is_new": record.is_new,
                "include_metadata": "true" if record.include_metadata else "false",
            }
        )
    return rows


def build_coverage_rows(rows: Sequence[CoverageRow]) -> list[dict[str, str]]:
    """Build serialisable coverage_report.tsv rows."""
    return [
        {
            "coverage_key": row.coverage_key,
            "required_count": str(row.required_count),
            "available_count": str(row.available_count),
            "selected_count": str(row.selected_count),
            "status": row.status,
            "details": row.details,
        }
        for row in rows
    ]


def write_outputs(
    outdir: Path,
    selected: Sequence[CandidateRecord],
    selection_reasons: dict[str, str],
    coverage_rows: Sequence[CoverageRow],
) -> None:
    """Write the generated medium-run inputs and reports."""
    run_acceptance_tests.write_csv(
        outdir / "sample_sheet.csv",
        run_acceptance_tests.SAMPLE_COLUMNS,
        build_sample_rows(selected),
    )
    run_acceptance_tests.write_tsv(
        outdir / "metadata.tsv",
        run_acceptance_tests.METADATA_COLUMNS,
        build_metadata_rows(selected),
    )
    run_acceptance_tests.write_tsv(
        outdir / "selection_report.tsv",
        (
            "selection_order",
            "accession",
            "phylum",
            "selection_reason",
            "role_tags",
            "is_new",
            "include_metadata",
        ),
        build_selection_rows(selected, selection_reasons),
    )
    run_acceptance_tests.write_tsv(
        outdir / "coverage_report.tsv",
        (
            "coverage_key",
            "required_count",
            "available_count",
            "selected_count",
            "status",
            "details",
        ),
        build_coverage_rows(coverage_rows),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Build medium real-run inputs from a candidate manifest."""
    args = parse_args(argv)
    try:
        candidates = parse_candidate_rows(args.candidate_tsv)
        selected, selection_reasons, coverage_rows = select_medium_cohort(
            candidates,
            args.sample_count,
        )
        write_outputs(args.outdir, selected, selection_reasons, coverage_rows)
    except BuildRealRunInputsError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
