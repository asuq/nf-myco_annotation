#!/usr/bin/env python3
"""Merge paired CheckM2 reports and assign gcode and QC status."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


LOGGER = logging.getLogger(__name__)
STRICT_DELTA_RULE = "strict_delta"
DELTA_THEN_ELEVEN_RULE = "delta_then_11"
GCODE_RULE_CHOICES = (STRICT_DELTA_RULE, DELTA_THEN_ELEVEN_RULE)

OUTPUT_COLUMNS = (
    "accession",
    "Completeness_gcode4",
    "Completeness_gcode11",
    "Contamination_gcode4",
    "Contamination_gcode11",
    "Coding_Density_gcode4",
    "Coding_Density_gcode11",
    "Average_Gene_Length_gcode4",
    "Average_Gene_Length_gcode11",
    "Total_Coding_Sequences_gcode4",
    "Total_Coding_Sequences_gcode11",
    "Gcode",
    "Low_quality",
    "checkm2_status",
    "warnings",
)
METRIC_ALIASES = {
    "Completeness": ("Completeness",),
    "Contamination": ("Contamination",),
    "Coding_Density": ("Coding_Density",),
    "Average_Gene_Length": ("Average_Gene_Length", "Avg_Gene_Length"),
    "Total_Coding_Sequences": ("Total_Coding_Sequences",),
}
SHARED_STAT_ALIASES = {
    "Name": ("Name",),
    "Genome_Size": ("Genome_Size",),
    "GC_Content": ("GC_Content",),
    "Contig_N50": ("Contig_N50", "N50"),
}
REQUIRED_SHARED_STATS = tuple(SHARED_STAT_ALIASES)


@dataclass(frozen=True)
class ParsedCheckM2Report:
    """Represent the parsed numeric metrics from one CheckM2 report."""

    metrics: dict[str, float]
    shared_stats: dict[str, float | str]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise paired CheckM2 reports for one sample."
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Original sample accession for the output row.",
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
    parser.add_argument(
        "--gcode-rule",
        choices=GCODE_RULE_CHOICES,
        default=STRICT_DELTA_RULE,
        help="Rule used to resolve gcode from paired CheckM2 reports.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def read_report_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV report into row dictionaries."""
    if not path.is_file():
        raise ValueError(f"Missing CheckM2 report: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"CheckM2 report is empty: {path}")
    if len(rows) != 1:
        raise ValueError(
            f"CheckM2 report must contain exactly one data row: {path}"
        )
    return rows


def first_present(row: dict[str, str], candidates: Sequence[str]) -> str:
    """Return the first present field from an alias list."""
    for candidate in candidates:
        if candidate in row:
            return row[candidate]
    raise ValueError(f"Missing required CheckM2 field aliases: {', '.join(candidates)}")


def parse_float(value: str, field_name: str, path: Path) -> float:
    """Parse a floating-point field or raise a descriptive error."""
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(
            f"Could not parse {field_name} value {value!r} in {path}."
        ) from error


def parse_report(path: Path) -> ParsedCheckM2Report:
    """Parse the required metrics from a single CheckM2 quality report."""
    row = read_report_rows(path)[0]

    metrics = {
        metric_name: parse_float(
            first_present(row, aliases),
            metric_name,
            path,
        )
        for metric_name, aliases in METRIC_ALIASES.items()
    }

    shared_stats: dict[str, float | str] = {}
    for stat_name, aliases in SHARED_STAT_ALIASES.items():
        value = None
        for alias in aliases:
            if alias in row:
                value = row[alias]
                break
        if value is None or not value.strip():
            continue
        if stat_name == "Name":
            shared_stats[stat_name] = value.strip()
        else:
            shared_stats[stat_name] = parse_float(value, stat_name, path)

    return ParsedCheckM2Report(metrics=metrics, shared_stats=shared_stats)


def format_metric(value: float | None) -> str:
    """Format numeric output fields or emit `NA`."""
    if value is None:
        return "NA"
    if math.isclose(value, round(value), rel_tol=0.0, abs_tol=1e-9):
        return str(int(round(value)))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def reports_have_consistent_shared_stats(
    report_four: ParsedCheckM2Report,
    report_eleven: ParsedCheckM2Report,
) -> bool:
    """Check whether shared assembly statistics match between two reports."""
    required_keys = set(REQUIRED_SHARED_STATS)
    shared_keys_four = set(report_four.shared_stats)
    shared_keys_eleven = set(report_eleven.shared_stats)
    if shared_keys_four != required_keys or shared_keys_eleven != required_keys:
        return False

    for key in REQUIRED_SHARED_STATS:
        left = report_four.shared_stats[key]
        right = report_eleven.shared_stats[key]
        if isinstance(left, str) or isinstance(right, str):
            if left != right:
                return False
            continue
        if not math.isclose(left, right, rel_tol=0.0, abs_tol=1e-9):
            return False
    return True


def empty_output_row(accession: str) -> dict[str, str]:
    """Create a default NA-filled output row."""
    row = {column: "NA" for column in OUTPUT_COLUMNS}
    row["accession"] = accession
    row["warnings"] = ""
    return row


def assign_low_quality(
    row: dict[str, str],
    report: ParsedCheckM2Report,
) -> None:
    """Populate Low_quality from the chosen report metrics."""
    low_quality_score = (
        report.metrics["Completeness"] - 5 * report.metrics["Contamination"]
    )
    row["Low_quality"] = "true" if low_quality_score <= 50 else "false"


def assign_gcode_from_valid_pair(
    row: dict[str, str],
    report_four: ParsedCheckM2Report,
    report_eleven: ParsedCheckM2Report,
    *,
    gcode_rule: str,
    warnings: list[str],
) -> None:
    """Assign gcode from a valid paired CheckM2 comparison."""
    completeness_four = report_four.metrics["Completeness"]
    completeness_eleven = report_eleven.metrics["Completeness"]

    if completeness_four - completeness_eleven > 10:
        row["Gcode"] = "4"
        assign_low_quality(row, report_four)
        return

    if completeness_eleven - completeness_four > 10:
        row["Gcode"] = "11"
        assign_low_quality(row, report_eleven)
        return

    if gcode_rule == DELTA_THEN_ELEVEN_RULE:
        row["Gcode"] = "11"
        assign_low_quality(row, report_eleven)
        return

    row["Gcode"] = "NA"
    row["Low_quality"] = "NA"
    warnings.append("gcode_na")


def build_output_row(
    accession: str,
    report_four: ParsedCheckM2Report | None,
    report_eleven: ParsedCheckM2Report | None,
    *,
    gcode_rule: str,
    warnings: list[str],
) -> dict[str, str]:
    """Create the final output row from two optional parsed reports."""
    row = empty_output_row(accession)

    if report_four is not None:
        row["Completeness_gcode4"] = format_metric(report_four.metrics["Completeness"])
        row["Contamination_gcode4"] = format_metric(report_four.metrics["Contamination"])
        row["Coding_Density_gcode4"] = format_metric(report_four.metrics["Coding_Density"])
        row["Average_Gene_Length_gcode4"] = format_metric(
            report_four.metrics["Average_Gene_Length"]
        )
        row["Total_Coding_Sequences_gcode4"] = format_metric(
            report_four.metrics["Total_Coding_Sequences"]
        )

    if report_eleven is not None:
        row["Completeness_gcode11"] = format_metric(
            report_eleven.metrics["Completeness"]
        )
        row["Contamination_gcode11"] = format_metric(
            report_eleven.metrics["Contamination"]
        )
        row["Coding_Density_gcode11"] = format_metric(
            report_eleven.metrics["Coding_Density"]
        )
        row["Average_Gene_Length_gcode11"] = format_metric(
            report_eleven.metrics["Average_Gene_Length"]
        )
        row["Total_Coding_Sequences_gcode11"] = format_metric(
            report_eleven.metrics["Total_Coding_Sequences"]
        )

    if report_four is not None and report_eleven is not None and reports_have_consistent_shared_stats(
        report_four,
        report_eleven,
    ):
        assign_gcode_from_valid_pair(
            row,
            report_four,
            report_eleven,
            gcode_rule=gcode_rule,
            warnings=warnings,
        )
    else:
        row["Gcode"] = "NA"
        row["Low_quality"] = "NA"

    failure_warnings = {
        "checkm2_gcode4_failed",
        "checkm2_gcode11_failed",
        "inconsistent_shared_stats",
    }
    row["checkm2_status"] = (
        "failed" if any(warning in failure_warnings for warning in warnings) else "done"
    )
    row["warnings"] = ";".join(dict.fromkeys(warnings))
    return row


def write_output(path: Path, row: dict[str, str]) -> None:
    """Write the final single-row TSV output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_COLUMNS), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def run_summary(
    accession: str,
    gcode4_report: Path,
    gcode11_report: Path,
    output: Path,
    *,
    gcode_rule: str,
) -> None:
    """Summarise paired CheckM2 reports into one per-sample TSV row."""
    warnings: list[str] = []

    try:
        report_four = parse_report(gcode4_report)
    except ValueError as error:
        LOGGER.warning(str(error))
        warnings.append("checkm2_gcode4_failed")
        report_four = None

    try:
        report_eleven = parse_report(gcode11_report)
    except ValueError as error:
        LOGGER.warning(str(error))
        warnings.append("checkm2_gcode11_failed")
        report_eleven = None

    if (
        report_four is not None
        and report_eleven is not None
        and not reports_have_consistent_shared_stats(report_four, report_eleven)
    ):
        warnings.append("inconsistent_shared_stats")

    row = build_output_row(
        accession,
        report_four,
        report_eleven,
        gcode_rule=gcode_rule,
        warnings=warnings,
    )
    write_output(output, row)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CheckM2 summarisation CLI."""
    args = parse_args(argv)
    configure_logging()
    run_summary(
        accession=args.accession,
        gcode4_report=args.gcode4_report,
        gcode11_report=args.gcode11_report,
        output=args.output,
        gcode_rule=args.gcode_rule,
    )
    LOGGER.info("Wrote CheckM2 summary for %s.", args.accession)
    return 0


if __name__ == "__main__":
    sys.exit(main())
