#!/usr/bin/env python3
"""Rescue ANI outputs from published per-sample pipeline results."""

from __future__ import annotations

import argparse
import csv
import logging
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import build_fastani_inputs
import build_master_table
import build_sample_status
import cluster_ani
import select_ani_representatives
import summarise_16s
import summarise_busco
import summarise_checkm2
import validate_inputs

LOGGER = logging.getLogger(__name__)

ASSUMED_SOURCE_COMMIT = "595edef"
CLUSTER_TABLE_COLUMNS = ("Accession", "Cluster_ID", "Matrix_Name")
PUBLISHED_STAGED_SOURCE = "published_staged"
GENOME_FASTA_FALLBACK_SOURCE = "genome_fasta_fallback"


class RescueAniError(RuntimeError):
    """Raised when ANI rescue cannot continue safely."""


@dataclass(frozen=True)
class ResolvedSample:
    """Store the resolved inputs for one rescued sample."""

    accession: str
    internal_id: str
    genome_fasta: str
    is_new: str
    assembly_level: str
    staged_fasta: Path
    staged_source: str


class SingleUseLineageAction(argparse.Action):
    """Accept one multi-value BUSCO lineage option and reject repeats."""

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[str] | None,
        option_string: str | None = None,
    ) -> None:
        """Store the provided BUSCO lineages and fail on repeated usage."""
        if getattr(namespace, self.dest, None) is not None:
            parser.error(f"{option_string} may only be supplied once.")
        if values is None:
            parser.error(f"{option_string} requires at least one lineage value.")
        if isinstance(values, str):
            resolved_values = [values]
        else:
            resolved_values = list(values)
        setattr(namespace, self.dest, resolved_values)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Rescue ANI outputs from a published nf-myco results directory."
    )
    parser.add_argument(
        "--source-outdir",
        required=True,
        type=Path,
        help="Published outdir from the failed pipeline run.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Original metadata table used for the failed run.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Fresh rescue outdir for reconstructed ANI artefacts.",
    )
    parser.add_argument(
        "--busco-lineage",
        action=SingleUseLineageAction,
        nargs="+",
        required=True,
        help=(
            "BUSCO lineages from the failed run in original order. Supply the option "
            "once with one or more values; the first lineage remains the primary "
            "ANI-scoring BUSCO column."
        ),
    )
    parser.add_argument(
        "--validated-samples",
        type=Path,
        help="Override source validated_samples.tsv. Defaults to source outdir tables.",
    )
    parser.add_argument(
        "--initial-status",
        type=Path,
        help="Override source validation-time sample_status.tsv seed.",
    )
    parser.add_argument(
        "--gcode-rule",
        choices=summarise_checkm2.GCODE_RULE_CHOICES,
        default=summarise_checkm2.STRICT_DELTA_RULE,
        help="CheckM2 gcode rule used for rescue. Default: strict_delta.",
    )
    parser.add_argument(
        "--ani-threshold",
        type=float,
        default=0.95,
        help="ANI threshold as a fraction in (0,1). Default: 0.95.",
    )
    parser.add_argument(
        "--ani-score-profile",
        choices=("default", "isolate", "mag"),
        default="default",
        help="ANI representative scoring profile. Default: default.",
    )
    parser.add_argument(
        "--fastani-binary",
        default="fastANI",
        help=(
            "FastANI executable or shell-style command prefix. "
            'Examples: "fastANI" or "singularity exec fastani.sif fastANI".'
        ),
    )
    parser.add_argument(
        "--seqtk-binary",
        default="seqtk",
        help=(
            "seqtk executable or shell-style command prefix used by the "
            'assembly-stats helper. Examples: "seqtk" or '
            '"singularity exec seqtk.sif seqtk".'
        ),
    )
    parser.add_argument(
        "--assembly-stats",
        type=Path,
        help=(
            "Override source assembly_stats.tsv. Defaults to "
            "source outdir cohort/assembly_stats/assembly_stats.tsv."
        ),
    )
    parser.add_argument(
        "--recalculate-assembly-stats",
        action="store_true",
        help=(
            "Recalculate assembly stats from staged FASTA files instead of "
            "reusing the published source table."
        ),
    )
    parser.add_argument(
        "--assembly-stats-jobs",
        type=int,
        default=1,
        help="Worker count for rescued assembly-stat calculation. Default: 1.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        force=True,
    )


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header list and row dictionaries."""
    if not path.is_file():
        raise RescueAniError(f"Missing TSV file: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RescueAniError(f"TSV file is empty or missing a header: {path}")
        header = [column.strip() for column in reader.fieldnames]
        rows = [{key: str(value or "").strip() for key, value in row.items()} for row in reader]
    return header, rows


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write a TSV file with a fixed header order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(header),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in header})


def read_single_row(path: Path) -> dict[str, str]:
    """Read a TSV expected to contain exactly one data row."""
    _header, rows = read_tsv_rows(path)
    if len(rows) != 1:
        raise RescueAniError(f"Expected exactly one row in {path}, found {len(rows)}.")
    return rows[0]


def copy_source_table(source: Path, destination: Path) -> Path:
    """Copy one source table into the rescue outdir."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def count_tsv_rows(path: Path) -> int:
    """Return the number of data rows in one TSV."""
    _header, rows = read_tsv_rows(path)
    return len(rows)


def format_accession_list(accessions: Sequence[str]) -> str:
    """Format accession-like identifiers for concise rescue messages."""
    if len(accessions) <= 5:
        return ", ".join(accessions)
    preview = ", ".join(accessions[:5])
    return f"{preview}, ... (+{len(accessions) - 5} more)"


def format_count_pairs(pairs: Sequence[tuple[str, int]]) -> str:
    """Format ordered key-count pairs for one log line."""
    if not pairs:
        return "none"
    return "; ".join(f"{key}={count}" for key, count in pairs)


def resolve_validated_samples_path(args: argparse.Namespace) -> Path:
    """Resolve the validated-samples source path."""
    return args.validated_samples or args.source_outdir / "tables" / "validated_samples.tsv"


def resolve_validation_warnings_path(source_outdir: Path) -> Path:
    """Resolve the validation-warnings source path."""
    return source_outdir / "tables" / "validation_warnings.tsv"


def resolve_source_assembly_stats_path(args: argparse.Namespace) -> Path:
    """Resolve the source assembly-stats table used by rescue."""
    return (
        args.assembly_stats
        or args.source_outdir / "cohort" / "assembly_stats" / "assembly_stats.tsv"
    )


def validate_outdirs(source_outdir: Path, rescue_outdir: Path) -> None:
    """Validate the source and destination outdir contract."""
    if not source_outdir.is_dir():
        raise RescueAniError(f"Source outdir is missing or not a directory: {source_outdir}")
    if source_outdir.resolve() == rescue_outdir.resolve():
        raise RescueAniError("--outdir must differ from --source-outdir for rescue runs.")


def load_validated_samples(path: Path) -> list[dict[str, str]]:
    """Load and validate the rescued sample manifest."""
    header, rows = read_tsv_rows(path)
    build_master_table.require_columns(
        header,
        build_master_table.VALIDATED_SAMPLE_REQUIRED_COLUMNS,
        "validated_samples.tsv",
    )
    if not rows:
        raise RescueAniError(f"validated_samples.tsv is empty: {path}")
    return rows


def resolve_staged_fasta(source_outdir: Path, sample_row: dict[str, str]) -> tuple[Path, str]:
    """Resolve the staged FASTA for one sample or fall back to genome_fasta."""
    accession = sample_row["accession"]
    internal_id = sample_row["internal_id"]
    staged_candidate = (
        source_outdir / "samples" / accession / "staged" / f"{internal_id}.fasta"
    )
    if staged_candidate.is_file():
        return staged_candidate.resolve(), PUBLISHED_STAGED_SOURCE

    fallback = Path(sample_row["genome_fasta"]).expanduser()
    if fallback.is_file():
        return fallback.resolve(), GENOME_FASTA_FALLBACK_SOURCE

    raise RescueAniError(
        "Could not resolve staged FASTA for accession "
        f"{accession!r}. Expected {staged_candidate} or existing genome_fasta {fallback}."
    )


def build_resolved_samples(
    source_outdir: Path,
    validated_rows: Sequence[dict[str, str]],
) -> list[ResolvedSample]:
    """Build the resolved sample list used by all rescue stages."""
    resolved: list[ResolvedSample] = []
    for row in validated_rows:
        staged_fasta, staged_source = resolve_staged_fasta(source_outdir, row)
        resolved.append(
            ResolvedSample(
                accession=row["accession"],
                internal_id=row["internal_id"],
                genome_fasta=row["genome_fasta"],
                is_new=row["is_new"],
                assembly_level=row["assembly_level"],
                staged_fasta=staged_fasta,
                staged_source=staged_source,
            )
        )
    return resolved


def log_startup_configuration(
    args: argparse.Namespace,
    *,
    validated_samples_source: Path,
    validation_warnings_source: Path,
    source_assembly_stats: Path,
) -> None:
    """Log the resolved rescue configuration before execution."""
    LOGGER.info(
        "Rescue configuration: source_outdir=%s rescue_outdir=%s metadata=%s.",
        args.source_outdir,
        args.outdir,
        args.metadata,
    )
    if args.initial_status is not None:
        LOGGER.info(
            "Rescue sources: validated_samples=%s initial_status_override=%s.",
            validated_samples_source,
            args.initial_status,
        )
    else:
        LOGGER.info(
            "Rescue sources: validated_samples=%s validation_warnings=%s initial_status=rebuild.",
            validated_samples_source,
            validation_warnings_source,
        )
    LOGGER.info(
        "Rescue options: busco_lineages=%s gcode_rule=%s ani_threshold=%s ani_score_profile=%s.",
        ",".join(args.busco_lineage),
        args.gcode_rule,
        args.ani_threshold,
        args.ani_score_profile,
    )
    LOGGER.info(
        "Rescue assembly stats: source=%s mode=%s.",
        source_assembly_stats,
        "recalculate" if args.recalculate_assembly_stats else "copy",
    )
    LOGGER.info("Rescue tools: fastani_binary=%s.", args.fastani_binary)
    LOGGER.info(
        "Rescue assembly-stats tools: seqtk_binary=%s assembly_stats_jobs=%d%s.",
        args.seqtk_binary,
        args.assembly_stats_jobs,
        "" if args.recalculate_assembly_stats else " (ignored unless --recalculate-assembly-stats)",
    )


def log_resolved_samples(samples: Sequence[ResolvedSample]) -> None:
    """Log the resolved staged FASTA origin for each sample."""
    fallback_count = 0
    for sample in samples:
        if sample.staged_source == GENOME_FASTA_FALLBACK_SOURCE:
            fallback_count += 1
        LOGGER.info(
            "Resolved sample %s: internal_id=%s staged_fasta=%s source=%s.",
            sample.accession,
            sample.internal_id,
            sample.staged_fasta,
            sample.staged_source,
        )
    LOGGER.info(
        "Resolved %d sample(s); genome_fasta_fallback=%d.",
        len(samples),
        fallback_count,
    )


def load_validation_warnings(
    path: Path,
) -> list[validate_inputs.ValidationWarning]:
    """Load published validation warnings when available."""
    if not path.is_file():
        LOGGER.warning(
            "Validation warnings table is missing at %s; rebuilding the initial-status seed without validation warnings.",
            path,
        )
        return []

    header, rows = read_tsv_rows(path)
    build_master_table.require_columns(
        header,
        validate_inputs.VALIDATION_WARNING_COLUMNS,
        "validation_warnings.tsv",
    )
    return [
        validate_inputs.ValidationWarning(
            accession=row["accession"],
            warning_code=row["warning_code"],
            message=row["message"],
        )
        for row in rows
    ]


def rebuild_initial_status_seed(
    *,
    validated_rows: Sequence[dict[str, str]],
    validation_warnings_path: Path,
    busco_lineages: Sequence[str],
    output_path: Path,
) -> Path:
    """Rebuild the validation-time sample-status seed for rescue."""
    sample_status_columns = validate_inputs.resolve_sample_status_columns(
        busco_lineages=busco_lineages,
    )
    warnings = load_validation_warnings(validation_warnings_path)
    records = [
        validate_inputs.SampleRecord(
            values=dict(row),
            internal_id=row["internal_id"],
            metadata_present=False,
        )
        for row in validated_rows
    ]
    rows = validate_inputs.build_initial_sample_status_rows(
        records=records,
        warnings=warnings,
        sample_status_columns=sample_status_columns,
    )
    write_tsv(output_path, sample_status_columns, rows)
    LOGGER.info(
        "Rebuilt rescue initial-status seed from validated_samples and validation_warnings: warnings_source=%s output=%s row_count=%d.",
        validation_warnings_path,
        output_path,
        len(rows),
    )
    return output_path


def prepare_initial_status_seed(
    *,
    args: argparse.Namespace,
    validated_rows: Sequence[dict[str, str]],
    validation_warnings_path: Path,
    output_path: Path,
) -> Path:
    """Write the initial-status seed that rescue will pass into final reporting."""
    if args.initial_status is not None:
        copied_path = copy_source_table(args.initial_status, output_path)
        LOGGER.info(
            "Using explicit rescue initial-status override: source=%s output=%s.",
            args.initial_status,
            copied_path,
        )
        return copied_path

    return rebuild_initial_status_seed(
        validated_rows=validated_rows,
        validation_warnings_path=validation_warnings_path,
        busco_lineages=args.busco_lineage,
        output_path=output_path,
    )


def recover_sixteen_s(
    source_outdir: Path,
    rescue_outdir: Path,
    samples: Sequence[ResolvedSample],
) -> Path:
    """Rebuild per-sample and combined 16S summaries from published Barrnap outputs."""
    rows: list[dict[str, str]] = []
    for sample in samples:
        source_dir = source_outdir / "samples" / sample.accession / "barrnap"
        output_dir = rescue_outdir / "samples" / sample.accession / "16s"
        summarise_16s.summarise_hits(
            accession=sample.accession,
            rrna_gff=source_dir / "rrna.gff",
            rrna_fasta=source_dir / "rrna.fa",
            outdir=output_dir,
        )
        status_path = output_dir / "16S_status.tsv"
        status_row = read_single_row(status_path)
        rows.append(status_row)
        LOGGER.info(
            "Recovered 16S for %s: barrnap_dir=%s status=%s output=%s.",
            sample.accession,
            source_dir,
            status_row.get("16S", "NA"),
            status_path,
        )

    combined_path = rescue_outdir / "tables" / "16s_statuses.tsv"
    write_tsv(combined_path, summarise_16s.STATUS_COLUMNS, rows)
    LOGGER.info("Wrote combined 16S summary to %s with %d row(s).", combined_path, len(rows))
    return combined_path


def recover_checkm2(
    source_outdir: Path,
    rescue_outdir: Path,
    samples: Sequence[ResolvedSample],
    *,
    gcode_rule: str,
) -> Path:
    """Rebuild per-sample and combined CheckM2 summaries from published reports."""
    rows: list[dict[str, str]] = []
    for sample in samples:
        output_dir = rescue_outdir / "samples" / sample.accession / "checkm2"
        output_path = output_dir / "checkm2_summary.tsv"
        summarise_checkm2.run_summary(
            accession=sample.accession,
            gcode4_report=(
                source_outdir
                / "samples"
                / sample.accession
                / "checkm2_gcode4"
                / "quality_report.tsv"
            ),
            gcode11_report=(
                source_outdir
                / "samples"
                / sample.accession
                / "checkm2_gcode11"
                / "quality_report.tsv"
            ),
            output=output_path,
            gcode_rule=gcode_rule,
        )
        summary_row = read_single_row(output_path)
        rows.append(summary_row)
        LOGGER.info(
            "Recovered CheckM2 for %s: gcode_rule=%s Gcode=%s Low_quality=%s output=%s.",
            sample.accession,
            gcode_rule,
            summary_row.get("Gcode", "NA"),
            summary_row.get("Low_quality", "NA"),
            output_path,
        )

    combined_path = rescue_outdir / "tables" / "checkm2_summaries.tsv"
    write_tsv(combined_path, summarise_checkm2.OUTPUT_COLUMNS, rows)
    LOGGER.info(
        "Wrote combined CheckM2 summary to %s with %d row(s).",
        combined_path,
        len(rows),
    )
    return combined_path


def recover_busco_lineage(
    source_outdir: Path,
    rescue_outdir: Path,
    samples: Sequence[ResolvedSample],
    *,
    lineage: str,
) -> Path:
    """Rebuild one per-lineage BUSCO summary table."""
    rows: list[dict[str, str]] = []
    busco_column = summarise_busco.busco_column_name(lineage)
    header = ("accession", "lineage", busco_column, "busco_status", "warnings")
    LOGGER.info("Rebuilding BUSCO lineage %s for %d sample(s).", lineage, len(samples))

    for sample in samples:
        summary_path = (
            source_outdir
            / "samples"
            / sample.accession
            / "busco"
            / lineage
            / "short_summary.json"
        )
        output_path = (
            rescue_outdir
            / "samples"
            / sample.accession
            / "busco"
            / lineage
            / f"busco_summary_{lineage}.tsv"
        )
        try:
            busco_summary = summarise_busco.parse_summary(summary_path)
            busco_status = "done"
            warnings = ""
        except ValueError as error:
            LOGGER.warning(
                "BUSCO rescue fell back to NA for %s (%s): %s",
                sample.accession,
                lineage,
                error,
            )
            busco_summary = "NA"
            busco_status = "failed"
            warnings = "busco_summary_failed"

        summarise_busco.write_output(
            accession=sample.accession,
            lineage=lineage,
            busco_summary=busco_summary,
            busco_status=busco_status,
            warnings=warnings,
            output=output_path,
        )
        summary_row = read_single_row(output_path)
        rows.append(summary_row)
        LOGGER.info(
            "Recovered BUSCO for %s: lineage=%s value=%s status=%s output=%s.",
            sample.accession,
            lineage,
            summary_row.get(busco_column, "NA"),
            summary_row.get("busco_status", "NA"),
            output_path,
        )

    combined_path = rescue_outdir / "tables" / f"busco_summary_{lineage}.tsv"
    write_tsv(combined_path, header, rows)
    done_count = sum(1 for row in rows if row.get("busco_status", "") == "done")
    failed_count = sum(1 for row in rows if row.get("busco_status", "") == "failed")
    LOGGER.info(
        "Wrote combined BUSCO summary for %s to %s with %d row(s); done=%d failed=%d.",
        lineage,
        combined_path,
        len(rows),
        done_count,
        failed_count,
    )
    return combined_path


def write_staged_manifest(
    rescue_outdir: Path,
    samples: Sequence[ResolvedSample],
) -> Path:
    """Write the staged-manifest TSV consumed by assembly-stats and ANI prep."""
    manifest_path = rescue_outdir / "cohort" / "staged_genomes.tsv"
    rows = [
        {
            "accession": sample.accession,
            "internal_id": sample.internal_id,
            "staged_filename": str(sample.staged_fasta),
        }
        for sample in samples
    ]
    write_tsv(manifest_path, ("accession", "internal_id", "staged_filename"), rows)
    LOGGER.info("Wrote staged manifest to %s with %d row(s).", manifest_path, len(rows))
    return manifest_path


def run_subprocess(
    command: Sequence[str | Path],
    *,
    cwd: Path,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run one subprocess command and capture text output."""
    try:
        return subprocess.run(
            [str(part) for part in command],
            cwd=cwd,
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError as error:
        raise RescueAniError(f"Command not found while running: {command[0]}") from error


def parse_command_prefix(command: str) -> list[str]:
    """Split one shell-style command prefix into argv tokens."""
    try:
        tokens = shlex.split(command)
    except ValueError as error:
        raise RescueAniError(f"Could not parse command prefix {command!r}.") from error
    if not tokens:
        raise RescueAniError("Command prefix must not be empty.")
    return tokens


def build_tool_wrapper_env(
    command: str,
    *,
    executable_name: str,
    wrapper_dir: Path,
) -> dict[str, str] | None:
    """Return a PATH-prefixed environment exposing one wrapped executable."""
    tokens = parse_command_prefix(command)
    if len(tokens) == 1 and tokens[0] == executable_name:
        return None

    wrapper_dir.mkdir(parents=True, exist_ok=True)
    wrapper_path = wrapper_dir / executable_name
    quoted_command = " ".join(shlex.quote(token) for token in tokens)
    wrapper_path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f'exec {quoted_command} "$@"\n',
        encoding="utf-8",
    )
    wrapper_path.chmod(0o755)

    env = os.environ.copy()
    current_path = env.get("PATH", "")
    env["PATH"] = f"{wrapper_dir}{os.pathsep}{current_path}" if current_path else str(wrapper_dir)
    return env


def validate_source_assembly_stats(
    path: Path,
    *,
    validated_accessions: Sequence[str],
) -> int:
    """Validate a source assembly-stats table before rescue reuses it."""
    validated_set = set(validated_accessions)
    try:
        assembly_stats_index = build_master_table.load_assembly_stats_index(
            path,
            validated_set,
        )
    except build_master_table.MasterTableError as error:
        raise RescueAniError(
            "Rescue source assembly-stats table is unavailable or invalid at "
            f"{path}: {error}. Use --recalculate-assembly-stats to rebuild it."
        ) from error

    missing_accessions = sorted(validated_set - set(assembly_stats_index))
    if missing_accessions:
        raise RescueAniError(
            "Rescue source assembly-stats table is incomplete at "
            f"{path}: missing validated accession rows for "
            f"{format_accession_list(missing_accessions)}. "
            "Use --recalculate-assembly-stats to rebuild it."
        )

    missing_metrics: list[str] = []
    for accession in sorted(validated_set):
        row = assembly_stats_index[accession]
        missing_columns = [
            column
            for column in build_master_table.ASSEMBLY_STATS_COLUMNS
            if build_master_table.is_missing(row.get(column, ""))
        ]
        if missing_columns:
            missing_metrics.append(f"{accession}[{','.join(missing_columns)}]")

    if missing_metrics:
        raise RescueAniError(
            "Rescue source assembly-stats table is incomplete at "
            f"{path}: missing required metric values for "
            f"{format_accession_list(missing_metrics)}. "
            "Use --recalculate-assembly-stats to rebuild it."
        )

    return len(assembly_stats_index)


def copy_source_assembly_stats(
    source_path: Path,
    output_path: Path,
    *,
    validated_accessions: Sequence[str],
) -> Path:
    """Validate and copy source assembly stats into the rescue outdir."""
    row_count = validate_source_assembly_stats(
        source_path,
        validated_accessions=validated_accessions,
    )
    copied_path = copy_source_table(source_path, output_path)
    LOGGER.info(
        "Copied rescue assembly stats from source table: source=%s output=%s row_count=%d.",
        source_path,
        copied_path,
        row_count,
    )
    return copied_path


def run_calculate_assembly_stats(
    staged_manifest: Path,
    output_path: Path,
    *,
    seqtk_binary: str,
    jobs: int,
) -> Path:
    """Run the existing assembly-stat helper over the resolved staged genomes."""
    script_path = Path(__file__).resolve().parent / "calculate_assembly_stats.sh"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info(
        "Running assembly stats: jobs=%d seqtk_binary=%s staged_manifest=%s output=%s.",
        jobs,
        seqtk_binary,
        staged_manifest,
        output_path,
    )
    with tempfile.TemporaryDirectory(prefix="seqtk_wrapper_", dir=output_path.parent) as tmpdir_name:
        env = build_tool_wrapper_env(
            seqtk_binary,
            executable_name="seqtk",
            wrapper_dir=Path(tmpdir_name),
        )
        result = run_subprocess(
            [
                script_path,
                "--staged-manifest",
                staged_manifest,
                "--jobs",
                str(jobs),
                "--output",
                output_path,
            ],
            cwd=Path(__file__).resolve().parents[1],
            env=env,
        )
    if result.returncode != 0:
        error_output = result.stderr.strip() or result.stdout.strip()
        raise RescueAniError(
            "calculate_assembly_stats.sh failed: "
            + (error_output or f"exit code {result.returncode}")
        )
    if not output_path.is_file():
        raise RescueAniError(f"Assembly stats output was not created: {output_path}")
    LOGGER.info("Wrote assembly stats to %s with %d row(s).", output_path, count_tsv_rows(output_path))
    return output_path


def prepare_rescue_assembly_stats(
    *,
    args: argparse.Namespace,
    staged_manifest: Path,
    samples: Sequence[ResolvedSample],
    source_assembly_stats: Path,
    output_path: Path,
) -> Path:
    """Prepare the rescued assembly-stats table by copy or recalculation."""
    validated_accessions = [sample.accession for sample in samples]
    if args.recalculate_assembly_stats:
        assembly_stats_path = run_calculate_assembly_stats(
            staged_manifest,
            output_path,
            seqtk_binary=args.seqtk_binary,
            jobs=args.assembly_stats_jobs,
        )
        LOGGER.info(
            "Rescue assembly stats ready: mode=recalculated output=%s row_count=%d.",
            assembly_stats_path,
            count_tsv_rows(assembly_stats_path),
        )
        return assembly_stats_path

    assembly_stats_path = copy_source_assembly_stats(
        source_assembly_stats,
        output_path,
        validated_accessions=validated_accessions,
    )
    LOGGER.info(
        "Rescue assembly stats ready: mode=copied output=%s row_count=%d.",
        assembly_stats_path,
        count_tsv_rows(assembly_stats_path),
    )
    return assembly_stats_path


def table_has_rows(path: Path) -> bool:
    """Return True when a TSV contains at least one data row."""
    _header, rows = read_tsv_rows(path)
    return bool(rows)


def normalise_fastani_outputs(fastani_dir: Path) -> None:
    """Normalise FastANI matrix output to the stable rescue filenames."""
    matrix_output = fastani_dir / "fastani.matrix"
    if matrix_output.exists():
        matrix_output.unlink()

    candidate_paths = [fastani_dir / "fastani.tsv.matrix"]
    candidate_paths.extend(
        sorted(
            path
            for path in fastani_dir.glob("*.matrix")
            if path.name != "fastani.matrix"
        )
    )
    candidate = next((path for path in candidate_paths if path.is_file()), None)
    if candidate is None:
        matrix_output.write_text("", encoding="utf-8")
        return
    shutil.copy2(candidate, matrix_output)


def summarise_ani_exclusions(
    path: Path,
) -> tuple[list[dict[str, str]], int, int, list[tuple[str, int]]]:
    """Return ANI inclusion counts and grouped exclusion reasons."""
    _header, rows = read_tsv_rows(path)
    included_count = 0
    excluded_count = 0
    reason_counts: dict[str, int] = {}
    for row in rows:
        if row.get("ani_included", "").strip() == "true":
            included_count += 1
            continue
        excluded_count += 1
        reason = row.get("ani_exclusion_reason", "").strip() or "unspecified"
        reason_counts[reason] = reason_counts.get(reason, 0) + 1
    return rows, included_count, excluded_count, sorted(reason_counts.items())


def log_ani_preparation(fastani_dir: Path) -> None:
    """Log ANI inclusion, exclusion, and reason summaries."""
    ani_metadata_path = fastani_dir / "ani_metadata.tsv"
    ani_exclusions_path = fastani_dir / "ani_exclusions.tsv"
    exclusion_rows, included_count, excluded_count, reason_counts = summarise_ani_exclusions(
        ani_exclusions_path
    )
    LOGGER.info(
        "Built ANI inputs under %s: ani_metadata_rows=%d ani_exclusion_rows=%d.",
        fastani_dir,
        count_tsv_rows(ani_metadata_path),
        len(exclusion_rows),
    )
    LOGGER.info(
        "ANI eligibility summary: included=%d excluded=%d.",
        included_count,
        excluded_count,
    )
    for row in exclusion_rows:
        reason = row.get("ani_exclusion_reason", "").strip() or "none"
        LOGGER.info(
            "ANI gating for %s: included=%s reason=%s.",
            row.get("accession", ""),
            row.get("ani_included", ""),
            reason,
        )
    LOGGER.info("ANI exclusion reasons: %s.", format_count_pairs(reason_counts))


def summarise_cluster_outputs(cluster_path: Path, representative_path: Path) -> tuple[int, int, int]:
    """Return cluster row, unique cluster, and representative counts."""
    _cluster_header, cluster_rows = read_tsv_rows(cluster_path)
    _representative_header, representative_rows = read_tsv_rows(representative_path)
    cluster_count = len(
        {row.get("Cluster_ID", "").strip() for row in cluster_rows if row.get("Cluster_ID", "").strip()}
    )
    return len(cluster_rows), cluster_count, len(representative_rows)


def log_cluster_outputs(cluster_path: Path, ani_summary_path: Path, representative_path: Path) -> None:
    """Log cluster and representative output sizes."""
    cluster_rows, cluster_count, representative_count = summarise_cluster_outputs(
        cluster_path,
        representative_path,
    )
    LOGGER.info(
        "ANI cluster outputs: cluster_path=%s ani_summary=%s cluster_rows=%d cluster_count=%d representative_count=%d.",
        cluster_path,
        ani_summary_path,
        cluster_rows,
        cluster_count,
        representative_count,
    )
    LOGGER.info(
        "ANI representatives ready: output=%s count=%d.",
        representative_path,
        representative_count,
    )


def run_fastani(
    *,
    fastani_binary: str,
    fastani_dir: Path,
    ani_threshold: float,
) -> Path:
    """Run FastANI over the rescued eligible genomes and normalise the matrix output."""
    _ = ani_threshold
    fastani_dir.mkdir(parents=True, exist_ok=True)
    path_list = fastani_dir / "fastani_paths.txt"
    fastani_tsv = fastani_dir / "fastani.tsv"
    fastani_log = fastani_dir / "fastani.log"
    matrix_output = fastani_dir / "fastani.matrix"
    fastani_command = parse_command_prefix(fastani_binary)
    LOGGER.info(
        "Running FastANI: command=%s cwd=%s path_list=%s.",
        shlex.join(fastani_command),
        fastani_dir,
        path_list,
    )

    result = run_subprocess(
        [
            *fastani_command,
            "--rl",
            path_list,
            "--ql",
            path_list,
            "--matrix",
            "-t",
            "1",
            "-o",
            fastani_tsv,
        ],
        cwd=fastani_dir,
    )
    fastani_log.write_text((result.stdout or "") + (result.stderr or ""), encoding="utf-8")
    if not fastani_tsv.exists():
        fastani_tsv.write_text("", encoding="utf-8")

    normalise_fastani_outputs(fastani_dir)

    if "Could not open " in fastani_log.read_text(encoding="utf-8"):
        raise RescueAniError(
            "FastANI could not open one or more rescued staged inputs. "
            f"See {fastani_log}."
        )
    if path_list.read_text(encoding="utf-8").strip() and not matrix_output.read_text(
        encoding="utf-8"
    ).strip():
        raise RescueAniError(
            "FastANI did not produce a matrix for the rescued eligible inputs. "
            f"See {fastani_log}."
        )
    LOGGER.info(
        "FastANI outputs ready: tsv=%s log=%s matrix=%s.",
        fastani_tsv,
        fastani_log,
        matrix_output,
    )
    return matrix_output


def write_empty_ani_outputs(cluster_dir: Path) -> tuple[Path, Path, Path]:
    """Write empty ANI outputs when no rescued samples are eligible."""
    cluster_path = cluster_dir / "cluster.tsv"
    ani_summary_path = cluster_dir / "ani_summary.tsv"
    ani_representatives_path = cluster_dir / "ani_representatives.tsv"
    write_tsv(cluster_path, CLUSTER_TABLE_COLUMNS, [])
    write_tsv(ani_summary_path, select_ani_representatives.ANI_SUMMARY_COLUMNS, [])
    write_tsv(
        ani_representatives_path,
        select_ani_representatives.ANI_REPRESENTATIVE_COLUMNS,
        [],
    )
    return cluster_path, ani_summary_path, ani_representatives_path


def run_cluster_and_representatives(
    *,
    fastani_dir: Path,
    cluster_dir: Path,
    ani_threshold: float,
    ani_score_profile: str,
) -> tuple[Path, Path, Path]:
    """Run clustering and representative selection from rescued ANI inputs."""
    cluster_dir.mkdir(parents=True, exist_ok=True)
    cluster_path = cluster_dir / "cluster.tsv"
    ani_summary_path = cluster_dir / "ani_summary.tsv"
    ani_representatives_path = cluster_dir / "ani_representatives.tsv"

    cluster_exit_code = cluster_ani.main(
        [
            "--ani-matrix",
            str(fastani_dir / "fastani.matrix"),
            "--ani-metadata",
            str(fastani_dir / "ani_metadata.tsv"),
            "--threshold",
            str(ani_threshold),
            "--outdir",
            str(cluster_dir),
        ]
    )
    if cluster_exit_code != 0:
        raise RescueAniError("cluster_ani.py failed during rescue.")

    representative_exit_code = select_ani_representatives.main(
        [
            "--ani-clusters",
            str(cluster_path),
            "--ani-metadata",
            str(fastani_dir / "ani_metadata.tsv"),
            "--ani-matrix",
            str(fastani_dir / "fastani.matrix"),
            "--ani-score-profile",
            ani_score_profile,
            "--ani-summary-output",
            str(ani_summary_path),
            "--ani-representatives-output",
            str(ani_representatives_path),
        ]
    )
    if representative_exit_code != 0:
        raise RescueAniError("select_ani_representatives.py failed during rescue.")
    return cluster_path, ani_summary_path, ani_representatives_path


def build_partial_master_table(
    *,
    validated_samples: Path,
    metadata: Path,
    busco_lineages: Sequence[str],
    checkm2_summary: Path,
    sixteen_s_summary: Path,
    busco_summaries: Sequence[Path],
    ani_summary: Path,
    assembly_stats: Path,
    output_path: Path,
) -> Path:
    """Build the partial rescued master table."""
    command = [
        "--validated-samples",
        str(validated_samples),
        "--metadata",
        str(metadata),
        "--checkm2",
        str(checkm2_summary),
        "--16s-status",
        str(sixteen_s_summary),
        "--ani",
        str(ani_summary),
        "--assembly-stats",
        str(assembly_stats),
        "--output",
        str(output_path),
    ]
    for lineage in busco_lineages:
        command.extend(["--busco-lineage", lineage])
    for busco_summary in busco_summaries:
        command.extend(["--busco", str(busco_summary)])

    exit_code = build_master_table.main(command)
    if exit_code != 0:
        raise RescueAniError("build_master_table.py failed during rescue.")
    LOGGER.info("Rescued master table ready: output=%s.", output_path)
    return output_path


def build_partial_sample_status(
    *,
    validated_samples: Path,
    initial_status: Path,
    metadata: Path,
    busco_lineages: Sequence[str],
    checkm2_summary: Path,
    sixteen_s_summary: Path,
    busco_summaries: Sequence[Path],
    ani_summary: Path,
    assembly_stats: Path,
    output_path: Path,
) -> Path:
    """Build the partial rescued sample-status table."""
    command = [
        "--validated-samples",
        str(validated_samples),
        "--initial-status",
        str(initial_status),
        "--metadata",
        str(metadata),
        "--checkm2",
        str(checkm2_summary),
        "--16s-status",
        str(sixteen_s_summary),
        "--ani",
        str(ani_summary),
        "--assembly-stats",
        str(assembly_stats),
        "--primary-busco-column",
        summarise_busco.primary_busco_column(busco_lineages),
        "--output",
        str(output_path),
    ]
    for lineage in busco_lineages:
        command.extend(["--busco-lineage", lineage])
    for busco_summary in busco_summaries:
        command.extend(["--busco", str(busco_summary)])

    exit_code = build_sample_status.main(command)
    if exit_code != 0:
        raise RescueAniError("build_sample_status.py failed during rescue.")
    LOGGER.info("Rescued sample status ready: output=%s.", output_path)
    return output_path


def write_rescue_provenance(
    *,
    rescue_outdir: Path,
    source_outdir: Path,
    metadata: Path,
    busco_lineages: Sequence[str],
    gcode_rule: str,
    ani_threshold: float,
    ani_score_profile: str,
    sample_count: int,
    ani_exclusions_path: Path,
    cluster_path: Path,
    representative_path: Path,
) -> Path:
    """Write a one-row rescue provenance table."""
    ani_exclusion_header, ani_exclusion_rows = read_tsv_rows(ani_exclusions_path)
    _ = ani_exclusion_header
    _cluster_header, cluster_rows = read_tsv_rows(cluster_path)
    _representative_header, representative_rows = read_tsv_rows(representative_path)

    included_samples = sum(
        1 for row in ani_exclusion_rows if row.get("ani_included", "") == "true"
    )
    excluded_samples = sum(
        1 for row in ani_exclusion_rows if row.get("ani_included", "") == "false"
    )
    provenance_path = rescue_outdir / "tables" / "rescue_provenance.tsv"
    header = (
        "source_outdir",
        "source_commit",
        "metadata_path",
        "busco_lineages",
        "gcode_rule",
        "ani_threshold",
        "ani_score_profile",
        "validated_sample_count",
        "ani_included_sample_count",
        "ani_excluded_sample_count",
        "cluster_row_count",
        "representative_row_count",
    )
    write_tsv(
        provenance_path,
        header,
        [
            {
                "source_outdir": str(source_outdir.resolve()),
                "source_commit": ASSUMED_SOURCE_COMMIT,
                "metadata_path": str(metadata.resolve()),
                "busco_lineages": ";".join(busco_lineages),
                "gcode_rule": gcode_rule,
                "ani_threshold": str(ani_threshold),
                "ani_score_profile": ani_score_profile,
                "validated_sample_count": str(sample_count),
                "ani_included_sample_count": str(included_samples),
                "ani_excluded_sample_count": str(excluded_samples),
                "cluster_row_count": str(len(cluster_rows)),
                "representative_row_count": str(len(representative_rows)),
            }
        ],
    )
    LOGGER.info(
        "Rescue provenance ready: output=%s validated=%d ani_included=%d ani_excluded=%d cluster_rows=%d representative_rows=%d.",
        provenance_path,
        sample_count,
        included_samples,
        excluded_samples,
        len(cluster_rows),
        len(representative_rows),
    )
    return provenance_path


def run_rescue(args: argparse.Namespace) -> None:
    """Run the full ANI rescue workflow."""
    validate_outdirs(args.source_outdir, args.outdir)
    validated_samples_source = resolve_validated_samples_path(args)
    validation_warnings_source = resolve_validation_warnings_path(args.source_outdir)
    source_assembly_stats = resolve_source_assembly_stats_path(args)
    log_startup_configuration(
        args,
        validated_samples_source=validated_samples_source,
        validation_warnings_source=validation_warnings_source,
        source_assembly_stats=source_assembly_stats,
    )
    if not args.metadata.is_file():
        raise RescueAniError(f"Metadata table does not exist: {args.metadata}")

    validated_rows = load_validated_samples(validated_samples_source)
    LOGGER.info(
        "Loaded %d validated sample(s) from %s.",
        len(validated_rows),
        validated_samples_source,
    )
    samples = build_resolved_samples(args.source_outdir, validated_rows)
    log_resolved_samples(samples)

    tables_dir = args.outdir / "tables"
    copied_validated_samples = copy_source_table(
        validated_samples_source, tables_dir / "validated_samples.tsv"
    )
    copied_initial_status = prepare_initial_status_seed(
        args=args,
        validated_rows=validated_rows,
        validation_warnings_path=validation_warnings_source,
        output_path=tables_dir / "source_sample_status.tsv",
    )

    LOGGER.info("Rebuilding per-sample 16S summaries for %d sample(s).", len(samples))
    sixteen_s_summary = recover_sixteen_s(args.source_outdir, args.outdir, samples)

    LOGGER.info("Rebuilding per-sample CheckM2 summaries.")
    checkm2_summary = recover_checkm2(
        args.source_outdir,
        args.outdir,
        samples,
        gcode_rule=args.gcode_rule,
    )

    LOGGER.info("Rebuilding per-lineage BUSCO summaries.")
    busco_summary_paths = [
        recover_busco_lineage(
            args.source_outdir,
            args.outdir,
            samples,
            lineage=lineage,
        )
        for lineage in args.busco_lineage
    ]

    LOGGER.info("Reconstructing staged manifest.")
    staged_manifest = write_staged_manifest(args.outdir, samples)
    assembly_stats = prepare_rescue_assembly_stats(
        args=args,
        staged_manifest=staged_manifest,
        samples=samples,
        source_assembly_stats=source_assembly_stats,
        output_path=args.outdir / "cohort" / "assembly_stats" / "assembly_stats.tsv",
    )

    LOGGER.info("Building rescued ANI eligibility and FastANI inputs.")
    fastani_dir = args.outdir / "cohort" / "fastani"
    build_fastani_inputs.run_build_fastani_inputs(
        validated_samples=copied_validated_samples,
        metadata=args.metadata,
        staged_manifest=staged_manifest,
        checkm2=checkm2_summary,
        sixteen_s_status=sixteen_s_summary,
        busco=busco_summary_paths,
        primary_busco_column=summarise_busco.primary_busco_column(args.busco_lineage),
        assembly_stats=assembly_stats,
        outdir=fastani_dir,
    )
    log_ani_preparation(fastani_dir)

    cluster_dir = args.outdir / "cohort" / "ani_clusters"
    if table_has_rows(fastani_dir / "ani_metadata.tsv"):
        LOGGER.info("Running FastANI, clustering, and representative selection.")
        run_fastani(
            fastani_binary=args.fastani_binary,
            fastani_dir=fastani_dir,
            ani_threshold=args.ani_threshold,
        )
        cluster_path, ani_summary, ani_representatives = run_cluster_and_representatives(
            fastani_dir=fastani_dir,
            cluster_dir=cluster_dir,
            ani_threshold=args.ani_threshold,
            ani_score_profile=args.ani_score_profile,
        )
    else:
        LOGGER.info(
            "No ANI-eligible samples remained after rescue gating; writing empty ANI outputs."
        )
        (fastani_dir / "fastani.tsv").write_text("", encoding="utf-8")
        (fastani_dir / "fastani.log").write_text("", encoding="utf-8")
        (fastani_dir / "fastani.matrix").write_text("", encoding="utf-8")
        cluster_path, ani_summary, ani_representatives = write_empty_ani_outputs(cluster_dir)
    log_cluster_outputs(cluster_path, ani_summary, ani_representatives)

    LOGGER.info("Building partial rescued final tables.")
    master_table_path = build_partial_master_table(
        validated_samples=copied_validated_samples,
        metadata=args.metadata,
        busco_lineages=args.busco_lineage,
        checkm2_summary=checkm2_summary,
        sixteen_s_summary=sixteen_s_summary,
        busco_summaries=busco_summary_paths,
        ani_summary=ani_summary,
        assembly_stats=assembly_stats,
        output_path=tables_dir / "master_table.tsv",
    )
    sample_status_path = build_partial_sample_status(
        validated_samples=copied_validated_samples,
        initial_status=copied_initial_status,
        metadata=args.metadata,
        busco_lineages=args.busco_lineage,
        checkm2_summary=checkm2_summary,
        sixteen_s_summary=sixteen_s_summary,
        busco_summaries=busco_summary_paths,
        ani_summary=ani_summary,
        assembly_stats=assembly_stats,
        output_path=tables_dir / "sample_status.tsv",
    )
    provenance_path = write_rescue_provenance(
        rescue_outdir=args.outdir,
        source_outdir=args.source_outdir,
        metadata=args.metadata,
        busco_lineages=args.busco_lineage,
        gcode_rule=args.gcode_rule,
        ani_threshold=args.ani_threshold,
        ani_score_profile=args.ani_score_profile,
        sample_count=len(samples),
        ani_exclusions_path=fastani_dir / "ani_exclusions.tsv",
        cluster_path=cluster_path,
        representative_path=ani_representatives,
    )
    LOGGER.info(
        "Rescue outputs ready: master_table=%s sample_status=%s provenance=%s.",
        master_table_path,
        sample_status_path,
        provenance_path,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the rescue CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_rescue(args)
    except (
        RescueAniError,
        build_fastani_inputs.FastAniInputError,
        build_master_table.MasterTableError,
        build_sample_status.SampleStatusError,
        validate_inputs.ValidationError,
        OSError,
    ) as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("ANI rescue completed under %s.", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
