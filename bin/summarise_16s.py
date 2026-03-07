#!/usr/bin/env python3
"""Summarise Barrnap 16S calls for one sample."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


LOGGER = logging.getLogger(__name__)

STATUS_COLUMNS = (
    "accession",
    "16S",
    "best_16S_header",
    "best_16S_length",
    "include_in_all_best_16S",
    "warnings",
)
TRUE_TOKENS = {"true", "t", "yes", "y", "1"}
FALSE_TOKENS = {"false", "f", "no", "n", "0"}


@dataclass(frozen=True)
class FastaRecord:
    """Represent one FASTA record."""

    header: str
    sequence: str


@dataclass(frozen=True)
class RrnaHit:
    """Represent one paired Barrnap GFF and FASTA hit."""

    order_index: int
    score: float
    is_16s: bool
    is_partial: bool
    record: FastaRecord


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise Barrnap GFF/FASTA outputs for one sample."
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Original sample accession for the output row.",
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
    parser.add_argument(
        "--is-atypical",
        default="false",
        help="Whether the sample is atypical and should be excluded from cohort 16S.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_boolean(value: str) -> bool:
    """Parse a boolean-like CLI argument."""
    token = value.strip().lower()
    if token in TRUE_TOKENS:
        return True
    if token in FALSE_TOKENS:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def parse_gff_hits(path: Path) -> list[tuple[bool, bool, float]]:
    """Parse Barrnap GFF rows into ordered hit metadata."""
    if not path.is_file():
        raise ValueError(f"Missing Barrnap GFF file: {path}")

    parsed_hits: list[tuple[bool, bool, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF line in {path}: {line}")
            attributes = fields[8].lower()
            is_16s = "16s" in attributes
            is_partial = "partial" in attributes
            score_value = fields[5]
            if score_value == ".":
                score = float("inf")
            else:
                try:
                    score = float(score_value)
                except ValueError as error:
                    raise ValueError(
                        f"Could not parse Barrnap score {score_value!r} in {path}."
                    ) from error
            parsed_hits.append((is_16s, is_partial, score))
    return parsed_hits


def parse_fasta_records(path: Path) -> list[FastaRecord]:
    """Parse a FASTA file into ordered records."""
    if not path.is_file():
        raise ValueError(f"Missing Barrnap FASTA file: {path}")

    records: list[FastaRecord] = []
    header: str | None = None
    sequence_lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(
                        FastaRecord(header=header, sequence="".join(sequence_lines))
                    )
                header = line[1:].strip()
                sequence_lines = []
                continue
            if header is None:
                raise ValueError(f"Malformed FASTA file without a header: {path}")
            sequence_lines.append(line.strip())

    if header is not None:
        records.append(FastaRecord(header=header, sequence="".join(sequence_lines)))
    return records


def pair_hits(
    gff_hits: Sequence[tuple[bool, bool, float]],
    fasta_records: Sequence[FastaRecord],
) -> list[RrnaHit]:
    """Pair GFF rows to FASTA records in file order."""
    if len(gff_hits) != len(fasta_records):
        raise ValueError(
            "Barrnap GFF and FASTA outputs contain different numbers of hits."
        )

    paired_hits: list[RrnaHit] = []
    for order_index, (gff_hit, record) in enumerate(
        zip(gff_hits, fasta_records, strict=True),
        start=1,
    ):
        gff_is_16s, gff_is_partial, score = gff_hit
        header_text = record.header.lower()
        is_16s = gff_is_16s or "16s" in header_text
        is_partial = gff_is_partial or "partial" in header_text
        paired_hits.append(
            RrnaHit(
                order_index=order_index,
                score=score,
                is_16s=is_16s,
                is_partial=is_partial,
                record=record,
            )
        )
    return paired_hits


def choose_best_hit(hits: Sequence[RrnaHit]) -> RrnaHit | None:
    """Choose the best 16S hit using the locked Barrnap tie-break rules."""
    candidate_hits = [hit for hit in hits if hit.is_16s]
    if not candidate_hits:
        return None

    intact_hits = [hit for hit in candidate_hits if not hit.is_partial]
    ranking_pool = intact_hits if intact_hits else candidate_hits
    return min(
        ranking_pool,
        key=lambda hit: (
            hit.score,
            -len(hit.record.sequence),
            hit.order_index,
        ),
    )


def write_best_fasta(path: Path, hit: RrnaHit | None) -> None:
    """Write the selected best-hit FASTA record or an empty file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if hit is None:
        path.write_text("", encoding="utf-8")
        return

    wrapped_sequence = [
        hit.record.sequence[index : index + 80]
        for index in range(0, len(hit.record.sequence), 80)
    ]
    path.write_text(
        f">{hit.record.header}\n" + "\n".join(wrapped_sequence) + "\n",
        encoding="utf-8",
    )


def write_status(path: Path, row: dict[str, str]) -> None:
    """Write the one-row 16S status TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(STATUS_COLUMNS), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def summarise_hits(
    accession: str,
    rrna_gff: Path,
    rrna_fasta: Path,
    outdir: Path,
    is_atypical: bool,
) -> None:
    """Summarise one sample's Barrnap outputs into 16S status files."""
    warnings: list[str] = []
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        hits = pair_hits(parse_gff_hits(rrna_gff), parse_fasta_records(rrna_fasta))
    except ValueError as error:
        LOGGER.warning(str(error))
        warnings.append("invalid_barrnap_output")
        row = {
            "accession": accession,
            "16S": "NA",
            "best_16S_header": "NA",
            "best_16S_length": "NA",
            "include_in_all_best_16S": "false",
            "warnings": ";".join(warnings),
        }
        write_best_fasta(outdir / "best_16S.fna", None)
        write_status(outdir / "16S_status.tsv", row)
        return

    sixteen_s_hits = [hit for hit in hits if hit.is_16s]
    best_hit = choose_best_hit(hits)

    if not sixteen_s_hits:
        status = "No"
    elif any(not hit.is_partial for hit in sixteen_s_hits):
        status = "Yes"
    else:
        status = "partial"

    include_in_all = best_hit is not None and not is_atypical
    row = {
        "accession": accession,
        "16S": status,
        "best_16S_header": best_hit.record.header if best_hit is not None else "NA",
        "best_16S_length": (
            str(len(best_hit.record.sequence)) if best_hit is not None else "NA"
        ),
        "include_in_all_best_16S": "true" if include_in_all else "false",
        "warnings": "",
    }

    write_best_fasta(outdir / "best_16S.fna", best_hit)
    write_status(outdir / "16S_status.tsv", row)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the 16S summarisation CLI."""
    args = parse_args(argv)
    configure_logging()
    summarise_hits(
        accession=args.accession,
        rrna_gff=args.rrna_gff,
        rrna_fasta=args.rrna_fasta,
        outdir=args.outdir,
        is_atypical=normalise_boolean(args.is_atypical),
    )
    LOGGER.info("Wrote 16S summary for %s.", args.accession)
    return 0


if __name__ == "__main__":
    sys.exit(main())
