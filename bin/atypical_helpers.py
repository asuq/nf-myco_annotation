#!/usr/bin/env python3
"""Shared helpers for interpreting NCBI atypical-warning metadata."""

from __future__ import annotations

from collections.abc import Callable, Mapping


ATYPICAL_WARNINGS_COLUMN = "Atypical_Warnings"
UNVERIFIED_SOURCE_ORGANISM_EXCEPTION = "unverified source organism"
MISSING_VALUE_TOKENS = {"", "na", "n/a", "null", "none"}


def is_missing_atypical_warning(value: str | None) -> bool:
    """Return True when an atypical-warning value should be treated as missing."""
    if value is None:
        return True
    return value.strip().lower() in MISSING_VALUE_TOKENS


def classify_atypical_warning(value: str | None) -> tuple[bool, bool]:
    """Classify one raw atypical-warning value into atypical and exception flags."""
    if is_missing_atypical_warning(value):
        return False, False

    lowered = value.casefold()
    return True, UNVERIFIED_SOURCE_ORGANISM_EXCEPTION in lowered


def detect_atypical_flags(
    metadata_row: Mapping[str, str],
    *,
    find_column_by_normalised_name: Callable[[tuple[str, ...], str], str | None],
) -> tuple[bool, bool]:
    """Read and classify the atypical-warning value from one metadata row."""
    atypical_column = find_column_by_normalised_name(
        tuple(metadata_row),
        ATYPICAL_WARNINGS_COLUMN,
    )
    if atypical_column is None:
        return False, False
    return classify_atypical_warning(metadata_row.get(atypical_column))
