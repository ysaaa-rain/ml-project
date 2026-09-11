"""Schema and sequence-level validation for promoter metadata."""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any

import pandas as pd


REQUIRED_COLUMNS = (
    "sequence_id",
    "species",
    "source_dataset",
    "sequence",
)

RECOMMENDED_COLUMNS = (
    "sigma_factor_type",
    "evidence_level",
)

DNA_ALPHABET = frozenset("ACGTN")


@dataclass
class ValidationReport:
    """Machine-readable validation summary."""

    row_count: int
    valid_row_count: int
    invalid_row_count: int
    missing_required_columns: list[str]
    missing_value_counts: dict[str, int]
    invalid_sequence_ids: list[str]
    invalid_characters: dict[str, list[str]]
    duplicate_sequence_id_count: int
    duplicate_sequence_count: int
    warnings: list[str]

    @property
    def is_valid(self) -> bool:
        return not self.missing_required_columns and self.invalid_row_count == 0

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["is_valid"] = self.is_valid
        return result


def normalize_sequence(value: Any) -> str:
    """Normalize a DNA sequence while preserving N as an unknown base.

    Whitespace and common FASTA separators are removed. RNA ``U`` is converted
    to ``T`` explicitly because the project consumes promoter DNA sequences.
    Other characters are retained so validation can report them precisely.
    """

    if pd.isna(value):
        return ""
    sequence = str(value).upper().replace("U", "T")
    return "".join(sequence.split())


def invalid_characters(sequence: str) -> list[str]:
    """Return sorted non-DNA characters found in a normalized sequence."""

    return sorted(set(sequence) - DNA_ALPHABET)


def validate_metadata(frame: pd.DataFrame) -> ValidationReport:
    """Validate the minimum metadata contract without mutating ``frame``."""

    missing_required = [column for column in REQUIRED_COLUMNS if column not in frame]
    missing_value_counts: dict[str, int] = {}
    for column in REQUIRED_COLUMNS:
        if column in frame:
            missing_value_counts[column] = int(frame[column].isna().sum())

    invalid_ids: list[str] = []
    invalid_chars: dict[str, list[str]] = {}
    invalid_rows: set[int] = set()
    if "sequence" in frame:
        for index, raw_sequence in frame["sequence"].items():
            sequence = normalize_sequence(raw_sequence)
            chars = invalid_characters(sequence)
            if not sequence or chars:
                invalid_rows.add(index)
                sequence_id = str(frame.loc[index, "sequence_id"]) if "sequence_id" in frame else str(index)
                invalid_ids.append(sequence_id)
                if chars:
                    invalid_chars[sequence_id] = chars

    duplicate_id_count = 0
    if "sequence_id" in frame:
        duplicate_id_count = int(frame["sequence_id"].duplicated(keep=False).sum())

    duplicate_sequence_count = 0
    if "sequence" in frame:
        normalized = frame["sequence"].map(normalize_sequence)
        duplicate_sequence_count = int(normalized.duplicated(keep=False).sum())

    warnings: list[str] = []
    for column in RECOMMENDED_COLUMNS:
        if column not in frame:
            warnings.append(f"recommended column missing: {column}")

    valid_rows = len(frame) - len(invalid_rows)
    return ValidationReport(
        row_count=len(frame),
        valid_row_count=valid_rows,
        invalid_row_count=len(invalid_rows),
        missing_required_columns=missing_required,
        missing_value_counts=missing_value_counts,
        invalid_sequence_ids=invalid_ids,
        invalid_characters=invalid_chars,
        duplicate_sequence_id_count=duplicate_id_count,
        duplicate_sequence_count=duplicate_sequence_count,
        warnings=warnings,
    )
