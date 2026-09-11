"""Reproducible preprocessing pipeline for promoter metadata."""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Iterable

import pandas as pd

from .schema import REQUIRED_COLUMNS, ValidationReport, normalize_sequence, validate_metadata


DEFAULT_UNKNOWN = "unknown"


def load_metadata(path: str | Path) -> pd.DataFrame:
    """Load a CSV/TSV metadata table based on its extension."""

    input_path = Path(path)
    if input_path.suffix.lower() in {".csv"}:
        frame = pd.read_csv(input_path)
    else:
        frame = pd.read_csv(input_path, sep="\t")
    return frame


def _require_columns(frame: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"metadata is missing required columns: {', '.join(missing)}")


def normalize_metadata(
    frame: pd.DataFrame,
    *,
    max_n_fraction: float = 0.1,
    validation_sources: Iterable[str] = (),
) -> tuple[pd.DataFrame, ValidationReport]:
    """Normalize, validate, deduplicate, and assign deterministic data splits.

    Rows with empty/invalid DNA, excessive ``N`` content, missing IDs, or
    missing source/species values are removed from the returned clean table.
    Every removal is represented in the returned validation report.
    """

    if not 0 <= max_n_fraction <= 1:
        raise ValueError("max_n_fraction must be between 0 and 1")
    _require_columns(frame)

    clean = frame.copy()
    for column in ("sequence_id", "species", "source_dataset"):
        clean[column] = clean[column].fillna("").astype(str).str.strip()
    clean["sequence"] = clean["sequence"].map(normalize_sequence)
    for column in ("sigma_factor_type", "evidence_level"):
        if column not in clean:
            clean[column] = DEFAULT_UNKNOWN
        clean[column] = clean[column].fillna(DEFAULT_UNKNOWN).astype(str).str.strip()
        clean.loc[clean[column].eq(""), column] = DEFAULT_UNKNOWN

    validation = validate_metadata(clean)
    valid_mask = (
        clean["sequence_id"].ne("")
        & clean["species"].ne("")
        & clean["source_dataset"].ne("")
        & clean["sequence"].ne("")
        & clean["sequence"].map(lambda sequence: not set(sequence) - set("ACGTN"))
        & clean["sequence"].map(lambda sequence: sequence.count("N") / len(sequence) <= max_n_fraction)
    )
    clean = clean.loc[valid_mask].copy()
    clean["sequence_length"] = clean["sequence"].str.len().astype(int)
    clean["gc_fraction"] = clean["sequence"].map(
        lambda sequence: (sequence.count("G") + sequence.count("C")) / len(sequence)
    )

    clean = clean.drop_duplicates(subset=["sequence_id"], keep="first")
    clean = clean.drop_duplicates(subset=["sequence"], keep="first")

    validation_source_set = {str(source) for source in validation_sources if str(source).strip()}
    clean["split"] = clean["source_dataset"].map(
        lambda source: "validation" if source in validation_source_set else "discovery"
    )
    clean = clean.sort_values("sequence_id").reset_index(drop=True)
    return clean, validation


def build_quality_report(
    original: pd.DataFrame,
    clean: pd.DataFrame,
    validation: ValidationReport,
    *,
    max_n_fraction: float,
) -> dict:
    """Build an auditable quality report for a preprocessing run."""

    return {
        "schema_version": "m1.0",
        "input_rows": int(len(original)),
        "output_rows": int(len(clean)),
        "removed_rows": int(len(original) - len(clean)),
        "max_n_fraction": max_n_fraction,
        "validation": validation.to_dict(),
        "output_counts": {
            "species": clean["species"].value_counts().to_dict() if not clean.empty else {},
            "source_dataset": clean["source_dataset"].value_counts().to_dict() if not clean.empty else {},
            "sigma_factor_type": clean["sigma_factor_type"].value_counts().to_dict() if not clean.empty else {},
            "split": clean["split"].value_counts().to_dict() if not clean.empty else {},
        },
        "length": {
            "min": int(clean["sequence_length"].min()) if not clean.empty else None,
            "median": float(clean["sequence_length"].median()) if not clean.empty else None,
            "max": int(clean["sequence_length"].max()) if not clean.empty else None,
        },
        "gc_fraction": {
            "min": float(clean["gc_fraction"].min()) if not clean.empty else None,
            "median": float(clean["gc_fraction"].median()) if not clean.empty else None,
            "max": float(clean["gc_fraction"].max()) if not clean.empty else None,
        },
    }


def shuffled_sequence(sequence: str, rng: random.Random) -> str:
    """Shuffle a sequence while preserving its length and base composition."""

    bases = list(sequence)
    rng.shuffle(bases)
    return "".join(bases)


def write_background_fasta(
    frame: pd.DataFrame,
    output_path: str | Path,
    *,
    replicates: int = 1,
    seed: int = 20260911,
) -> None:
    """Write sequence-shuffled background controls with deterministic IDs."""

    if replicates < 1:
        raise ValueError("replicates must be at least 1")
    rng = random.Random(seed)
    rows = []
    for _, row in frame.iterrows():
        for replicate in range(1, replicates + 1):
            rows.append(
                (
                    f"{row['sequence_id']}__shuffle{replicate}",
                    shuffled_sequence(row["sequence"], rng),
                )
            )
    write_fasta(rows, output_path)


def write_fasta(records: Iterable[tuple[str, str]], output_path: str | Path) -> None:
    """Write ``(identifier, sequence)`` records in a wrapped FASTA format."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(f"{sequence[start:start + 80]}\n")


def write_quality_report(report: dict, output_path: str | Path) -> None:
    """Write a stable, UTF-8 JSON quality report."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
