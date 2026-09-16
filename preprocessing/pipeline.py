"""Reproducible preprocessing pipeline for promoter metadata."""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Iterable

import pandas as pd

from .homology import cluster_sequences, homology_aware_split
from .schema import REQUIRED_COLUMNS, ValidationReport, normalize_sequence, validate_metadata
from .shuffle import SHUFFLE_METHODS, shuffle_sequence
from .windowing import apply_tss_window


DEFAULT_UNKNOWN = "unknown"

# Sentinel meaning "no TSS windowing requested". ``None`` is a valid explicit
# choice (keep the native sequence), so it cannot double as the sentinel.
_NO_WINDOW = object()


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
    tss_window: tuple[int, int] | None | object = _NO_WINDOW,
    tss_column: str = "tss_position",
    strand_column: str = "strand",
    homology_clustering: bool = False,
    identity_threshold: float = 0.90,
    kmer_size: int = 15,
) -> tuple[pd.DataFrame, ValidationReport]:
    """Normalize, validate, deduplicate, and assign deterministic data splits.

    Rows with empty/invalid DNA, excessive ``N`` content, missing IDs, or
    missing source/species values are removed from the returned clean table.
    Every removal is represented in the returned validation report.

    Pipeline order matters and is deliberate:

    1. **TSS windowing** (optional) is applied *before* validation, because a
       window that extends past the sequence ends must be recorded as a dropped
       record, not silently truncated.
    2. Quality filtering and exact deduplication run on the resulting sequences,
       so identical TSS windows are deduplicated even when the source records
       differed.
    3. **Homology clustering** (optional) runs *after* deduplication and forces
       whole clusters onto one side of the discovery/validation split, so
       near-identical sequences cannot leak across the boundary.

    Pass ``tss_window=(upstream, downstream)`` to anchor every sequence to its
    transcription start site and normalise strand orientation. The annotation
    columns named by ``tss_column`` and ``strand_column`` must then be present.
    """

    if not 0 <= max_n_fraction <= 1:
        raise ValueError("max_n_fraction must be between 0 and 1")
    _require_columns(frame)

    tss_windowing_report: dict | None = None
    source_frame = frame
    if tss_window is not _NO_WINDOW and tss_window is not None:
        upstream, downstream = tss_window
        missing = [
            column
            for column in (tss_column, strand_column, "sequence_id")
            if column not in source_frame
        ]
        if missing:
            raise ValueError(
                "TSS windowing requires the annotation columns: " + ", ".join(missing)
            )
        source_frame, tss_windowing_report = apply_tss_window(
            source_frame,
            upstream=int(upstream),
            downstream=int(downstream),
            tss_column=tss_column,
            strand_column=strand_column,
        )

    clean = source_frame.copy()
    for column in ("sequence_id", "species", "source_dataset"):
        clean[column] = clean[column].fillna("").astype(str).str.strip()
    clean["sequence"] = clean["sequence"].map(normalize_sequence)
    for column in ("sigma_factor_type", "evidence_level"):
        if column not in clean:
            clean[column] = DEFAULT_UNKNOWN
        clean[column] = clean[column].fillna(DEFAULT_UNKNOWN).astype(str).str.strip()
        clean.loc[clean[column].eq(""), column] = DEFAULT_UNKNOWN

    validation = validate_metadata(clean)
    missing_identity_mask = (
        clean["sequence_id"].eq("")
        | clean["species"].eq("")
        | clean["source_dataset"].eq("")
    )
    empty_sequence_mask = clean["sequence"].eq("")
    invalid_character_mask = clean["sequence"].map(lambda sequence: bool(set(sequence) - set("ACGTN")))
    high_n_mask = clean["sequence"].map(
        lambda sequence: bool(sequence) and sequence.count("N") / len(sequence) > max_n_fraction
    )
    valid_mask = (
        clean["sequence_id"].ne("")
        & clean["species"].ne("")
        & clean["source_dataset"].ne("")
        & clean["sequence"].ne("")
        & clean["sequence"].map(lambda sequence: not set(sequence) - set("ACGTN"))
        & clean["sequence"].map(
            lambda sequence: bool(sequence)
            and sequence.count("N") / len(sequence) <= max_n_fraction
        )
    )
    clean = clean.loc[valid_mask].copy()
    quality_filter_counts = {
        "missing_identity": int(missing_identity_mask.sum()),
        "empty_sequence": int(empty_sequence_mask.sum()),
        "invalid_characters": int(invalid_character_mask.sum()),
        "high_n_fraction": int(high_n_mask.sum()),
    }
    duplicate_id_mask = clean["sequence_id"].duplicated(keep="first")
    quality_filter_counts["duplicate_sequence_id"] = int(duplicate_id_mask.sum())
    clean["sequence_length"] = clean["sequence"].str.len().astype(int)
    clean["gc_fraction"] = clean["sequence"].map(
        lambda sequence: (sequence.count("G") + sequence.count("C")) / len(sequence)
    )

    clean = clean.drop_duplicates(subset=["sequence_id"], keep="first")
    duplicate_sequence_mask = clean["sequence"].duplicated(keep="first")
    quality_filter_counts["duplicate_sequence"] = int(duplicate_sequence_mask.sum())
    clean = clean.drop_duplicates(subset=["sequence"], keep="first")

    validation_source_set = {str(source) for source in validation_sources if str(source).strip()}
    clean["split"] = clean["source_dataset"].map(
        lambda source: "validation" if source in validation_source_set else "discovery"
    )

    homology_report: dict | None = None
    split_report: dict | None = None
    if homology_clustering and not clean.empty:
        clean, homology_report = cluster_sequences(
            clean,
            identity_threshold=identity_threshold,
            kmer_size=kmer_size,
        )
        clean, split_report = homology_aware_split(clean)

    clean = clean.sort_values("sequence_id").reset_index(drop=True)
    clean.attrs["quality_filter_counts"] = quality_filter_counts
    clean.attrs["tss_windowing"] = tss_windowing_report
    clean.attrs["homology_clustering"] = homology_report
    clean.attrs["homology_split"] = split_report
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
        "schema_version": "m1.1",
        "input_rows": int(len(original)),
        "output_rows": int(len(clean)),
        "removed_rows": int(len(original) - len(clean)),
        "removed_by_reason": clean.attrs.get("quality_filter_counts", {}),
        "max_n_fraction": max_n_fraction,
        "validation": validation.to_dict(),
        "tss_windowing": clean.attrs.get("tss_windowing"),
        "homology_clustering": clean.attrs.get("homology_clustering"),
        "homology_split": clean.attrs.get("homology_split"),
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
    """Shuffle a sequence while preserving its length and base composition.

    Retained for backward compatibility; it is the N0 mononucleotide control.
    Use :func:`preprocessing.shuffle.shuffle_sequence` to select between the N0
    and N1 controls.
    """

    bases = list(sequence)
    rng.shuffle(bases)
    return "".join(bases)


def write_background_fasta(
    frame: pd.DataFrame,
    output_path: str | Path,
    *,
    replicates: int = 1,
    seed: int = 20260911,
    method: str = "mononucleotide",
) -> None:
    """Write shuffled background controls with deterministic IDs.

    ``method`` selects the negative control: ``mononucleotide`` is the N0
    control that only preserves single-base composition, ``dinucleotide`` is the
    stricter N1 control that also preserves dinucleotide composition.
    """

    if replicates < 1:
        raise ValueError("replicates must be at least 1")
    if method not in SHUFFLE_METHODS:
        raise ValueError(f"method must be one of {SHUFFLE_METHODS}")
    rng = random.Random(seed)
    rows = []
    for _, row in frame.iterrows():
        for replicate in range(1, replicates + 1):
            rows.append(
                (
                    f"{row['sequence_id']}__{method}{replicate}",
                    shuffle_sequence(row["sequence"], rng, method=method),
                )
            )
    write_fasta(rows, output_path)


def write_background_controls(
    frame: pd.DataFrame,
    output_dir: str | Path,
    *,
    replicates: int = 1,
    seed: int = 20260911,
    methods: Iterable[str] = ("mononucleotide", "dinucleotide"),
) -> dict:
    """Write one background FASTA per shuffle method.

    Both controls share the same seed so that the only difference between them
    is the shuffle algorithm, not the random stream position.
    """

    output_dir = Path(output_dir)
    written: dict[str, str] = {}
    for method in methods:
        if method not in SHUFFLE_METHODS:
            raise ValueError(f"method must be one of {SHUFFLE_METHODS}")
        path = output_dir / f"background_{method}.fasta"
        write_background_fasta(frame, path, replicates=replicates, seed=seed, method=method)
        written[method] = str(path.resolve())
    return written


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
