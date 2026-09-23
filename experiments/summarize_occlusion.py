"""Summarize fixed-rule L2 occlusion windows in TSS-relative coordinates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from preprocessing.provenance import sha256_file


def run(
    input_path: str | Path,
    windows_path: str | Path,
    output_dir: str | Path,
    *,
    top_k: int = 3,
) -> dict:
    """Keep the top-k windows per sequence using a fixed deterministic rule."""

    input_path = Path(input_path).resolve()
    windows_path = Path(windows_path).resolve()
    output_dir = Path(output_dir).resolve()
    if top_k < 1:
        raise ValueError("top_k must be positive")
    metadata = pd.read_csv(input_path, sep="\t")
    windows = pd.read_csv(windows_path, sep="\t")
    required_metadata = {"sequence_id", "sigma_factor_type", "tss_offset_in_window"}
    missing_metadata = sorted(required_metadata - set(metadata.columns))
    if missing_metadata:
        raise ValueError(f"L2 input is missing columns: {', '.join(missing_metadata)}")
    required_windows = {"sequence_id", "start", "end", "importance", "original_fragment"}
    missing_windows = sorted(required_windows - set(windows.columns))
    if missing_windows:
        raise ValueError(f"occlusion table is missing columns: {', '.join(missing_windows)}")
    if metadata["sequence_id"].astype(str).duplicated().any():
        raise ValueError("L2 input contains duplicate sequence_id values")
    if windows.empty:
        raise ValueError("occlusion table is empty")

    metadata = metadata.copy()
    metadata["sigma_factor_type"] = metadata["sigma_factor_type"].fillna("unknown").astype(str)
    metadata.loc[metadata["sigma_factor_type"].str.strip().eq(""), "sigma_factor_type"] = "unknown"
    metadata["tss_offset_in_window"] = pd.to_numeric(metadata["tss_offset_in_window"], errors="coerce")
    if metadata["tss_offset_in_window"].isna().any():
        raise ValueError("L2 input contains missing or non-numeric tss_offset_in_window")

    metadata_columns = ["sequence_id", "sigma_factor_type", "tss_offset_in_window"]
    for column in ("species", "source_dataset", "split"):
        if column in metadata:
            metadata_columns.append(column)
    frame = windows.merge(
        metadata[metadata_columns],
        on="sequence_id",
        how="left",
        validate="many_to_one",
    )
    if frame["sigma_factor_type"].isna().any():
        raise ValueError("occlusion table contains sequence_id values absent from the L2 input")
    for column in ("start", "end", "importance"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    frame["relative_start"] = frame["start"] - frame["tss_offset_in_window"]
    frame["relative_end"] = frame["end"] - frame["tss_offset_in_window"]
    frame = frame.sort_values(
        ["sequence_id", "importance", "start"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    frame["rank_within_sequence"] = frame.groupby("sequence_id", sort=False).cumcount() + 1
    top = frame.loc[frame["rank_within_sequence"] <= top_k].copy()
    if top["sequence_id"].nunique() != metadata["sequence_id"].nunique():
        missing_ids = sorted(
            set(metadata["sequence_id"].astype(str)) - set(top["sequence_id"].astype(str))
        )
        raise ValueError(f"some selected sequences have no top windows: {missing_ids[:5]}")

    output_dir.mkdir(parents=True, exist_ok=True)
    top_columns = [
        "sequence_id",
        "sigma_factor_type",
        "start",
        "end",
        "relative_start",
        "relative_end",
        "importance",
        "rank_within_sequence",
        "original_fragment",
    ]
    top[top_columns].to_csv(output_dir / "top_occlusion_windows.tsv", sep="\t", index=False)

    group_summary = (
        top.groupby("sigma_factor_type", sort=True)
        .agg(
            sequence_count=("sequence_id", "nunique"),
            top_window_count=("sequence_id", "size"),
            importance_median=("importance", "median"),
            importance_q1=("importance", lambda values: values.quantile(0.25)),
            importance_q3=("importance", lambda values: values.quantile(0.75)),
            relative_start_median=("relative_start", "median"),
            relative_start_q1=("relative_start", lambda values: values.quantile(0.25)),
            relative_start_q3=("relative_start", lambda values: values.quantile(0.75)),
        )
        .reset_index()
    )
    group_summary.to_csv(output_dir / "occlusion_group_summary.tsv", sep="\t", index=False)

    importance = top["importance"]
    result = {
        "schema_version": "occlusion-summary-1.0",
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "windows": str(windows_path),
        "windows_sha256": sha256_file(windows_path),
        "sequence_count": int(metadata["sequence_id"].nunique()),
        "window_count": int(len(frame)),
        "top_window_count": int(len(top)),
        "top_k_per_sequence": top_k,
        "group_summary": group_summary.to_dict(orient="records"),
        "top_importance_median": float(importance.median()),
        "top_importance_q1": float(importance.quantile(0.25)),
        "top_importance_q3": float(importance.quantile(0.75)),
        "selection_rule": f"For each fixed representative sequence, retain the top {top_k} windows by cosine-distance change; ties are resolved by the smaller window start.",
        "boundary": "Exploratory source-internal L2 evidence. No resampling stability, independent validation, or biological significance is claimed.",
    }
    (output_dir / "occlusion_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--windows", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--top-k", type=int, default=3)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.input, args.windows, args.output_dir, top_k=args.top_k),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
