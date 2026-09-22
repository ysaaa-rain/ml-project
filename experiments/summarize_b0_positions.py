"""Summarize B0 hit positions in a TSS-aligned input table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def run(input_path: str | Path, hits_path: str | Path, output_dir: str | Path) -> dict:
    input_path = Path(input_path)
    hits_path = Path(hits_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(input_path, sep="\t")
    hits = pd.read_csv(hits_path, sep="\t")
    required = {"sequence_id", "tss_offset_in_window"}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(f"TSS-aligned input is missing columns: {', '.join(missing)}")
    if hits.empty:
        raise ValueError("B0 hit table is empty")
    offsets = metadata[["sequence_id", "tss_offset_in_window"]]
    frame = hits.merge(offsets, on="sequence_id", how="left", validate="many_to_one")
    if frame["tss_offset_in_window"].isna().any():
        raise ValueError("some B0 hits have no TSS offset in the input table")
    frame["relative_start"] = frame["start"] - frame["tss_offset_in_window"]
    summary = (
        frame.groupby(["motif_name", "strand"], sort=True)
        .agg(
            hit_count=("sequence_id", "size"),
            sequence_count=("sequence_id", "nunique"),
            relative_start_median=("relative_start", "median"),
            relative_start_q1=("relative_start", lambda values: values.quantile(0.25)),
            relative_start_q3=("relative_start", lambda values: values.quantile(0.75)),
            relative_start_min=("relative_start", "min"),
            relative_start_max=("relative_start", "max"),
        )
        .reset_index()
    )
    summary.to_csv(output_dir / "position_summary.tsv", sep="\t", index=False)

    figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, sharey=True)
    for axis, strand in zip(axes, ("+", "-")):
        subset = frame.loc[frame["strand"].eq(strand)]
        for motif_name, motif_frame in subset.groupby("motif_name", sort=True):
            axis.hist(
                motif_frame["relative_start"],
                bins=range(-60, 17),
                alpha=0.55,
                label=motif_name,
            )
        axis.set_title(f"B0 exact-consensus hits, scan strand {strand}")
        axis.set_ylabel("hit count")
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend(fontsize=8)
    axes[-1].set_xlabel("hit start relative to TSS (bp)")
    figure.suptitle("Exploratory B0 positions in TSS-aligned RegulonDB subset")
    figure.tight_layout()
    figure.savefig(output_dir / "position_distribution.png", dpi=160)
    plt.close(figure)

    positive = summary.loc[summary["strand"].eq("+")]
    result = {
        "schema_version": "b0-position-summary-1.0",
        "input": str(input_path.resolve()),
        "hits": str(hits_path.resolve()),
        "hit_count": int(len(frame)),
        "position_summary": summary.to_dict(orient="records"),
        "positive_strand_medians": {
            str(row.motif_name): float(row.relative_start_median)
            for row in positive.itertuples()
        },
        "boundary": "Exploratory exact-consensus positions; not FIMO p-values or q-values.",
    }
    (output_dir / "position_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--hits", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    result = run(args.input, args.hits, args.output_dir)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
