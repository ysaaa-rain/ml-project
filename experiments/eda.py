"""M2 exploratory analysis for the standardized promoter table."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def _group_summary(frame: pd.DataFrame) -> pd.DataFrame:
    group_columns = ["species", "source_dataset", "sigma_factor_type", "split"]
    summary = (
        frame.groupby(group_columns, dropna=False)
        .agg(
            sequence_count=("sequence_id", "count"),
            length_min=("sequence_length", "min"),
            length_median=("sequence_length", "median"),
            length_max=("sequence_length", "max"),
            gc_fraction_mean=("gc_fraction", "mean"),
            gc_fraction_median=("gc_fraction", "median"),
        )
        .reset_index()
        .sort_values(group_columns)
    )
    return summary


def _position_composition(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for position in range(int(frame["sequence_length"].max())):
        bases = frame["sequence"].str[position].dropna()
        counts = bases.value_counts()
        total = int(counts.sum())
        rows.append(
            {
                "position_1based": position + 1,
                "n_sequences": total,
                "a_fraction": counts.get("A", 0) / total if total else 0.0,
                "c_fraction": counts.get("C", 0) / total if total else 0.0,
                "g_fraction": counts.get("G", 0) / total if total else 0.0,
                "t_fraction": counts.get("T", 0) / total if total else 0.0,
                "n_fraction": counts.get("N", 0) / total if total else 0.0,
            }
        )
    return pd.DataFrame(rows)


def _save_plots(frame: pd.DataFrame, output_dir: Path) -> list[str]:
    plot_paths: list[str] = []

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].hist(frame["sequence_length"], bins=min(20, max(1, frame["sequence_length"].nunique())))
    axes[0].set_title("Promoter length")
    axes[0].set_xlabel("Length (bp)")
    axes[0].set_ylabel("Count")
    axes[1].hist(frame["gc_fraction"], bins=10, range=(0, 1))
    axes[1].set_title("GC fraction")
    axes[1].set_xlabel("GC fraction")
    axes[1].set_ylabel("Count")
    path = output_dir / "length_gc_distribution.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    plot_paths.append(path.name)

    counts = frame["sigma_factor_type"].value_counts().sort_values()
    fig, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
    counts.plot.barh(ax=axis, color="#4472C4")
    axis.set_title("Sequences by sigma factor type")
    axis.set_xlabel("Count")
    axis.set_ylabel("Sigma factor type")
    path = output_dir / "sigma_factor_counts.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    plot_paths.append(path.name)
    return plot_paths


def run_eda(input_path: str | Path, output_dir: str | Path) -> dict:
    """Run deterministic descriptive analysis on a cleaned TSV table."""

    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(input_path, sep="\t")
    required = {"sequence", "sequence_id", "sequence_length", "gc_fraction"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"clean table is missing EDA columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError("cannot run EDA on an empty clean table")

    group_summary = _group_summary(frame)
    position_composition = _position_composition(frame)
    group_summary.to_csv(output_dir / "group_summary.tsv", sep="\t", index=False)
    position_composition.to_csv(output_dir / "position_composition.tsv", sep="\t", index=False)
    plot_paths = _save_plots(frame, output_dir)
    summary = {
        "input": str(input_path.resolve()),
        "rows": int(len(frame)),
        "species": sorted(frame["species"].dropna().unique().tolist()),
        "source_datasets": sorted(frame["source_dataset"].dropna().unique().tolist()),
        "sigma_factor_types": sorted(frame["sigma_factor_type"].dropna().unique().tolist()),
        "plots": plot_paths,
    }
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Cleaned promoter TSV")
    parser.add_argument("--output-dir", required=True, help="EDA output directory")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary = run_eda(args.input, args.output_dir)
    print(f"rows={summary['rows']}")
    print(f"output_dir={Path(args.output_dir).resolve()}")


if __name__ == "__main__":
    main()
