"""启动子序列质量、组成与已知元件的探索性分析。"""

from __future__ import annotations

import argparse
import math
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
    """Summarize base composition and plug-in information at each position.

    Information is reported relative to a uniform A/C/G/T background (2-H).
    It is a descriptive positional statistic, not an evolutionary-conservation
    estimate and not a bias-corrected motif score.
    """

    rows = []
    tss_offset = None
    if "tss_offset_in_window" in frame and frame["tss_offset_in_window"].notna().any():
        offsets = frame["tss_offset_in_window"].dropna().astype(int).unique()
        if len(offsets) != 1:
            raise ValueError("EDA expects one shared tss_offset_in_window for positional plots")
        tss_offset = int(offsets[0])

    for position in range(int(frame["sequence_length"].max())):
        bases = frame["sequence"].str[position].dropna()
        counts = bases.value_counts()
        total = int(counts.sum())
        canonical_total = sum(int(counts.get(base, 0)) for base in "ACGT")
        fractions = {
            base: int(counts.get(base, 0)) / canonical_total if canonical_total else 0.0
            for base in "ACGT"
        }
        entropy = -sum(value * math.log2(value) for value in fractions.values() if value > 0)
        consensus_base = max("ACGT", key=lambda base: (fractions[base], -"ACGT".index(base)))
        rows.append(
            {
                "position_1based": position + 1,
                "relative_position_bp": position - tss_offset if tss_offset is not None else None,
                "n_sequences": total,
                "canonical_base_count": canonical_total,
                "a_fraction": fractions["A"],
                "c_fraction": fractions["C"],
                "g_fraction": fractions["G"],
                "t_fraction": fractions["T"],
                "n_fraction": counts.get("N", 0) / total if total else 0.0,
                "shannon_entropy_bits": entropy,
                "information_bits_vs_uniform": 2.0 - entropy if canonical_total else 0.0,
                "consensus_base": consensus_base if canonical_total else "N",
                "consensus_fraction": fractions[consensus_base] if canonical_total else 0.0,
            }
        )
    return pd.DataFrame(rows)


def _position_composition_by_sigma(frame: pd.DataFrame) -> pd.DataFrame:
    if "sigma_factor_type" not in frame:
        raise ValueError("clean table is missing EDA column: sigma_factor_type")
    pieces = []
    labels = frame["sigma_factor_type"].fillna("unknown").astype(str)
    for label in sorted(labels.unique()):
        summary = _position_composition(frame.loc[labels.eq(label)])
        summary.insert(0, "sigma_factor_type", label)
        pieces.append(summary)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def _missingness(frame: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "column": column,
                "row_count": len(frame),
                "missing_count": int(frame[column].isna().sum()),
                "missing_fraction": float(frame[column].isna().mean()),
            }
            for column in frame.columns
        ]
    )


def _save_plots(
    frame: pd.DataFrame,
    position_composition: pd.DataFrame,
    position_by_sigma: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
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

    x = position_composition["relative_position_bp"].where(
        position_composition["relative_position_bp"].notna(),
        position_composition["position_1based"],
    )
    fig, axis = plt.subplots(figsize=(9, 4), constrained_layout=True)
    axis.plot(x, position_composition["information_bits_vs_uniform"], color="#4472C4")
    axis.set_title("Positional information relative to uniform A/C/G/T")
    axis.set_xlabel("Position relative to TSS (bp)" if position_composition["relative_position_bp"].notna().any() else "Position (1-based)")
    axis.set_ylabel("Information (bits)")
    path = output_dir / "position_information_overall.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    plot_paths.append(path.name)

    fig, axis = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for label, subset in position_by_sigma.groupby("sigma_factor_type", sort=True):
        x = subset["relative_position_bp"].where(
            subset["relative_position_bp"].notna(), subset["position_1based"]
        )
        axis.plot(x, subset["information_bits_vs_uniform"], label=label, linewidth=1.2)
    axis.set_title("Positional information by sigma group")
    axis.set_xlabel("Position relative to TSS (bp)" if position_by_sigma["relative_position_bp"].notna().any() else "Position (1-based)")
    axis.set_ylabel("Information (bits)")
    axis.legend(fontsize=8, ncol=2)
    path = output_dir / "position_information_by_sigma.png"
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
    position_by_sigma = _position_composition_by_sigma(frame)
    group_summary.to_csv(output_dir / "group_summary.tsv", sep="\t", index=False)
    position_composition.to_csv(output_dir / "position_composition.tsv", sep="\t", index=False)
    position_by_sigma.to_csv(output_dir / "position_composition_by_sigma.tsv", sep="\t", index=False)
    missingness = _missingness(frame)
    missingness.to_csv(output_dir / "missingness.tsv", sep="\t", index=False)
    plot_paths = _save_plots(frame, position_composition, position_by_sigma, output_dir)
    summary = {
        "input": str(input_path.resolve()),
        "rows": int(len(frame)),
        "species": sorted(frame["species"].dropna().unique().tolist()),
        "source_datasets": sorted(frame["source_dataset"].dropna().unique().tolist()),
        "sigma_factor_types": sorted(frame["sigma_factor_type"].dropna().unique().tolist()),
        "position_count": int(len(position_composition)),
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
