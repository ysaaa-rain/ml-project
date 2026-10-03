"""Create descriptive M2 summaries for the six TJU Pan promoter tables."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from experiments.prepare_tjupan_motif_inputs import REQUIRED_COLUMNS, SPECIES
from preprocessing.provenance import sha256_file


def _load_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = sorted(set(REQUIRED_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing columns: {missing}")
    if frame.empty:
        raise ValueError(f"{path} is empty")
    frame["seq"] = frame["seq"].str.strip().str.upper()
    frame["label"] = frame["label"].str.strip()
    if not frame["label"].isin(["0", "1"]).all():
        raise ValueError(f"{path} contains labels other than 0 or 1")
    invalid = ~frame["seq"].str.fullmatch(r"[ACGT]+")
    if invalid.any():
        raise ValueError(f"{path} contains {int(invalid.sum())} non-ACGT or empty sequences")
    frame["sequence_length"] = frame["seq"].str.len()
    frame["gc_fraction"] = frame["seq"].map(
        lambda sequence: (sequence.count("G") + sequence.count("C")) / len(sequence)
    )
    return frame


def _information_rows(slug: str, positives: pd.DataFrame) -> list[dict]:
    rows = []
    width = int(positives["sequence_length"].max())
    for position in range(width):
        counts = Counter(sequence[position] for sequence in positives["seq"] if len(sequence) > position)
        total = sum(counts.values())
        fractions = {base: counts[base] / total if total else 0.0 for base in "ACGT"}
        entropy = -sum(value * math.log2(value) for value in fractions.values() if value > 0)
        rows.append(
            {
                "species": slug,
                "position_1based": position + 1,
                "n_sequences": total,
                "a_fraction": fractions["A"],
                "c_fraction": fractions["C"],
                "g_fraction": fractions["G"],
                "t_fraction": fractions["T"],
                "shannon_entropy_bits": entropy,
                "information_bits_vs_uniform": 2.0 - entropy,
            }
        )
    return rows


def run(dataset_root: Path, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    composition_rows = []
    position_rows = []
    overlap_rows = []
    source_hashes = {}
    split_hashes = {}

    for slug, directory in SPECIES.items():
        species_dir = dataset_root / directory
        dataset_path = species_dir / "Dataset.csv"
        frame = _load_table(dataset_path)
        source_hashes[str(dataset_path)] = sha256_file(dataset_path)
        positives = frame.loc[frame["label"].eq("1")]
        positive_path = species_dir / "positive_samples.csv"
        positive_file = _load_table(positive_path)
        source_hashes[str(positive_path)] = sha256_file(positive_path)
        if set(positives["seq_id"]) != set(positive_file["seq_id"]):
            raise ValueError(f"{directory}: positive_samples.csv IDs differ from label=1")

        for label, label_name in (("1", "positive"), ("0", "negative")):
            subset = frame.loc[frame["label"].eq(label)]
            summary_rows.append(
                {
                    "species": slug,
                    "class": label_name,
                    "records": len(subset),
                    "unique_sequences": subset["seq"].nunique(),
                    "duplicate_sequence_records": len(subset) - subset["seq"].nunique(),
                    "length_min": int(subset["sequence_length"].min()),
                    "length_median": float(subset["sequence_length"].median()),
                    "length_max": int(subset["sequence_length"].max()),
                    "gc_fraction_mean": float(subset["gc_fraction"].mean()),
                    "gc_fraction_median": float(subset["gc_fraction"].median()),
                }
            )
            for sequence_type, count in subset["seq_type"].value_counts().items():
                composition_rows.append(
                    {"species": slug, "class": label_name, "seq_type": sequence_type, "records": int(count)}
                )

        position_rows.extend(_information_rows(slug, positives))
        sequence_splits: dict[str, set[str]] = defaultdict(set)
        for split in ("train", "dev", "test"):
            split_path = species_dir / f"{split}.csv"
            split_frame = _load_table(split_path)
            split_hashes[str(split_path)] = sha256_file(split_path)
            for sequence in split_frame["seq"]:
                sequence_splits[sequence].add(split)
        split_overlaps = [splits for splits in sequence_splits.values() if len(splits) > 1]
        overlap_rows.append(
            {
                "species": slug,
                "exact_sequences_in_multiple_splits": len(split_overlaps),
                "split_rule": "train/dev/test are partitions of one source dataset, not independent datasets",
            }
        )

    summaries = pd.DataFrame(summary_rows)
    composition = pd.DataFrame(composition_rows)
    positions = pd.DataFrame(position_rows)
    overlaps = pd.DataFrame(overlap_rows)
    outputs = {
        "species_class_summary.tsv": summaries,
        "sequence_type_counts.tsv": composition,
        "positive_position_composition.tsv": positions,
        "split_overlap_summary.tsv": overlaps,
    }
    for filename, frame in outputs.items():
        frame.to_csv(output_dir / filename, sep="\t", index=False)

    _plot_class_counts(summaries, output_dir / "class_counts_by_species.png")
    _plot_gc(summaries, output_dir / "gc_distribution_by_species.png")
    _plot_position_information(positions, output_dir / "positive_position_information.png")

    output_hashes = {
        path.name: sha256_file(path)
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path.name != "run_manifest.json"
    }
    manifest = {
        "schema_version": "tjupan-m2-eda-1.0",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "source": "TJU course cloud drive; reg_and_gen/Datasets",
        "input_files": source_hashes | split_hashes,
        "species_count": len(SPECIES),
        "record_count": int(summaries["records"].sum()),
        "positive_record_count": int(summaries.loc[summaries["class"].eq("positive"), "records"].sum()),
        "negative_record_count": int(summaries.loc[summaries["class"].eq("negative"), "records"].sum()),
        "coordinate_interpretation": "sequence-index only; strand orientation has not been independently verified",
        "sigma_annotation": "not present in the TJU Pan CSV files",
        "position_statistic": "descriptive 2-H relative to uniform A/C/G/T; not bias-corrected and not an evolutionary conservation score",
        "outputs": output_hashes,
    }
    manifest_path = output_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def _plot_class_counts(summary: pd.DataFrame, path: Path) -> None:
    pivot = summary.pivot(index="species", columns="class", values="records").fillna(0)
    pivot = pivot.reindex(columns=["positive", "negative"])
    pivot.plot.bar(figsize=(9, 4), color=["#4472C4", "#ED7D31"])
    plt.ylabel("Record count")
    plt.xlabel("Species")
    plt.title("TJU Pan promoter and control records")
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def _plot_gc(summary: pd.DataFrame, path: Path) -> None:
    fig, axis = plt.subplots(figsize=(9, 4), constrained_layout=True)
    species = list(summary["species"].drop_duplicates())
    for offset, class_name, color in ((-0.16, "positive", "#4472C4"), (0.16, "negative", "#ED7D31")):
        values = [summary.loc[(summary["species"].eq(item)) & (summary["class"].eq(class_name)), "gc_fraction_mean"].iloc[0] for item in species]
        axis.scatter([i + offset for i in range(len(species))], values, label=class_name, color=color)
    axis.set_xticks(range(len(species)), species, rotation=25, ha="right")
    axis.set_ylabel("Mean GC fraction")
    axis.set_title("Mean GC fraction by species and class")
    axis.legend()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_position_information(positions: pd.DataFrame, path: Path) -> None:
    fig, axis = plt.subplots(figsize=(9, 4), constrained_layout=True)
    for species, subset in positions.groupby("species", sort=True):
        axis.plot(subset["position_1based"], subset["information_bits_vs_uniform"], label=species, linewidth=1.1)
    axis.set_xlabel("Position in sequence (1-based; orientation unverified)")
    axis.set_ylabel("Descriptive information (bits)")
    axis.set_title("Positive-class positional composition")
    axis.legend(fontsize=8, ncol=2)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", default="data/raw/tju_pan_promoter/reg_and_gen/Datasets")
    parser.add_argument("--output-dir", default="results/eda/tjupan_m2_20261003")
    args = parser.parse_args()
    manifest = run(Path(args.dataset_root), Path(args.output_dir))
    print(f"records={manifest['record_count']} positives={manifest['positive_record_count']} negatives={manifest['negative_record_count']}")
    print(f"output_dir={Path(args.output_dir).resolve()}")


if __name__ == "__main__":
    main()
