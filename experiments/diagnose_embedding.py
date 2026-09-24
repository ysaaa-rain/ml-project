"""Run descriptive bias and stability diagnostics for an embedding run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def _safe_spearman(left: pd.Series, right: pd.Series) -> float | None:
    values = pd.concat([left, right], axis=1).apply(pd.to_numeric, errors="coerce").dropna()
    if len(values) < 2 or values.iloc[:, 0].nunique() < 2 or values.iloc[:, 1].nunique() < 2:
        return None
    result = spearmanr(values.iloc[:, 0], values.iloc[:, 1])
    return float(result.statistic)


def _categorical_alignment(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for column in ("sigma_factor_type", "species", "source_dataset", "split", "evidence_level"):
        if column not in frame:
            continue
        labels = frame[column].fillna("unknown").astype(str)
        if labels.nunique() < 2:
            continue
        rows.append(
            {
                "label": column,
                "n_labels": int(labels.nunique()),
                "adjusted_rand_index": float(adjusted_rand_score(labels, frame["cluster"])),
                "normalized_mutual_info": float(normalized_mutual_info_score(labels, frame["cluster"])),
            }
        )
    return pd.DataFrame(
        rows,
        columns=["label", "n_labels", "adjusted_rand_index", "normalized_mutual_info"],
    )


def _numeric_correlations(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for feature in ("sequence_length", "gc_fraction"):
        if feature not in frame:
            continue
        for coordinate in ("pca_1", "pca_2"):
            if coordinate not in frame:
                continue
            left = pd.to_numeric(frame[feature], errors="coerce")
            right = pd.to_numeric(frame[coordinate], errors="coerce")
            valid = pd.concat([left, right], axis=1).dropna()
            pearson = None
            if len(valid) >= 2 and valid.iloc[:, 0].nunique() >= 2 and valid.iloc[:, 1].nunique() >= 2:
                pearson = float(np.corrcoef(valid.iloc[:, 0], valid.iloc[:, 1])[0, 1])
            rows.append(
                {
                    "feature": feature,
                    "coordinate": coordinate,
                    "n": int(len(valid)),
                    "pearson": pearson,
                    "spearman": _safe_spearman(left, right),
                }
            )
    return pd.DataFrame(rows, columns=["feature", "coordinate", "n", "pearson", "spearman"])


def _nearest_neighbor_summary(frame: pd.DataFrame, neighbors: pd.DataFrame) -> dict:
    if neighbors.empty:
        return {"row_count": 0}
    lookup = frame.set_index("sequence_id")
    rows = neighbors.copy()
    source = lookup["source_dataset"].astype(str)
    species = lookup["species"].astype(str)
    sigma = lookup["sigma_factor_type"].fillna("unknown").astype(str)
    cluster = lookup["cluster"]
    rows["same_source"] = [source.get(a) == source.get(b) for a, b in zip(rows.sequence_id, rows.neighbor_id)]
    rows["same_species"] = [species.get(a) == species.get(b) for a, b in zip(rows.sequence_id, rows.neighbor_id)]
    rows["same_cluster"] = [cluster.get(a) == cluster.get(b) for a, b in zip(rows.sequence_id, rows.neighbor_id)]
    known_sigma = [
        sigma.get(a) not in {None, "", "unknown"} and sigma.get(b) not in {None, "", "unknown"}
        for a, b in zip(rows.sequence_id, rows.neighbor_id)
    ]
    rows["known_sigma_pair"] = known_sigma
    rows["same_known_sigma"] = [
        bool(known and sigma.get(a) == sigma.get(b))
        for known, a, b in zip(known_sigma, rows.sequence_id, rows.neighbor_id)
    ]
    summary = {
        "row_count": int(len(rows)),
        "mean_cosine_distance": float(rows["cosine_distance"].mean()),
        "same_source_rate": float(rows["same_source"].mean()),
        "same_species_rate": float(rows["same_species"].mean()),
        "same_cluster_rate": float(rows["same_cluster"].mean()),
        "known_sigma_pair_count": int(rows["known_sigma_pair"].sum()),
        "same_known_sigma_rate": (
            float(rows.loc[rows["known_sigma_pair"], "same_known_sigma"].mean())
            if rows["known_sigma_pair"].any()
            else None
        ),
    }
    return summary


def run(embedding_dir: str | Path, output_dir: str | Path) -> dict:
    embedding_dir = Path(embedding_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    clusters_path = embedding_dir / "embedding_clusters.tsv"
    neighbors_path = embedding_dir / "nearest_neighbors.tsv"
    frame = pd.read_csv(clusters_path, sep="\t")
    neighbors = pd.read_csv(neighbors_path, sep="\t")
    required = {"sequence_id", "cluster", "sequence_length", "gc_fraction", "pca_1", "pca_2"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"embedding_clusters.tsv is missing columns: {', '.join(missing)}")
    if frame["sequence_id"].duplicated().any():
        raise ValueError("embedding_clusters.tsv contains duplicate sequence_id values")

    cluster_summary = (
        frame.groupby("cluster", sort=True)
        .agg(
            sequence_count=("sequence_id", "size"),
            length_mean=("sequence_length", "mean"),
            length_median=("sequence_length", "median"),
            gc_mean=("gc_fraction", "mean"),
            gc_median=("gc_fraction", "median"),
        )
        .reset_index()
    )
    cluster_summary.to_csv(output_dir / "cluster_summary.tsv", sep="\t", index=False)
    alignment = _categorical_alignment(frame)
    alignment.to_csv(output_dir / "categorical_alignment.tsv", sep="\t", index=False)
    correlations = _numeric_correlations(frame)
    correlations.to_csv(output_dir / "numeric_correlations.tsv", sep="\t", index=False)
    neighbor_summary = _nearest_neighbor_summary(frame, neighbors)

    summary = {
        "schema_version": "embedding-diagnostics-1.0",
        "sequence_count": int(len(frame)),
        "cluster_count": int(frame["cluster"].nunique()),
        "cluster_summary": cluster_summary.to_dict(orient="records"),
        "categorical_alignment": alignment.to_dict(orient="records"),
        "numeric_correlations": correlations.replace({np.nan: None}).to_dict(orient="records"),
        "nearest_neighbor": neighbor_summary,
        "boundary": "Descriptive diagnostics only; this run does not establish cross-source reproducibility or biological significance.",
    }
    (output_dir / "diagnostics_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embedding-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    summary = run(args.embedding_dir, args.output_dir)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
