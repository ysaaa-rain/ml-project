"""Permutation diagnostic for sigma-label agreement in a fixed kNN graph."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(embedding_dir: str | Path, output: str | Path, *, permutations: int = 2000, seed: int = 20260928) -> dict:
    if permutations < 1:
        raise ValueError("permutations must be positive")
    embedding_dir = Path(embedding_dir)
    clusters_path = embedding_dir / "embedding_clusters.tsv"
    neighbors_path = embedding_dir / "nearest_neighbors.tsv"
    frame = pd.read_csv(clusters_path, sep="\t")
    neighbors = pd.read_csv(neighbors_path, sep="\t")
    required = {"sequence_id", "sigma_factor_type", "gc_fraction"}
    if not required <= set(frame.columns):
        raise ValueError(f"embedding_clusters.tsv missing {sorted(required - set(frame.columns))}")
    if not {"sequence_id", "neighbor_id"} <= set(neighbors.columns):
        raise ValueError("nearest_neighbors.tsv missing sequence_id or neighbor_id")
    if frame.sequence_id.isna().any() or frame.sequence_id.duplicated().any():
        raise ValueError("sequence_id must be unique and nonempty")
    ids = frame.sequence_id.astype(str)
    lookup = pd.Series(np.arange(len(frame)), index=ids)
    if not neighbors.sequence_id.isin(ids).all() or not neighbors.neighbor_id.isin(ids).all():
        raise ValueError("nearest_neighbors.tsv references an unknown sequence_id")
    source = lookup.loc[neighbors.sequence_id].to_numpy()
    target = lookup.loc[neighbors.neighbor_id].to_numpy()
    labels = frame.sigma_factor_type.fillna("unknown").astype(str).to_numpy()
    known = (labels != "unknown") & (labels != "")
    valid_edge = known[source] & known[target]
    source, target = source[valid_edge], target[valid_edge]
    if len(source) == 0 or known.sum() < 2:
        raise ValueError("no neighbor edges connect two known sigma labels")
    gc = pd.to_numeric(frame.gc_fraction, errors="coerce").to_numpy()
    if not np.isfinite(gc[known]).all():
        raise ValueError("known sigma records require finite gc_fraction")
    # Duplicate directed kNN edges are retained; labels are permuted per sequence,
    # so shared vertices and reciprocal edges remain dependent in the null model.
    observed = float(np.mean(labels[source] == labels[target]))
    rng = np.random.default_rng(seed)
    known_idx = np.flatnonzero(known)
    bins = np.asarray(pd.qcut(gc[known], q=min(10, len(known_idx)), labels=False, duplicates="drop"))
    groups = [known_idx[bins == value] for value in np.unique(bins) if pd.notna(value)]
    if not groups:
        groups = [known_idx]
    null_global = np.empty(permutations)
    null_gc = np.empty(permutations)
    for iteration in range(permutations):
        shuffled = labels.copy()
        shuffled[known_idx] = rng.permutation(labels[known_idx])
        null_global[iteration] = np.mean(shuffled[source] == shuffled[target])
        shuffled = labels.copy()
        for group in groups:
            shuffled[group] = rng.permutation(labels[group])
        null_gc[iteration] = np.mean(shuffled[source] == shuffled[target])

    def summarize(null: np.ndarray) -> dict:
        return {
            "mean": float(null.mean()),
            "percentile_2_5": float(np.quantile(null, 0.025)),
            "percentile_97_5": float(np.quantile(null, 0.975)),
            "observed_minus_null_mean": float(observed - null.mean()),
            "p_one_sided_enrichment": float((1 + np.count_nonzero(null >= observed)) / (permutations + 1)),
        }

    result = {
        "schema_version": "sigma-neighbor-permutation-1.0",
        "sequence_count": int(len(frame)),
        "known_sigma_sequence_count": int(known.sum()),
        "known_sigma_directed_edge_count": int(len(source)),
        "observed_same_sigma_rate": observed,
        "permutations": permutations,
        "seed": seed,
        "gc_strata_count": len(groups),
        "global_label_permutation": summarize(null_global),
        "gc_stratified_label_permutation": summarize(null_gc),
        "input_sha256": {clusters_path.name: _sha256(clusters_path), neighbors_path.name: _sha256(neighbors_path)},
        "boundary": "Fixed-graph within-source label diagnostic; not cluster stability, independent validation, or a biological conclusion.",
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embedding-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--permutations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260928)
    args = parser.parse_args()
    print(json.dumps(run(args.embedding_dir, args.output, permutations=args.permutations, seed=args.seed), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
