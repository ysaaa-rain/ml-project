"""Run the L0-L2 DNA embedding pipeline from a reproducible YAML config."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from models.dna_embedding import build_embedder
from models.embedding_analysis import (
    cluster_embeddings,
    label_alignment,
    nearest_neighbor_table,
    occlusion_importance,
    pca_coordinates,
    windows_to_frame,
)
from preprocessing.provenance import sha256_file


REQUIRED_COLUMNS = ("sequence_id", "sequence", "species", "source_dataset", "split")


def _resolve(base_dir: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (base_dir / path).resolve()


def _validate_input(frame: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"clean embedding input is missing columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError("embedding input is empty")
    invalid = frame["sequence"].astype(str).map(lambda value: bool(set(value.upper()) - set("ACGTN")))
    if invalid.any():
        raise ValueError("embedding input contains non-DNA characters; run preprocessing first")
    if frame["sequence_id"].astype(str).duplicated().any():
        raise ValueError("embedding input contains duplicate sequence_id values")


def run(config_path: str | Path) -> dict:
    config_path = Path(config_path).resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    base_dir = config_path.parent.parent if config_path.parent.name == "configs" else config_path.parent
    input_path = _resolve(base_dir, str(config["input"]))
    output_dir = _resolve(base_dir, str(config["output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(input_path, sep="\t")
    _validate_input(frame)
    embedder = build_embedder(config.get("embedding", {}))
    sequence_ids = frame["sequence_id"].astype(str).tolist()
    sequences = frame["sequence"].astype(str).str.upper().tolist()
    embeddings = embedder.embed(sequences)
    if embeddings.ndim != 2 or embeddings.shape[0] != len(frame):
        raise RuntimeError("embedding backend returned an invalid matrix shape")

    np.save(output_dir / "sequence_embeddings.npy", embeddings.astype(np.float32))
    embeddings_sha256 = sha256_file(output_dir / "sequence_embeddings.npy")
    index_columns = [
        column
        for column in (
            "sequence_id",
            "species",
            "source_dataset",
            "sigma_factor_type",
            "evidence_level",
            "split",
            "sequence_length",
            "gc_fraction",
        )
        if column in frame
    ]
    frame[index_columns].to_csv(output_dir / "embedding_index.tsv", sep="\t", index=False)

    analysis = config.get("analysis", {})
    seed = int(analysis.get("seed", 20260913))
    clusters, silhouette = cluster_embeddings(
        embeddings,
        n_clusters=int(analysis.get("n_clusters", 3)),
        seed=seed,
    )
    coordinates = pca_coordinates(embeddings)
    cluster_frame = frame[index_columns].copy()
    cluster_frame["cluster"] = clusters
    cluster_frame["pca_1"] = coordinates[:, 0]
    cluster_frame["pca_2"] = coordinates[:, 1]
    cluster_frame.to_csv(output_dir / "embedding_clusters.tsv", sep="\t", index=False)
    neighbors = nearest_neighbor_table(
        embeddings,
        sequence_ids,
        neighbors=int(analysis.get("neighbors", 5)),
    )
    neighbors.to_csv(output_dir / "nearest_neighbors.tsv", sep="\t", index=False)

    alignment = {}
    for column in ("sigma_factor_type", "species", "source_dataset", "split"):
        if column in frame:
            alignment[column] = label_alignment(clusters, frame[column].fillna("unknown").astype(str))

    occlusion_config = config.get("occlusion", {})
    windows = []
    if occlusion_config.get("enabled", False):
        maximum = min(int(occlusion_config.get("max_sequences", 10)), len(frame))
        for row in frame.iloc[:maximum].itertuples():
            windows.extend(
                occlusion_importance(
                    embedder,
                    sequence_id=str(row.sequence_id),
                    sequence=str(row.sequence),
                    window_size=int(occlusion_config.get("window_size", 12)),
                    stride=int(occlusion_config.get("stride", 3)),
                    mask_base=str(occlusion_config.get("mask_base", "N")),
                )
            )
    windows_to_frame(windows).to_csv(output_dir / "occlusion_windows.tsv", sep="\t", index=False)

    scientific_use = bool(embedder.metadata.get("scientific_use"))
    summary = {
        "schema_version": "embedding-mainline-1.0",
        "evidence_level": "model_inference_pending_biological_validation" if scientific_use else "smoke_only",
        "sequence_count": len(frame),
        "embedding_shape": list(embeddings.shape),
        "silhouette": silhouette,
        "label_alignment_diagnostic": alignment,
        "occlusion_window_count": len(windows),
        "warning": None if scientific_use else "Smoke backend output is not scientific evidence.",
    }
    (output_dir / "embedding_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "config": str(config_path),
        "config_sha256": sha256_file(config_path),
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "embedding": embedder.metadata,
        "embeddings_sha256": embeddings_sha256,
        "analysis": {"seed": seed, **analysis},
        "summary": summary,
        "outputs": {
            "embeddings": "sequence_embeddings.npy",
            "index": "embedding_index.tsv",
            "clusters": "embedding_clusters.tsv",
            "neighbors": "nearest_neighbors.tsv",
            "occlusion": "occlusion_windows.tsv",
            "summary": "embedding_summary.json",
        },
    }
    (output_dir / "embedding_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML configuration file")
    args = parser.parse_args()
    manifest = run(args.config)
    summary = manifest["summary"]
    print(f"evidence_level={summary['evidence_level']}")
    print(f"sequence_count={summary['sequence_count']}")
    print(f"embedding_shape={summary['embedding_shape']}")


if __name__ == "__main__":
    main()
