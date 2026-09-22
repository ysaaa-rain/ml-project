import json

import pandas as pd

from experiments.diagnose_embedding import run


def test_diagnose_embedding_writes_bias_and_neighbor_summaries(tmp_path):
    embedding_dir = tmp_path / "embedding"
    output_dir = tmp_path / "diagnostics"
    embedding_dir.mkdir()
    clusters = pd.DataFrame(
        {
            "sequence_id": ["a", "b", "c", "d"],
            "species": ["ecoli", "ecoli", "bs", "bs"],
            "source_dataset": ["rdb", "rdb", "dbtbs", "dbtbs"],
            "sigma_factor_type": ["Sigma70", "unknown", "SigA", "SigA"],
            "split": ["discovery", "discovery", "validation", "validation"],
            "sequence_length": [81, 81, 42, 44],
            "gc_fraction": [0.4, 0.45, 0.3, 0.35],
            "cluster": [0, 0, 1, 1],
            "pca_1": [-1.0, -0.8, 0.8, 1.0],
            "pca_2": [0.1, 0.2, -0.1, -0.2],
        }
    )
    clusters.to_csv(embedding_dir / "embedding_clusters.tsv", sep="\t", index=False)
    neighbors = pd.DataFrame(
        {
            "sequence_id": ["a", "c", "b", "d"],
            "rank": [1, 1, 1, 1],
            "neighbor_id": ["b", "d", "a", "c"],
            "cosine_distance": [0.1, 0.2, 0.1, 0.2],
        }
    )
    neighbors.to_csv(embedding_dir / "nearest_neighbors.tsv", sep="\t", index=False)

    summary = run(embedding_dir, output_dir)
    assert summary["sequence_count"] == 4
    assert summary["cluster_count"] == 2
    assert summary["nearest_neighbor"]["same_source_rate"] == 1.0
    assert (output_dir / "cluster_summary.tsv").exists()
    payload = json.loads((output_dir / "diagnostics_summary.json").read_text(encoding="utf-8"))
    assert payload["schema_version"] == "embedding-diagnostics-1.0"
