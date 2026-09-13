import json

import numpy as np
import pandas as pd
import yaml

from experiments.run_embedding import run
from models.dna_embedding import HashKmerSmokeEmbedder, build_embedder
from models.embedding_analysis import (
    cluster_embeddings,
    nearest_neighbor_table,
    occlusion_importance,
)


def test_hash_backend_is_explicit_and_deterministic():
    embedder = HashKmerSmokeEmbedder(dimension=16, kmer_size=3, seed=7)
    first = embedder.embed(["ACGTACGT", "TTTTAAAA"])
    second = embedder.embed(["ACGTACGT", "TTTTAAAA"])
    assert first.shape == (2, 16)
    assert np.allclose(first, second)
    assert embedder.metadata["scientific_use"] is False


def test_hash_backend_requires_smoke_gate():
    try:
        build_embedder({"backend": "hash-smoke"})
    except ValueError as exc:
        assert "allow_test_backend" in str(exc)
    else:
        raise AssertionError("hash-smoke must be explicitly enabled")


def test_embedding_analysis_excludes_self_and_supports_occlusion():
    embedder = HashKmerSmokeEmbedder(dimension=24, kmer_size=3, seed=11)
    sequence_ids = ["a", "b", "c"]
    matrix = embedder.embed(["AAAACCCC", "AAAACCCA", "GGGGTTTT"])
    clusters, _ = cluster_embeddings(matrix, n_clusters=2, seed=11)
    neighbors = nearest_neighbor_table(matrix, sequence_ids, neighbors=1)
    windows = occlusion_importance(
        embedder,
        sequence_id="a",
        sequence="AAAACCCC",
        window_size=4,
        stride=2,
    )
    assert len(clusters) == 3
    assert len(neighbors) == 3
    assert not (neighbors["sequence_id"] == neighbors["neighbor_id"]).any()
    assert windows
    assert all(window.end > window.start for window in windows)


def test_embedding_pipeline_writes_traceable_smoke_outputs(tmp_path):
    input_path = tmp_path / "promoters_clean.tsv"
    output_dir = tmp_path / "embedding_output"
    frame = pd.DataFrame(
        {
            "sequence_id": ["s1", "s2", "s3", "s4"],
            "species": ["a", "a", "b", "b"],
            "source_dataset": ["discovery", "discovery", "validation", "validation"],
            "sigma_factor_type": ["sigma70", "sigma70", "sigmaB", "sigmaB"],
            "split": ["discovery", "discovery", "validation", "validation"],
            "sequence": ["AAAACCCCGGGG", "AAAACCCCGGGA", "TTTTGGGGCCCC", "TTTTGGGGCCCA"],
        }
    )
    frame.to_csv(input_path, sep="\t", index=False)
    config = {
        "input": str(input_path),
        "output_dir": str(output_dir),
        "embedding": {
            "backend": "hash-smoke",
            "allow_test_backend": True,
            "dimension": 16,
            "kmer_size": 3,
            "seed": 17,
        },
        "analysis": {"n_clusters": 2, "neighbors": 1, "seed": 17},
        "occlusion": {"enabled": True, "max_sequences": 1, "window_size": 4, "stride": 2},
    }
    config_path = tmp_path / "embedding.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    manifest = run(config_path)
    summary = json.loads((output_dir / "embedding_summary.json").read_text(encoding="utf-8"))
    assert manifest["summary"]["evidence_level"] == "smoke_only"
    assert summary["embedding_shape"] == [4, 16]
    assert (output_dir / "sequence_embeddings.npy").exists()
    assert (output_dir / "nearest_neighbors.tsv").exists()
    assert (output_dir / "occlusion_windows.tsv").exists()
