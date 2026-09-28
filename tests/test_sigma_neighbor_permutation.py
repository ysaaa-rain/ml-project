import pandas as pd
import pytest

from experiments.sigma_neighbor_permutation import run


def test_permutation_is_reproducible_and_keeps_known_edge_count(tmp_path):
    source = tmp_path / "embedding"
    source.mkdir()
    pd.DataFrame({
        "sequence_id": ["a", "b", "c", "d", "e"],
        "sigma_factor_type": ["S1", "S1", "S2", "S2", "unknown"],
        "gc_fraction": [0.4, 0.4, 0.6, 0.6, 0.5],
    }).to_csv(source / "embedding_clusters.tsv", sep="\t", index=False)
    pd.DataFrame({
        "sequence_id": ["a", "b", "c", "d", "a"],
        "neighbor_id": ["b", "a", "d", "c", "e"],
    }).to_csv(source / "nearest_neighbors.tsv", sep="\t", index=False)
    first = run(source, tmp_path / "one.json", permutations=40, seed=4)
    second = run(source, tmp_path / "two.json", permutations=40, seed=4)
    assert first == second
    assert first["known_sigma_directed_edge_count"] == 4
    assert first["observed_same_sigma_rate"] == 1.0
    assert first["global_label_permutation"]["p_one_sided_enrichment"] >= 1 / 41


def test_permutation_rejects_unknown_neighbor(tmp_path):
    source = tmp_path / "embedding"
    source.mkdir()
    pd.DataFrame({"sequence_id": ["a", "b"], "sigma_factor_type": ["S1", "S2"], "gc_fraction": [0.5, 0.5]}).to_csv(
        source / "embedding_clusters.tsv", sep="\t", index=False
    )
    pd.DataFrame({"sequence_id": ["a"], "neighbor_id": ["missing"]}).to_csv(
        source / "nearest_neighbors.tsv", sep="\t", index=False
    )
    with pytest.raises(ValueError, match="unknown sequence_id"):
        run(source, tmp_path / "out.json", permutations=10)
