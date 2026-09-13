"""Tests for near-homolog clustering and leakage-free split assignment."""

import random
from pathlib import Path

import pandas as pd
import pytest

from preprocessing.homology import (
    cluster_sequences,
    homology_aware_split,
    kmer_sketch,
    pair_identity,
    write_cluster_report,
)


def _random_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _mutate(sequence: str, *, positions: int, seed: int) -> str:
    """Change exactly ``positions`` bases, always to a different base."""

    rng = random.Random(seed)
    bases = list(sequence)
    chosen = rng.sample(range(len(bases)), positions)
    for index in chosen:
        bases[index] = rng.choice([base for base in "ACGT" if base != bases[index]])
    return "".join(bases)


def test_pair_identity_is_one_for_identical_sequences():
    sequence = _random_dna(100, 1)
    assert pair_identity(sequence, sequence) == 1.0


def test_pair_identity_tolerates_a_truncated_end():
    sequence = _random_dna(120, 2)
    # Same sequence with 4 bases missing at the end, as two records describing
    # the same promoter often differ. Identity is scored over the longer length.
    truncated = sequence[:-4]
    assert pair_identity(sequence, truncated) == pytest.approx(116 / 120)


def test_kmer_sketch_is_strand_invariant():
    sequence = _random_dna(60, 3)
    complement = str.maketrans("ACGT", "TGCA")
    reverse_complement = sequence.translate(complement)[::-1]
    assert kmer_sketch(sequence, 15) == kmer_sketch(reverse_complement, 15)


def test_cluster_groups_near_duplicates_and_keeps_distinct_sequences_separate():
    base = _random_dna(200, 10)
    near = _mutate(base, positions=5, seed=11)      # 97.5% identity
    near2 = _mutate(base, positions=8, seed=12)     # 96% identity
    far = _random_dna(200, 13)                      # unrelated

    frame = pd.DataFrame(
        {
            "sequence_id": ["base", "near", "near2", "far"],
            "species": ["E", "E", "E", "B"],
            "source_dataset": ["regulondb", "regulondb", "regulondb", "dbtbs"],
            "sequence": [base, near, near2, far],
        }
    )
    clustered, summary = cluster_sequences(frame, identity_threshold=0.90)

    groups = clustered.groupby("cluster_id")["sequence_id"].apply(set).to_dict()
    assert len(groups) == 2, groups
    assert {"base", "near", "near2"} in groups.values()
    assert {"far"} in groups.values()

    assert summary["cluster_count"] == 2
    assert summary["multi_member_cluster_count"] == 3
    assert summary["largest_cluster_size"] == 3
    assert summary["homologous_pair_count"] >= 3


def test_clustering_is_deterministic_and_uses_stable_representatives():
    base = _random_dna(150, 20)
    near = _mutate(base, positions=4, seed=21)
    frame = pd.DataFrame(
        {
            "sequence_id": ["zzz", "aaa"],
            "species": ["E", "E"],
            "source_dataset": ["regulondb", "regulondb"],
            "sequence": [base, near],
        }
    )
    first, first_summary = cluster_sequences(frame, identity_threshold=0.90)
    second, second_summary = cluster_sequences(frame, identity_threshold=0.90)
    assert first_summary == second_summary
    assert list(first["cluster_id"]) == list(second["cluster_id"])
    # The representative is the smallest sequence_id in the cluster.
    assert set(first["cluster_id"]) == {"aaa"}


def test_homology_aware_split_removes_cross_source_leakage():
    base = _random_dna(200, 30)
    near = _mutate(base, positions=5, seed=31)
    unrelated = _random_dna(200, 32)

    frame = pd.DataFrame(
        {
            "sequence_id": ["a_reg", "a_dbt", "b_reg", "b_dbt"],
            "species": ["E", "E", "B", "B"],
            "source_dataset": ["regulondb", "dbtbs", "regulondb", "dbtbs"],
            # a_reg and a_dbt are near-identical but come from different sources:
            # this is exactly the leakage case the clustering must catch.
            "sequence": [base, near, unrelated, _mutate(unrelated, positions=5, seed=33)],
            "split": ["discovery", "validation", "discovery", "validation"],
        }
    )
    clustered, _ = cluster_sequences(frame, identity_threshold=0.90)
    assert clustered.groupby("cluster_id")["source_dataset"].nunique().max() == 2

    split_frame, report = homology_aware_split(clustered)

    # No cluster may straddle the two sides any more.
    assert split_frame.groupby("cluster_id")["split"].nunique().max() == 1
    assert report["leakage_violations_after"] == 0
    assert report["clusters_spanning_sides_before"] == 2
    assert report["clusters_reassigned"] == 2

    # The reassigned clusters must be listed for review.
    reassigned_ids = {entry["cluster_id"] for entry in report["reassigned_clusters"]}
    assert reassigned_ids == set(split_frame["cluster_id"])


def test_homology_aware_split_leaves_clean_clusters_untouched():
    frame = pd.DataFrame(
        {
            "sequence_id": ["x", "y"],
            "sequence": [_random_dna(120, 40), _random_dna(120, 41)],
            "source_dataset": ["regulondb", "dbtbs"],
            "split": ["discovery", "validation"],
            "cluster_id": ["x", "y"],
        }
    )
    split_frame, report = homology_aware_split(frame)
    assert list(split_frame["split"]) == ["discovery", "validation"]
    assert report["clusters_reassigned"] == 0
    assert report["leakage_violations_after"] == 0


def test_cluster_report_lists_multi_source_clusters(tmp_path):
    base = _random_dna(200, 50)
    near = _mutate(base, positions=4, seed=51)
    frame = pd.DataFrame(
        {
            "sequence_id": ["s1", "s2"],
            "species": ["E", "B"],
            "source_dataset": ["regulondb", "dbtbs"],
            "sequence": [base, near],
            "split": ["discovery", "validation"],
        }
    )
    clustered, _ = cluster_sequences(frame, identity_threshold=0.90)
    info = write_cluster_report(clustered, tmp_path / "clusters.tsv")

    assert info["cluster_count"] == 1
    assert info["clusters_spanning_multiple_sources"] == 1
    text = (tmp_path / "clusters.tsv").read_text(encoding="utf-8")
    # Sources are listed sorted, so the cluster spanning both is visible.
    assert "dbtbs,regulondb" in text


def test_cluster_sequences_requires_sequence_column():
    frame = pd.DataFrame({"sequence_id": ["a"], "species": ["E"], "source_dataset": ["s"]})
    with pytest.raises(ValueError):
        cluster_sequences(frame)


def test_identity_threshold_is_respected():
    base = _random_dna(200, 60)
    # 20 mutations in 200 bp: identity is exactly 0.90.
    borderline = _mutate(base, positions=20, seed=61)
    assert pair_identity(base, borderline) == pytest.approx(0.90)

    frame = pd.DataFrame(
        {
            "sequence_id": ["a", "b"],
            "species": ["E", "E"],
            "source_dataset": ["regulondb", "regulondb"],
            "sequence": [base, borderline],
        }
    )
    # The threshold is inclusive, so 0.90 joins at 0.90 but not at 0.91.
    joined, joined_summary = cluster_sequences(frame, identity_threshold=0.90)
    strict, strict_summary = cluster_sequences(frame, identity_threshold=0.91)
    loose, loose_summary = cluster_sequences(frame, identity_threshold=0.80)

    assert strict_summary["cluster_count"] == 2
    assert joined_summary["cluster_count"] == 1
    assert loose_summary["cluster_count"] == 1
    assert strict["cluster_id"].nunique() == 2
    assert joined["cluster_id"].nunique() == 1


def test_clustering_through_normalize_metadata_prevents_leakage():
    from preprocessing.pipeline import normalize_metadata

    base = _random_dna(200, 70)
    near = _mutate(base, positions=5, seed=71)
    frame = pd.DataFrame(
        {
            "sequence_id": ["p1", "p2"],
            "species": ["E", "B"],
            "source_dataset": ["regulondb", "dbtbs"],
            "sequence": [base, near],
        }
    )
    clean, _ = normalize_metadata(
        frame,
        validation_sources=["dbtbs"],
        homology_clustering=True,
        identity_threshold=0.90,
    )
    # Default splitting would put these on opposite sides; clustering must not.
    assert clean["split"].nunique() == 1
    report = clean.attrs["homology_split"]
    assert report["leakage_violations_after"] == 0
    assert report["clusters_reassigned"] == 1


def test_clustering_disabled_by_default_keeps_legacy_behaviour():
    from preprocessing.pipeline import normalize_metadata

    frame = pd.DataFrame(
        {
            "sequence_id": ["p1", "p2"],
            "species": ["E", "B"],
            "source_dataset": ["regulondb", "dbtbs"],
            "sequence": [_random_dna(120, 80), _random_dna(120, 81)],
        }
    )
    clean, _ = normalize_metadata(frame, validation_sources=["dbtbs"])
    assert set(clean["split"]) == {"discovery", "validation"}
    assert clean.attrs.get("homology_clustering") is None
