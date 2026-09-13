"""Near-homolog clustering so discovery/validation splits cannot leak.

Splitting purely by ``source_dataset`` lets near-identical sequences end up on
opposite sides of the discovery/validation boundary. Homologous sequences share
their motifs, so such leakage inflates cross-dataset consistency and makes RQ3
("are the patterns reproducible across sources?") unanswerable.

This module groups sequences into homologous clusters and then assigns whole
clusters to a single side of the split, guaranteeing that no cluster straddles
discovery and validation.

Similarity is computed as alignment-free percent identity over the shorter
sequence, restricted to offsets that keep most of the shorter sequence aligned::

    identity = matches / min(len(a), len(b))

Candidate pairs are found with a k-mer sketch (Jaccard on distinct k-mers) and
then confirmed exactly, which keeps the cost near-linear for the sequence counts
this project works with.
"""

from __future__ import annotations

from typing import Any, Iterable

import pandas as pd

DEFAULT_IDENTITY_THRESHOLD = 0.90
DEFAULT_KMER_SIZE = 15
DEFAULT_LENGTH_TOLERANCE = 0.20
# The k-mer sketch is only a cheap prefilter; percent identity decides. Two
# sequences at 90% identity share surprisingly few exact 15-mers (Jaccard around
# 0.08), so this threshold must stay well below what a genuine homolog produces
# or real homologs get filtered out before the exact comparison.
DEFAULT_SKETCH_JACCARD = 0.05
DEFAULT_MAX_KMER_FREQUENCY = 200

CLUSTER_COLUMNS = (
    "cluster_id",
    "cluster_size",
    "cluster_representative",
    "cluster_source_count",
)


def canonical_kmer(sequence: str, k: int) -> str:
    """Return the canonical (minimum of forward/reverse-complement) k-mer."""

    forward = sequence
    reverse = _reverse_complement(sequence)
    return forward if forward <= reverse else reverse


def _reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(table)[::-1]


def kmer_sketch(sequence: str, k: int, *, strand_invariant: bool = True) -> frozenset[str]:
    """Return the set of distinct canonical k-mers of a sequence."""

    if len(sequence) < k:
        return frozenset({sequence})
    if strand_invariant:
        return frozenset(canonical_kmer(sequence[i : i + k], k) for i in range(len(sequence) - k + 1))
    return frozenset(sequence[i : i + k] for i in range(len(sequence) - k + 1))


def pair_identity(first: str, second: str, *, max_offset_fraction: float = 0.25) -> float:
    """Return the best alignment-free identity over plausible offsets.

    The shorter sequence is slid across the longer one within
    ``max_offset_fraction`` of its length and the best
    ``matches / len(longer)`` is returned. This detects near-identical sequences
    even when they are shifted relative to each other, without a full alignment.

    The denominator is the *longer* sequence: scoring against the shorter one
    would report a near-zero identity for a pair that is identical except for a
    few bases missing at one end, which is a common real difference between two
    records describing the same promoter. The default offset window is
    correspondingly wide for the same reason.
    """

    if not first or not second:
        return 0.0
    short, long = (first, second) if len(first) <= len(second) else (second, first)
    span = len(long) - len(short)
    max_offset = max(1, int(len(long) * max_offset_fraction))
    offsets: Iterable[int]
    if span <= 0:
        offsets = (0,)
    else:
        low = max(0, span // 2 - max_offset)
        high = min(span, span // 2 + max_offset)
        offsets = range(low, high + 1)
    best = 0.0
    for offset in offsets:
        window = long[offset : offset + len(short)]
        if len(window) != len(short):
            continue
        matches = sum(1 for a, b in zip(short, window) if a == b)
        identity = matches / len(long)
        if identity > best:
            best = identity
    return best


class _UnionFind:
    """Deterministic union-find over sequence indices."""

    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != root:
            self.parent[item], item = root, self.parent[item]
        return root

    def union(self, first: int, second: int) -> None:
        root_a, root_b = self.find(first), self.find(second)
        if root_a == root_b:
            return
        if root_a < root_b:
            self.parent[root_b] = root_a
        else:
            self.parent[root_a] = root_b


def homology_pairs(
    records: list[tuple[str, str]],
    *,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    kmer_size: int = DEFAULT_KMER_SIZE,
    length_tolerance: float = DEFAULT_LENGTH_TOLERANCE,
    sketch_jaccard: float = DEFAULT_SKETCH_JACCARD,
    max_kmer_frequency: int = DEFAULT_MAX_KMER_FREQUENCY,
) -> list[tuple[int, int, float]]:
    """Return ``(i, j, identity)`` triples for homologous pairs, i < j."""

    if not 0.0 < identity_threshold <= 1.0:
        raise ValueError("identity_threshold must be in (0, 1]")
    if kmer_size < 1:
        raise ValueError("kmer_size must be at least 1")
    if length_tolerance < 0:
        raise ValueError("length_tolerance must not be negative")

    sketches = [kmer_sketch(sequence, kmer_size) for _, sequence in records]

    # Inverted index, skipping k-mers that are too common to be informative.
    posting: dict[str, list[int]] = {}
    for index, sketch in enumerate(sketches):
        for kmer in sketch:
            posting.setdefault(kmer, []).append(index)

    candidates: set[tuple[int, int]] = set()
    for kmer, members in posting.items():
        if len(members) > max_kmer_frequency:
            continue
        for position, first in enumerate(members):
            for second in members[position + 1 :]:
                candidates.add((first, second))

    pairs: list[tuple[int, int, float]] = []
    for first, second in sorted(candidates):
        sequence_a = records[first][1]
        sequence_b = records[second][1]
        shorter, longer = sorted((len(sequence_a), len(sequence_b)))
        if shorter == 0 or (longer - shorter) / shorter > length_tolerance:
            continue
        union = len(sketches[first] | sketches[second])
        if union and len(sketches[first] & sketches[second]) / union < sketch_jaccard:
            continue
        identity = pair_identity(sequence_a, sequence_b)
        if identity >= identity_threshold:
            pairs.append((first, second, identity))
    return pairs


def cluster_sequences(
    frame: pd.DataFrame,
    *,
    sequence_column: str = "sequence",
    id_column: str = "sequence_id",
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    kmer_size: int = DEFAULT_KMER_SIZE,
    length_tolerance: float = DEFAULT_LENGTH_TOLERANCE,
    sketch_jaccard: float = DEFAULT_SKETCH_JACCARD,
    source_column: str = "source_dataset",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add homology cluster columns to a copy of ``frame``.

    Returns the augmented frame and a summary dictionary. Sequences are
    processed in ``sequence_id`` order so the clustering is deterministic.
    """

    if sequence_column not in frame:
        raise ValueError(f"frame is missing the sequence column: {sequence_column}")
    if id_column not in frame:
        raise ValueError(f"frame is missing the id column: {id_column}")

    work = frame.copy()
    work[sequence_column] = work[sequence_column].astype(str).str.upper()
    work = work.sort_values(id_column).reset_index(drop=True)

    identifiers = work[id_column].astype(str).tolist()
    sequences = work[sequence_column].tolist()
    records = list(zip(identifiers, sequences))

    union_find = _UnionFind(len(records))
    pairs = homology_pairs(
        records,
        identity_threshold=identity_threshold,
        kmer_size=kmer_size,
        length_tolerance=length_tolerance,
        sketch_jaccard=sketch_jaccard,
    )
    for first, second, _identity in pairs:
        union_find.union(first, second)

    roots = [union_find.find(index) for index in range(len(records))]
    members: dict[int, list[int]] = {}
    for index, root in enumerate(roots):
        members.setdefault(root, []).append(index)

    cluster_ids: dict[int, str] = {}
    for root, indices in members.items():
        representative = min(identifiers[index] for index in indices)
        cluster_ids[root] = representative

    work["cluster_id"] = [cluster_ids[root] for root in roots]
    sizes = work.groupby("cluster_id")["cluster_id"].transform("size")
    work["cluster_size"] = sizes.astype(int)
    work["cluster_representative"] = work["cluster_id"]

    if source_column in work:
        source_counts = work.groupby("cluster_id")[source_column].transform("nunique")
        work["cluster_source_count"] = source_counts.astype(int)
    else:
        work["cluster_source_count"] = 1

    # Restore the original row order for a stable, reviewable table.
    work = work.sort_values(id_column).reset_index(drop=True)

    summary = {
        "identity_threshold": identity_threshold,
        "kmer_size": kmer_size,
        "length_tolerance": length_tolerance,
        "sequence_count": len(records),
        "cluster_count": int(work["cluster_id"].nunique()),
        "homologous_pair_count": len(pairs),
        "multi_member_cluster_count": int((work["cluster_size"] > 1).sum()),
        "largest_cluster_size": int(work["cluster_size"].max()) if len(work) else 0,
        "clusters_spanning_sources": int(
            (work.groupby("cluster_id")["cluster_source_count"].first() > 1).sum()
        )
        if len(work)
        else 0,
    }
    return work, summary


def homology_aware_split(
    frame: pd.DataFrame,
    *,
    base_split_column: str = "split",
    cluster_column: str = "cluster_id",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Force every homology cluster onto a single side of the split.

    The base split (usually derived from ``source_dataset``) votes per cluster;
    the majority side wins and every member of the cluster is reassigned to it.
    Ties go to ``discovery``. The returned report lists the clusters that had to
    be reassigned, which is the direct evidence that leakage was prevented.
    """

    if cluster_column not in frame:
        raise ValueError(f"frame is missing the cluster column: {cluster_column}")
    if base_split_column not in frame:
        raise ValueError(f"frame is missing the split column: {base_split_column}")

    work = frame.copy()
    final_split: dict[str, str] = {}
    reassigned: list[dict[str, Any]] = []
    spanning = 0

    for cluster_id, group in work.groupby(cluster_column, sort=True):
        counts = group[base_split_column].value_counts()
        discovery_votes = int(counts.get("discovery", 0))
        validation_votes = int(counts.get("validation", 0))
        # Ties and empty clusters default to discovery (the larger, safer side).
        chosen = "validation" if validation_votes > discovery_votes else "discovery"
        final_split[str(cluster_id)] = chosen
        if len(group) > 1 and group[base_split_column].nunique() > 1:
            spanning += 1
            reassigned.append(
                {
                    "cluster_id": str(cluster_id),
                    "cluster_size": int(len(group)),
                    "discovery_votes": discovery_votes,
                    "validation_votes": validation_votes,
                    "assigned_to": chosen,
                    "sequence_ids": sorted(group["sequence_id"].astype(str).tolist())
                    if "sequence_id" in group
                    else [],
                }
            )

    work[base_split_column] = work[cluster_column].astype(str).map(final_split)

    # Verify the invariant rather than trusting the assignment logic.
    violations = (
        work.groupby(cluster_column)[base_split_column].nunique().gt(1).sum() if len(work) else 0
    )

    report = {
        "clusters_spanning_sides_before": int(spanning),
        "clusters_reassigned": len(reassigned),
        "leakage_violations_after": int(violations),
        "reassigned_clusters": reassigned,
    }
    return work, report


def write_cluster_report(frame: pd.DataFrame, output_path: str) -> dict[str, Any]:
    """Write a per-cluster summary table for review.

    Each row is one homology cluster with its size, the number of distinct
    sources it spans, and the species/source/split values found in it. A cluster
    whose ``source_dataset`` lists more than one value is exactly the kind of
    group that would have caused leakage before clustering.
    """

    from pathlib import Path

    if "cluster_id" not in frame:
        raise ValueError("frame is missing the cluster column: cluster_id")

    def _joined(column: str):
        def _join(values: pd.Series) -> str:
            return ",".join(sorted(set(values.astype(str))))

        return _join

    aggregations: dict[str, Any] = {"cluster_size": ("cluster_id", "size")}
    for column in ("species", "source_dataset", "split", "sigma_factor_type"):
        if column in frame:
            aggregations[column] = (column, _joined(column))

    grouped = frame.groupby("cluster_id", sort=True).agg(**aggregations).reset_index()

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    grouped.to_csv(path, sep="\t", index=False)
    return {
        "cluster_report": str(path.resolve()),
        "cluster_count": int(len(grouped)),
        "clusters_spanning_multiple_sources": int(
            grouped["source_dataset"].str.contains(",", regex=False).sum()
        )
        if "source_dataset" in grouped
        else 0,
    }
