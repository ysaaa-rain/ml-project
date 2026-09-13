"""Analysis utilities for the DNA embedding mainline."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import NearestNeighbors

from .dna_embedding import DNAEmbedder


def _normalized(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float32)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / np.where(norms == 0, 1.0, norms)


def pca_coordinates(embeddings: np.ndarray) -> np.ndarray:
    """Return deterministic two-dimensional PCA coordinates."""

    embeddings = np.asarray(embeddings, dtype=np.float32)
    if len(embeddings) == 0:
        return np.empty((0, 2), dtype=np.float32)
    if len(embeddings) == 1:
        return np.zeros((1, 2), dtype=np.float32)
    components = min(2, embeddings.shape[0], embeddings.shape[1])
    reduced = PCA(n_components=components, svd_solver="full").fit_transform(embeddings)
    if components == 1:
        reduced = np.column_stack([reduced[:, 0], np.zeros(len(reduced))])
    return reduced.astype(np.float32)


def cluster_embeddings(
    embeddings: np.ndarray,
    *,
    n_clusters: int,
    seed: int,
) -> tuple[np.ndarray, float | None]:
    """Cluster embeddings and report silhouette when it is well-defined."""

    embeddings = _normalized(embeddings)
    if len(embeddings) < 2:
        return np.zeros(len(embeddings), dtype=int), None
    n_clusters = max(2, min(n_clusters, len(embeddings) - 1))
    labels = KMeans(n_clusters=n_clusters, random_state=seed, n_init=20).fit_predict(embeddings)
    score = None
    if 1 < len(set(labels)) < len(labels):
        score = float(silhouette_score(embeddings, labels, metric="cosine"))
    return labels, score


def nearest_neighbor_table(
    embeddings: np.ndarray,
    sequence_ids: Sequence[str],
    *,
    neighbors: int,
) -> pd.DataFrame:
    """Return cosine nearest neighbors while excluding each sequence itself."""

    if len(embeddings) != len(sequence_ids):
        raise ValueError("embedding row count must match sequence_ids")
    if len(sequence_ids) < 2:
        return pd.DataFrame(columns=["sequence_id", "rank", "neighbor_id", "cosine_distance"])
    k = min(max(1, neighbors), len(sequence_ids) - 1)
    normalized = _normalized(embeddings)
    model = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(normalized)
    distances, indices = model.kneighbors(normalized)
    rows = []
    for row_index, sequence_id in enumerate(sequence_ids):
        rank = 0
        for distance, neighbor_index in zip(distances[row_index], indices[row_index]):
            if neighbor_index == row_index:
                continue
            rank += 1
            rows.append(
                {
                    "sequence_id": sequence_id,
                    "rank": rank,
                    "neighbor_id": sequence_ids[neighbor_index],
                    "cosine_distance": float(distance),
                }
            )
            if rank == k:
                break
    return pd.DataFrame(rows)


def label_alignment(cluster_labels: Sequence[int], labels: Sequence[str]) -> dict:
    """Diagnostic alignment only; labels do not turn the project into classification."""

    labels = [str(value) for value in labels]
    if len(cluster_labels) != len(labels) or len(set(labels)) < 2:
        return {"adjusted_rand_index": None, "normalized_mutual_info": None}
    return {
        "adjusted_rand_index": float(adjusted_rand_score(labels, cluster_labels)),
        "normalized_mutual_info": float(normalized_mutual_info_score(labels, cluster_labels)),
    }


@dataclass(frozen=True)
class OcclusionWindow:
    sequence_id: str
    start: int
    end: int
    importance: float
    original_fragment: str


def occlusion_importance(
    embedder: DNAEmbedder,
    *,
    sequence_id: str,
    sequence: str,
    window_size: int,
    stride: int,
    mask_base: str = "N",
) -> list[OcclusionWindow]:
    """Measure local importance as cosine distance after masking a window."""

    if window_size < 1 or stride < 1:
        raise ValueError("window_size and stride must be positive")
    if len(mask_base) != 1:
        raise ValueError("mask_base must be one character")
    starts = list(range(0, max(1, len(sequence) - window_size + 1), stride))
    final_start = max(0, len(sequence) - window_size)
    if final_start not in starts:
        starts.append(final_start)
    masked = []
    boundaries = []
    for start in starts:
        end = min(len(sequence), start + window_size)
        masked.append(sequence[:start] + mask_base * (end - start) + sequence[end:])
        boundaries.append((start, end))
    vectors = _normalized(embedder.embed([sequence, *masked]))
    original = vectors[0]
    windows = []
    for (start, end), vector in zip(boundaries, vectors[1:]):
        windows.append(
            OcclusionWindow(
                sequence_id=sequence_id,
                start=start,
                end=end,
                importance=float(1.0 - np.dot(original, vector)),
                original_fragment=sequence[start:end],
            )
        )
    return sorted(windows, key=lambda item: (-item.importance, item.start))


def windows_to_frame(windows: Sequence[OcclusionWindow]) -> pd.DataFrame:
    columns = ["sequence_id", "start", "end", "importance", "original_fragment"]
    return pd.DataFrame([asdict(window) for window in windows], columns=columns)
