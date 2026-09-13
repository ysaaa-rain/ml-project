"""Sequence shuffling controls for promoter motif analysis.

Two negative controls are used by this project:

``mononucleotide`` (N0)
    Preserves the length and the single-base composition of each sequence.
    Because it destroys dinucleotide frequencies, real sequences usually look
    far more motif-like than their shuffled counterparts, so a significant
    result against this background is weak evidence.

``dinucleotide`` (N1)
    Preserves the length, the single-base composition *and* the dinucleotide
    composition, following the Altschul-Erikson shuffle. This is the stricter
    control required before an E-value may be treated as trustworthy.

Both shuffles are deterministic for a given seed so that a run can be repeated
exactly.
"""

from __future__ import annotations

import random

# Kept local so this module stays independent of the schema definition.
BASES = "ACGT"

SHUFFLE_METHODS = ("mononucleotide", "dinucleotide")


def shuffled_mononucleotide(sequence: str, rng: random.Random) -> str:
    """Shuffle a sequence while preserving its length and base composition."""

    bases = list(sequence)
    rng.shuffle(bases)
    return "".join(bases)


def _dinucleotide_alphabet(sequence: str) -> tuple[list[str], dict[str, int]]:
    """Return the distinct bases and their index in first-appearance order."""

    order: dict[str, int] = {}
    for base in sequence:
        if base not in order:
            order[base] = len(order)
    return list(order), order


def shuffled_dinucleotide(sequence: str, rng: random.Random) -> str:
    """Shuffle a sequence while preserving its dinucleotide composition.

    Implements the Altschul-Erikson shuffle: the sequence is treated as an
    Eulerian trail over a multigraph whose vertices are bases and whose edges
    are the observed adjacent base pairs. A random Eulerian trail over the same
    edge multiset yields a sequence with identical mononucleotide and
    dinucleotide counts.

    The input is returned unchanged when it cannot be shuffled meaningfully:
    an empty sequence, a single base, or any character outside ``ACGT``. The
    first base of the input is fixed by construction (the trail must start at
    the vertex with out-degree exceeding in-degree).
    """

    if len(sequence) < 2:
        return sequence
    if set(sequence) - set(BASES):
        return sequence

    alphabet, _ = _dinucleotide_alphabet(sequence)
    edge_count = len(sequence) - 1

    # Adjacency as lists of edge indices, so parallel edges stay distinguishable.
    outgoing: dict[str, list[int]] = {base: [] for base in alphabet}
    head_of: list[str] = []
    for start, end in zip(sequence, sequence[1:]):
        head_of.append(end)
        outgoing[start].append(len(head_of) - 1)

    # Hierholzer's algorithm with randomised edge choice at each step.
    stack: list[str] = [sequence[0]]
    trail: list[int] = []
    while stack:
        vertex = stack[-1]
        available = outgoing[vertex]
        if available:
            index = rng.randrange(len(available))
            edge = available[index]
            available[index] = available[-1]
            available.pop()
            stack.append(head_of[edge])
            trail.append(edge)
        else:
            stack.pop()

    if len(trail) != edge_count:
        # The graph is not traversable from this start (cannot happen for a
        # sequence drawn from ACGT, but never return a corrupt result).
        return sequence

    # ``trail`` is already the forward Eulerian trail from ``sequence[0]``:
    # appending each chosen edge head reconstructs a valid sequence. Reversing
    # the trail here would break the directed edge adjacency (AAC -> ACA), so
    # it is deliberately not reversed.
    bases = [sequence[0]]
    bases.extend(head_of[edge] for edge in trail)
    shuffled = "".join(bases)
    if len(shuffled) != len(sequence):
        return sequence
    if dinucleotide_counts(shuffled) != dinucleotide_counts(sequence):
        return sequence
    return shuffled


def shuffle_sequence(sequence: str, rng: random.Random, *, method: str = "mononucleotide") -> str:
    """Dispatch to the requested shuffle method."""

    if method == "mononucleotide":
        return shuffled_mononucleotide(sequence, rng)
    if method == "dinucleotide":
        return shuffled_dinucleotide(sequence, rng)
    raise ValueError(f"unknown shuffle method: {method!r}; expected one of {SHUFFLE_METHODS}")


def dinucleotide_counts(sequence: str) -> dict[str, int]:
    """Return the dinucleotide counts of a sequence (helper for verification)."""

    counts: dict[str, int] = {}
    for start, end in zip(sequence, sequence[1:]):
        pair = start + end
        counts[pair] = counts.get(pair, 0) + 1
    return counts
