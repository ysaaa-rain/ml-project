from collections import Counter
import random

from experiments.validate_regulondb_cross_dataset import (
    _candidate_index,
    exact_mcnemar_p,
    _similarity_category,
    dinucleotide_shuffle,
    reverse_complement,
)


def test_dinucleotide_shuffle_preserves_sequence_and_edge_counts():
    sequence = "ACGTTGCA" * 10 + "A"
    shuffled = dinucleotide_shuffle(sequence, random.Random(7))
    assert len(shuffled) == len(sequence) == 81
    assert set(shuffled) <= set("ACGT")
    assert Counter(zip(shuffled, shuffled[1:])) == Counter(zip(sequence, sequence[1:]))


def test_cross_dataset_overlap_catches_exact_reverse_and_ninety_percent_identity():
    sequence = "ACGTTGCA" * 10 + "A"
    index = _candidate_index([sequence])
    assert _similarity_category(sequence, [sequence], index) == "exact_same_orientation"
    assert _similarity_category(reverse_complement(sequence), [sequence], index) == "exact_reverse_complement"
    chars = list(sequence)
    for position in range(8):
        chars[position] = next(base for base in "ACGT" if base != chars[position])
    assert _similarity_category("".join(chars), [sequence], index) == "high_similarity_ge_90pct"


def test_exact_mcnemar_uses_paired_discordant_sequences():
    assert exact_mcnemar_p(0, 0) == 1.0
    assert exact_mcnemar_p(3, 3) == 1.0
    assert exact_mcnemar_p(10, 0) == exact_mcnemar_p(0, 10)
    assert exact_mcnemar_p(10, 0) < 0.01
