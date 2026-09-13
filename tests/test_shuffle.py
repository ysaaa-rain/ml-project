"""Tests for the N0/N1 sequence shuffling controls."""

import itertools
import random
from collections import Counter

import pandas as pd
import pytest

from preprocessing.pipeline import write_background_fasta
from preprocessing.shuffle import (
    SHUFFLE_METHODS,
    dinucleotide_counts,
    shuffle_sequence,
    shuffled_dinucleotide,
    shuffled_mononucleotide,
)


def test_mononucleotide_shuffle_preserves_length_and_composition():
    sequence = "AACCGGTTACGT"
    rng = random.Random(1)
    result = shuffled_mononucleotide(sequence, rng)
    assert len(result) == len(sequence)
    assert Counter(result) == Counter(sequence)


def test_dinucleotide_shuffle_preserves_mono_and_dinucleotide_counts():
    sequence = "AATTAAACGTACGT"
    rng = random.Random(7)
    for _ in range(50):
        result = shuffled_dinucleotide(sequence, rng)
        assert len(result) == len(sequence)
        assert Counter(result) == Counter(sequence)
        assert dinucleotide_counts(result) == dinucleotide_counts(sequence)


def test_dinucleotide_shuffle_holds_for_every_short_sequence():
    # Exhaustive over all sequences up to length 6: the shuffle must never
    # change the dinucleotide composition.
    for length in range(2, 7):
        for bases in itertools.product("ACGT", repeat=length):
            sequence = "".join(bases)
            rng = random.Random(99)
            result = shuffled_dinucleotide(sequence, rng)
            assert dinucleotide_counts(result) == dinucleotide_counts(sequence), sequence


def test_mononucleotide_shuffle_does_break_dinucleotides():
    # This is the reason N1 exists: the N0 control is usually too easy to beat.
    sequence = "AAAAACGTACGTACGTACGT"
    rng = random.Random(3)
    broke = 0
    for _ in range(20):
        if dinucleotide_counts(shuffled_mononucleotide(sequence, rng)) != dinucleotide_counts(sequence):
            broke += 1
    assert broke > 0


def test_shuffle_is_deterministic_for_a_seed():
    sequence = "AACCGGTTACGTACGT"
    assert shuffled_dinucleotide(sequence, random.Random(11)) == shuffled_dinucleotide(
        sequence, random.Random(11)
    )
    assert shuffled_mononucleotide(sequence, random.Random(11)) == shuffled_mononucleotide(
        sequence, random.Random(11)
    )


def test_dinucleotide_shuffle_returns_short_or_invalid_sequences_unchanged():
    rng = random.Random(5)
    assert shuffled_dinucleotide("", rng) == ""
    assert shuffled_dinucleotide("A", rng) == "A"
    assert shuffled_dinucleotide("ACGTN", rng) == "ACGTN"


def test_shuffle_sequence_dispatches_and_rejects_unknown_method():
    rng = random.Random(2)
    assert shuffle_sequence("ACGTACGT", rng, method="mononucleotide")
    assert shuffle_sequence("ACGTACGT", rng, method="dinucleotide")
    with pytest.raises(ValueError):
        shuffle_sequence("ACGTACGT", rng, method="dinucleotide_shuffle")


def test_background_fasta_records_the_shuffle_method(tmp_path):
    frame = pd.DataFrame(
        {"sequence_id": ["a", "b"], "sequence": ["AACCGGTTACGT", "TTGGCCAATTGG"]}
    )
    mono = tmp_path / "mono.fasta"
    di = tmp_path / "di.fasta"
    write_background_fasta(frame, mono, replicates=1, seed=4, method="mononucleotide")
    write_background_fasta(frame, di, replicates=1, seed=4, method="dinucleotide")

    mono_text = mono.read_text(encoding="utf-8")
    di_text = di.read_text(encoding="utf-8")

    # Identifiers name the control so the two backgrounds cannot be confused.
    assert ">a__mononucleotide1" in mono_text
    assert ">a__dinucleotide1" in di_text
    assert "__shuffle" not in mono_text

    # The dinucleotide control must preserve dinucleotide composition.
    di_sequence = di_text.splitlines()[1]
    assert dinucleotide_counts(di_sequence) == dinucleotide_counts("AACCGGTTACGT")


def test_background_fasta_rejects_unknown_method(tmp_path):
    frame = pd.DataFrame({"sequence_id": ["a"], "sequence": ["ACGTACGT"]})
    with pytest.raises(ValueError):
        write_background_fasta(frame, tmp_path / "x.fasta", method="nope")


def test_shuffle_methods_constant_lists_both_controls():
    assert set(SHUFFLE_METHODS) == {"mononucleotide", "dinucleotide"}
