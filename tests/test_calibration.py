"""Tests for the B0 threshold calibration and its use by the scan."""

import json
import random
from pathlib import Path

import pandas as pd
import pytest

from experiments.calibrate_b0 import (
    best_score_per_sequence,
    build_length_background,
    calibrate_per_length,
    lookup_threshold,
    threshold_grid,
)
from experiments.known_element_scan import (
    known_motifs,
    load_calibration,
    perfect_match_threshold,
    run_scan,
)


def _random_records(count: int, length: int, seed: int):
    rng = random.Random(seed)
    return [
        (f"s{index}", "".join(rng.choice("ACGT") for _ in range(length)))
        for index in range(count)
    ]


def test_threshold_grid_contains_the_exact_perfect_match_score():
    motifs = known_motifs()
    perfect = max(motif.max_score for motif in motifs)
    grid = threshold_grid(motifs)

    # The perfect-match score is not a round number; the grid must contain it
    # exactly, otherwise a genuinely perfect hit fails the comparison.
    assert perfect in grid
    assert any(abs(value - perfect) < 1e-12 for value in grid)
    # And it must not contain a value just above it that would reject it.
    just_below = [value for value in grid if value < perfect]
    assert max(just_below) < perfect
    assert min(value for value in grid if value >= perfect) == perfect


def test_threshold_grid_is_sorted_and_covers_the_range():
    grid = threshold_grid(known_motifs())
    assert grid == sorted(grid)
    assert grid[0] == 0.0
    assert grid[-1] == max(motif.max_score for motif in known_motifs())


def test_perfect_match_threshold_matches_a_consensus_hit():
    motifs = known_motifs()
    pwm = motifs[0]
    from models.pwm import scan_pwm

    hits = scan_pwm(pwm, "x", "TTGACA", min_score=-100.0, both_strands=False)
    assert max(hit.score for hit in hits) == pytest.approx(perfect_match_threshold())


def test_best_score_per_sequence_picks_the_maximum():
    records = [("a", "TTGACAGGGG"), ("b", "GGGGGGGGGG")]
    best = best_score_per_sequence(known_motifs(), records)
    # Every sequence gets an entry, including one with no positive-scoring window.
    assert set(best) == {"a", "b"}
    assert best["a"] >= perfect_match_threshold() - 1e-9
    assert best["b"] < best["a"]


def test_build_length_background_reaches_the_requested_size():
    groups = {30: ["ACGT" * 7 + "AC", "TTTT" * 7 + "TT"]}
    background = build_length_background(
        groups, method="dinucleotide", min_sequences=500, max_replicates=5000, seed=1
    )
    assert len(background[30]) >= 500
    # Composition is preserved for every generated record.
    from preprocessing.shuffle import dinucleotide_counts

    for _, sequence in background[30][:20]:
        assert len(sequence) == 30
        assert dinucleotide_counts(sequence) in (
            dinucleotide_counts("ACGT" * 7 + "AC"),
            dinucleotide_counts("TTTT" * 7 + "TT"),
        )


def test_calibrate_per_length_sets_a_threshold_at_or_below_the_perfect_score():
    # Short sequences contain very few windows, so a spurious perfect consensus
    # match is rare and a low false positive rate is reachable.
    groups = {40: [sequence for _, sequence in _random_records(30, 40, 5)]}
    table, scores = calibrate_per_length(
        groups,
        known_motifs(),
        target_false_positive_rate=0.05,
        min_background_per_length=400,
        max_replicates=5000,
        seed=2,
    )
    entry = table[40]
    perfect = max(motif.max_score for motif in known_motifs())
    assert entry["threshold"] is not None
    assert entry["threshold"] <= perfect
    assert entry["target_met"] is True
    assert entry["achieved_false_positive_rate"] <= 0.05
    assert entry["background_sequences"] >= 400
    assert 40 in scores


def test_calibrate_per_length_reports_an_unreachable_target_honestly():
    # A 200 bp background contains ~195 windows per sequence, so the chance of an
    # incidental perfect consensus match is itself far above 5%. No threshold on
    # the grid can reach the target, and the record must say so instead of
    # presenting the fallback as if it had been calibrated.
    groups = {200: [sequence for _, sequence in _random_records(30, 200, 6)]}
    table, _ = calibrate_per_length(
        groups,
        known_motifs(),
        target_false_positive_rate=0.05,
        min_background_per_length=400,
        max_replicates=5000,
        seed=7,
    )
    entry = table[200]
    perfect = max(motif.max_score for motif in known_motifs())
    assert entry["target_met"] is False
    assert entry["threshold"] == pytest.approx(perfect)
    assert entry["requires_perfect_match"] is True
    assert entry["achieved_false_positive_rate"] > 0.05


def test_calibrate_per_length_reports_when_the_target_cannot_be_met():
    # A sequence consisting of the consensus repeated gives a background that
    # almost always contains a perfect match, so no threshold can reach a low
    # false positive rate. The tool must say so instead of pretending otherwise.
    consensus = "TTGACA"
    groups = {len(consensus) * 4: [consensus * 4]}
    table, _ = calibrate_per_length(
        groups,
        known_motifs(),
        target_false_positive_rate=0.01,
        min_background_per_length=200,
        max_replicates=200,
        seed=3,
    )
    entry = table[len(consensus) * 4]
    assert entry["target_met"] is False
    assert entry["requires_perfect_match"] is True
    assert entry["threshold"] == pytest.approx(max(motif.max_score for motif in known_motifs()))


def test_lookup_threshold_prefers_the_unrounded_value():
    perfect = max(motif.max_score for motif in known_motifs())
    calibration = {
        "per_length_thresholds": {
            "60": {
                "threshold": round(perfect, 4),          # 10.355, above the real max
                "threshold_full_precision": perfect,     # 10.3549502...
            }
        },
        "perfect_match_score": round(perfect, 4),
        "perfect_match_score_full_precision": perfect,
    }
    threshold, source = lookup_threshold(calibration, 60)
    assert source == "exact_length"
    # Using the rounded value would reject a perfect match.
    assert threshold == pytest.approx(perfect)
    assert threshold < round(perfect, 4)


def test_lookup_threshold_does_not_extend_a_calibration_to_other_lengths():
    perfect = max(motif.max_score for motif in known_motifs())
    calibration = {
        # A permissive threshold, but only justified for 30 bp sequences.
        "per_length_thresholds": {"30": {"threshold": 6.0, "threshold_full_precision": 6.0}},
        "perfect_match_score_full_precision": perfect,
    }
    threshold, source = lookup_threshold(calibration, 62)
    # Applying the 30 bp threshold here would be unsafe: the false positive rate
    # depends on the window count, so an uncalibrated length must stay strict.
    assert source == "uncalibrated_length_perfect_match"
    assert threshold == pytest.approx(perfect)


def test_lookup_threshold_falls_back_to_perfect_match_without_a_table():
    perfect = max(motif.max_score for motif in known_motifs())
    threshold, source = lookup_threshold({"perfect_match_score_full_precision": perfect}, 80)
    assert source == "perfect_match_fallback"
    assert threshold == pytest.approx(perfect)


def test_load_calibration_rejects_a_foreign_file(tmp_path):
    path = tmp_path / "not_a_calibration.json"
    path.write_text(json.dumps({"hello": "world"}), encoding="utf-8")
    with pytest.raises(ValueError):
        load_calibration(path)
    with pytest.raises(FileNotFoundError):
        load_calibration(tmp_path / "missing.json")


def test_scan_without_calibration_reports_only_perfect_matches(tmp_path):
    frame = pd.DataFrame(
        {
            "sequence_id": ["perfect", "one_mismatch"],
            "sequence": ["TTGACA" + "G" * 24, "TTGACT" + "G" * 24],
        }
    )
    source = tmp_path / "promoters.tsv"
    frame.to_csv(source, sep="\t", index=False)

    summary = run_scan(source, tmp_path / "out")
    # Only the zero-mismatch record may pass the default threshold.
    assert summary["hit_count"] == 1
    hits = pd.read_csv(tmp_path / "out" / "known_element_hits.tsv", sep="\t")
    assert list(hits["sequence_id"]) == ["perfect"]
    assert summary["thresholds"]["perfect_match_default"]["threshold"] == pytest.approx(
        perfect_match_threshold()
    )


def test_scan_rejects_near_consensus_hits_that_the_old_default_accepted(tmp_path):
    """Regression guard for the old min_score=2.0 default.

    A single mismatch scores about 6.5, far above 2.0, so the previous default
    reported it. That is exactly the false-positive behaviour being fixed.
    """

    frame = pd.DataFrame(
        {"sequence_id": ["one_mismatch"], "sequence": ["TTGACT" + "G" * 24]}
    )
    source = tmp_path / "promoters.tsv"
    frame.to_csv(source, sep="\t", index=False)

    assert run_scan(source, tmp_path / "strict")["hit_count"] == 0
    override = run_scan(source, tmp_path / "legacy", min_score=2.0)
    assert override["hit_count"] == 1


def test_scan_applies_the_per_length_calibrated_threshold(tmp_path):
    perfect = max(motif.max_score for motif in known_motifs())
    calibration = {
        "record_type": "b0_threshold_calibration",
        "perfect_match_score_full_precision": perfect,
        "per_length_thresholds": {
            # A deliberately permissive threshold for the 30 bp record only.
            "30": {"threshold": 6.0, "threshold_full_precision": 6.0},
        },
    }
    calibration_path = tmp_path / "b0_calibration.json"
    calibration_path.write_text(json.dumps(calibration), encoding="utf-8")

    frame = pd.DataFrame(
        {
            "sequence_id": ["len30_mismatch", "len60_mismatch"],
            "sequence": ["TTGACT" + "G" * 24, "TTGACT" + "G" * 54],
        }
    )
    source = tmp_path / "promoters.tsv"
    frame.to_csv(source, sep="\t", index=False)

    summary = run_scan(source, tmp_path / "out", calibration_path=calibration_path)
    # The 30 bp record is judged with threshold 6.0 and passes; the 60 bp record
    # has no entry so it falls back to the perfect-match score and fails.
    hits = pd.read_csv(tmp_path / "out" / "known_element_hits.tsv", sep="\t")
    assert list(hits["sequence_id"]) == ["len30_mismatch"]
    assert "exact_length" in summary["thresholds"]
    assert summary["calibration"] == str(calibration_path.resolve())
