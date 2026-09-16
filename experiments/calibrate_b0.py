"""Calibrate and pre-register the B0 known-element scan threshold.

The B0 baseline previously used ``min_score = 2.0``. Because the scoring scale
is roughly bimodal, a threshold of 2.0 is not a "weak" filter: it can admit
almost every real and shuffled sequence. The exact effect is recalculated from
the input whenever this script runs and written to the calibration record; no
fixed demo hit count is embedded in the source.

The residual background hits are perfect 6-mer matches of the consensus, which
occur by chance with a probability that depends on how many windows a sequence
contains. A single absolute threshold therefore cannot control the false
positive rate across sequence lengths: a 30 bp window has far more chance of
containing a spurious consensus match than a 100 bp one.

This script therefore calibrates **per sequence length**:

1. it scans the real sequences to find which lengths are present;
2. it builds a large dinucleotide-shuffled background for each length;
3. for each length it searches a score grid for the lowest threshold whose
   background false positive rate is at or below the target;
4. it writes the resulting threshold table plus the decision record, so the
   threshold is documented *before* biological results are inspected.

The output is a calibration record, not a biological result. Because the score
scale is discrete, the achievable false positive rate jumps in steps; the record
reports the rate that was actually achieved rather than only the target.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from models.pwm import scan_many
from preprocessing.provenance import sha256_file
from preprocessing.shuffle import SHUFFLE_METHODS, shuffle_sequence

# Pre-registered fallback used when no calibration record is available for a
# sequence length: only a perfect (zero-mismatch) consensus match is accepted.
PERFECT_MATCH_FRACTION = 1.0

DEFAULT_TARGET_FALSE_POSITIVE_RATE = 0.05
DEFAULT_MIN_BACKGROUND_PER_LENGTH = 2000
DEFAULT_RESEED_REPLICATES = 5000
CALIBRATION_FILENAME = "b0_calibration.json"


def build_length_background(
    groups: dict[int, list[str]],
    *,
    method: str,
    min_sequences: int,
    max_replicates: int,
    seed: int,
) -> dict[int, list[tuple[str, str]]]:
    """Build a background with at least ``min_sequences`` records per length.

    One shuffled copy of each real sequence of that length forms a cycle; the
    cycle is repeated (each repetition resampling the same composition) until
    the pool reaches ``min_sequences``. Short length groups are therefore not
    calibrated from a handful of samples, which matters because the false
    positive rate is estimated from how often a background sequence reaches a
    given score.

    ``max_replicates`` is only a safety valve against a pathological input; it
    is deliberately generous so that it does not silently shrink the pool.
    """

    rng = random.Random(seed)
    background: dict[int, list[tuple[str, str]]] = {}
    for length in sorted(groups):
        sequences = groups[length]
        if not sequences:
            continue
        pool: list[tuple[str, str]] = []
        replicate = 0
        while len(pool) < min_sequences and replicate < max_replicates:
            for index, sequence in enumerate(sequences):
                pool.append(
                    (
                        f"len{length}_rep{replicate}_seq{index}",
                        shuffle_sequence(sequence, rng, method=method),
                    )
                )
            replicate += 1
        background[length] = pool
    return background


def best_score_per_sequence(
    motifs,
    records: Iterable[tuple[str, str]],
    *,
    both_strands: bool = True,
) -> dict[str, float]:
    """Return the best PWM score found in each sequence.

    Every input sequence gets an entry, using negative infinity when no window
    scores above zero. Dropping those sequences would understate the background
    denominator and therefore overstate the false positive rate.
    """

    materialized = list(records)
    best: dict[str, float] = {identifier: float("-inf") for identifier, _ in materialized}
    hits = scan_many(motifs, materialized, min_score=-100.0, both_strands=both_strands)
    for hit in hits:
        if hit.score > best[hit.sequence_id]:
            best[hit.sequence_id] = hit.score
    return best


def threshold_grid(motifs, *, mid_step: float = 0.25, top_step: float = 0.005) -> list[float]:
    """Candidate thresholds covering the whole score range.

    Two details matter here and both are load-bearing:

    * the grid includes the *exact* theoretical maximum score, because the
      perfect-match score (about 10.35495 for a 6 bp consensus) is not a round
      number, and a grid built by repeated addition accumulates error and can
      land just above it;
    * thresholds are never rounded before being compared with scores. Comparing
      against a rounded-up value such as 10.355 would reject genuinely perfect
      matches while also reporting a false positive rate of zero.

    The candidates are therefore constructed explicitly rather than stepped:
    evenly spaced values across the sub-optimal range, the neighbourhood of the
    perfect-match score sampled at ``top_step``, and the perfect-match score
    itself.
    """

    perfect = max(motif.max_score for motif in motifs)
    values: set[float] = {perfect}

    # Evenly spaced candidates across the range where mismatched matches live.
    index = 0
    while index * mid_step < perfect:
        values.add(index * mid_step)
        index += 1

    # Fine sampling in the band that decides between "one mismatch" and
    # "perfect", which is where the achievable false positive rate changes.
    top_low = perfect - 0.5
    index = 0
    while top_low + index * top_step < perfect:
        values.add(top_low + index * top_step)
        index += 1
    return sorted(values)


def calibrate_per_length(
    groups: dict[int, list[str]],
    motifs,
    *,
    target_false_positive_rate: float = DEFAULT_TARGET_FALSE_POSITIVE_RATE,
    shuffle_method: str = "dinucleotide",
    min_background_per_length: int = DEFAULT_MIN_BACKGROUND_PER_LENGTH,
    max_replicates: int = DEFAULT_RESEED_REPLICATES,
    seed: int = 20260911,
) -> tuple[dict[int, dict[str, Any]], dict[int, list[float]]]:
    """Find the lowest threshold meeting the target rate for each length.

    Returns the per-length threshold table and the raw background best-score
    lists (used for the score distribution plot).
    """

    if shuffle_method not in SHUFFLE_METHODS:
        raise ValueError(f"shuffle_method must be one of {SHUFFLE_METHODS}")
    if not 0.0 < target_false_positive_rate < 1.0:
        raise ValueError("target_false_positive_rate must be in (0, 1)")

    background = build_length_background(
        groups,
        method=shuffle_method,
        min_sequences=min_background_per_length,
        max_replicates=max_replicates,
        seed=seed,
    )
    grid = threshold_grid(motifs)

    table: dict[int, dict[str, Any]] = {}
    raw_scores: dict[int, list[float]] = {}
    for length in sorted(background):
        records = background[length]
        best = best_score_per_sequence(motifs, records)
        scores = sorted(best.values(), reverse=True)
        raw_scores[length] = scores
        total = len(records)

        # Walk the grid upward and stop at the first threshold that is already
        # sustainable. Because the background best-score distribution is highly
        # discrete (zero mismatches score ~10.35, one mismatch only ~6.5), the
        # achieved rate usually drops to zero in a single step; the smallest
        # qualifying threshold is therefore the best precision/recall
        # trade-off available.
        chosen: float | None = None
        achieved = 0
        achieved_rate = 1.0
        target_met = False
        for threshold in grid:
            above = sum(1 for score in scores if score >= threshold)
            rate = above / total if total else 0.0
            if rate <= target_false_positive_rate:
                chosen = threshold
                achieved = above
                achieved_rate = rate
                target_met = True
                break
        if not target_met:
            # No grid point reaches the target on this input. Fall back to the
            # perfect-match score and record that the target was not met, rather
            # than reporting a threshold as if it were justified.
            chosen = max(motif.max_score for motif in motifs)
            achieved = sum(1 for score in scores if score >= chosen)
            achieved_rate = achieved / total if total else 0.0

        highest_score = scores[0] if scores else 0.0
        perfect_score = max(motif.max_score for motif in motifs)
        table[length] = {
            "sequence_length": length,
            "background_sequences": total,
            "target_false_positive_rate": target_false_positive_rate,
            "target_met": target_met,
            "threshold": chosen,
            "threshold_full_precision": chosen,
            "threshold_rounded": round(chosen, 4) if chosen is not None else None,
            "background_sequences_above_threshold": achieved,
            "achieved_false_positive_rate": round(achieved_rate, 5),
            "max_background_best_score": round(highest_score, 4),
            "perfect_match_score": round(perfect_score, 4),
            "requires_perfect_match": chosen is not None and chosen >= perfect_score,
            "next_higher_background_score": round(
                max((score for score in scores if chosen is not None and score < chosen), default=0.0),
                4,
            ),
        }
    return table, raw_scores


def lookup_threshold(
    calibration: dict[str, Any],
    sequence_length: int,
) -> tuple[float, str]:
    """Return ``(threshold, source)`` for a sequence length.

    Only an exact length match uses a calibrated threshold. Falling back to the
    nearest calibrated length would apply a threshold derived for one window
    count to a sequence with a different number of windows, and because the
    false positive rate depends strongly on window count that is not safe: a
    permissive threshold calibrated for short sequences would leak into longer
    ones. Uncalibrated lengths therefore use the strict perfect-match rule.
    """

    table = calibration.get("per_length_thresholds", {})

    def _threshold(entry: dict[str, Any]) -> float | None:
        # Prefer the unrounded value: a rounded threshold can sit just above the
        # perfect-match score and silently reject genuine hits.
        for key in ("threshold_full_precision", "threshold"):
            value = entry.get(key)
            if value is not None:
                return float(value)
        return None

    key = str(sequence_length)
    if key in table:
        value = _threshold(table[key])
        if value is not None:
            return value, "exact_length"

    perfect = calibration.get("perfect_match_score_full_precision")
    if perfect is None:
        perfect = calibration.get("perfect_match_score", 0.0)
    source = "uncalibrated_length_perfect_match" if table else "perfect_match_fallback"
    return float(perfect or 0.0), source


def _plot_calibration(
    real_by_length: dict[int, list[float]],
    background_scores: dict[int, list[float]],
    table: dict[int, dict[str, Any]],
    output_path: Path,
    perfect_score: float,
) -> None:
    lengths = sorted(table)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)

    axis = axes[0]
    for length in lengths:
        scores = background_scores.get(length, [])
        if not scores:
            continue
        axis.hist(
            scores,
            bins=[index * 0.5 for index in range(0, int(perfect_score / 0.5) + 3)],
            alpha=0.5,
            label=f"background {length}bp (n={len(scores)})",
        )
    axis.set_yscale("log")
    axis.set_xlabel("best PWM score per sequence")
    axis.set_ylabel("sequence count (log scale)")
    axis.set_title("Background best-score distribution by length")
    axis.legend(fontsize=7)

    axis = axes[1]
    for length in lengths:
        scores = real_by_length.get(length, [])
        if scores:
            axis.scatter(
                [length] * len(scores),
                scores,
                s=18,
                alpha=0.8,
                label="real" if length == lengths[0] else None,
                color="#4472C4",
            )
        entry = table[length]
        if entry["threshold"] is not None:
            axis.scatter(
                [length],
                [entry["threshold"]],
                marker="_",
                s=400,
                color="#d93025",
                label="calibrated threshold" if length == lengths[0] else None,
            )
    axis.axhline(perfect_score, color="#188038", linestyle=":", linewidth=1.2, label="perfect match")
    axis.set_xlabel("sequence length (bp)")
    axis.set_ylabel("best PWM score per sequence")
    axis.set_title("Calibrated threshold vs real sequence scores")
    axis.legend(fontsize=7)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=160)
    plt.close(figure)


def calibrate(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    shuffle_method: str = "dinucleotide",
    seed: int = 20260911,
    target_false_positive_rate: float = DEFAULT_TARGET_FALSE_POSITIVE_RATE,
    min_background_per_length: int = DEFAULT_MIN_BACKGROUND_PER_LENGTH,
    max_replicates: int = DEFAULT_RESEED_REPLICATES,
) -> dict:
    """Run the full calibration and write the decision record."""

    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(input_path, sep="\t")
    for column in ("sequence_id", "sequence"):
        if column not in frame:
            raise ValueError(f"clean table is missing required column: {column}")

    real_records = list(zip(frame["sequence_id"].astype(str), frame["sequence"].astype(str)))
    groups: dict[int, list[str]] = {}
    for _, sequence in real_records:
        groups.setdefault(len(sequence), []).append(sequence)

    from experiments.known_element_scan import known_motifs

    motifs = known_motifs()
    perfect_score = round(max(motif.max_score for motif in motifs), 4)

    table, background_scores = calibrate_per_length(
        groups,
        motifs,
        target_false_positive_rate=target_false_positive_rate,
        shuffle_method=shuffle_method,
        min_background_per_length=min_background_per_length,
        max_replicates=max_replicates,
        seed=seed,
    )

    real_best = best_score_per_sequence(motifs, real_records)
    real_by_length: dict[int, list[float]] = {}
    for (identifier, sequence), score in zip(real_records, [real_best.get(i) for i, _ in real_records]):
        if score is not None:
            real_by_length.setdefault(len(sequence), []).append(score)

    _plot_calibration(
        real_by_length,
        background_scores,
        table,
        output_dir / "b0_score_distribution.png",
        perfect_score,
    )

    threshold_rows = pd.DataFrame(
        [
            {
                "sequence_length": entry["sequence_length"],
                "threshold": entry["threshold"],
                "background_sequences": entry["background_sequences"],
                "background_sequences_above_threshold": entry["background_sequences_above_threshold"],
                "achieved_false_positive_rate": entry["achieved_false_positive_rate"],
                "requires_perfect_match": entry["requires_perfect_match"],
                "max_background_best_score": entry["max_background_best_score"],
            }
            for entry in table.values()
        ]
    )
    threshold_rows.to_csv(output_dir / "b0_threshold_table.tsv", sep="\t", index=False)

    # Compare against the previously hard-coded default of 2.0, as the concrete
    # evidence for why the change was needed.
    legacy_hits = scan_many(motifs, real_records, min_score=2.0, both_strands=True)
    legacy_sequences = len({hit.sequence_id for hit in legacy_hits})

    targets_met = [entry["target_met"] for entry in table.values()]
    if all(targets_met):
        verdict = (
            "every sequence length has a threshold meeting the target false positive "
            "rate; thresholds are ready for use"
        )
    elif any(targets_met):
        unmet = sorted(length for length, entry in table.items() if not entry["target_met"])
        verdict = (
            "some sequence lengths could not reach the target false positive rate "
            f"(lengths: {unmet}); those lengths fall back to the perfect-match score, "
            "so sensitivity there is limited to exact consensus matches"
        )
    else:
        verdict = (
            "no sequence length can reach the target false positive rate, not even with "
            "the perfect-match score; the input is unusable for known-element recovery "
            "calibration and the resulting thresholds must not be presented as "
            "calibrated. This is expected for the small synthetic demo fixture, whose "
            "sequences are built from exact consensus elements"
        )

    record = {
        "record_type": "b0_threshold_calibration",
        "record_version": "2.0",
        "note": (
            "Thresholds are locked before biological results are inspected. "
            "Changing them requires a new calibration record."
        ),
        "verdict": verdict,
        "target_met_for_all_lengths": all(targets_met) if targets_met else False,
        "input": str(input_path.resolve()),
        "input_sha256": sha256_file(input_path),
        "sequence_count": len(real_records),
        "sequence_lengths": sorted(groups),
        "seed": seed,
        "shuffle_method": shuffle_method,
        "per_length_background_size": min_background_per_length,
        "target_false_positive_rate": target_false_positive_rate,
        "perfect_match_score": perfect_score,
        "perfect_match_score_full_precision": max(motif.max_score for motif in motifs),
        "perfect_match_fraction": PERFECT_MATCH_FRACTION,
        "motifs": [
            {
                "name": motif.name,
                "width": motif.width,
                "max_score": round(motif.max_score, 4),
                "reverse_complement_consensus_score": round(motif.background_max_score(), 4),
            }
            for motif in motifs
        ],
        "per_length_thresholds": {
            str(length): entry for length, entry in sorted(table.items())
        },
        "previous_default_threshold": 2.0,
        "previous_default_real_hits": len(legacy_hits),
        "previous_default_real_sequences_with_hit": legacy_sequences,
        "previous_default_sequence_fraction": round(legacy_sequences / len(real_records), 4)
        if real_records
        else None,
        "artifacts": {
            "threshold_table": str((output_dir / "b0_threshold_table.tsv").resolve()),
            "score_distribution_plot": str((output_dir / "b0_score_distribution.png").resolve()),
            "record": str((output_dir / CALIBRATION_FILENAME).resolve()),
        },
    }

    (output_dir / CALIBRATION_FILENAME).write_text(
        json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return record


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Cleaned promoter TSV")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--shuffle-method",
        choices=SHUFFLE_METHODS,
        default="dinucleotide",
        help="Background control to calibrate against (default: the stricter N1 shuffle)",
    )
    parser.add_argument("--seed", type=int, default=20260911)
    parser.add_argument(
        "--target-false-positive-rate",
        type=float,
        default=DEFAULT_TARGET_FALSE_POSITIVE_RATE,
        help="Maximum acceptable fraction of background sequences carrying a hit",
    )
    parser.add_argument(
        "--background-per-length",
        type=int,
        default=DEFAULT_MIN_BACKGROUND_PER_LENGTH,
        help="Minimum number of background sequences generated per sequence length",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    record = calibrate(
        args.input,
        args.output_dir,
        shuffle_method=args.shuffle_method,
        seed=args.seed,
        target_false_positive_rate=args.target_false_positive_rate,
        min_background_per_length=args.background_per_length,
    )
    print(f"sequence_count={record['sequence_count']}")
    print(f"sequence_lengths={record['sequence_lengths']}")
    print(f"perfect_match_score={record['perfect_match_score']}")
    for length, entry in record["per_length_thresholds"].items():
        print(
            f"  length={length}bp threshold={entry['threshold']} "
            f"achieved_fpr={entry['achieved_false_positive_rate']} "
            f"(background n={entry['background_sequences']})"
        )
    print(
        "previous_default_threshold=2.0 "
        f"matched {record['previous_default_sequence_fraction']} of real sequences"
    )
    print(f"output_dir={Path(args.output_dir).resolve()}")


if __name__ == "__main__":
    main()
