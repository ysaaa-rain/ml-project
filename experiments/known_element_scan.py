"""B0 scan for known bacterial promoter core-element consensus sequences.

Thresholds are calibrated per sequence length by ``experiments.calibrate_b0``.
The reason is that a single absolute score cannot control the false positive
rate across lengths: the consensus elements are only 6 bp, so a 30 bp window has
a much higher chance of containing a spurious perfect match than a 100 bp one.

When no calibration record is available the scan accepts only a perfect
(zero-mismatch) consensus match, which is the conservative default recorded in
``experiments.calibrate_b0.PERFECT_MATCH_FRACTION``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from models.pwm import PWM, consensus_pwm, scan_many


# Consensus sequences of the two core promoter elements.
MINUS_35_CONSENSUS = "TTGACA"
MINUS_10_CONSENSUS = "TATAAT"

# Fallback threshold as a fraction of the theoretical maximum PWM score, used
# only when no calibration record is supplied. 1.0 means "perfect match only".
# Calibrated runs override this per sequence length; see experiments/calibrate_b0.py.
CALIBRATED_MIN_SCORE_FRACTION = 1.0

DEFAULT_CALIBRATION_FILENAME = "b0_calibration.json"


def known_motifs() -> tuple[PWM, PWM]:
    """Return the B0 known-element PWMs (-35 box, -10 box)."""

    return (
        consensus_pwm("minus_35_box", MINUS_35_CONSENSUS),
        consensus_pwm("minus_10_box", MINUS_10_CONSENSUS),
    )


def perfect_match_threshold() -> float:
    """Return the score of a zero-mismatch consensus match."""

    return min(motif.max_score for motif in known_motifs())


def calibrated_min_score(fraction: float = CALIBRATED_MIN_SCORE_FRACTION) -> float:
    """Return the fallback threshold for the given fraction of the maximum score."""

    return perfect_match_threshold() * fraction


# Retained under the original name so existing callers keep working.
KNOWN_MOTIFS = known_motifs()


def load_calibration(path: str | Path) -> dict:
    """Load a calibration record written by ``experiments.calibrate_b0``."""

    calibration_path = Path(path)
    if not calibration_path.is_file():
        raise FileNotFoundError(f"calibration record not found: {calibration_path}")
    record = json.loads(calibration_path.read_text(encoding="utf-8"))
    if not isinstance(record, dict) or "per_length_thresholds" not in record:
        raise ValueError(f"not a B0 calibration record: {calibration_path}")
    return record


def resolve_threshold(calibration: dict | None, sequence_length: int) -> tuple[float, str]:
    """Return ``(threshold, source)`` for one sequence length."""

    if calibration is None:
        return perfect_match_threshold(), "perfect_match_default"
    from experiments.calibrate_b0 import lookup_threshold

    return lookup_threshold(calibration, sequence_length)


def run_scan(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    min_score: float | None = None,
    calibration_path: str | Path | None = None,
) -> dict:
    """Run the B0 scan and write the hit and summary tables.

    Precedence: an explicit ``min_score`` wins, then ``calibration_path``, then
    the perfect-match default.
    """

    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(input_path, sep="\t")
    records = list(zip(frame["sequence_id"].astype(str), frame["sequence"].astype(str)))
    length_by_id = {identifier: len(sequence) for identifier, sequence in records}

    calibration: dict | None = None
    if min_score is None and calibration_path is not None:
        calibration = load_calibration(calibration_path)

    # Scan down to the theoretical minimum so that a per-length threshold can be
    # applied to the same hit set; the minimum score of a 6 bp consensus PWM is
    # well above -100, so this floor keeps every candidate window without
    # discarding anything a length-specific threshold could have accepted.
    hits = scan_many(KNOWN_MOTIFS, records, min_score=-100.0, both_strands=True)

    threshold_usage: dict[str, dict] = {}
    kept = []
    for hit in hits:
        if min_score is not None:
            threshold, source = min_score, "explicit_override"
        else:
            threshold, source = resolve_threshold(calibration, length_by_id.get(hit.sequence_id, 0))
        usage = threshold_usage.setdefault(
            source, {"threshold": threshold, "hit_count": 0, "sequences": set()}
        )
        usage["hit_count"] += 1
        usage["sequences"].add(hit.sequence_id)
        if hit.score >= threshold:
            kept.append(hit)

    hits_frame = pd.DataFrame([hit.__dict__ for hit in kept])
    if hits_frame.empty:
        hits_frame = pd.DataFrame(
            columns=["motif_name", "sequence_id", "start", "end", "strand", "score", "p_value_proxy", "matched_sequence"]
        )
    hits_frame.to_csv(output_dir / "known_element_hits.tsv", sep="\t", index=False)
    summary = (
        hits_frame.groupby("motif_name", as_index=False)
        .agg(sequence_count=("sequence_id", "nunique"), hit_count=("sequence_id", "count"))
        if not hits_frame.empty
        else pd.DataFrame(columns=["motif_name", "sequence_count", "hit_count"])
    )
    summary.to_csv(output_dir / "known_element_summary.tsv", sep="\t", index=False)

    for usage in threshold_usage.values():
        usage["sequence_count"] = len(usage.pop("sequences"))
    return {
        "input": str(input_path.resolve()),
        "sequence_count": len(records),
        "hit_count": len(hits_frame),
        "thresholds": threshold_usage,
        "calibration": str(Path(calibration_path).resolve()) if calibration_path else None,
        "min_score": min_score,
        "min_score_fraction": CALIBRATED_MIN_SCORE_FRACTION if min_score is None else None,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Cleaned promoter TSV")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--calibration",
        default=None,
        help=(
            "Calibration record from experiments.calibrate_b0. Strongly "
            "recommended: without it only perfect consensus matches are reported."
        ),
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=None,
        help=(
            "Override all thresholds with one absolute score. Overriding "
            "invalidates comparison with previously reported B0 results."
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary = run_scan(
        args.input,
        args.output_dir,
        min_score=args.min_score,
        calibration_path=args.calibration,
    )
    print(f"sequence_count={summary['sequence_count']}")
    print(f"hit_count={summary['hit_count']}")
    for source, usage in sorted(summary["thresholds"].items()):
        print(
            f"  threshold_source={source} threshold={usage['threshold']:.4f} "
            f"hits={usage['hit_count']}"
        )
    if summary["calibration"]:
        print(f"calibration={summary['calibration']}")
    else:
        print("calibration=none (perfect-match default; run experiments.calibrate_b0 first)")


if __name__ == "__main__":
    main()
