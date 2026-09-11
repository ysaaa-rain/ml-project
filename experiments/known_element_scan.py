"""B0 scan for known bacterial promoter core-element consensus sequences."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from models.pwm import consensus_pwm, scan_many


KNOWN_MOTIFS = (
    consensus_pwm("minus_35_box", "TTGACA"),
    consensus_pwm("minus_10_box", "TATAAT"),
)


def run_scan(input_path: str | Path, output_dir: str | Path, *, min_score: float = 2.0) -> dict:
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(input_path, sep="\t")
    records = list(zip(frame["sequence_id"].astype(str), frame["sequence"].astype(str)))
    hits = scan_many(KNOWN_MOTIFS, records, min_score=min_score, both_strands=True)
    hits_frame = pd.DataFrame([hit.__dict__ for hit in hits])
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
    return {
        "input": str(input_path.resolve()),
        "sequence_count": len(records),
        "hit_count": len(hits_frame),
        "min_score": min_score,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Cleaned promoter TSV")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-score", type=float, default=2.0)
    args = parser.parse_args()
    summary = run_scan(args.input, args.output_dir, min_score=args.min_score)
    print(f"sequence_count={summary['sequence_count']}")
    print(f"hit_count={summary['hit_count']}")


if __name__ == "__main__":
    main()
