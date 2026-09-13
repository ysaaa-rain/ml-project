"""Command-line entry point for the M1 promoter data build."""

from __future__ import annotations

import argparse
from pathlib import Path

from .pipeline import (
    build_quality_report,
    load_metadata,
    normalize_metadata,
    write_background_fasta,
    write_fasta,
    write_quality_report,
)
from .shuffle import SHUFFLE_METHODS


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input CSV/TSV metadata file")
    parser.add_argument("--output-dir", required=True, help="Directory for cleaned data and reports")
    parser.add_argument("--max-n-fraction", type=float, default=0.1)
    parser.add_argument("--validation-source", action="append", default=[])
    parser.add_argument("--background-replicates", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260911)
    parser.add_argument(
        "--background-method",
        choices=SHUFFLE_METHODS,
        default="mononucleotide",
        help="N0 mononucleotide or N1 dinucleotide-preserving shuffle",
    )
    parser.add_argument(
        "--tss-window",
        nargs=2,
        type=int,
        metavar=("UPSTREAM", "DOWNSTREAM"),
        default=None,
        help=(
            "Extract a TSS-anchored window and normalise strand orientation. "
            "Requires tss_position and strand columns in the input."
        ),
    )
    parser.add_argument(
        "--homology-clustering",
        action="store_true",
        help="Cluster near-homologs and force each cluster onto one side of the split",
    )
    parser.add_argument("--identity-threshold", type=float, default=0.90)
    parser.add_argument("--kmer-size", type=int, default=15)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir)
    original = load_metadata(args.input)
    clean, validation = normalize_metadata(
        original,
        max_n_fraction=args.max_n_fraction,
        validation_sources=args.validation_source,
        tss_window=tuple(args.tss_window) if args.tss_window else None,
        homology_clustering=args.homology_clustering,
        identity_threshold=args.identity_threshold,
        kmer_size=args.kmer_size,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    clean.to_csv(output_dir / "promoters_clean.tsv", sep="\t", index=False)
    write_fasta(
        ((row.sequence_id, row.sequence) for row in clean.itertuples()),
        output_dir / "promoters_clean.fasta",
    )
    write_background_fasta(
        clean,
        output_dir / "background_shuffled.fasta",
        replicates=args.background_replicates,
        seed=args.seed,
        method=args.background_method,
    )
    report = build_quality_report(
        original,
        clean,
        validation,
        max_n_fraction=args.max_n_fraction,
    )
    report["seed"] = args.seed
    report["validation_sources"] = args.validation_source
    report["background_method"] = args.background_method
    write_quality_report(report, output_dir / "quality_report.json")
    print(f"input_rows={len(original)}")
    print(f"output_rows={len(clean)}")
    print(f"output_dir={output_dir.resolve()}")


if __name__ == "__main__":
    main()
