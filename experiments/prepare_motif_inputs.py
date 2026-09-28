"""Prepare audited, source-internal FASTA inputs for motif discovery."""

from __future__ import annotations

import argparse
import json
import random
import re
from pathlib import Path

import pandas as pd

from preprocessing.pipeline import write_fasta
from preprocessing.provenance import sha256_file
from preprocessing.shuffle import dinucleotide_counts, shuffle_sequence


def run(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    source: str = "regulondb",
    upstream: int = 60,
    downstream: int = 20,
    min_group_size: int = 50,
    seed: int = 20260922,
) -> dict:
    input_path = Path(input_path).resolve()
    output_dir = Path(output_dir).resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"output directory is not empty: {output_dir}")
    if min_group_size < 1 or upstream < 0 or downstream < 0:
        raise ValueError("group size must be positive and window sizes nonnegative")
    frame = pd.read_csv(input_path, sep="\t", dtype={"sequence_id": str}).sort_values("sequence_id")
    required = {
        "sequence_id", "sequence", "source_dataset", "split", "sigma_factor_type",
        "tss_position", "strand", "window_upstream", "window_downstream",
        "tss_offset_in_window", "window_strand",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"input is missing required columns: {', '.join(missing)}")
    if frame.empty or frame.sequence_id.isna().any() or frame.sequence_id.duplicated().any():
        raise ValueError("sequence_id must be nonempty and unique")
    if not frame.source_dataset.eq(source).all() or not frame.split.eq("discovery").all():
        raise ValueError("motif discovery input must contain only the specified source and discovery split")
    if frame.tss_position.isna().any() or frame.strand.isna().any():
        raise ValueError("motif discovery input requires audited TSS and strand fields")
    if not frame.window_upstream.eq(upstream).all() or not frame.window_downstream.eq(downstream).all():
        raise ValueError("input window does not match the locked TSS coordinates")
    if not frame.tss_offset_in_window.eq(upstream).all():
        raise ValueError("TSS offset does not match the locked window")
    expected_strand = frame.strand.map({"+": 1, "-": -1})
    if expected_strand.isna().any() or not frame.window_strand.eq(expected_strand).all():
        raise ValueError("window strand does not match the audited source strand")
    expected_length = upstream + downstream + 1
    sequences = frame.sequence.astype(str).str.upper()
    if not sequences.str.fullmatch(f"[ACGT]{{{expected_length}}}").all():
        raise ValueError("all sequences must have the locked length and contain only A/C/G/T")

    labels = frame.sigma_factor_type.fillna("unknown").astype(str)
    known = sorted(label for label in labels.unique() if label != "unknown")
    if any(not re.fullmatch(r"[A-Za-z0-9_-]+", label) for label in known):
        raise ValueError("sigma labels must be safe ASCII file names")
    included = [label for label in known if int((labels == label).sum()) >= min_group_size]

    ids = frame.sequence_id.astype(str).tolist()
    real = sequences.tolist()
    rng = random.Random(seed)
    control = [shuffle_sequence(sequence, rng, method="dinucleotide") for sequence in real]
    if any(dinucleotide_counts(a) != dinucleotide_counts(b) for a, b in zip(real, control)):
        raise ValueError("N1 background did not preserve dinucleotide counts")

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, dict] = {}

    def save(name: str, indices: list[int]) -> None:
        primary = output_dir / f"{name}.fasta"
        background = output_dir / f"{name}_n1.fasta"
        write_fasta(((ids[i], real[i]) for i in indices), primary)
        write_fasta(((f"{ids[i]}__dinucleotide1", control[i]) for i in indices), background)
        outputs[name] = {
            "sequence_count": len(indices),
            "primary_fasta": primary.name,
            "primary_sha256": sha256_file(primary),
            "n1_fasta": background.name,
            "n1_sha256": sha256_file(background),
            "n1_changed_count": sum(real[i] != control[i] for i in indices),
        }

    save("all", list(range(len(ids))))
    for label in included:
        save(label, [i for i, value in enumerate(labels) if value == label])
    manifest = {
        "schema_version": "motif-inputs-1.0",
        "status": "inputs_prepared_no_motif_discovery_or_fimo_run",
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "source_dataset": source,
        "split": "discovery",
        "coordinate_window": [-upstream, downstream],
        "sequence_length": expected_length,
        "control": "one dinucleotide-preserving shuffle per primary sequence",
        "seed": seed,
        "minimum_group_size": min_group_size,
        "excluded_groups": {label: int((labels == label).sum()) for label in known if label not in included},
        "unknown_sigma_count": int((labels == "unknown").sum()),
        "outputs": outputs,
    }
    (output_dir / "motif_input_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--source", default="regulondb")
    parser.add_argument("--upstream", type=int, default=60)
    parser.add_argument("--downstream", type=int, default=20)
    parser.add_argument("--min-group-size", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260922)
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir, source=args.source, upstream=args.upstream,
                         downstream=args.downstream, min_group_size=args.min_group_size, seed=args.seed),
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
