"""Build the auditable promoter dataset consumed by the DNA embedding mainline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from experiments.eda import run_eda
from experiments.preflight import collect_environment
from preprocessing.pipeline import (
    build_quality_report,
    load_metadata,
    normalize_metadata,
    write_background_fasta,
    write_fasta,
    write_quality_report,
)
from preprocessing.provenance import sha256_file
from preprocessing.shuffle import SHUFFLE_METHODS


def _parse_tss_window(config: dict) -> tuple[int, int] | None:
    """Read the optional ``tss_window`` config entry as (upstream, downstream)."""

    raw = config.get("tss_window")
    if raw is None:
        return None
    if isinstance(raw, dict):
        upstream = raw.get("upstream")
        downstream = raw.get("downstream")
    elif isinstance(raw, (list, tuple)) and len(raw) == 2:
        upstream, downstream = raw
    else:
        raise ValueError(
            "tss_window must be {upstream: int, downstream: int} or a two-element list"
        )
    if upstream is None or downstream is None:
        raise ValueError("tss_window requires both upstream and downstream")
    upstream, downstream = int(upstream), int(downstream)
    if upstream < 0 or downstream < 0:
        raise ValueError("tss_window values must not be negative")
    return upstream, downstream


def run(config_path: str | Path) -> dict:
    config_path = Path(config_path).resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    base_dir = config_path.parent.parent if config_path.parent.name == "configs" else config_path.parent
    input_path = (base_dir / config["input"]).resolve()
    output_dir = (base_dir / config["output_dir"]).resolve()
    original = load_metadata(input_path)
    clean, validation = normalize_metadata(
        original,
        max_n_fraction=float(config.get("max_n_fraction", 0.1)),
        validation_sources=config.get("validation_sources", []),
        tss_window=_parse_tss_window(config),
        tss_column=str(config.get("tss_column", "tss_position")),
        strand_column=str(config.get("strand_column", "strand")),
        homology_clustering=bool(config.get("homology_clustering", False)),
        identity_threshold=float(config.get("identity_threshold", 0.90)),
        kmer_size=int(config.get("kmer_size", 15)),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    clean.to_csv(output_dir / "promoters_clean.tsv", sep="\t", index=False)
    write_fasta(
        ((row.sequence_id, row.sequence) for row in clean.itertuples()),
        output_dir / "promoters_clean.fasta",
    )
    seed = int(config.get("seed", 20260911))
    replicates = int(config.get("background_replicates", 1))

    background_methods = tuple(config.get("background_methods", SHUFFLE_METHODS))
    for method in background_methods:
        if method not in SHUFFLE_METHODS:
            raise ValueError(f"unknown background method: {method!r}")

    # Both negative controls are written from the same seed, so the only
    # difference between them is the shuffle algorithm itself. The N0 file keeps
    # its historical name so existing steps keep working.
    write_background_fasta(
        clean,
        output_dir / "background_shuffled.fasta",
        replicates=replicates,
        seed=seed,
        method="mononucleotide",
    )
    if "dinucleotide" in background_methods:
        write_background_fasta(
            clean,
            output_dir / "background_dinucleotide.fasta",
            replicates=replicates,
            seed=seed,
            method="dinucleotide",
        )

    quality_report = build_quality_report(
        original,
        clean,
        validation,
        max_n_fraction=float(config.get("max_n_fraction", 0.1)),
    )
    quality_report.update(
        {
            "config": str(config_path),
            "config_sha256": sha256_file(config_path),
            "input": str(input_path),
            "input_sha256": sha256_file(input_path),
            "seed": seed,
            "validation_sources": config.get("validation_sources", []),
            "background_methods": list(background_methods),
        }
    )
    write_quality_report(quality_report, output_dir / "quality_report.json")
    eda_summary = run_eda(output_dir / "promoters_clean.tsv", output_dir)
    manifest = {
        "quality_report": quality_report,
        "eda": eda_summary,
        "environment": collect_environment(),
    }
    (output_dir / "m1_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML configuration file")
    args = parser.parse_args()
    manifest = run(args.config)
    print(f"input_rows={manifest['quality_report']['input_rows']}")
    print(f"output_rows={manifest['quality_report']['output_rows']}")
    print(f"output_dir={Path(manifest['eda']['input']).parent}")


if __name__ == "__main__":
    main()
