"""Run the complete M1 data build and M2 EDA smoke workflow."""

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
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    clean.to_csv(output_dir / "promoters_clean.tsv", sep="\t", index=False)
    write_fasta(
        ((row.sequence_id, row.sequence) for row in clean.itertuples()),
        output_dir / "promoters_clean.fasta",
    )
    seed = int(config.get("seed", 20260911))
    write_background_fasta(
        clean,
        output_dir / "background_shuffled.fasta",
        replicates=int(config.get("background_replicates", 1)),
        seed=seed,
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
