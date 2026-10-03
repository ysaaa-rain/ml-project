"""Validate the six TJU Pan promoter tables and prepare motif FASTA inputs."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from preprocessing.provenance import sha256_file


SPECIES = {
    "bacillus_subtilis": "Bacillus_subtilis",
    "baumannii": "Baumanii",
    "bradyrhizobium": "Bradyrhizobium",
    "diphtheria": "Diphtheria",
    "escherichia_coli": "Escherichia_coli",
    "staphylococcus": "Staphylococcus",
}
REQUIRED_COLUMNS = ("seq_id", "seq_type", "seq", "label")
DNA_RE = re.compile(r"[ACGT]+")
SEQUENCE_RE = re.compile(r"[ACGT]{81}")


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not set(REQUIRED_COLUMNS).issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns {REQUIRED_COLUMNS}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} is empty")
    return rows


def _validate_rows(
    rows: list[dict[str, str]], path: Path, *, require_81bp: bool = True
) -> list[dict[str, str]]:
    normalized: list[dict[str, str]] = []
    seen_ids: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        seq_id = row["seq_id"].strip()
        seq = row["seq"].strip().upper()
        label = row["label"].strip()
        if not seq_id or seq_id in seen_ids:
            raise ValueError(f"{path}:{line_number}: seq_id must be nonempty and unique: {seq_id!r}")
        if label not in {"0", "1"}:
            raise ValueError(f"{path}:{line_number}: label must be 0 or 1, got {label!r}")
        if not DNA_RE.fullmatch(seq):
            raise ValueError(f"{path}:{line_number}: sequence must contain only A/C/G/T")
        if require_81bp and not SEQUENCE_RE.fullmatch(seq):
            raise ValueError(f"{path}:{line_number}: expected an 81 bp A/C/G/T sequence")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", seq_id):
            raise ValueError(f"{path}:{line_number}: seq_id is unsafe for FASTA: {seq_id!r}")
        seen_ids.add(seq_id)
        normalized.append(
            {
                "seq_id": seq_id,
                "seq_type": row["seq_type"].strip(),
                "seq": seq,
                "label": label,
            }
        )
    return normalized


def _write_fasta(rows: list[dict[str, str]], path: Path, *, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for row in rows:
            handle.write(f">{row['seq_id']}|label={label}\n{row['seq']}\n")


def prepare(dataset_root: Path, output_dir: Path, report_path: Path) -> dict[str, Any]:
    report: dict[str, Any] = {
        "schema_version": "tjupan-motif-inputs-1.0",
        "status": "inputs_validated_and_written",
        "source": "TJU course cloud drive; reg_and_gen/Datasets",
        "upstream_reference": "https://github.com/KevinHZS/PromLoop",
        "coordinate_claim": "81 bp, relative to TSS [-60,+20], as stated by upstream README",
        "orientation_status": "not independently verified; CSV has no strand field",
        "sigma_annotation": "not provided by these CSV files",
        "species": {},
    }

    for slug, directory in SPECIES.items():
        species_dir = dataset_root / directory
        dataset_path = species_dir / "Dataset.csv"
        positive_path = species_dir / "positive_samples.csv"
        dataset = _validate_rows(_read_rows(dataset_path), dataset_path, require_81bp=False)
        positives_file = _validate_rows(_read_rows(positive_path), positive_path)
        positives_dataset = [row for row in dataset if row["label"] == "1"]
        negatives = [row for row in dataset if row["label"] == "0"]
        valid_negatives = [row for row in negatives if SEQUENCE_RE.fullmatch(row["seq"])]
        excluded_controls = [
            {"seq_id": row["seq_id"], "length": len(row["seq"])}
            for row in negatives
            if not SEQUENCE_RE.fullmatch(row["seq"])
        ]

        expected_by_id = {row["seq_id"]: (row["seq"], row["label"]) for row in positives_dataset}
        provided_by_id = {row["seq_id"]: (row["seq"], row["label"]) for row in positives_file}
        if expected_by_id != provided_by_id:
            raise ValueError(f"{directory}: positive_samples.csv differs from Dataset.csv label=1 rows")
        if any(row["label"] != "1" for row in positives_file):
            raise ValueError(f"{directory}: positive_samples.csv contains non-positive labels")

        positive_target = output_dir / f"{slug}.positive.fasta"
        control_target = output_dir / f"{slug}.control.fasta"
        _write_fasta(positives_file, positive_target, label="1")
        _write_fasta(valid_negatives, control_target, label="0")

        split_rows: dict[str, list[dict[str, str]]] = {}
        split_ids: dict[str, set[str]] = {}
        for split in ("train", "dev", "test"):
            split_path = species_dir / f"{split}.csv"
            split_rows[split] = _validate_rows(
                _read_rows(split_path), split_path, require_81bp=False
            )
            split_ids[split] = {row["seq_id"] for row in split_rows[split]}
        if set.union(*split_ids.values()) != {row["seq_id"] for row in dataset}:
            raise ValueError(f"{directory}: train/dev/test IDs do not exactly cover Dataset.csv")
        if sum(map(len, split_rows.values())) != len(dataset):
            raise ValueError(f"{directory}: train/dev/test contain duplicate IDs across splits")

        seq_split_membership: dict[str, set[str]] = defaultdict(set)
        for split, rows in split_rows.items():
            for row in rows:
                seq_split_membership[row["seq"]].add(split)
        sequence_leakage = {
            seq: sorted(splits)
            for seq, splits in seq_split_membership.items()
            if len(splits) > 1
        }

        positive_seq_counts = Counter(row["seq"] for row in positives_file)
        report["species"][slug] = {
            "source_directory": directory,
            "dataset_csv": {
                "path": str(dataset_path),
                "sha256": sha256_file(dataset_path),
                "records": len(dataset),
                "positives": len(positives_dataset),
                "negatives": len(negatives),
                "invalid_length_records": [
                    {"seq_id": row["seq_id"], "length": len(row["seq"])}
                    for row in dataset
                    if not SEQUENCE_RE.fullmatch(row["seq"])
                ],
                "seq_type_values": dict(Counter(row["seq_type"] for row in dataset)),
            },
            "positive_samples_csv": {
                "path": str(positive_path),
                "sha256": sha256_file(positive_path),
                "records": len(positives_file),
                "unique_sequences": len(positive_seq_counts),
                "duplicate_sequence_records": sum(n - 1 for n in positive_seq_counts.values()),
                "fasta": str(positive_target),
                "fasta_sha256": sha256_file(positive_target),
            },
            "negative_control": {
                "records": len(valid_negatives),
                "excluded_invalid_length_records": excluded_controls,
                "unique_sequences": len({row["seq"] for row in valid_negatives}),
                "fasta": str(control_target),
                "fasta_sha256": sha256_file(control_target),
            },
            "split_records": {split: len(rows) for split, rows in split_rows.items()},
            "cross_split_exact_sequence_count": len(sequence_leakage),
            "positive_file_matches_label_1": True,
        }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-root",
        default="data/raw/tju_pan_promoter/reg_and_gen/Datasets",
        help="Directory containing the six species subdirectories",
    )
    parser.add_argument("--output-dir", default="tmp/tjupan_motif_inputs_20261003")
    parser.add_argument("--report", default="results/motif/tjupan_m3_20261003/input_manifest.json")
    args = parser.parse_args()
    report = prepare(Path(args.dataset_root), Path(args.output_dir), Path(args.report))
    for slug, item in report["species"].items():
        print(
            f"{slug}: records={item['dataset_csv']['records']} "
            f"positive={item['positive_samples_csv']['records']} "
            f"negative={item['negative_control']['records']} "
            f"duplicates={item['positive_samples_csv']['duplicate_sequence_records']} "
            f"cross_split_exact={item['cross_split_exact_sequence_count']}"
        )
    print(f"manifest={args.report}")


if __name__ == "__main__":
    main()
