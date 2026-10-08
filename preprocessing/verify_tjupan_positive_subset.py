"""Strictly audit whether TJU positive_samples.csv equals Dataset.csv label=1.

Read-only. Checks six species, every row and every field, unique IDs and
frozen raw-file SHA256. No data processing or motif experiment is performed.

Run:
    python -m preprocessing.verify_tjupan_positive_subset \
      --report tmp/tju_positive_subset_verification.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FOLDERS = (
    "Bacillus_subtilis",
    "Baumanii",
    "Bradyrhizobium",
    "Diphtheria",
    "Escherichia_coli",
    "Staphylococcus",
)
REQUIRED = {"seq_id", "seq_type", "seq", "label"}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames or []
        if len(headers) != len(set(headers)) or not REQUIRED.issubset(headers):
            raise ValueError(f"{path}: invalid, repeated or missing header: {headers}")
        rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError(f"{path}: row has more columns than header")
    if any(any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"{path}: row has missing columns")
    return headers, rows


def verify_pair(dataset_path: Path, positive_path: Path) -> dict:
    dheads, dataset = load_rows(dataset_path)
    pheads, positives = load_rows(positive_path)
    errors: list[str] = []
    if set(dheads) != set(pheads):
        errors.append(f"headers differ: dataset={dheads} positives={pheads}")

    def index(rows: list[dict[str, str]], where: str) -> dict[str, dict[str, str]]:
        ids = [row["seq_id"] for row in rows]
        duplicate = sorted(k for k, n in Counter(ids).items() if n > 1)
        if duplicate:
            errors.append(f"{where} duplicate seq_id: {duplicate[:20]} ({len(duplicate)} total)")
        if any(not x for x in ids):
            errors.append(f"{where} has empty seq_id")
        return {row["seq_id"]: row for row in rows}

    all_by_id = index(dataset, "Dataset.csv")
    pos_by_id = index(positives, "positive_samples.csv")
    reference = {key: row for key, row in all_by_id.items() if row["label"] == "1"}
    non_positive = sorted(key for key, row in pos_by_id.items() if row["label"] != "1")
    missing = sorted(set(reference) - set(pos_by_id))
    unexpected = sorted(set(pos_by_id) - set(reference))
    differing = sorted(
        key for key in set(reference) & set(pos_by_id)
        if reference[key] != pos_by_id[key]
    )
    if non_positive:
        errors.append(f"positive_samples has non-label=1 rows: {non_positive[:20]}")
    if missing:
        errors.append(f"positive_samples missing Dataset label=1 IDs: {missing[:20]}")
    if unexpected:
        errors.append(f"positive_samples has extra IDs: {unexpected[:20]}")
    if differing:
        errors.append(f"field-level mismatch IDs: {differing[:20]}")
    if len(reference) != len([row for row in dataset if row["label"] == "1"]):
        errors.append("Dataset label=1 count differs from unique-ID count")
    if len(pos_by_id) != len(positives):
        errors.append("positive_samples contains duplicate IDs")
    return {
        "dataset_rows": len(dataset),
        "dataset_positive_rows": sum(row["label"] == "1" for row in dataset),
        "positive_samples_rows": len(positives),
        "dataset_unique_ids": len(all_by_id),
        "positive_unique_ids": len(pos_by_id),
        "missing_positive_ids": len(missing),
        "unexpected_positive_ids": len(unexpected),
        "different_records": len(differing),
        "non_positive_records": len(non_positive),
        "errors": errors,
        "identical_positive_rows_by_id": not errors,
    }


def run(root: Path = ROOT, manifest_path: Path | None = None) -> dict:
    root = root.resolve()
    manifest_path = manifest_path or root / "data/processed/pr01_02_data_v2/manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    hashes = manifest.get("raw_sha256", {})
    result = {
        "status": "passed",
        "scope": "All six TJU species; exact row equality by unique ID, including every column",
        "frozen_manifest": str(manifest_path),
        "species": {},
    }
    for species in FOLDERS:
        base = root / "data/raw/tju_pan_promoter/reg_and_gen/Datasets" / species
        dpath, ppath = base / "Dataset.csv", base / "positive_samples.csv"
        detail = verify_pair(dpath, ppath)
        for path in (dpath, ppath):
            relative = path.relative_to(root).as_posix()
            expected = hashes.get(relative)
            actual = file_sha256(path)
            detail[path.name + "_sha256"] = actual
            if expected is None or actual != expected:
                detail["errors"].append(f"raw SHA256 mismatch/missing: {relative}")
        detail["passed"] = not detail["errors"]
        result["species"][species] = detail
        if not detail["passed"]:
            result["status"] = "failed"
    result["totals"] = {
        "dataset_positive_rows": sum(d["dataset_positive_rows"] for d in result["species"].values()),
        "positive_samples_rows": sum(d["positive_samples_rows"] for d in result["species"].values()),
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--report", type=Path, help="Optional JSON evidence output path")
    args = parser.parse_args()
    try:
        result = run(root=args.root, manifest_path=args.manifest)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"verification failed: {exc}", file=sys.stderr)
        return 2
    payload = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
    print(payload)
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
