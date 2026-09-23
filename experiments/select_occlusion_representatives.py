"""Select deterministic, sigma-balanced representatives for L2 occlusion."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from preprocessing.provenance import sha256_file


def run(
    input_path: str | Path,
    output_path: str | Path,
    *,
    group_column: str = "sigma_factor_type",
    per_group: int = 4,
) -> dict:
    """Select the lexicographically first rows in every metadata group."""

    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    if per_group < 1:
        raise ValueError("per_group must be positive")
    frame = pd.read_csv(input_path, sep="\t")
    required = {"sequence_id", "sequence", group_column}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"representative input is missing columns: {', '.join(missing)}")
    if frame["sequence_id"].astype(str).duplicated().any():
        raise ValueError("representative input contains duplicate sequence_id values")

    frame = frame.copy()
    frame[group_column] = frame[group_column].fillna("unknown").astype(str)
    frame.loc[frame[group_column].str.strip().eq(""), group_column] = "unknown"
    counts = frame[group_column].value_counts(sort=False).sort_index()
    insufficient = counts[counts < per_group]
    if not insufficient.empty:
        details = ", ".join(f"{name}={int(count)}" for name, count in insufficient.items())
        raise ValueError(f"some groups have fewer than {per_group} rows: {details}")

    selected_parts = []
    for group_name in sorted(counts.index.astype(str).tolist()):
        group = frame.loc[frame[group_column].eq(group_name)].sort_values(
            "sequence_id", kind="mergesort"
        )
        selected_parts.append(group.head(per_group))
    selected = pd.concat(selected_parts, ignore_index=True)
    selected = selected.sort_values([group_column, "sequence_id"], kind="mergesort").reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(output_path, sep="\t", index=False)

    manifest = {
        "schema_version": "l2-representative-selection-1.0",
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "output": str(output_path),
        "output_sha256": sha256_file(output_path),
        "group_column": group_column,
        "per_group": per_group,
        "input_group_counts": {str(name): int(count) for name, count in counts.items()},
        "selected_group_counts": {
            str(name): int(count)
            for name, count in selected[group_column].value_counts(sort=False).sort_index().items()
        },
        "sequence_count": int(len(selected)),
        "selection_rule": "Within each sigma group, select the lexicographically first sequence_id values.",
        "boundary": "Deterministic source-internal representatives for exploratory L2 occlusion; not a validation sample.",
    }
    manifest_path = output_path.with_name(f"{output_path.stem}.selection.json")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--group-column", default="sigma_factor_type")
    parser.add_argument("--per-group", type=int, default=4)
    args = parser.parse_args()
    print(
        json.dumps(
            run(
                args.input,
                args.output,
                group_column=args.group_column,
                per_group=args.per_group,
            ),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
