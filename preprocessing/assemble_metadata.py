"""Assemble already-audited source tables into one promoter metadata table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .pipeline import load_metadata
from .provenance import sha256_file
from .schema import REQUIRED_COLUMNS


def assemble_source_tables(
    input_paths: Iterable[str | Path],
    *,
    output_path: str | Path,
    report_path: str | Path,
) -> dict[str, Any]:
    """Concatenate canonical source tables without silently repairing records."""

    paths = [Path(path).resolve() for path in input_paths]
    if not paths:
        raise ValueError("at least one input table is required")
    frames: list[pd.DataFrame] = []
    source_reports: list[dict[str, Any]] = []
    for path in paths:
        frame = load_metadata(path)
        missing = sorted(set(REQUIRED_COLUMNS) - set(frame.columns))
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        if frame["source_dataset"].isna().any() or frame["source_dataset"].astype(str).str.strip().eq("").any():
            raise ValueError(f"{path} contains missing source_dataset values")
        if frame["sequence_id"].isna().any() or frame["sequence_id"].astype(str).str.strip().eq("").any():
            raise ValueError(f"{path} contains missing sequence_id values")
        frames.append(frame)
        source_reports.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "rows": int(len(frame)),
                "sources": sorted(frame["source_dataset"].astype(str).unique().tolist()),
                "species": sorted(frame["species"].dropna().astype(str).unique().tolist()),
                "sigma_counts": frame["sigma_factor_type"].fillna("unknown").astype(str).value_counts().to_dict()
                if "sigma_factor_type" in frame
                else {},
            }
        )

    merged = pd.concat(frames, ignore_index=True, sort=False)
    duplicate_ids = merged["sequence_id"].astype(str).duplicated(keep=False)
    if duplicate_ids.any():
        examples = sorted(merged.loc[duplicate_ids, "sequence_id"].astype(str).unique().tolist())[:10]
        raise ValueError(
            "sequence_id must be globally unique; duplicate examples: " + ", ".join(examples)
        )
    merged = merged.sort_values(["source_dataset", "sequence_id"]).reset_index(drop=True)
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, sep="\t", index=False)
    report = {
        "input_tables": source_reports,
        "output_path": str(output),
        "output_sha256": sha256_file(output),
        "output_rows": int(len(merged)),
        "species_counts": merged["species"].fillna("unknown").astype(str).value_counts().to_dict(),
        "source_counts": merged["source_dataset"].astype(str).value_counts().to_dict(),
        "sigma_counts": merged["sigma_factor_type"].fillna("unknown").astype(str).value_counts().to_dict()
        if "sigma_factor_type" in merged
        else {},
        "missing_counts": {
            column: int(merged[column].isna().sum())
            for column in merged.columns
            if merged[column].isna().any()
        },
    }
    report_file = Path(report_path).resolve()
    report_file.parent.mkdir(parents=True, exist_ok=True)
    report_file.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, help="Canonical source TSV; repeat for each source")
    parser.add_argument("--output", required=True, help="Combined metadata TSV")
    parser.add_argument("--report", required=True, help="Assembly manifest JSON")
    args = parser.parse_args()
    report = assemble_source_tables(args.input, output_path=args.output, report_path=args.report)
    print(f"output_rows={report['output_rows']}")
    print(f"output_sha256={report['output_sha256']}")
    print(f"output={report['output_path']}")


if __name__ == "__main__":
    main()
