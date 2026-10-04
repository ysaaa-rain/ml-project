"""Test whether TJU Pan E. coli motifs associate with a second strength dataset."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import itertools
import json
import math
from pathlib import Path
import subprocess
import time
from typing import Any, Iterable

import numpy as np
from scipy.stats import mannwhitneyu, rankdata, spearmanr, t as student_t

from preprocessing.provenance import sha256_file


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STRENGTH = ROOT / "data/raw/tju_pan_promoter/strenth/data/supplementary/E_coli.txt"
DEFAULT_TJUPAN_FASTA = ROOT / "tmp/tjupan_fimo_inputs_20261003/escherichia_coli.combined.fasta"
DEFAULT_MOTIF_FILE = ROOT / "results/motif/tjupan_m3_20261003/summary/matrices/escherichia_coli.meme"
DEFAULT_MOTIF_MANIFEST = ROOT / "results/motif/tjupan_m3_20261003/summary/summary_manifest.json"
DEFAULT_FIMO = ROOT / "tmp/meme-suite-5.5.9/bin/fimo"
DEFAULT_RUN_DIR = ROOT / "tmp/experiment5_strength/tjupan_ecoli_20261003"
DEFAULT_RESULT_DIR = ROOT / "results/motif/tjupan_m3_20261003/strength"
DNA = set("ACGT")


def _portable_command(command: list[str]) -> list[str]:
    """Avoid machine-specific absolute paths in committed run manifests."""
    normalized: list[str] = []
    for item in command:
        path = Path(item)
        if path.is_absolute():
            try:
                item = path.resolve().relative_to(ROOT.resolve()).as_posix()
            except ValueError:
                pass
        normalized.append(item)
    return normalized
HIT_THRESHOLDS = (0.05, 0.01, 0.001)


def _read_strength_table(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"strength file has no header: {path}")
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        if not {"strength", "promoter"} <= set(reader.fieldnames):
            raise ValueError(f"expected strength and promoter columns, got {reader.fieldnames}")
        rows: list[dict[str, Any]] = []
        for row_number, raw in enumerate(reader, start=2):
            sequence = (raw.get("promoter") or "").strip().upper()
            try:
                strength = float((raw.get("strength") or "").strip())
            except ValueError as error:
                raise ValueError(f"invalid strength at row {row_number}") from error
            if len(sequence) != 50 or not set(sequence) <= DNA:
                raise ValueError(f"expected a 50-nt A/C/G/T promoter at row {row_number}")
            if not math.isfinite(strength) or strength <= 0:
                raise ValueError(f"strength must be finite and positive at row {row_number}")
            rows.append({"sequence": sequence, "strength": strength})
    if not rows:
        raise ValueError("strength table is empty")
    if len({row["sequence"] for row in rows}) != len(rows):
        raise ValueError("strength table contains duplicate sequences; resolve sample semantics first")
    return rows


def _read_81bp_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if header is None:
            return
        sequence = "".join(chunks).upper()
        if len(sequence) != 81 or not set(sequence) <= DNA:
            raise ValueError(f"expected an 81-nt A/C/G/T sequence in {path}")
        if "|label=1" in header:
            label = "positive"
        elif "|label=0" in header:
            label = "control"
        else:
            raise ValueError(f"missing label=1/0 in TJU Pan FASTA header: {header}")
        records.append((header, sequence, label))

    for line_number, raw_line in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header, chunks = line[1:], []
        elif header is None:
            raise ValueError(f"sequence before FASTA header at {path}:{line_number}")
        else:
            chunks.append(line)
    finish()
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    return records


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _find_overlaps(
    strength_rows: list[dict[str, Any]], tjupan_records: list[tuple[str, str, str]]
) -> tuple[set[str], dict[str, int]]:
    windows: dict[str, set[str]] = {"positive": set(), "control": set()}
    for _, sequence, label in tjupan_records:
        for start in range(len(sequence) - 50 + 1):
            windows[label].add(sequence[start : start + 50])
    matches: dict[str, set[str]] = {}
    for label, candidates in windows.items():
        candidates_with_rc = candidates | {_reverse_complement(item) for item in candidates}
        matches[label] = {row["sequence"] for row in strength_rows} & candidates_with_rc
    excluded = set().union(*matches.values())
    stats = {
        "tju_pan_positive_sequences": len({seq for _, seq, label in tjupan_records if label == "positive"}),
        "tju_pan_control_sequences": len({seq for _, seq, label in tjupan_records if label == "control"}),
        "strength_sequences_matching_positive_any_50nt_window_or_reverse_complement": len(matches["positive"]),
        "strength_sequences_matching_control_any_50nt_window_or_reverse_complement": len(matches["control"]),
        "strength_sequences_matching_both_labels": len(matches["positive"] & matches["control"]),
        "strength_sequences_excluded_for_overlap": len(excluded),
    }
    return excluded, stats


def _write_scan_inputs(
    retained: list[dict[str, Any]], fasta_path: Path, background_path: Path
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    indexed: list[dict[str, Any]] = []
    base_counts = {base: 0 for base in "ACGT"}
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for index, row in enumerate(retained, start=1):
            sequence_id = f"strength_{index:06d}"
            sequence = row["sequence"]
            row_record = {
                "sequence_id": sequence_id,
                "sequence": sequence,
                "strength": row["strength"],
                "log10_strength": math.log10(row["strength"]),
                "gc_fraction": (sequence.count("G") + sequence.count("C")) / len(sequence),
            }
            indexed.append(row_record)
            handle.write(f">{sequence_id}\n{sequence}\n")
            for base in sequence:
                base_counts[base] += 1
    total = sum(base_counts.values())
    frequencies = {base: base_counts[base] / total for base in "ACGT"}
    background_path.parent.mkdir(parents=True, exist_ok=True)
    background_path.write_text(
        "".join(f"{base} {frequencies[base]:.12g}\n" for base in "ACGT"), encoding="ascii"
    )
    return {"base_counts": base_counts, "frequencies": frequencies}, indexed


def _motif_ids(motif_file: Path) -> list[str]:
    ids = [line.split()[1] for line in motif_file.read_text(encoding="utf-8").splitlines() if line.startswith("MOTIF ")]
    if len(ids) != 10 or len(ids) != len(set(ids)):
        raise ValueError(f"expected 10 unique E. coli motifs, got {len(ids)}")
    return ids


def _bh(p_values: Iterable[float]) -> list[float]:
    values = [float(value) for value in p_values]
    if any(not 0 <= value <= 1 for value in values):
        raise ValueError("p-values must be in [0, 1]")
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.0] * len(values)
    running = 1.0
    count = len(values)
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original_index] * count / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def _partial_spearman(x: np.ndarray, y: np.ndarray, covariate: np.ndarray) -> tuple[float, float]:
    if len(x) != len(y) or len(x) != len(covariate) or len(x) < 4:
        raise ValueError("partial Spearman requires aligned arrays with at least four samples")
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(covariate)
    design = np.column_stack((np.ones(len(rz)), rz))
    residual_x = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    residual_y = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    if np.std(residual_x) == 0 or np.std(residual_y) == 0:
        return 0.0, 1.0
    rho = float(np.corrcoef(residual_x, residual_y)[0, 1])
    rho = min(1.0, max(-1.0, rho))
    degrees = len(x) - 3
    if abs(rho) == 1:
        return rho, 0.0
    statistic = abs(rho) * math.sqrt(degrees / max(1e-300, 1 - rho * rho))
    return rho, float(2 * student_t.sf(statistic, degrees))


def _consume_fimo_stdout(
    lines: Iterable[str], motif_ids: list[str], sequence_ids: set[str]
) -> tuple[dict[str, dict[str, dict[str, Any]]], int]:
    iterator = iter(lines)
    header_line = next((line for line in iterator if line.strip() and not line.startswith("#")), None)
    if header_line is None:
        raise ValueError("FIMO emitted no TSV header")
    reader = csv.DictReader(itertools.chain([header_line], iterator), delimiter="\t")
    required = {"motif_id", "sequence_name", "score", "p-value"}
    if not reader.fieldnames or not required <= set(reader.fieldnames):
        raise ValueError(f"unexpected FIMO text columns: {reader.fieldnames}")
    state = {
        motif_id: {
            sequence_id: {"best_score": -math.inf, "min_p_value": 1.0, **{f"hits_p_le_{threshold:g}": 0 for threshold in HIT_THRESHOLDS}}
            for sequence_id in sequence_ids
        }
        for motif_id in motif_ids
    }
    site_count = 0
    for row in reader:
        motif_id, sequence_id = row["motif_id"], row["sequence_name"]
        if motif_id not in state or sequence_id not in sequence_ids:
            raise ValueError(f"FIMO output contains unknown motif or sequence: {motif_id}, {sequence_id}")
        try:
            score, p_value = float(row["score"]), float(row["p-value"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid FIMO score/p-value row: {row}") from error
        if not math.isfinite(score) or not 0 <= p_value <= 1:
            raise ValueError(f"invalid FIMO score/p-value row: {row}")
        item = state[motif_id][sequence_id]
        item["best_score"] = max(item["best_score"], score)
        item["min_p_value"] = min(item["min_p_value"], p_value)
        for threshold in HIT_THRESHOLDS:
            if p_value <= threshold:
                item[f"hits_p_le_{threshold:g}"] += 1
        site_count += 1
    return state, site_count


def _build_association_tables(
    fimo_state: dict[str, dict[str, dict[str, Any]]], indexed_rows: list[dict[str, Any]], motif_ids: list[str]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    sequence_ids = [row["sequence_id"] for row in indexed_rows]
    by_id = {row["sequence_id"]: row for row in indexed_rows}
    log_strength = np.asarray([by_id[sequence_id]["log10_strength"] for sequence_id in sequence_ids])
    gc_fraction = np.asarray([by_id[sequence_id]["gc_fraction"] for sequence_id in sequence_ids])
    main: list[dict[str, Any]] = []
    threshold_rows: list[dict[str, Any]] = []
    for motif_id in motif_ids:
        results = [fimo_state[motif_id][sequence_id] for sequence_id in sequence_ids]
        scores = np.asarray([item["best_score"] if math.isfinite(item["best_score"]) else -1e6 for item in results])
        spearman = spearmanr(scores, log_strength)
        partial_rho, partial_p = _partial_spearman(scores, log_strength, gc_fraction)
        row: dict[str, Any] = {
            "motif_id": motif_id,
            "n_sequences": len(sequence_ids),
            "spearman_best_site_score_rho": float(spearman.statistic),
            "spearman_best_site_score_p": float(spearman.pvalue),
            "gc_adjusted_partial_spearman_rho": partial_rho,
            "gc_adjusted_partial_spearman_p": partial_p,
            "median_best_site_score": float(np.median(scores)),
            "median_log10_strength": float(np.median(log_strength)),
        }
        for threshold in HIT_THRESHOLDS:
            key = f"hits_p_le_{threshold:g}"
            present = np.asarray([item[key] > 0 for item in results])
            group_hit, group_no_hit = log_strength[present], log_strength[~present]
            if len(group_hit) and len(group_no_hit):
                test = mannwhitneyu(group_hit, group_no_hit, alternative="two-sided", method="asymptotic")
                delta_log = float(np.median(group_hit) - np.median(group_no_hit))
                effect = float(2 * test.statistic / (len(group_hit) * len(group_no_hit)) - 1)
                p_value = float(test.pvalue)
                hit_median = float(np.median(group_hit))
                no_hit_median = float(np.median(group_no_hit))
            else:
                delta_log, effect, p_value, hit_median, no_hit_median = 0.0, 0.0, 1.0, math.nan, math.nan
            threshold_rows.append(
                {
                    "motif_id": motif_id,
                    "site_p_threshold": threshold,
                    "n_with_hit": int(present.sum()),
                    "n_without_hit": int((~present).sum()),
                    "median_log10_strength_with_hit": hit_median,
                    "median_log10_strength_without_hit": no_hit_median,
                    "median_difference_log10": delta_log,
                    "median_fold_difference": 10**delta_log if abs(delta_log) < 300 else math.inf,
                    "cliffs_delta": effect,
                    "mann_whitney_p": p_value,
                }
            )
        main.append(row)
    for field, p_field in (
        ("spearman_best_site_score_q_bh", "spearman_best_site_score_p"),
        ("gc_adjusted_partial_spearman_q_bh", "gc_adjusted_partial_spearman_p"),
    ):
        adjusted = _bh(row[p_field] for row in main)
        for row, q_value in zip(main, adjusted, strict=True):
            row[field] = q_value
    adjusted_threshold = _bh(row["mann_whitney_p"] for row in threshold_rows)
    for row, q_value in zip(threshold_rows, adjusted_threshold, strict=True):
        row["mann_whitney_q_bh_30_tests"] = q_value
    return main, threshold_rows


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(
    strength_file: Path = DEFAULT_STRENGTH,
    tjupan_fasta: Path = DEFAULT_TJUPAN_FASTA,
    motif_file: Path = DEFAULT_MOTIF_FILE,
    motif_manifest: Path = DEFAULT_MOTIF_MANIFEST,
    fimo: Path = DEFAULT_FIMO,
    run_dir: Path = DEFAULT_RUN_DIR,
    result_dir: Path = DEFAULT_RESULT_DIR,
) -> dict[str, Any]:
    for path in (strength_file, tjupan_fasta, motif_file, motif_manifest, fimo):
        if not path.is_file():
            raise FileNotFoundError(path)
    source_motif_manifest = json.loads(motif_manifest.read_text(encoding="utf-8"))
    expected_matrix = source_motif_manifest.get("output_sha256", {}).get("matrices/escherichia_coli.meme", {}).get("sha256")
    if source_motif_manifest.get("status") != "complete" or not expected_matrix:
        raise ValueError("complete E. coli MEME summary is required")
    if sha256_file(motif_file) != expected_matrix:
        raise ValueError("E. coli MEME matrix hash differs from its summary manifest")
    if result_dir.exists() and any(result_dir.iterdir()):
        raise FileExistsError(f"result directory is not empty; use a new --result-dir: {result_dir}")
    if run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(f"run directory is not empty; use a new --run-dir: {run_dir}")

    version = subprocess.run([str(fimo), "--version"], capture_output=True, text=True, check=True, timeout=15)
    fimo_version = (version.stdout or version.stderr).strip().splitlines()[0]
    if fimo_version != "5.5.9":
        raise ValueError(f"expected FIMO 5.5.9, got {fimo_version}")

    strength_rows = _read_strength_table(strength_file)
    tjupan_records = _read_81bp_fasta(tjupan_fasta)
    excluded, overlap_stats = _find_overlaps(strength_rows, tjupan_records)
    retained = [row for row in strength_rows if row["sequence"] not in excluded]
    if not retained:
        raise ValueError("no strength sequences remain after overlap removal")
    fasta_path, background_path = run_dir / "strength_disjoint.fasta", run_dir / "target.background"
    background, indexed_rows = _write_scan_inputs(retained, fasta_path, background_path)
    motif_ids = _motif_ids(motif_file)
    result_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)
    # _write_scan_inputs creates the run directory before the collision check, so this is guarded above.
    stream_audit_path = run_dir / "fimo_stream_audit.txt"
    stderr_path = run_dir / "fimo_stderr.log"
    command = [
        str(fimo), "--text", "--skip-matched-sequence", "--no-qvalue", "--thresh", "1",
        "--motif-pseudo", "0.1", "--bfile", str(background_path), "--max-stored-scores", "15000000",
        str(motif_file), str(fasta_path),
    ]
    started_at = datetime.now().astimezone().isoformat(timespec="seconds")
    started_clock = time.monotonic()
    with stderr_path.open("w", encoding="utf-8") as stderr_handle:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr_handle, text=True, bufsize=1)
        if process.stdout is None:
            raise RuntimeError("failed to open FIMO standard output")
        # Keep only per-sequence summaries in memory. The full raw site table is not retained.
        site_state, scanned_sites = _consume_fimo_stdout(process.stdout, motif_ids, {row["sequence_id"] for row in indexed_rows})
        return_code = process.wait()
    stream_audit_path.write_text(f"streamed_site_rows={scanned_sites}\n", encoding="utf-8")
    elapsed = time.monotonic() - started_clock
    if return_code != 0:
        tail = "\n".join(stderr_path.read_text(encoding="utf-8").splitlines()[-30:])
        raise RuntimeError(f"FIMO exited with status {return_code}; see {stderr_path}\n{tail}")

    main_rows, threshold_rows = _build_association_tables(site_state, indexed_rows, motif_ids)
    main_path = result_dir / "motif_strength_associations.tsv"
    threshold_path = result_dir / "motif_strength_hit_threshold_sensitivity.tsv"
    _write_tsv(main_path, main_rows)
    _write_tsv(threshold_path, threshold_rows)
    manifest_path = result_dir / "motif_strength_manifest.json"
    manifest: dict[str, Any] = {
        "schema_version": "pr01-02-exp5-strength-1.0",
        "status": "complete",
        "purpose": "关联 TJU Pan E. coli de novo motif 得分与另一份 dRNA-seq 派生启动子数值标签，不训练强度预测模型",
        "source_description": "TJU cloud shared strength files; E. coli K-12 MG1655 natural promoter set described in mSystems 2025; numeric labels are the supplied dRNA-seq-derived values",
        "source_files": {
            "strength_table": {"path": str(strength_file.relative_to(ROOT)), "sha256": sha256_file(strength_file), "bytes": strength_file.stat().st_size},
            "tjupan_ecoli_fasta": {"path": str(tjupan_fasta.relative_to(ROOT)), "sha256": sha256_file(tjupan_fasta), "bytes": tjupan_fasta.stat().st_size},
            "meme_matrix": {"path": str(motif_file.relative_to(ROOT)), "sha256": sha256_file(motif_file), "bytes": motif_file.stat().st_size},
            "meme_summary_manifest": {"path": str(motif_manifest.relative_to(ROOT)), "sha256": sha256_file(motif_manifest)},
        },
        "fimo": {"version": fimo_version, "path": str(fimo.relative_to(ROOT))},
        "counts": {
            "strength_records": len(strength_rows),
            "unique_strength_sequences": len(strength_rows),
            "retained_after_conservative_overlap_exclusion": len(retained),
            "fimo_sites_streamed": scanned_sites,
            **overlap_stats,
        },
        "parameters": {
            "sequence_window_overlap_rule": "remove any 50-mer equal to a 50-nt substring of any TJU Pan E. coli 81-nt positive/control sequence, in either orientation; windows are checked at all 32 offsets",
            "motif_count": len(motif_ids),
            "scan_strands": "both strands according to MEME matrix; reverse-complement scanning enabled",
            "fimo_site_p_threshold_for_streaming": 1.0,
            "motif_pseudocount": 0.1,
            "background": "zero-order frequencies estimated from retained strength-set sequences",
            "primary_test": "Spearman correlation between each promoter's maximum FIMO site score across both strands/positions and log10 supplied dRNA-seq value; BH across 10 motifs",
            "composition_sensitivity": "partial Spearman correlation on ranked variables, adjusting for promoter GC fraction; BH across 10 motifs",
            "hit_sensitivity": "Mann-Whitney U tests on log10 values, sequence-positive if any site p<=0.05, 0.01, or 0.001; BH across all 30 motif-threshold tests",
            "raw_site_output": "streamed and aggregated; individual sites, IDs, and DNA sequences are not retained in report-facing results",
        },
        "background": background,
        "strength_distribution": {
            "minimum": min(row["strength"] for row in retained),
            "median": float(np.median([row["strength"] for row in retained])),
            "maximum": max(row["strength"] for row in retained),
            "log10_median": float(np.median([row["log10_strength"] for row in indexed_rows])),
        },
        "execution": {"started_at": started_at, "elapsed_seconds": round(elapsed, 3), "return_code": return_code, "command": _portable_command(command)},
        "outputs": {},
    }
    for path in (main_path, threshold_path):
        manifest["outputs"][path.name] = {"sha256": sha256_file(path), "bytes": path.stat().st_size, "rows": len(main_rows) if path == main_path else len(threshold_rows)}
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strength-file", type=Path, default=DEFAULT_STRENGTH)
    parser.add_argument("--tjupan-fasta", type=Path, default=DEFAULT_TJUPAN_FASTA)
    parser.add_argument("--motif-file", type=Path, default=DEFAULT_MOTIF_FILE)
    parser.add_argument("--motif-manifest", type=Path, default=DEFAULT_MOTIF_MANIFEST)
    parser.add_argument("--fimo", type=Path, default=DEFAULT_FIMO)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    args = parser.parse_args()
    manifest = run(**vars(args))
    print(f"Experiment 5 status: {manifest['status']}")
    print(f"Strength records retained: {manifest['counts']['retained_after_conservative_overlap_exclusion']}")
    print(f"FIMO site rows streamed: {manifest['counts']['fimo_sites_streamed']}")
    print(f"Results: {args.result_dir}")


if __name__ == "__main__":
    main()
