"""Measure known-reference PWM hits across raw-p and sequence-coordinate thresholds."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import time
from typing import Any

import numpy as np

from experiments.prepare_tjupan_motif_inputs import SPECIES
from experiments.summarize_tjupan_fimo import benjamini_hochberg, fisher_exact_two_sided
from experiments.summarize_tjupan_reference_overlap import REFERENCES, _read_fasta, _write_tsv
from preprocessing.provenance import sha256_file


def _read_raw_hits(path: Path) -> list[dict[str, Any]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter="\t")
    required = {"motif_id", "sequence_name", "start", "stop", "strand", "score", "p-value"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise ValueError(f"FIMO output is missing columns: {path}")
    hits: list[dict[str, Any]] = []
    for row in reader:
        hit = {
            "motif_id": row["motif_id"],
            "sequence_name": row["sequence_name"],
            "start": int(row["start"]),
            "stop": int(row["stop"]),
            "strand": row["strand"],
            "score": float(row["score"]),
            "p_value": float(row["p-value"]),
        }
        if hit["motif_id"] not in REFERENCES or hit["strand"] not in {"+", "-"}:
            raise ValueError(f"unexpected FIMO hit: {row}")
        if not 1 <= hit["start"] <= hit["stop"] <= 81 or not 0 <= hit["p_value"] <= 1:
            raise ValueError(f"invalid FIMO hit: {row}")
        hits.append(hit)
    return hits


def _read_fasta_with_gc(path: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    identifier: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if identifier is None:
            return
        sequence = "".join(chunks).upper()
        if len(sequence) != 81 or set(sequence) - set("ACGT"):
            raise ValueError(f"expected labeled 81-bp DNA sequence in {path}: {identifier}")
        label = identifier.rsplit("|label=", 1)[-1]
        if label not in {"0", "1"} or identifier in records:
            raise ValueError(f"invalid or repeated FASTA record: {identifier}")
        records[identifier] = {
            "label": label,
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "gc_count": sequence.count("G") + sequence.count("C"),
        }

    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            identifier, chunks = line[1:], []
        elif identifier is None:
            raise ValueError(f"sequence before FASTA header in {path}")
        else:
            chunks.append(line)
    finish()
    return records


def summarize_hits(
    species: str,
    sequences: dict[str, dict[str, str]],
    hits: list[dict[str, Any]],
    thresholds: tuple[float, ...],
) -> list[dict[str, Any]]:
    by_motif: dict[str, list[dict[str, Any]]] = {motif: [] for motif in REFERENCES}
    for hit in hits:
        if hit["sequence_name"] not in sequences:
            raise ValueError(f"unknown sequence in FIMO output: {hit['sequence_name']}")
        by_motif[hit["motif_id"]].append(hit)

    unique_by_label = {
        label: {record["sequence_sha256"] for record in sequences.values() if record["label"] == label}
        for label in ("1", "0")
    }
    shared = unique_by_label["1"] & unique_by_label["0"]
    rows: list[dict[str, Any]] = []
    p_groups: dict[float, list[dict[str, Any]]] = {threshold: [] for threshold in thresholds}

    for motif in REFERENCES:
        motif_hits = by_motif[motif]
        for threshold in thresholds:
            site_count = sum(hit["p_value"] <= threshold for hit in motif_hits)
            best_by_label: dict[str, dict[str, dict[str, Any]]] = {"1": {}, "0": {}}
            for hit in motif_hits:
                if hit["p_value"] > threshold:
                    continue
                record = sequences[hit["sequence_name"]]
                label = record["label"]
                sequence_hash = record["sequence_sha256"]
                current = best_by_label[label].get(sequence_hash)
                if current is None or (hit["p_value"], hit["start"], hit["strand"]) < (
                    current["p_value"], current["start"], current["strand"]
                ):
                    best_by_label[label][sequence_hash] = hit

            all_hit_ids = {label: set(best_by_label[label]) for label in ("1", "0")}
            pos_hit, ctl_hit = all_hit_ids["1"] - shared, all_hit_ids["0"] - shared
            pos_total, ctl_total = unique_by_label["1"] - shared, unique_by_label["0"] - shared
            p_value = fisher_exact_two_sided(
                len(pos_hit), len(pos_total) - len(pos_hit),
                len(ctl_hit), len(ctl_total) - len(ctl_hit),
            )

            def position_summary(label: str) -> tuple[float | None, float | None]:
                selected = list(best_by_label[label].values())
                midpoints = [(hit["start"] + hit["stop"]) / 2 for hit in selected]
                plus_fraction = sum(hit["strand"] == "+" for hit in selected) / len(selected) if selected else None
                return (statistics.median(midpoints) if midpoints else None, plus_fraction)

            pos_median, pos_plus = position_summary("1")
            ctl_median, ctl_plus = position_summary("0")
            row = {
                "species": species,
                "reference_motif": motif,
                "site_p_cutoff": threshold,
                "raw_site_count": site_count,
                "positive_unique_hit_sequences": len(all_hit_ids["1"]),
                "positive_unique_sequences": len(unique_by_label["1"]),
                "positive_coverage": len(all_hit_ids["1"]) / len(unique_by_label["1"]),
                "control_unique_hit_sequences": len(all_hit_ids["0"]),
                "control_unique_sequences": len(unique_by_label["0"]),
                "control_coverage": len(all_hit_ids["0"]) / len(unique_by_label["0"]),
                "coverage_difference": (
                    len(all_hit_ids["1"]) / len(unique_by_label["1"])
                    - len(all_hit_ids["0"]) / len(unique_by_label["0"])
                ),
                "cross_label_shared_sequences_excluded_from_fisher": len(shared),
                "fisher_positive_exclusive_n": len(pos_total),
                "fisher_control_exclusive_n": len(ctl_total),
                "fisher_p_value": p_value,
                "fisher_bh_q_24_tests": None,
                "positive_best_hit_median_midpoint_1based": pos_median,
                "positive_best_hit_plus_strand_fraction": pos_plus,
                "control_best_hit_median_midpoint_1based": ctl_median,
                "control_best_hit_plus_strand_fraction": ctl_plus,
            }
            rows.append(row)
            p_groups[threshold].append(row)

    for group in p_groups.values():
        adjusted = benjamini_hochberg([row["fisher_p_value"] for row in group])
        for row, q_value in zip(group, adjusted, strict=True):
            row["fisher_bh_q_24_tests"] = q_value
    return rows


def summarize_position_bins(
    species: str,
    sequences: dict[str, dict[str, str]],
    hits: list[dict[str, Any]],
    thresholds: tuple[float, ...],
) -> list[dict[str, Any]]:
    """Count unique sequences with any qualifying site by input coordinate/strand."""
    output: list[dict[str, Any]] = []
    for motif in REFERENCES:
        motif_hits = [hit for hit in hits if hit["motif_id"] == motif]
        for threshold in thresholds:
            counts: dict[tuple[str, int, str], set[str]] = {}
            for hit in motif_hits:
                if hit["p_value"] > threshold:
                    continue
                sequence = sequences[hit["sequence_name"]]
                key = (sequence["label"], hit["start"], hit["strand"])
                counts.setdefault(key, set()).add(sequence["sequence_sha256"])
            for (label, start, strand), sequence_hashes in sorted(counts.items()):
                output.append({
                    "species": species,
                    "reference_motif": motif,
                    "site_p_cutoff": threshold,
                    "label": label,
                    "start_1based": start,
                    "strand": strand,
                    "unique_sequences_with_site": len(sequence_hashes),
                })
    return output


def summarize_ecoli_expected_windows(
    sequences: dict[str, dict[str, Any]],
    hits: list[dict[str, Any]],
    thresholds: tuple[float, ...],
    permutations: int = 10_000,
) -> list[dict[str, Any]]:
    """Test prespecified sigma70-like input windows separately on both strands."""
    windows = {"Ecoli_sigma70_minus10": (49, 54), "Ecoli_sigma70_minus35": (24, 29)}
    if permutations < 1:
        raise ValueError("permutations must be positive")
    rows: list[dict[str, Any]] = []
    permutation_p_values: list[float] = []
    for threshold in thresholds:
        threshold_rows: list[dict[str, Any]] = []
        for motif, (start_min, start_max) in windows.items():
            motif_hits = [hit for hit in hits if hit["motif_id"] == motif]
            for strand in ("+", "-"):
                hit_hashes = {"1": set(), "0": set()}
                for hit in motif_hits:
                    if hit["p_value"] <= threshold and start_min <= hit["start"] <= start_max and hit["strand"] == strand:
                        record = sequences[hit["sequence_name"]]
                        hit_hashes[record["label"]].add(record["sequence_sha256"])
                total_hashes = {
                    label: {r["sequence_sha256"] for r in sequences.values() if r["label"] == label}
                    for label in ("1", "0")
                }
                shared = total_hashes["1"] & total_hashes["0"]
                pos_total, ctl_total = total_hashes["1"] - shared, total_hashes["0"] - shared
                pos_hit, ctl_hit = hit_hashes["1"] - shared, hit_hashes["0"] - shared
                by_sequence_hash: dict[str, tuple[str, int]] = {}
                for record in sequences.values():
                    by_sequence_hash.setdefault(
                        record["sequence_sha256"], (record["label"], record["gc_count"])
                    )
                group_values: dict[int, list[tuple[int, int]]] = {}
                for sequence_hash, (label, gc_count) in by_sequence_hash.items():
                    if sequence_hash in shared:
                        continue
                    group_values.setdefault(gc_count, []).append((
                        int(label == "1"), int(sequence_hash in hit_hashes[label])
                    ))
                observed_residual = 0.0
                variance_total = 0.0
                group_parameters: list[tuple[int, int, int, float]] = []
                for group in group_values.values():
                    n_total = len(group)
                    n_positive = sum(is_positive for is_positive, _ in group)
                    n_hit = sum(is_hit for _, is_hit in group)
                    positive_hits_in_group = sum(is_hit for is_positive, is_hit in group if is_positive)
                    expected = n_positive * n_hit / n_total
                    variance = (
                        n_positive * (n_total - n_positive) * n_hit * (n_total - n_hit)
                        / (n_total * n_total * (n_total - 1))
                        if n_total > 1 else 0.0
                    )
                    observed_residual += positive_hits_in_group - expected
                    variance_total += variance
                    if variance > 0:
                        group_parameters.append((n_hit, n_total - n_hit, n_positive, expected))
                observed_z = observed_residual / np.sqrt(variance_total) if variance_total > 0 else 0.0
                seed_code = sum(ord(char) for char in motif + strand) + int(threshold * 1_000_000)
                rng = np.random.default_rng(20261003 + seed_code)
                simulated_residuals = np.zeros(permutations, dtype=float)
                for n_hit, n_nonhit, n_positive, expected in group_parameters:
                    simulated_residuals += rng.hypergeometric(
                        n_hit, n_nonhit, n_positive, size=permutations
                    ) - expected
                perm_p = (
                    (int(np.count_nonzero(np.abs(simulated_residuals) >= abs(observed_residual))) + 1)
                    / (permutations + 1)
                    if variance_total > 0 else 1.0
                )
                permutation_p_values.append(perm_p)
                threshold_rows.append({
                    "species": "escherichia_coli",
                    "reference_motif": motif,
                    "site_p_cutoff": threshold,
                    "input_start_window_1based": f"{start_min}-{start_max}",
                    "strand": strand,
                    "positive_unique_sequences_hit": len(pos_hit),
                    "positive_unique_sequences": len(pos_total),
                    "control_unique_sequences_hit": len(ctl_hit),
                    "control_unique_sequences": len(ctl_total),
                    "fisher_p_value": fisher_exact_two_sided(
                        len(pos_hit), len(pos_total) - len(pos_hit),
                        len(ctl_hit), len(ctl_total) - len(ctl_hit),
                    ),
                    "fisher_bh_q_4_tests": None,
                    "gc_stratified_permutation_count": permutations,
                    "gc_stratified_z_score": observed_z,
                    "gc_stratified_permutation_p_value": perm_p,
                    "gc_stratified_bh_q_12_tests": None,
                })
        adjusted = benjamini_hochberg([row["fisher_p_value"] for row in threshold_rows])
        for row, q_value in zip(threshold_rows, adjusted, strict=True):
            row["fisher_bh_q_4_tests"] = q_value
            rows.append(row)
    permutation_q = benjamini_hochberg(permutation_p_values)
    for row, q_value in zip(rows, permutation_q, strict=True):
        row["gc_stratified_bh_q_12_tests"] = q_value
    return rows


def run(
    *,
    reference_meme: Path,
    fimo_manifest_path: Path,
    fimo: Path,
    work_dir: Path,
    output_dir: Path,
    thresholds: tuple[float, ...] = (0.05, 0.01, 0.001),
) -> dict[str, Any]:
    if not thresholds or tuple(sorted(set(thresholds), reverse=True)) != thresholds:
        raise ValueError("thresholds must be unique and ordered from largest to smallest")
    if any(not 0 < value <= 1 for value in thresholds):
        raise ValueError("site p-value thresholds must be in (0, 1]")
    source = json.loads(fimo_manifest_path.read_text(encoding="utf-8"))
    if source.get("status") != "complete" or set(source.get("results", {})) != set(SPECIES):
        raise ValueError("a complete six-species FIMO input manifest is required")
    version_result = subprocess.run([str(fimo), "--version"], capture_output=True, text=True, check=True, timeout=15)
    version = (version_result.stdout or version_result.stderr).strip().splitlines()[0]
    if version != "5.5.9":
        raise ValueError(f"expected FIMO 5.5.9, got {version}")

    all_rows: list[dict[str, Any]] = []
    position_rows: list[dict[str, Any]] = []
    ecoli_window_rows: list[dict[str, Any]] = []
    species_runs: dict[str, Any] = {}
    work_dir.mkdir(parents=True, exist_ok=True)
    for slug in SPECIES:
        source_record = source["results"][slug]
        fasta = Path(source_record["combined_fasta"])
        background = Path(source_record["background_file"])
        for path, key in ((fasta, "combined_fasta_sha256"), (background, "background_file_sha256")):
            if not path.is_file() or sha256_file(path) != source_record[key]:
                raise ValueError(f"input file missing or changed: {path}")
        tsv_path = work_dir / f"{slug}.p{thresholds[0]:g}.tsv"
        log_path = work_dir / f"{slug}.stderr.log"
        command = [
            str(fimo), "--text", "--skip-matched-sequence", "--thresh", f"{thresholds[0]:g}",
            "--bfile", str(background), str(reference_meme), str(fasta),
        ]
        started = time.monotonic()
        with tsv_path.open("w", encoding="utf-8") as output:
            process = subprocess.run(command, stdout=output, stderr=subprocess.PIPE, text=True, check=False)
        log_path.write_text(process.stderr, encoding="utf-8")
        if process.returncode != 0:
            raise RuntimeError(f"FIMO failed for {slug}; see {log_path}")
        sequences = _read_fasta_with_gc(fasta)
        hits = _read_raw_hits(tsv_path)
        if len(sequences) != source_record["positive_count"] + source_record["control_count"]:
            raise ValueError(f"unexpected sequence count for {slug}")
        rows = summarize_hits(slug, sequences, hits, thresholds)
        all_rows.extend(rows)
        position_rows.extend(summarize_position_bins(slug, sequences, hits, thresholds))
        if slug == "escherichia_coli":
            ecoli_window_rows = summarize_ecoli_expected_windows(sequences, hits, thresholds)
        species_runs[slug] = {
            "combined_fasta_sha256": sha256_file(fasta),
            "background_sha256": sha256_file(background),
            "raw_pvalue_tsv_sha256": sha256_file(tsv_path),
            "raw_pvalue_tsv_bytes": tsv_path.stat().st_size,
            "raw_site_count_at_max_threshold": len(hits),
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "command": command,
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "pvalue_sensitivity.tsv"
    _write_tsv(summary_path, all_rows)
    position_path = output_dir / "pvalue_position_distribution.tsv"
    _write_tsv(position_path, position_rows)
    window_path = output_dir / "ecoli_sigma70_expected_window_tests.tsv"
    _write_tsv(window_path, ecoli_window_rows)
    manifest = {
        "schema_version": "reference-pwm-sensitivity-1.0",
        "status": "complete",
        "tool": "MEME Suite FIMO",
        "version": version,
        "reference_meme_sha256": sha256_file(reference_meme),
        "source_fimo_manifest_sha256": sha256_file(fimo_manifest_path),
        "parameters": {
            "max_site_p_value": thresholds[0],
            "reported_site_p_cutoffs": list(thresholds),
            "background": "same species-specific zero-order model estimated from label=0 controls",
            "strands": "both strands",
            "sequence_unit": "unique exact sequences for coverage; sequences shared across labels excluded from Fisher tests",
            "multiple_testing": "Benjamini-Hochberg over 24 species-by-reference Fisher tests separately at each site p cutoff",
            "position": "Input FASTA coordinates; expected windows test sigma70-like -10 starts 49-54 and -35 starts 24-29 under the stated [-60,+20] orientation",
            "expected_window_null": "10,000 Monte Carlo label permutations within exact GC-count strata; Benjamini-Hochberg over all 12 motif/strand/cutoff tests",
        },
        "species_runs": species_runs,
        "summary_sha256": sha256_file(summary_path),
        "position_distribution_sha256": sha256_file(position_path),
        "ecoli_expected_window_tests_sha256": sha256_file(window_path),
        "summary_rows": len(all_rows),
        "raw_sites_at_max_cutoff": sum(run["raw_site_count_at_max_threshold"] for run in species_runs.values()),
    }
    manifest_path = output_dir / "pvalue_sensitivity_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-meme", type=Path, default=Path("baselines/known_promoter_elements.meme"))
    parser.add_argument("--fimo-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo/fimo_run_manifest.json"))
    parser.add_argument("--fimo", type=Path, default=Path("tmp/meme-suite-5.5.9/bin/fimo"))
    parser.add_argument("--work-dir", type=Path, default=Path("tmp/fimo_sensitivity_reference_20261003"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/reference_elements"))
    args = parser.parse_args()
    manifest = run(
        reference_meme=args.reference_meme,
        fimo_manifest_path=args.fimo_manifest,
        fimo=args.fimo,
        work_dir=args.work_dir,
        output_dir=args.output_dir,
    )
    print(f"rows={manifest['summary_rows']}; raw_sites={manifest['raw_sites_at_max_cutoff']}")


if __name__ == "__main__":
    main()
