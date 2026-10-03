"""Create sequence-free summaries from completed TJU Pan FIMO runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES
from preprocessing.provenance import sha256_file


def _read_fasta(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    current: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if current is None:
            return
        if current in records:
            raise ValueError(f"duplicate FASTA identifier: {current}")
        sequence = "".join(chunks).upper()
        label = current.rsplit("|label=", 1)[-1]
        if label not in {"0", "1"} or len(sequence) != 81 or set(sequence) - set("ACGT"):
            raise ValueError(f"invalid sequence record in {path}: {current}")
        records[current] = {"label": label, "sequence": sequence}

    for raw_line in path.read_text(encoding="ascii").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            current = line[1:]
            chunks = []
        elif current is None:
            raise ValueError(f"sequence before FASTA header in {path}")
        else:
            chunks.append(line)
    finish()
    return records


def _read_hits(path: Path, species: str) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        data_lines = [line for line in handle if line.strip() and not line.startswith("#")]
    if not data_lines:
        return []
    reader = csv.DictReader(data_lines, delimiter="\t")
    required = {"motif_id", "sequence_name", "start", "stop", "strand", "p-value", "q-value"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise ValueError(f"FIMO TSV is missing required columns: {path}: {reader.fieldnames}")
    rows: list[dict[str, Any]] = []
    for row in reader:
        motif_id = row["motif_id"]
        prefix = f"{species}_"
        if not motif_id.startswith(prefix):
            raise ValueError(f"unexpected motif ID in {path}: {motif_id}")
        local_motif_id = motif_id[len(prefix):]
        start, stop = int(row["start"]), int(row["stop"])
        p_value, q_value = float(row["p-value"]), float(row["q-value"])
        if start < 1 or stop < start or stop > 81 or not 0 <= p_value <= 1 or not 0 <= q_value <= 1:
            raise ValueError(f"invalid FIMO hit coordinates/statistics: {row}")
        if row["strand"] not in {"+", "-"}:
            raise ValueError(f"invalid FIMO strand: {row['strand']}")
        rows.append({
            "motif_id": local_motif_id,
            "sequence_name": row["sequence_name"],
            "start": start,
            "stop": stop,
            "strand": row["strand"],
            "p_value": p_value,
            "q_value": q_value,
        })
    return rows


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Return a two-sided Fisher exact p-value for [[a,b],[c,d]]."""
    row1, row2 = a + b, c + d
    col1, total = a + c, a + b + c + d
    lower, upper = max(0, row1 + col1 - total), min(row1, col1)

    def log_probability(x: int) -> float:
        return _log_comb(col1, x) + _log_comb(total - col1, row1 - x) - _log_comb(total, row1)

    observed = log_probability(a)
    probabilities = [log_probability(x) for x in range(lower, upper + 1)]
    selected = [value for value in probabilities if value <= observed + 1e-12]
    maximum = max(selected)
    return min(1.0, math.exp(maximum) * sum(math.exp(value - maximum) for value in selected))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in p_values):
        raise ValueError("p-values must be finite and in [0, 1]")
    count = len(p_values)
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * count / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def summarize(
    run_manifest_path: Path,
    input_manifest_path: Path,
    input_dir: Path,
    meme_summary_path: Path,
    run_root: Path,
    output_dir: Path,
    q_threshold: float = 0.05,
) -> dict[str, Any]:
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    input_manifest = json.loads(input_manifest_path.read_text(encoding="utf-8"))
    if run_manifest.get("status") != "complete" or run_manifest.get("fimo_version") != "5.5.9":
        raise ValueError("all six species must have successful FIMO 5.5.9 runs")
    if run_manifest.get("parameters", {}).get("qvalue_threshold") != q_threshold:
        raise ValueError("summary q-value threshold differs from FIMO run manifest")
    with meme_summary_path.open(encoding="utf-8", newline="") as handle:
        motif_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(motif_rows) != 60:
        raise ValueError(f"expected 60 species-specific MEME motifs, found {len(motif_rows)}")

    output_dir.mkdir(parents=True, exist_ok=True)
    result_rows: list[dict[str, Any]] = []
    position_counts: dict[tuple[str, str, str, int, str], int] = {}
    source_hashes: dict[str, str] = {}
    for slug in SPECIES:
        run_record = run_manifest.get("results", {}).get(slug)
        if not run_record or run_record.get("status") != "success":
            raise ValueError(f"missing successful FIMO result for {slug}")
        if len(run_record.get("motif_pi0_estimates", {})) != 10:
            raise ValueError(f"FIMO per-motif q-value estimation audit is incomplete for {slug}")
        fasta_path = input_dir / f"{slug}.combined.fasta"
        fimo_tsv = run_root / slug / "fimo.tsv"
        if not fasta_path.is_file() or sha256_file(fasta_path) != run_record.get("combined_fasta_sha256"):
            raise ValueError(f"combined FIMO input missing or hash mismatch: {fasta_path}")
        if not fimo_tsv.is_file() or sha256_file(fimo_tsv) != run_record.get("fimo_tsv_sha256"):
            raise ValueError(f"FIMO TSV missing or hash mismatch: {fimo_tsv}")
        sequences = _read_fasta(fasta_path)
        positive_records = {name: item["sequence"] for name, item in sequences.items() if item["label"] == "1"}
        control_records = {name: item["sequence"] for name, item in sequences.items() if item["label"] == "0"}
        positive_unique = set(positive_records.values())
        control_unique = set(control_records.values())
        shared_sequences = positive_unique & control_unique
        if len(positive_records) != run_record["positive_count"] or len(control_records) != run_record["control_count"]:
            raise ValueError(f"input record counts differ from run manifest for {slug}")
        hits = _read_hits(fimo_tsv, slug)
        for hit in hits:
            if hit["sequence_name"] not in sequences:
                raise ValueError(f"FIMO hit references unknown sequence: {hit['sequence_name']}")

        expected_motifs = {row["motif_id"]: row for row in motif_rows if row["species"] == slug}
        if len(expected_motifs) != 10:
            raise ValueError(f"expected 10 MEME motifs for {slug}, found {len(expected_motifs)}")
        hit_names: dict[str, set[str]] = {motif_id: set() for motif_id in expected_motifs}
        site_counts: dict[str, int] = {motif_id: 0 for motif_id in expected_motifs}
        positive_site_counts: dict[str, int] = {motif_id: 0 for motif_id in expected_motifs}
        control_site_counts: dict[str, int] = {motif_id: 0 for motif_id in expected_motifs}
        positive_hit_sequences: dict[str, set[str]] = {motif_id: set() for motif_id in expected_motifs}
        control_hit_sequences: dict[str, set[str]] = {motif_id: set() for motif_id in expected_motifs}

        for hit in hits:
            motif_id = hit["motif_id"]
            if motif_id not in expected_motifs:
                raise ValueError(f"FIMO output contains unexpected motif {motif_id} for {slug}")
            label = sequences[hit["sequence_name"]]["label"]
            sequence = sequences[hit["sequence_name"]]["sequence"]
            site_counts[motif_id] += 1
            hit_names[motif_id].add(hit["sequence_name"])
            (positive_site_counts if label == "1" else control_site_counts)[motif_id] += 1
            (positive_hit_sequences if label == "1" else control_hit_sequences)[motif_id].add(sequence)
            key = (slug, motif_id, label, hit["start"], hit["strand"])
            position_counts[key] = position_counts.get(key, 0) + 1

        for motif_id, motif in expected_motifs.items():
            positive_hits = positive_hit_sequences[motif_id]
            control_hits = control_hit_sequences[motif_id]
            positive_test_hits = positive_hits - shared_sequences
            control_test_hits = control_hits - shared_sequences
            positive_test_total = positive_unique - shared_sequences
            control_test_total = control_unique - shared_sequences
            a = len(positive_test_hits)
            b = len(positive_test_total) - a
            c = len(control_test_hits)
            d = len(control_test_total) - c
            p_value = fisher_exact_two_sided(a, b, c, d)
            positive_coverage = len(positive_hits) / len(positive_unique)
            control_coverage = len(control_hits) / len(control_unique)
            result_rows.append({
                "species": slug,
                "motif_id": motif_id,
                "consensus": motif["consensus"],
                "width_bp": int(motif["width_bp"]),
                "meme_e_value": float(motif["e_value"]),
                "fimo_q_threshold": q_threshold,
                "fimo_significant_site_count": site_counts[motif_id],
                "fimo_hit_record_count": len(hit_names[motif_id]),
                "positive_hit_site_count": positive_site_counts[motif_id],
                "control_hit_site_count": control_site_counts[motif_id],
                "positive_hit_unique_sequences": len(positive_hits),
                "positive_unique_sequences": len(positive_unique),
                "positive_unique_sequence_coverage": positive_coverage,
                "control_hit_unique_sequences": len(control_hits),
                "control_unique_sequences": len(control_unique),
                "control_unique_sequence_coverage": control_coverage,
                "coverage_difference_positive_minus_control": positive_coverage - control_coverage,
                "cross_label_overlap_sequences_excluded_from_fisher": len(shared_sequences),
                "fisher_two_sided_p_value": p_value,
                "fisher_bh_q_value_60_tests": None,
            })

        source_hashes[f"{slug}/combined_fasta"] = run_record["combined_fasta_sha256"]
        source_hashes[f"{slug}/control_background"] = run_record["background_file_sha256"]
        source_hashes[f"{slug}/motif_matrix"] = run_record["matrix_sha256"]
        source_hashes[f"{slug}/fimo_tsv"] = run_record["fimo_tsv_sha256"]

    if len(result_rows) != 60:
        raise ValueError(f"expected 60 species/motif summaries, found {len(result_rows)}")
    adjusted = benjamini_hochberg([row["fisher_two_sided_p_value"] for row in result_rows])
    for row, q_value in zip(result_rows, adjusted, strict=True):
        row["fisher_bh_q_value_60_tests"] = q_value

    summary_path = output_dir / "motif_hit_summary.tsv"
    summary_columns = list(result_rows[0])
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result_rows)

    position_path = output_dir / "position_distribution.tsv"
    with position_path.open("w", encoding="utf-8", newline="") as handle:
        columns = ["species", "motif_id", "label", "start_1based_input", "strand", "significant_site_count"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(position_counts):
            slug, motif_id, label, start, strand = key
            writer.writerow({
                "species": slug,
                "motif_id": motif_id,
                "label": label,
                "start_1based_input": start,
                "strand": strand,
                "significant_site_count": position_counts[key],
            })

    outputs = [summary_path, position_path]
    manifest = {
        "purpose": "FIMO 显著位点的去序列汇总、唯一序列覆盖与正负类差异探索；不含样本 ID 或 DNA 序列",
        "status": "complete",
        "dataset": input_manifest.get("source"),
        "fimo_version": run_manifest["fimo_version"],
        "parameters": run_manifest["parameters"],
        "input_manifest_sha256": sha256_file(input_manifest_path),
        "fimo_run_manifest_sha256": sha256_file(run_manifest_path),
        "meme_summary_sha256": sha256_file(meme_summary_path),
        "species": list(SPECIES),
        "motif_count": len(result_rows),
        "fimo_significant_site_count": sum(row["fimo_significant_site_count"] for row in result_rows),
        "motifs_with_fisher_bh_q_le_0_05": sum(row["fisher_bh_q_value_60_tests"] <= 0.05 for row in result_rows),
        "fisher_test_scope": "unique exact sequences; exact sequences shared across positive/control labels excluded from the test; Benjamini-Hochberg correction across all 60 species/motif pairs",
        "fimo_qvalue_scope": "FIMO computes q-values separately for each motif across that motif's tested positions; it does not pool the 10 motifs. MEME Suite 5.5.9 logs one pi0 estimate from a uniformly sampled p-value reservoir per motif.",
        "data_scope_limit": "Motifs were discovered from the same positive sequences later scanned by FIMO; these are exploratory localization/association results, not independent validation. Negative controls are source dataset label=0 records, not necessarily experimentally verified non-promoters.",
        "coordinate_limit": "Coordinates are 1-based in input FASTA; strand/orientation is not validated, so no TSS-relative conversion is made.",
        "source_sha256": source_hashes,
        "output_sha256": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in outputs
        },
    }
    manifest_path = output_dir / "fimo_summary_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo/fimo_run_manifest.json"))
    parser.add_argument("--input-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/input_manifest.json"))
    parser.add_argument("--input-dir", type=Path, default=Path("tmp/tjupan_fimo_inputs_20261003"))
    parser.add_argument("--meme-summary", type=Path, default=Path("results/motif/tjupan_m3_20261003/summary/motif_summary.tsv"))
    parser.add_argument("--run-root", type=Path, default=Path("tmp/fimo_runs/tjupan_promloop_20261003"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo"))
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args()
    manifest = summarize(
        args.run_manifest, args.input_manifest, args.input_dir, args.meme_summary,
        args.run_root, args.output_dir, args.q_threshold,
    )
    print(
        f"FIMO summary: {manifest['fimo_significant_site_count']} significant sites; "
        f"{manifest['motifs_with_fisher_bh_q_le_0_05']} of 60 motif/species pairs have Fisher BH q≤0.05"
    )


if __name__ == "__main__":
    main()
