"""Cross-scan RegulonDB sigma-group MEME motifs and summarize group enrichment."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from datetime import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import subprocess
from typing import Any

from experiments.render_meme_logos import parse_meme_xml
from experiments.summarize_tjupan_fimo import benjamini_hochberg, fisher_exact_two_sided
from experiments.validate_regulondb_cross_dataset import exact_mcnemar_p
from preprocessing.provenance import sha256_file


ROOT = Path(__file__).resolve().parents[1]
GROUPS = ("Sigma70", "Sigma24", "Sigma32", "Sigma38", "Sigma28", "Sigma54")
THRESHOLDS = (0.05, 0.01, 0.001)
DEFAULT_INPUT_DIR = ROOT / "tmp/motif_inputs_20260930"
DEFAULT_INPUT_MANIFEST = DEFAULT_INPUT_DIR / "motif_input_manifest.json"
DEFAULT_MATRIX_DIR = ROOT / "results/motif/m3_20260930/matrices"
DEFAULT_MOTIF_RUN_DIR = ROOT / "tmp/meme_runs/regulondb_tss_20260930"
DEFAULT_COMPARISON_MANIFEST = ROOT / "results/motif/m3_20260930/comparison_manifest.json"
DEFAULT_FIMO = ROOT / "tmp/meme-suite-5.5.9/bin/fimo"
DEFAULT_WORK_DIR = ROOT / "tmp/regulondb_sigma_fimo_20261004"
DEFAULT_RESULT_DIR = ROOT / "results/motif/regulondb_sigma_fimo_20261004"
DNA = set("ACGT")


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if header is None:
            return
        sequence = "".join(chunks).upper()
        if len(sequence) != 81 or set(sequence) - DNA:
            raise ValueError(f"expected an 81-nt A/C/G/T sequence in {path}: {header}")
        records.append((header, sequence))

    for line_number, raw in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header, chunks = line[1:].split()[0], []
        elif header is None:
            raise ValueError(f"sequence before FASTA header at {path}:{line_number}")
        else:
            chunks.append(line)
    finish()
    if len({identifier for identifier, _ in records}) != len(records):
        raise ValueError(f"duplicate FASTA identifiers in {path}")
    return records


def _portable_path(path: Path) -> str:
    """Return repository-relative paths in committed manifests when possible."""
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _read_group_inputs(
    input_dir: Path, input_manifest_path: Path
) -> tuple[list[tuple[str, str]], dict[str, dict[str, str]], dict[str, Any]]:
    manifest = json.loads(input_manifest_path.read_text(encoding="utf-8"))
    if manifest.get("source_dataset") != "regulondb" or manifest.get("split") != "discovery":
        raise ValueError("input manifest is not the expected RegulonDB discovery set")
    fasta_records: list[tuple[str, str]] = []
    meta: dict[str, dict[str, str]] = {}
    positive_sets: dict[str, set[str]] = {}
    raw_sha256: dict[str, str] = {}
    counts: dict[str, dict[str, int]] = {}

    for group in GROUPS:
        record = manifest.get("outputs", {}).get(group)
        if not record:
            raise ValueError(f"missing {group} in input manifest")
        positive_path = input_dir / record["primary_fasta"]
        control_path = input_dir / record["n1_fasta"]
        if sha256_file(positive_path) != record["primary_sha256"] or sha256_file(control_path) != record["n1_sha256"]:
            raise ValueError(f"{group} FASTA hash differs from input manifest")
        positives, controls = _read_fasta(positive_path), _read_fasta(control_path)
        if len(positives) != record["sequence_count"] or len(controls) != record["sequence_count"]:
            raise ValueError(f"{group} primary/control counts do not match input manifest")
        positive_ids = {identifier for identifier, _ in positives}
        control_ids = {identifier.removesuffix("__dinucleotide1") for identifier, _ in controls}
        if positive_ids != control_ids:
            raise ValueError(f"{group} N1 controls do not map one-to-one to primary IDs")
        positive_sequences = [sequence for _, sequence in positives]
        if len(set(positive_sequences)) != len(positive_sequences):
            raise ValueError(f"{group} contains repeated positive sequences; deduplication policy required")
        positive_sets[group] = set(positive_sequences)
        counts[group] = {"positive_records": len(positives), "control_records": len(controls)}
        for label, records in (("positive", positives), ("control", controls)):
            for identifier, sequence in records:
                normalized_id = identifier.removesuffix("__dinucleotide1")
                sequence_name = f"{group}|{label}|{normalized_id}"
                if sequence_name in meta:
                    raise ValueError(f"combined FIMO sequence name collision: {sequence_name}")
                meta[sequence_name] = {"group": group, "label": label, "pair_id": normalized_id}
                fasta_records.append((sequence_name, sequence))
        raw_sha256[f"{group}_positive_fasta"] = record["primary_sha256"]
        raw_sha256[f"{group}_N1_fasta"] = record["n1_sha256"]

    overlap_audit: list[dict[str, Any]] = []
    for group_a, group_b in itertools.combinations(GROUPS, 2):
        a, b = positive_sets[group_a], positive_sets[group_b]
        b_rc = {_reverse_complement(sequence) for sequence in b}
        overlap_audit.append({
            "group_a": group_a,
            "group_b": group_b,
            "same_orientation_exact_overlap": len(a & b),
            "reverse_complement_exact_overlap": len(a & b_rc),
        })
        if a & b or a & b_rc:
            raise ValueError(f"exact or reverse-complement sequence overlap across {group_a}/{group_b}")
    return fasta_records, meta, {
        "groups": counts,
        "overlap_audit": overlap_audit,
        "input_sha256": raw_sha256,
        "input_manifest_sha256": sha256_file(input_manifest_path),
    }


def _matrix_metadata(matrix_path: Path, xml_path: Path) -> dict[str, dict[str, Any]]:
    xml_motifs = {motif["id"]: motif for motif in parse_meme_xml(xml_path)}
    matrix_ids: dict[str, str] = {}
    for line in matrix_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MOTIF "):
            parts = line.split(maxsplit=2)
            if len(parts) != 3:
                raise ValueError(f"expected MEME matrix ID and alternate ID in {matrix_path}: {line}")
            matrix_ids[parts[1]] = parts[2]
    if len(matrix_ids) != 10 or set(matrix_ids.values()) != set(xml_motifs):
        raise ValueError(f"matrix/XML motif identifiers differ in {matrix_path}")
    return {
        alt_id: {
            "id": alt_id,
            "matrix_id": matrix_id,
            "consensus": xml_motifs[alt_id]["consensus"],
            "width": xml_motifs[alt_id]["width"],
            "e_value": xml_motifs[alt_id]["e_value"],
            "ordinal": int(alt_id.rsplit("-", 1)[-1]),
        }
        for matrix_id, alt_id in matrix_ids.items()
    }


def _parse_fimo(path: Path, matrix_metadata: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        data_lines = [line for line in handle if line.strip() and not line.startswith("#")]
    if not data_lines:
        return []
    reader = csv.DictReader(data_lines, delimiter="\t")
    required = {"motif_id", "sequence_name", "start", "stop", "strand", "score", "p-value", "q-value"}
    if not reader.fieldnames or not required <= set(reader.fieldnames):
        raise ValueError(f"unexpected FIMO TSV columns in {path}: {reader.fieldnames}")
    by_alias = {item["id"]: item for item in matrix_metadata.values()}
    by_matrix_id = {item["matrix_id"]: item for item in matrix_metadata.values()}
    hits: list[dict[str, Any]] = []
    for row in reader:
        raw_id = row["motif_id"]
        motif = by_matrix_id.get(raw_id) or by_alias.get(raw_id)
        if motif is None:
            raise ValueError(f"unmapped FIMO motif ID {raw_id!r} in {path}")
        start, stop = int(row["start"]), int(row["stop"])
        p_value, q_value = float(row["p-value"]), float(row["q-value"])
        if start < 1 or stop < start or stop > 81 or not 0 <= p_value <= 1 or not 0 <= q_value <= 1:
            raise ValueError(f"invalid FIMO site row: {row}")
        if row["strand"] not in {"+", "-"} or not math.isfinite(float(row["score"])):
            raise ValueError(f"invalid FIMO site row: {row}")
        hits.append({
            "motif_id": motif["id"],
            "start": start,
            "stop": stop,
            "strand": row["strand"],
            "score": float(row["score"]),
            "p_value": p_value,
            "q_value": q_value,
            "sequence_name": row["sequence_name"],
        })
    return hits


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty result table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _nearest_pair_gap(first: dict[str, Any], second: dict[str, Any]) -> tuple[int, str]:
    if first["stop"] < second["start"]:
        return second["start"] - first["stop"] - 1, "a_before_b"
    if second["stop"] < first["start"]:
        return first["start"] - second["stop"] - 1, "b_before_a"
    overlap = min(first["stop"], second["stop"]) - max(first["start"], second["start"]) + 1
    order = "a_before_b" if first["start"] <= second["start"] else "b_before_a"
    return -overlap, order


def _paired_discordance(
    positive_hits: set[str], control_hits: set[str], sequence_meta: dict[str, dict[str, str]]
) -> tuple[int, int]:
    positive_pairs = {sequence_meta[name]["pair_id"] for name in positive_hits}
    control_pairs = {sequence_meta[name]["pair_id"] for name in control_hits}
    return len(positive_pairs - control_pairs), len(control_pairs - positive_pairs)


def _pair_spacing_rows(
    group: str,
    motif_metadata: dict[str, dict[str, Any]],
    hits_by_motif_sequence: dict[str, dict[str, list[dict[str, Any]]]],
    positive_ids: set[str],
    sequence_meta: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    motifs = sorted(motif_metadata, key=lambda motif: motif_metadata[motif]["ordinal"])
    rows: list[dict[str, Any]] = []

    def nearest_for_sequence(motif_a: str, motif_b: str, sequence_name: str, threshold: float) -> tuple[int, str] | None:
        sites_a = [hit for hit in hits_by_motif_sequence.get(motif_a, {}).get(sequence_name, []) if hit["p_value"] <= threshold]
        sites_b = [hit for hit in hits_by_motif_sequence.get(motif_b, {}).get(sequence_name, []) if hit["p_value"] <= threshold]
        candidates = [_nearest_pair_gap(a, b) for a in sites_a for b in sites_b]
        return min(candidates, key=lambda item: (abs(item[0]), item[0], item[1])) if candidates else None

    for threshold in THRESHOLDS:
        for motif_a, motif_b in itertools.combinations(motifs, 2):
            positive_records: list[tuple[int, str]] = []
            control_records: list[tuple[int, str]] = []
            positive_pair_ids: set[str] = set()
            control_pair_ids: set[str] = set()
            for sequence_name in positive_ids:
                paired_control = f"{group}|control|{sequence_meta[sequence_name]['pair_id']}"
                positive_gap = nearest_for_sequence(motif_a, motif_b, sequence_name, threshold)
                control_gap = nearest_for_sequence(motif_a, motif_b, paired_control, threshold)
                if positive_gap is not None:
                    positive_records.append(positive_gap)
                    positive_pair_ids.add(sequence_meta[sequence_name]["pair_id"])
                if control_gap is not None:
                    control_records.append(control_gap)
                    control_pair_ids.add(sequence_meta[paired_control]["pair_id"])
            positive_only = len(positive_pair_ids - control_pair_ids)
            control_only = len(control_pair_ids - positive_pair_ids)
            positive_gaps = [gap for gap, _ in positive_records]
            control_gaps = [gap for gap, _ in control_records]
            rows.append({
                "sigma_group": group,
                "motif_a": motif_a,
                "motif_a_consensus": motif_metadata[motif_a]["consensus"],
                "motif_b": motif_b,
                "motif_b_consensus": motif_metadata[motif_b]["consensus"],
                "site_p_threshold": threshold,
                "positive_unique_sequences": len(positive_ids),
                "positive_sequences_with_both_motifs": len(positive_records),
                "positive_cohit_rate": len(positive_records) / len(positive_ids),
                "N1_control_sequences": len(positive_ids),
                "N1_sequences_with_both_motifs": len(control_records),
                "N1_cohit_rate": len(control_records) / len(positive_ids),
                "cohit_rate_difference": len(positive_records) / len(positive_ids) - len(control_records) / len(positive_ids),
                "paired_positive_only": positive_only,
                "paired_N1_only": control_only,
                "paired_mcnemar_p": exact_mcnemar_p(positive_only, control_only),
                "paired_mcnemar_bh_q": None,
                "positive_median_nearest_signed_gap_bp": statistics.median(positive_gaps) if positive_gaps else "",
                "N1_median_nearest_signed_gap_bp": statistics.median(control_gaps) if control_gaps else "",
                "positive_overlap_sequences": sum(gap < 0 for gap in positive_gaps),
                "positive_motif_a_left_sequences": sum(order == "a_before_b" for _, order in positive_records),
                "positive_motif_b_left_sequences": sum(order == "b_before_a" for _, order in positive_records),
            })
    return rows


def run(
    *,
    input_dir: Path = DEFAULT_INPUT_DIR,
    input_manifest_path: Path = DEFAULT_INPUT_MANIFEST,
    matrix_dir: Path = DEFAULT_MATRIX_DIR,
    meme_run_dir: Path = DEFAULT_MOTIF_RUN_DIR,
    comparison_manifest_path: Path = DEFAULT_COMPARISON_MANIFEST,
    fimo: Path = DEFAULT_FIMO,
    work_dir: Path = DEFAULT_WORK_DIR,
    result_dir: Path = DEFAULT_RESULT_DIR,
) -> dict[str, Any]:
    if not fimo.is_file():
        raise FileNotFoundError(f"FIMO executable not found: {fimo}")
    if result_dir.exists() and any(result_dir.iterdir()):
        raise FileExistsError(f"result directory is not empty: {result_dir}")
    if work_dir.exists() and any(work_dir.iterdir()):
        raise FileExistsError(f"work directory is not empty: {work_dir}")
    input_records, sequence_meta, input_audit = _read_group_inputs(input_dir, input_manifest_path)
    comparison_manifest = json.loads(comparison_manifest_path.read_text(encoding="utf-8"))
    input_records_path = work_dir / "all_sigma_primary_and_N1.fasta"
    background_path = work_dir / "pooled_N1_background.txt"
    work_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    _write_fasta(input_records, input_records_path)

    background_counts = Counter(base for identifier, sequence in input_records if sequence_meta[identifier]["label"] == "control" for base in sequence)
    total_bases = sum(background_counts.values())
    background_frequencies = {base: background_counts[base] / total_bases for base in "ACGT"}
    background_path.write_text(
        "".join(f"{base} {background_frequencies[base]:.12g}\n" for base in "ACGT"), encoding="ascii"
    )
    backgrounds = {"counts": dict(background_counts), "frequencies": background_frequencies}

    ids_by_group_class: dict[tuple[str, str], set[str]] = defaultdict(set)
    for name, item in sequence_meta.items():
        ids_by_group_class[(item["group"], item["label"])].add(name)

    matrix_metadata: dict[str, dict[str, dict[str, Any]]] = {}
    matrix_hashes: dict[str, str] = {}
    input_audit["input_sha256"]["combined_fasta"] = sha256_file(input_records_path)
    input_audit["input_sha256"]["pooled_N1_background"] = sha256_file(background_path)
    for group in GROUPS:
        matrix_path = matrix_dir / f"{group}.meme"
        xml_path = meme_run_dir / group / "meme.xml"
        matrix_key = f"matrices/{group}.meme"
        expected_matrix_hash = comparison_manifest["output_sha256"].get(matrix_key)
        actual_matrix_hash = sha256_file(matrix_path)
        if expected_matrix_hash != actual_matrix_hash:
            raise ValueError(f"{group} motif matrix hash differs from comparison manifest")
        matrix_metadata[group] = _matrix_metadata(matrix_path, xml_path)
        matrix_hashes[matrix_key] = actual_matrix_hash

    coverage_index: dict[tuple[str, str, str, str, float], set[str]] = defaultdict(set)
    site_count_index: Counter[tuple[str, str, str, str, float]] = Counter()
    q_coverage_index: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    site_rows_by_motif_sequence: dict[str, dict[str, dict[str, list[dict[str, Any]]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    run_records: dict[str, Any] = {}
    fimo_version = subprocess.run([str(fimo), "--version"], capture_output=True, text=True, check=True).stdout.strip()

    for source_group in GROUPS:
        output_dir = work_dir / f"fimo_{source_group}"
        command = [
            str(fimo), "--oc", str(output_dir), "--bfile", str(background_path),
            "--thresh", "0.05", "--motif-pseudo", "0.1", "--max-stored-scores", "15000000",
            str(matrix_dir / f"{source_group}.meme"), str(input_records_path),
        ]
        process = subprocess.run(command, capture_output=True, text=True, check=False, timeout=900)
        log_path = work_dir / f"fimo_{source_group}.stderr.log"
        log_path.write_text(process.stderr, encoding="utf-8")
        if process.returncode != 0:
            raise RuntimeError(f"FIMO failed for {source_group}; see {log_path}")
        hits = _parse_fimo(output_dir / "fimo.tsv", matrix_metadata[source_group])
        run_records[source_group] = {
            "return_code": process.returncode,
            "command": [_portable_path(Path(part)) if part.startswith("/") else part for part in command],
            "site_rows_p_le_0.05": len(hits),
            "fimo_tsv_sha256": sha256_file(output_dir / "fimo.tsv"),
            "stderr_sha256": sha256_file(log_path),
        }
        valid_names = set(sequence_meta)
        for hit in hits:
            name = hit["sequence_name"]
            if name not in valid_names:
                raise ValueError(f"FIMO returned unknown sequence name: {name}")
            target = sequence_meta[name]
            motif_id = hit["motif_id"]
            site_rows_by_motif_sequence[source_group][motif_id][name].append(hit)
            for threshold in THRESHOLDS:
                if hit["p_value"] <= threshold:
                    key = (source_group, motif_id, target["group"], target["label"], threshold)
                    coverage_index[key].add(name)
                    site_count_index[key] += 1
            if hit["q_value"] <= 0.05:
                q_coverage_index[(source_group, motif_id, target["group"], target["label"])].add(name)

    coverage_rows: list[dict[str, Any]] = []
    specificity_rows: list[dict[str, Any]] = []
    enrichment_test_rows: list[dict[str, Any]] = []
    specificity_test_rows: list[dict[str, Any]] = []
    for source_group in GROUPS:
        for motif_id, motif in sorted(matrix_metadata[source_group].items(), key=lambda item: item[1]["ordinal"]):
            for target_group in GROUPS:
                positive_ids = ids_by_group_class[(target_group, "positive")]
                control_ids = ids_by_group_class[(target_group, "control")]
                row: dict[str, Any] = {
                    "motif_source_sigma": source_group,
                    "motif_id": motif_id,
                    "consensus": motif["consensus"],
                    "discovery_e_value": motif["e_value"],
                    "width_bp": motif["width"],
                    "target_sigma": target_group,
                    "positive_sequences": len(positive_ids),
                    "N1_control_sequences": len(control_ids),
                    "positive_sequences_FIMO_site_q_le_0.05": len(q_coverage_index[(source_group, motif_id, target_group, "positive")]),
                    "N1_sequences_FIMO_site_q_le_0.05": len(q_coverage_index[(source_group, motif_id, target_group, "control")]),
                }
                for threshold in THRESHOLDS:
                    pos_key = (source_group, motif_id, target_group, "positive", threshold)
                    nul_key = (source_group, motif_id, target_group, "control", threshold)
                    pos_hits, nul_hits = coverage_index[pos_key], coverage_index[nul_key]
                    pos_only, nul_only = _paired_discordance(pos_hits, nul_hits, sequence_meta)
                    mcnemar_p = exact_mcnemar_p(pos_only, nul_only)
                    row.update({
                        f"positive_sites_p_le_{threshold:g}": site_count_index[pos_key],
                        f"positive_hit_sequences_p_le_{threshold:g}": len(pos_hits),
                        f"positive_hit_rate_p_le_{threshold:g}": len(pos_hits) / len(positive_ids),
                        f"N1_sites_p_le_{threshold:g}": site_count_index[nul_key],
                        f"N1_hit_sequences_p_le_{threshold:g}": len(nul_hits),
                        f"N1_hit_rate_p_le_{threshold:g}": len(nul_hits) / len(control_ids),
                        f"paired_mcnemar_p_p_le_{threshold:g}": mcnemar_p,
                        f"paired_positive_only_p_le_{threshold:g}": pos_only,
                        f"paired_N1_only_p_le_{threshold:g}": nul_only,
                    })
                    enrichment_test_rows.append({
                        "row_index": len(coverage_rows),
                        "source_group": source_group,
                        "motif_id": motif_id,
                        "target_group": target_group,
                        "threshold": threshold,
                        "p_value": mcnemar_p,
                    })
                row["positive_sequences_FIMO_site_q_le_0.05_rate"] = row["positive_sequences_FIMO_site_q_le_0.05"] / len(positive_ids)
                row["N1_sequences_FIMO_site_q_le_0.05_rate"] = row["N1_sequences_FIMO_site_q_le_0.05"] / len(control_ids)
                coverage_rows.append(row)

                for threshold in THRESHOLDS:
                    target_hits = coverage_index[(source_group, motif_id, target_group, "positive", threshold)]
                    other_ids = set().union(*(ids_by_group_class[(other_group, "positive")] for other_group in GROUPS if other_group != target_group))
                    other_hits = set().union(*(coverage_index[(source_group, motif_id, other_group, "positive", threshold)] for other_group in GROUPS if other_group != target_group))
                    specificity_p = fisher_exact_two_sided(
                        len(target_hits), len(positive_ids) - len(target_hits),
                        len(other_hits), len(other_ids) - len(other_hits),
                    )
                    specificity_rows.append({
                        "motif_source_sigma": source_group,
                        "motif_id": motif_id,
                        "consensus": motif["consensus"],
                        "target_sigma": target_group,
                        "site_p_threshold": threshold,
                        "target_hit_sequences": len(target_hits),
                        "target_sequences": len(positive_ids),
                        "target_hit_rate": len(target_hits) / len(positive_ids),
                        "other_sigma_hit_sequences": len(other_hits),
                        "other_sigma_sequences": len(other_ids),
                        "other_sigma_hit_rate": len(other_hits) / len(other_ids),
                        "target_vs_other_fisher_p": specificity_p,
                        "target_vs_other_bh_q": None,
                        "target_is_motif_discovery_group": target_group == source_group,
                    })

    enrichment_q_by_family: dict[float, list[float]] = {}
    for threshold in THRESHOLDS:
        tests = [item for item in enrichment_test_rows if item["threshold"] == threshold]
        enrichment_q_by_family[threshold] = benjamini_hochberg([item["p_value"] for item in tests])
        for item, q_value in zip(tests, enrichment_q_by_family[threshold], strict=True):
            row = coverage_rows[item["row_index"]]
            row[f"paired_mcnemar_bh_q_p_le_{threshold:g}"] = q_value

    for threshold in THRESHOLDS:
        indices = [index for index, row in enumerate(specificity_rows) if row["site_p_threshold"] == threshold]
        adjusted = benjamini_hochberg([specificity_rows[index]["target_vs_other_fisher_p"] for index in indices])
        for index, q_value in zip(indices, adjusted, strict=True):
            specificity_rows[index]["target_vs_other_bh_q"] = q_value

    spacing_rows: list[dict[str, Any]] = []
    for group in GROUPS:
        own_positive_ids = ids_by_group_class[(group, "positive")]
        spacing_rows.extend(_pair_spacing_rows(group, matrix_metadata[group], site_rows_by_motif_sequence[group], own_positive_ids, sequence_meta))
    for threshold in THRESHOLDS:
        indices = [index for index, row in enumerate(spacing_rows) if row["site_p_threshold"] == threshold]
        adjusted = benjamini_hochberg([spacing_rows[index]["paired_mcnemar_p"] for index in indices])
        for index, q_value in zip(indices, adjusted, strict=True):
            spacing_rows[index]["paired_mcnemar_bh_q"] = q_value

    coverage_path = result_dir / "sigma_motif_coverage.tsv"
    specificity_path = result_dir / "sigma_motif_group_specificity.tsv"
    spacing_path = result_dir / "sigma_motif_pair_spacing.tsv"
    _write_tsv(coverage_path, coverage_rows)
    _write_tsv(specificity_path, specificity_rows)
    _write_tsv(spacing_path, spacing_rows)

    outputs = [coverage_path, specificity_path, spacing_path]
    result = {
        "schema_version": "pr01-02-regulondb-sigma-fimo-1.0",
        "status": "complete",
        "purpose": "以 RegulonDB 六个 sigma 组 MEME 矩阵交叉扫描六组启动子，比较组内对照富集、跨组命中频率及组内 motif 间距；属于辅助分析",
        "dataset": "RegulonDB E. coli K-12 TSS-aligned discovery subset",
        "groups": list(GROUPS),
        "input_counts": input_audit["groups"],
        "cross_group_exact_sequence_overlap": input_audit["overlap_audit"],
        "input_sha256": input_audit["input_sha256"],
        "input_manifest_sha256": input_audit["input_manifest_sha256"],
        "comparison_manifest_sha256": sha256_file(comparison_manifest_path),
        "matrix_sha256": matrix_hashes,
        "tool": {"fimo_version": fimo_version, "fimo_executable": _portable_path(fimo)},
        "parameters": {
            "all_target_sequences": "每个 sigma 组的正类与一一匹配 N1 对照；六组在每个矩阵任务中同时扫描",
            "background": "六组 N1 对照 pooled zero-order A/C/G/T frequencies",
            "strands": "both strands",
            "site_output_cutoff": 0.05,
            "reported_site_p_thresholds": list(THRESHOLDS),
            "site_q_threshold": 0.05,
            "motif_pseudocount": 0.1,
            "max_stored_scores": 15000000,
            "paired_control_test": "exact McNemar on positive/N1 sequence-paired hit presence; BH separately within each raw-p threshold across 360 source-motif by target-group tests",
            "group_specificity_test": "two-sided Fisher comparing each target sigma group's positive hit rate against the other five groups; BH separately within each raw-p threshold across 360 source-motif by target-group tests",
            "spacing": "for each discovery group's own 10 motifs, compare paired positive/N1 sequence co-occurrence using exact McNemar and BH separately across 270 group-motif pairs per raw-p threshold; report nearest signed-gap summaries in both classes",
            "selection_limit": "motifs and FIMO target sequences overlap for the source sigma group; group-origin comparisons are exploratory and require holdout replication",
            "position_limit": "RegulonDB inputs are TSS-oriented; site coordinates are 1-81 in input direction. This auxiliary result does not replace six-species TJU Pan analysis.",
        },
        "background": backgrounds,
        "fimo_runs": run_records,
        "rows": {
            "sigma_motif_coverage.tsv": len(coverage_rows),
            "sigma_motif_group_specificity.tsv": len(specificity_rows),
            "sigma_motif_pair_spacing.tsv": len(spacing_rows),
            "pair_spacing_bh_significant_by_threshold": {
                str(threshold): sum(1 for row in spacing_rows if row["site_p_threshold"] == threshold and row["paired_mcnemar_bh_q"] <= 0.05)
                for threshold in THRESHOLDS
            },
        },
        "output_sha256": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in outputs
        },
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    (result_dir / "sigma_fimo_manifest.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return result


def _completion_message(result: dict[str, Any]) -> str:
    summary_rows = sum(value for value in result["rows"].values() if isinstance(value, int))
    site_rows = sum(item["site_rows_p_le_0.05"] for item in result["fimo_runs"].values())
    return f"Completed cross-sigma FIMO: {summary_rows} summary rows, {site_rows} site rows across six scans"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--input-manifest", type=Path, default=DEFAULT_INPUT_MANIFEST)
    parser.add_argument("--matrix-dir", type=Path, default=DEFAULT_MATRIX_DIR)
    parser.add_argument("--meme-run-dir", type=Path, default=DEFAULT_MOTIF_RUN_DIR)
    parser.add_argument("--comparison-manifest", type=Path, default=DEFAULT_COMPARISON_MANIFEST)
    parser.add_argument("--fimo", type=Path, default=DEFAULT_FIMO)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    args = parser.parse_args()
    result = run(
        input_dir=args.input_dir,
        input_manifest_path=args.input_manifest,
        matrix_dir=args.matrix_dir,
        meme_run_dir=args.meme_run_dir,
        comparison_manifest_path=args.comparison_manifest,
        fimo=args.fimo,
        work_dir=args.work_dir,
        result_dir=args.result_dir,
    )
    print(_completion_message(result))


if __name__ == "__main__":
    main()
