"""Test de novo motif-pair combinations against matched dinucleotide shuffles."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from itertools import combinations
import math
from pathlib import Path
import random
import statistics
import subprocess
from typing import Any

from experiments.analyze_tjupan_motif_pairs import site_relation
from experiments.prepare_tjupan_motif_inputs import SPECIES
from experiments.summarize_tjupan_fimo import benjamini_hochberg, fisher_exact_two_sided
from experiments.validate_regulondb_cross_dataset import dinucleotide_shuffle, exact_mcnemar_p
from preprocessing.provenance import sha256_file


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if identifier is None:
            return
        sequence = "".join(chunks).upper()
        if len(sequence) != 81 or set(sequence) - set("ACGT"):
            raise ValueError(f"expected an 81-bp A/C/G/T sequence: {identifier}")
        records.append((identifier, sequence))

    for line_number, raw in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            identifier, chunks = line[1:], []
        elif identifier is None:
            raise ValueError(f"sequence before FASTA header at {path}:{line_number}")
        else:
            chunks.append(line)
    finish()
    return records


def unique_sequences(records: list[tuple[str, str]]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for _, sequence in records:
        if sequence not in seen:
            seen.add(sequence)
            result.append(sequence)
    return result


def summarize_pair_presence(
    positive_a: set[int],
    positive_b: set[int],
    null_a: set[int],
    null_b: set[int],
    sample_count: int,
) -> dict[str, Any]:
    all_indexes = positive_a | positive_b | null_a | null_b
    if all_indexes and min(all_indexes) < 0:
        raise ValueError("sequence indexes must be nonnegative")
    if any(index >= sample_count for index in positive_a | positive_b | null_a | null_b):
        raise ValueError("sequence index exceeds paired sample count")
    positive_both = positive_a & positive_b
    null_both = null_a & null_b
    discordant_positive = len(positive_both - null_both)
    discordant_null = len(null_both - positive_both)
    fisher_p = fisher_exact_two_sided(
        len(positive_both), len(positive_a - positive_b),
        len(positive_b - positive_a), sample_count - len(positive_a | positive_b),
    )
    return {
        "positive_a_hits": len(positive_a),
        "positive_b_hits": len(positive_b),
        "positive_pair_cooccurrence": len(positive_both),
        "positive_pair_cooccurrence_rate": len(positive_both) / sample_count if sample_count else 0.0,
        "null_a_hits": len(null_a),
        "null_b_hits": len(null_b),
        "null_pair_cooccurrence": len(null_both),
        "null_pair_cooccurrence_rate": len(null_both) / sample_count if sample_count else 0.0,
        "positive_only_pair_count": discordant_positive,
        "null_only_pair_count": discordant_null,
        "paired_mcnemar_p_value": exact_mcnemar_p(discordant_positive, discordant_null),
        "positive_pair_fisher_independence_p_value": fisher_p,
    }


def _parse_fimo_tsv(path: Path) -> list[dict[str, Any]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter="\t")
    required = {"motif_id", "sequence_name", "start", "stop", "strand", "p-value"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise ValueError(f"invalid FIMO text output: {path}")
    hits = []
    for row in reader:
        hit = {
            "motif_id": row["motif_id"],
            "sequence_name": row["sequence_name"],
            "start": int(row["start"]),
            "stop": int(row["stop"]),
            "strand": row["strand"],
            "p_value": float(row["p-value"]),
        }
        if not 1 <= hit["start"] <= hit["stop"] <= 81 or hit["strand"] not in {"+", "-"}:
            raise ValueError(f"invalid FIMO site: {row}")
        hits.append(hit)
    return hits


def _nearest_site_pair(
    sites_a: list[dict[str, Any]], sites_b: list[dict[str, Any]],
) -> tuple[str, int] | None:
    if not sites_a or not sites_b:
        return None
    choices: list[tuple[int, int, str, int]] = []
    for first in sites_a:
        for second in sites_b:
            relation, gap = site_relation(first, second)
            separation = 0 if relation == "overlapping" else gap
            choices.append((separation, abs(gap), relation, gap))
    _, _, relation, gap = min(choices, key=lambda row: (row[0], row[1], row[2]))
    return relation, gap


def _write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        if not rows:
            raise ValueError(f"columns are required for an empty TSV: {path}")
        columns = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _site_sets(
    hits: list[dict[str, Any]],
    sequence_ids: list[str],
    motif_ids: list[str],
    threshold: float,
) -> tuple[dict[str, dict[str, list[dict[str, Any]]]], dict[str, set[int]]]:
    index_by_id = {identifier: i for i, identifier in enumerate(sequence_ids)}
    sites: dict[str, dict[str, list[dict[str, Any]]]] = {
        identifier: {motif: [] for motif in motif_ids} for identifier in sequence_ids
    }
    present: dict[str, set[int]] = {motif: set() for motif in motif_ids}
    for hit in hits:
        if hit["sequence_name"] not in index_by_id or hit["motif_id"] not in present:
            raise ValueError(f"unexpected sequence or motif in FIMO output: {hit}")
        if hit["p_value"] <= threshold:
            sites[hit["sequence_name"]][hit["motif_id"]].append(hit)
            present[hit["motif_id"]].add(index_by_id[hit["sequence_name"]])
    return sites, present


def run(
    *,
    input_manifest_path: Path,
    fimo_manifest_path: Path,
    matrix_dir: Path,
    fimo: Path,
    work_dir: Path,
    output_dir: Path,
    thresholds: tuple[float, ...] = (0.05, 0.01, 0.001),
    replicates: int = 3,
    seed: int = 20261003,
) -> dict[str, Any]:
    if not thresholds or tuple(sorted(set(thresholds), reverse=True)) != thresholds:
        raise ValueError("thresholds must be unique and ordered from largest to smallest")
    if any(not 0 < threshold <= 1 for threshold in thresholds) or replicates < 1:
        raise ValueError("thresholds must be in (0,1] and replicates must be positive")
    input_manifest = json.loads(input_manifest_path.read_text(encoding="utf-8"))
    fimo_manifest = json.loads(fimo_manifest_path.read_text(encoding="utf-8"))
    if set(input_manifest.get("species", {})) != set(SPECIES) or set(fimo_manifest.get("results", {})) != set(SPECIES):
        raise ValueError("complete six-species input and FIMO manifests are required")
    if fimo_manifest.get("input_manifest_sha256") != sha256_file(input_manifest_path):
        raise ValueError("input manifest does not match the de novo FIMO run")
    version_run = subprocess.run([str(fimo), "--version"], capture_output=True, text=True, check=True, timeout=15)
    version = (version_run.stdout or version_run.stderr).strip().splitlines()[0]
    if version != "5.5.9":
        raise ValueError(f"expected FIMO 5.5.9, got {version}")

    work_dir.mkdir(parents=True, exist_ok=True)
    raw_records: dict[str, list[str]] = {}
    unique_positives: dict[str, list[str]] = {}
    motif_ids_by_species: dict[str, list[str]] = {}
    matrix_hashes: dict[str, str] = {}
    for slug in SPECIES:
        positive_fasta = Path(input_manifest["species"][slug]["positive_samples_csv"]["fasta"])
        if not positive_fasta.is_file():
            raise FileNotFoundError(positive_fasta)
        records = read_fasta(positive_fasta)
        raw_records[slug] = [sequence for _, sequence in records]
        unique_positives[slug] = unique_sequences(records)
        if not unique_positives[slug]:
            raise ValueError(f"no positive sequences for {slug}")
        matrix_path = matrix_dir / f"{slug}.meme"
        matrix_hashes[slug] = sha256_file(matrix_path)
        motif_ids = [line.split()[1] for line in matrix_path.read_text(encoding="utf-8").splitlines() if line.startswith("MOTIF ")]
        if len(motif_ids) != 10:
            raise ValueError(f"expected 10 motifs in {matrix_path}")
        motif_ids_by_species[slug] = motif_ids

    species_rep_runs: dict[str, list[dict[str, Any]]] = {slug: [] for slug in SPECIES}
    pair_rows: list[dict[str, Any]] = []
    spacing_rows: list[dict[str, Any]] = []
    independence_tests: dict[float, list[dict[str, Any]]] = {threshold: [] for threshold in thresholds}
    paired_tests: dict[float, list[dict[str, Any]]] = {threshold: [] for threshold in thresholds}

    for slug in SPECIES:
        positives = unique_positives[slug]
        motif_ids = motif_ids_by_species[slug]
        fasta_path = work_dir / f"{slug}.combination_inputs.fasta"
        background_path = work_dir / f"{slug}.shuffle.background"
        base_counts = Counter(base for sequence in positives for base in sequence)
        base_total = sum(base_counts.values())
        background_path.write_text(
            "".join(f"{base} {base_counts[base] / base_total:.12g}\n" for base in "ACGT"),
            encoding="ascii",
        )
        for replicate in range(1, replicates + 1):
            rng = random.Random(seed + (list(SPECIES).index(slug) * 100) + replicate - 1)
            null_sequences = [dinucleotide_shuffle(sequence, rng) for sequence in positives]
            if len(set(null_sequences)) != len(null_sequences) or set(positives) & set(null_sequences):
                raise ValueError(f"{slug} replicate {replicate}: duplicate or cross-class shuffle sequence")
            sequence_ids: list[str] = []
            fasta_records: list[tuple[str, str]] = []
            for i, sequence in enumerate(positives, start=1):
                identifier = f"positive_{i:06d}|label=1"
                sequence_ids.append(identifier)
                fasta_records.append((identifier, sequence))
            for i, sequence in enumerate(null_sequences, start=1):
                identifier = f"null_{i:06d}|label=0"
                sequence_ids.append(identifier)
                fasta_records.append((identifier, sequence))
            with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
                for identifier, sequence in fasta_records:
                    handle.write(f">{identifier}\n{sequence}\n")
            matrix_path = matrix_dir / f"{slug}.meme"
            raw_path = work_dir / f"{slug}.replicate_{replicate:02d}.fimo.tsv"
            command = [
                str(fimo), "--text", "--skip-matched-sequence", "--thresh", f"{thresholds[0]:g}",
                "--bfile", str(background_path), str(matrix_path), str(fasta_path),
            ]
            started = subprocess.run(command, capture_output=True, text=True, check=False, timeout=600)
            if started.returncode != 0:
                raise RuntimeError(f"FIMO failed for {slug} replicate {replicate}: {started.stderr[-3000:]}")
            raw_path.write_text(started.stdout, encoding="utf-8")
            hits = _parse_fimo_tsv(raw_path)
            positive_ids = sequence_ids[:len(positives)]
            null_ids = sequence_ids[len(positives):]
            paired_indexes = set(range(len(positives)))
            for threshold in thresholds:
                sites, present = _site_sets(hits, sequence_ids, motif_ids, threshold)
                for motif_a, motif_b in combinations(motif_ids, 2):
                    positive_a = present[motif_a] & paired_indexes
                    positive_b = present[motif_b] & paired_indexes
                    null_a = {index - len(positives) for index in present[motif_a] if index >= len(positives)}
                    null_b = {index - len(positives) for index in present[motif_b] if index >= len(positives)}
                    paired = summarize_pair_presence(positive_a, positive_b, null_a, null_b, len(positives))
                    contingency = {
                        "species": slug,
                        "shuffle_replicate": replicate,
                        "site_p_cutoff": threshold,
                        "motif_a": motif_a,
                        "motif_b": motif_b,
                        "positive_unique_sequences": len(positives),
                        **paired,
                        "positive_pair_independence_bh_q_270_tests": None,
                        "paired_shuffle_bh_q_270_tests": None,
                        "paired_shuffle_worst_replicate_p_value": None,
                        "paired_shuffle_worst_replicate_bh_q_270_tests": None,
                    }
                    pair_rows.append(contingency)
                    if replicate == 1:
                        independence_tests[threshold].append(contingency)
                    paired_tests[threshold].append(contingency)

                    def nearest_for_class(ids: list[str], class_sites: dict[str, dict[str, list[dict[str, Any]]]]) -> tuple[list[int], Counter[str]]:
                        gaps: list[int] = []
                        orders: Counter[str] = Counter()
                        for identifier in ids:
                            nearest = _nearest_site_pair(class_sites[identifier][motif_a], class_sites[identifier][motif_b])
                            if nearest is not None:
                                relation, gap = nearest
                                gaps.append(gap)
                                orders[relation] += 1
                        return gaps, orders

                    positive_gaps, positive_orders = nearest_for_class(positive_ids, sites)
                    null_gaps, null_orders = nearest_for_class(null_ids, sites)
                    spacing_rows.append({
                        "species": slug,
                        "shuffle_replicate": replicate,
                        "site_p_cutoff": threshold,
                        "motif_a": motif_a,
                        "motif_b": motif_b,
                        "positive_pair_sequences": len(positive_gaps),
                        "positive_signed_gap_median_bp": statistics.median(positive_gaps) if positive_gaps else "",
                        "positive_signed_gap_q1_bp": statistics.quantiles(positive_gaps, n=4, method="inclusive")[0] if len(positive_gaps) > 1 else (positive_gaps[0] if positive_gaps else ""),
                        "positive_signed_gap_q3_bp": statistics.quantiles(positive_gaps, n=4, method="inclusive")[2] if len(positive_gaps) > 1 else (positive_gaps[0] if positive_gaps else ""),
                        "positive_A_before_B": positive_orders["A_before_B"],
                        "positive_B_before_A": positive_orders["B_before_A"],
                        "positive_overlapping": positive_orders["overlapping"],
                        "null_pair_sequences": len(null_gaps),
                        "null_signed_gap_median_bp": statistics.median(null_gaps) if null_gaps else "",
                        "null_A_before_B": null_orders["A_before_B"],
                        "null_B_before_A": null_orders["B_before_A"],
                        "null_overlapping": null_orders["overlapping"],
                        "coordinate_note": "signed gaps and order are relative to input-string coordinates; orientation is unverified",
                    })
            species_rep_runs[slug].append({
                "replicate": replicate,
                "seed": seed + (list(SPECIES).index(slug) * 100) + replicate - 1,
                "positive_records_before_dedup": len(raw_records[slug]),
                "positive_unique_sequences": len(positives),
                "shuffle_unique_sequences": len(set(null_sequences)),
                "positive_fasta_sha256": sha256_file(fasta_path),
                "background_sha256": sha256_file(background_path),
                "raw_fimo_tsv_sha256": sha256_file(raw_path),
                "command": command,
            })

    for threshold in thresholds:
        independence_rows = independence_tests[threshold]
        independence_q = benjamini_hochberg([
            row["positive_pair_fisher_independence_p_value"] for row in independence_rows
        ])
        independence_q_by_pair = {
            (row["species"], row["motif_a"], row["motif_b"]): q_value
            for row, q_value in zip(independence_rows, independence_q, strict=True)
        }

        replicate_q_by_pair: dict[tuple[str, str, str, int], float] = {}
        paired_rows = paired_tests[threshold]
        for replicate in range(1, replicates + 1):
            replicate_rows = [row for row in paired_rows if row["shuffle_replicate"] == replicate]
            replicate_q = benjamini_hochberg([row["paired_mcnemar_p_value"] for row in replicate_rows])
            replicate_q_by_pair.update({
                (row["species"], row["motif_a"], row["motif_b"], replicate): q_value
                for row, q_value in zip(replicate_rows, replicate_q, strict=True)
            })

        paired_by_pair: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        for row in paired_rows:
            key = (row["species"], row["motif_a"], row["motif_b"])
            paired_by_pair[key].append(row["paired_mcnemar_p_value"])
        pair_keys = sorted(paired_by_pair)
        worst_p_values = [max(paired_by_pair[key]) for key in pair_keys]
        worst_q_values = benjamini_hochberg(worst_p_values)
        worst_q_by_pair = dict(zip(pair_keys, worst_q_values, strict=True))

        for row in paired_rows:
            key = (row["species"], row["motif_a"], row["motif_b"])
            row["positive_pair_independence_bh_q_270_tests"] = independence_q_by_pair[key]
            row["paired_shuffle_bh_q_270_tests"] = replicate_q_by_pair[
                (*key, row["shuffle_replicate"])
            ]
            row["paired_shuffle_worst_replicate_p_value"] = max(paired_by_pair[key])
            row["paired_shuffle_worst_replicate_bh_q_270_tests"] = worst_q_by_pair[key]

    output_dir.mkdir(parents=True, exist_ok=True)
    pairs_path = output_dir / "motif_combination_tests.tsv"
    spacing_path = output_dir / "motif_combination_spacing.tsv"
    summary_rows: list[dict[str, Any]] = []
    for threshold in thresholds:
        for slug in SPECIES:
            species_rows = [row for row in paired_tests[threshold] if row["species"] == slug]
            unique_pair_rows = [row for row in species_rows if row["shuffle_replicate"] == 1]
            robust_rows = [
                row for row in unique_pair_rows
                if row["paired_shuffle_worst_replicate_bh_q_270_tests"] <= 0.05
            ]
            replicate_significant = {
                str(replicate): sum(
                    row["species"] == slug
                    and row["shuffle_replicate"] == replicate
                    and row["paired_shuffle_bh_q_270_tests"] <= 0.05
                    for row in paired_tests[threshold]
                )
                for replicate in range(1, replicates + 1)
            }
            deltas_by_pair: list[float] = []
            positive_excess = null_excess = 0
            for row in robust_rows:
                same_pair = [
                    candidate for candidate in species_rows
                    if candidate["motif_a"] == row["motif_a"]
                    and candidate["motif_b"] == row["motif_b"]
                ]
                mean_delta = statistics.mean(
                    candidate["positive_pair_cooccurrence"] - candidate["null_pair_cooccurrence"]
                    for candidate in same_pair
                )
                deltas_by_pair.append(100 * mean_delta / row["positive_unique_sequences"])
                if mean_delta > 0:
                    positive_excess += 1
                elif mean_delta < 0:
                    null_excess += 1
            summary_rows.append({
                "species": slug,
                "site_p_cutoff": threshold,
                "unique_motif_pairs": len(unique_pair_rows),
                "positive_pair_fisher_bh_q_le_0_05": sum(
                    row["positive_pair_independence_bh_q_270_tests"] <= 0.05
                    for row in unique_pair_rows
                ),
                **{
                    f"matched_shuffle_seed_{replicate}_bh_q_le_0_05": replicate_significant[str(replicate)]
                    for replicate in range(1, replicates + 1)
                },
                "matched_shuffle_worst_replicate_bh_q_le_0_05": len(robust_rows),
                "robust_pairs_positive_excess": positive_excess,
                "robust_pairs_null_excess": null_excess,
                "median_robust_pair_coverage_difference_percentage_points": (
                    statistics.median(deltas_by_pair) if deltas_by_pair else ""
                ),
            })

    _write_tsv(pairs_path, pair_rows)
    _write_tsv(spacing_path, spacing_rows)
    summary_path = output_dir / "motif_combination_summary.tsv"
    _write_tsv(summary_path, summary_rows)
    manifest = {
        "schema_version": "tjupan-motif-combination-test-1.0",
        "status": "complete",
        "purpose": "sequence-level motif combination and spacing comparisons against paired dinucleotide-preserving nulls",
        "tool": "MEME Suite FIMO",
        "version": version,
        "dataset": "TJU Pan six-species positive promoter sequences",
        "motif_count": 60,
        "species": list(SPECIES),
        "motif_pairs_per_species": 45,
        "replicates_per_species": replicates,
        "seed_rule": "base seed + 100 * species index + replicate index - 1",
        "site_p_thresholds": list(thresholds),
        "sequence_unit": "unique exact positive sequences; one dinucleotide-preserving shuffle paired per positive sequence and replicate",
        "background_model": "per-species mononucleotide frequencies; dinucleotide shuffle preserves each positive sequence's dinucleotide and mononucleotide counts",
        "independence_test": "two-sided Fisher exact test of motif A x motif B presence within true positive sequences; BH across 270 species x motif-pair tests per site-p cutoff",
        "matched_null_test": "paired exact McNemar test comparing motif-pair co-occurrence for each true sequence and its dinucleotide shuffle; BH across 270 species x motif-pair tests separately per replicate and per cutoff; robust q-value uses each pair's worst (largest) p-value across all shuffle replicates, then BH across 270 unique pairs",
        "spacing": "one nearest site pair per co-occurring sequence; positive gap counts intervening bases, negative gap counts overlap; order uses input-string orientation",
        "interpretation_limit": "motifs were discovered on the same positives; this is an exploratory combination analysis and not independent validation; no TSS-relative orientation is assumed",
        "input_sha256": {
            "input_manifest": sha256_file(input_manifest_path),
            "fimo_run_manifest": sha256_file(fimo_manifest_path),
            "motif_matrices": matrix_hashes,
        },
        "paired_rows": len(pair_rows),
        "spacing_rows": len(spacing_rows),
        "positive_pair_independence_tests_bh_q_le_0_05_by_cutoff": {
            f"{threshold:g}": sum(row["positive_pair_independence_bh_q_270_tests"] <= 0.05 for row in independence_tests[threshold])
            for threshold in thresholds
        },
        "matched_shuffle_replicate_bh_q_le_0_05_by_cutoff": {
            f"{threshold:g}": {
                str(replicate): sum(
                    row["paired_shuffle_bh_q_270_tests"] <= 0.05
                    for row in paired_tests[threshold] if row["shuffle_replicate"] == replicate
                )
                for replicate in range(1, replicates + 1)
            }
            for threshold in thresholds
        },
        "matched_shuffle_worst_replicate_bh_q_le_0_05_by_cutoff": {
            f"{threshold:g}": sum(
                row["paired_shuffle_worst_replicate_bh_q_270_tests"] <= 0.05
                for row in paired_tests[threshold] if row["shuffle_replicate"] == 1
            )
            for threshold in thresholds
        },
        "output_sha256": {path.name: sha256_file(path) for path in (pairs_path, spacing_path, summary_path)},
        "species_runs": species_rep_runs,
    }
    manifest_path = output_dir / "motif_combination_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/input_manifest.json"))
    parser.add_argument("--fimo-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo/fimo_run_manifest.json"))
    parser.add_argument("--matrix-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/summary/matrices"))
    parser.add_argument("--fimo", type=Path, default=Path("tmp/meme-suite-5.5.9/bin/fimo"))
    parser.add_argument("--work-dir", type=Path, default=Path("tmp/tjupan_motif_combinations_20261003"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/combinations"))
    parser.add_argument("--replicates", type=int, default=3)
    args = parser.parse_args()
    manifest = run(
        input_manifest_path=args.input_manifest,
        fimo_manifest_path=args.fimo_manifest,
        matrix_dir=args.matrix_dir,
        fimo=args.fimo,
        work_dir=args.work_dir,
        output_dir=args.output_dir,
        replicates=args.replicates,
    )
    print("independence BH q<=0.05:", manifest["positive_pair_independence_tests_bh_q_le_0_05_by_cutoff"])
    print("paired-shuffle replicate BH q<=0.05:", manifest["matched_shuffle_replicate_bh_q_le_0_05_by_cutoff"])
    print("paired-shuffle worst-replicate BH q<=0.05:", manifest["matched_shuffle_worst_replicate_bh_q_le_0_05_by_cutoff"])


if __name__ == "__main__":
    main()
