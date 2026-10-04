"""Test E. coli TJU Pan motifs on sequence-disjoint RegulonDB promoters."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import subprocess
from typing import Any

from experiments.summarize_tjupan_fimo import benjamini_hochberg, fisher_exact_two_sided
from preprocessing.provenance import sha256_file


DNA = "ACGT"
BLOCKS = tuple((start, start + 9) for start in range(0, 81, 9))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def dinucleotide_shuffle(sequence: str, rng: random.Random) -> str:
    if len(sequence) != 81 or set(sequence) - set(DNA):
        raise ValueError("expected an 81-bp A/C/G/T sequence")
    original_counts = Counter(zip(sequence, sequence[1:]))
    for _ in range(64):
        edges: dict[str, list[str]] = defaultdict(list)
        for left, right in zip(sequence, sequence[1:]):
            edges[left].append(right)
        for destinations in edges.values():
            rng.shuffle(destinations)
        stack, path = [sequence[0]], []
        while stack:
            node = stack[-1]
            if edges[node]:
                stack.append(edges[node].pop())
            else:
                path.append(stack.pop())
        shuffled = "".join(reversed(path))
        if len(shuffled) == len(sequence) and Counter(zip(shuffled, shuffled[1:])) == original_counts:
            if shuffled != sequence:
                return shuffled
    return sequence


def _validate_sequence(sequence: str, source: str) -> str:
    sequence = sequence.strip().upper()
    if len(sequence) != 81 or set(sequence) - set(DNA):
        raise ValueError(f"{source}: expected 81 A/C/G/T bases")
    return sequence


def _candidate_index(sequences: list[str]) -> dict[str, set[int]]:
    index: dict[str, set[int]] = defaultdict(set)
    for i, sequence in enumerate(sequences):
        for start, stop in BLOCKS:
            index[sequence[start:stop]].add(i)
    return index


def _similarity_category(sequence: str, tju_sequences: list[str], index: dict[str, set[int]]) -> str | None:
    if len(sequence) != 81:
        raise ValueError("similarity comparison requires 81-bp sequences")
    reverse = reverse_complement(sequence)
    candidates: set[int] = set()
    for query in (sequence, reverse):
        for start, stop in BLOCKS:
            candidates.update(index.get(query[start:stop], ()))
    best = 81
    for i in candidates:
        target = tju_sequences[i]
        direct = sum(a == b for a, b in zip(sequence, target))
        rc = sum(a == b for a, b in zip(reverse, target))
        best = min(best, 81 - min(81, max(direct, rc)))
        if best == 0:
            return "exact_same_orientation" if sequence == target else "exact_reverse_complement"
    return "high_similarity_ge_90pct" if best <= 8 else None


def _read_csv_sequences(path: Path) -> tuple[list[str], list[str]]:
    positives: list[str] = []
    controls: list[str] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"seq", "label"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"invalid TJU Pan CSV columns: {path}")
        for row in reader:
            sequence = row["seq"].strip().upper()
            if row["label"] == "1":
                positives.append(_validate_sequence(sequence, str(path)))
            elif row["label"] == "0" and len(sequence) == 81 and not set(sequence) - set(DNA):
                controls.append(sequence)
    if not positives:
        raise ValueError("no TJU Pan positive sequences found")
    return positives, controls


def _read_regulondb(path: Path) -> tuple[list[str], dict[str, Any]]:
    records: dict[str, dict[str, str]] = {}
    row_count = 0
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence", "sequence_length", "source_dataset", "tss_position", "strand"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"invalid RegulonDB columns: {path}")
        for row in reader:
            row_count += 1
            if row["source_dataset"] != "regulondb" or int(row["sequence_length"]) != 81 or float(row["tss_position"]) != 61:
                continue
            sequence = _validate_sequence(row["sequence"], str(path))
            records.setdefault(sequence, row)
    return list(records), {"source_rows": row_count, "unique_81bp_tss61_sequences": len(records)}


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")


def _parse_fimo(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    required = {"motif_id", "sequence_name", "start", "stop", "strand", "p-value", "q-value"}
    if rows and not required.issubset(rows[0]):
        raise ValueError(f"unexpected FIMO output columns: {path}")
    hits = []
    for row in rows:
        hits.append({
            "motif_id": row["motif_id"],
            "sequence_name": row["sequence_name"],
            "start": int(row["start"]),
            "stop": int(row["stop"]),
            "strand": row["strand"],
            "p_value": float(row["p-value"]),
            "q_value": float(row["q-value"]),
        })
    return hits


def _write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        if not rows:
            raise ValueError(f"columns required for empty table: {path}")
        columns = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def exact_mcnemar_p(discordant_positive: int, discordant_null: int) -> float:
    """Two-sided exact McNemar test from paired discordant outcomes."""
    if min(discordant_positive, discordant_null) < 0:
        raise ValueError("discordant counts must be nonnegative")
    total = discordant_positive + discordant_null
    if total == 0:
        return 1.0
    tail_end = min(discordant_positive, discordant_null)
    log_terms = [
        math.lgamma(total + 1) - math.lgamma(k + 1) - math.lgamma(total - k + 1) - total * math.log(2)
        for k in range(tail_end + 1)
    ]
    peak = max(log_terms)
    lower_tail = math.exp(peak) * sum(math.exp(value - peak) for value in log_terms)
    return min(1.0, 2.0 * lower_tail)


def run(
    *,
    tju_csv: Path,
    regulondb_tsv: Path,
    motif_meme: Path,
    fimo: Path,
    work_dir: Path,
    output_dir: Path,
    seed: int = 20261003,
    shuffle_replicates: int = 5,
    q_threshold: float = 0.05,
) -> dict[str, Any]:
    if shuffle_replicates < 1:
        raise ValueError("shuffle_replicates must be positive")
    tju_positive_records, tju_controls = _read_csv_sequences(tju_csv)
    tju_positives = sorted(set(tju_positive_records))
    regulondb_sequences, regdb_counts = _read_regulondb(regulondb_tsv)
    index = _candidate_index(tju_positives)
    category_counts: Counter[str] = Counter()
    independent: list[str] = []
    for sequence in regulondb_sequences:
        category = _similarity_category(sequence, tju_positives, index)
        if category is None:
            independent.append(sequence)
        else:
            category_counts[category] += 1
    if not independent:
        raise ValueError("all RegulonDB sequences overlap the TJU Pan discovery sequences")

    work_dir.mkdir(parents=True, exist_ok=True)
    version_run = subprocess.run([str(fimo), "--version"], capture_output=True, text=True, check=True, timeout=15)
    version = (version_run.stdout or version_run.stderr).strip().splitlines()[0]
    if version != "5.5.9":
        raise ValueError(f"expected FIMO 5.5.9, got {version}")
    motif_ids = [
        line.split()[1]
        for line in motif_meme.read_text(encoding="utf-8").splitlines()
        if line.startswith("MOTIF ")
    ]
    if len(motif_ids) != 10:
        raise ValueError(f"expected 10 E. coli MEME motifs, found {len(motif_ids)}")
    summary_rows: list[dict[str, Any]] = []
    criteria = (("raw_p", 0.05), ("raw_p", 0.01), ("raw_p", 0.001), ("site_q", q_threshold))
    fimo_replicates: list[dict[str, Any]] = []
    for replicate in range(1, shuffle_replicates + 1):
        rng = random.Random(seed + replicate - 1)
        shuffled: list[str] = []
        unchanged = 0
        for sequence in independent:
            control = dinucleotide_shuffle(sequence, rng)
            unchanged += control == sequence
            shuffled.append(control)
        sequences = {"1": independent, "0": shuffled}
        shared = set(sequences["1"]) & set(sequences["0"])
        if len(set(shuffled)) != len(shuffled) or shared:
            raise ValueError("a shuffle replicate contains duplicate or cross-label sequences; choose another seed")
        id_sequence: dict[str, str] = {}
        fasta_records: list[tuple[str, str]] = []
        label_by_id: dict[str, str] = {}
        for label in ("1", "0"):
            for i, sequence in enumerate(sequences[label], start=1):
                identifier = f"{label}_{i:06d}|label={label}"
                fasta_records.append((identifier, sequence))
                id_sequence[identifier] = sequence
                label_by_id[identifier] = label
        input_fasta = work_dir / f"regulondb_independent_with_dinucleotide_null_r{replicate:02d}.fasta"
        _write_fasta(fasta_records, input_fasta)
        counts = Counter(base for sequence in shuffled for base in sequence)
        total = sum(counts.values())
        background = work_dir / f"dinucleotide_null_background_r{replicate:02d}.txt"
        background.write_text("".join(f"{base} {counts[base] / total:.12g}\n" for base in "ACGT"), encoding="ascii")
        run_key = hashlib.sha256((sha256_file(input_fasta) + sha256_file(motif_meme) + sha256_file(background)).encode()).hexdigest()[:12]
        fimo_dir = work_dir / f"fimo_{run_key}"
        if fimo_dir.exists():
            raise FileExistsError(f"FIMO output directory already exists; choose a new work directory: {fimo_dir}")
        command = [
            str(fimo), "--oc", str(fimo_dir), "--bfile", str(background), "--thresh", "0.05",
            "--max-stored-scores", "1000000", str(motif_meme), str(input_fasta),
        ]
        process = subprocess.run(command, capture_output=True, text=True, check=False, timeout=600)
        stderr_path = work_dir / f"fimo_{run_key}.stderr.log"
        stderr_path.write_text(process.stderr, encoding="utf-8")
        if process.returncode != 0:
            raise RuntimeError(f"FIMO failed; see {stderr_path}")
        hits = _parse_fimo(fimo_dir / "fimo.tsv")
        expected_ids = {identifier for identifier, _ in fasta_records}
        if any(hit["sequence_name"] not in expected_ids for hit in hits):
            raise ValueError("FIMO output contains an unknown sequence ID")
        by_motif: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for hit in hits:
            by_motif[hit["motif_id"]].append(hit)
        if set(by_motif) - set(motif_ids):
            raise ValueError("FIMO output contains an unknown motif ID")
        for motif in motif_ids:
            by_motif.setdefault(motif, [])
        replicate_rows: list[dict[str, Any]] = []
        by_criterion: dict[str, list[dict[str, Any]]] = {}
        for motif in motif_ids:
            motif_hits = by_motif[motif]
            total_unique = {label: {hashlib.sha256(s.encode("ascii")).hexdigest() for s in sequences[label]} for label in ("1", "0")}
            shared_hashes = {hashlib.sha256(s.encode("ascii")).hexdigest() for s in shared}
            for criterion, cutoff in criteria:
                selected = [hit for hit in motif_hits if (hit["p_value"] <= cutoff if criterion == "raw_p" else hit["q_value"] <= cutoff)]
                hit_ids = {
                    label: {
                        hashlib.sha256(id_sequence[hit["sequence_name"]].encode("ascii")).hexdigest()
                        for hit in selected if label_by_id[hit["sequence_name"]] == label
                    }
                    for label in ("1", "0")
                }
                positive_test, control_test = total_unique["1"] - shared_hashes, total_unique["0"] - shared_hashes
                positive_hit, control_hit = hit_ids["1"] - shared_hashes, hit_ids["0"] - shared_hashes
                p_value = fisher_exact_two_sided(
                    len(positive_hit), len(positive_test) - len(positive_hit),
                    len(control_hit), len(control_test) - len(control_hit),
                )
                positive_hits = [hit for hit in selected if label_by_id[hit["sequence_name"]] == "1"]
                positive_pair_keys = {
                    hit["sequence_name"].split("_", 1)[1].split("|", 1)[0]
                    for hit in selected if label_by_id[hit["sequence_name"]] == "1"
                }
                null_pair_keys = {
                    hit["sequence_name"].split("_", 1)[1].split("|", 1)[0]
                    for hit in selected if label_by_id[hit["sequence_name"]] == "0"
                }
                discordant_positive = len(positive_pair_keys - null_pair_keys)
                discordant_null = len(null_pair_keys - positive_pair_keys)
                row = {
                    "shuffle_replicate": replicate,
                    "motif_id": motif,
                    "hit_criterion": criterion,
                    "site_cutoff": cutoff,
                    "positive_sites": len(positive_hits),
                    "positive_unique_sequences_hit": len(hit_ids["1"]),
                    "positive_unique_sequences": len(total_unique["1"]),
                    "positive_hit_coverage": len(hit_ids["1"]) / len(total_unique["1"]),
                    "dinucleotide_null_sites": sum(1 for hit in selected if label_by_id[hit["sequence_name"]] == "0"),
                    "dinucleotide_null_unique_sequences_hit": len(hit_ids["0"]),
                    "dinucleotide_null_unique_sequences": len(total_unique["0"]),
                    "dinucleotide_null_hit_coverage": len(hit_ids["0"]) / len(total_unique["0"]),
                    "coverage_difference": len(hit_ids["1"]) / len(total_unique["1"]) - len(hit_ids["0"]) / len(total_unique["0"]),
                    "fisher_two_sided_p_value": p_value,
                    "fisher_bh_q_10_motifs_within_replicate_criterion": None,
                    "paired_positive_only_hits": discordant_positive,
                    "paired_null_only_hits": discordant_null,
                    "paired_exact_mcnemar_p_value": exact_mcnemar_p(discordant_positive, discordant_null),
                    "paired_mcnemar_bh_q_10_motifs_within_replicate_criterion": None,
                    "positive_hit_median_start_1based": statistics.median(hit["start"] for hit in positive_hits) if positive_hits else "",
                    "positive_hit_plus_strand_fraction": sum(hit["strand"] == "+" for hit in positive_hits) / len(positive_hits) if positive_hits else "",
                }
                replicate_rows.append(row)
                by_criterion.setdefault(f"{criterion}:{cutoff}", []).append(row)
        for rows in by_criterion.values():
            adjusted = benjamini_hochberg([row["fisher_two_sided_p_value"] for row in rows])
            paired_adjusted = benjamini_hochberg([row["paired_exact_mcnemar_p_value"] for row in rows])
            for row, q_value, paired_q in zip(rows, adjusted, paired_adjusted, strict=True):
                row["fisher_bh_q_10_motifs_within_replicate_criterion"] = q_value
                row["paired_mcnemar_bh_q_10_motifs_within_replicate_criterion"] = paired_q
        summary_rows.extend(replicate_rows)
        fimo_replicates.append({
            "replicate": replicate,
            "seed": seed + replicate - 1,
            "unchanged_shuffles": unchanged,
            "positive_control_shared_sequences": len(shared),
            "positive_unique_sequences": len(set(independent)),
            "dinucleotide_null_unique_sequences": len(set(shuffled)),
            "input_fasta_sha256": sha256_file(input_fasta),
            "background_sha256": sha256_file(background),
            "fimo_tsv_sha256": sha256_file(fimo_dir / "fimo.tsv"),
            "run_key": run_key,
            "command": command,
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "regulondb_external_motif_validation.tsv"
    _write_tsv(summary_path, summary_rows)
    overlap_path = output_dir / "dataset_overlap_summary.tsv"
    overlap_rows = [
        {"category": category, "unique_regulondb_sequences": count}
        for category, count in sorted(category_counts.items())
    ]
    overlap_rows.append({"category": "independent_sequences_retained", "unique_regulondb_sequences": len(independent)})
    _write_tsv(overlap_path, overlap_rows, ["category", "unique_regulondb_sequences"])
    manifest = {
        "schema_version": "regulondb-independent-motif-validation-1.0",
        "status": "complete",
        "purpose": "sequence-disjoint external check of TJU Pan E. coli MEME motifs",
        "tju_pan_positive_records": len(tju_positive_records),
        "tju_pan_unique_positive_sequences": len(tju_positives),
        "regulondb": regdb_counts,
        "overlap_filter": {
            "exact_same_orientation": category_counts["exact_same_orientation"],
            "exact_reverse_complement": category_counts["exact_reverse_complement"],
            "ungapped_hamming_identity_at_least_73_of_81": category_counts["high_similarity_ge_90pct"],
            "threshold": "at least 73/81 identical positions (>=90%), direct or reverse-complement orientation",
            "independent_sequences_retained": len(independent),
        },
        "validation_null": {
            "method": "five independently seeded Euler-trail dinucleotide-preserving shuffles per independent RegulonDB sequence; each replicate is tested separately",
            "seed": seed,
            "replicates": shuffle_replicates,
        },
        "fimo": {
            "version": version,
            "site_p_threshold_for_q_calculation": 0.05,
            "reported_site_q_threshold": q_threshold,
            "q_scope": "per motif across all positive and dinucleotide-shuffled sequences scanned",
            "reported_sensitivity": "FIMO site q<=0.05 and raw site p<=0.05, 0.01, 0.001; paired exact McNemar and unpaired Fisher tests are BH-adjusted across 10 motifs within each replicate and hit criterion",
            "background": "mononucleotide frequencies of the dinucleotide-preserving shuffled sequences",
            "replicate_runs": fimo_replicates,
        },
        "source_sha256": {
            "tju_pan_ecoli_csv": sha256_file(tju_csv),
            "regulondb_cleaned_promoters": sha256_file(regulondb_tsv),
            "ecoli_motif_matrix": sha256_file(motif_meme),
        },
        "summary_sha256": sha256_file(summary_path),
        "overlap_summary_sha256": sha256_file(overlap_path),
        "interpretation": "FIMO hits on sequence-disjoint positives above a matched dinucleotide-preserving null support motif recurrence; they do not alone establish functional activity or sigma specificity.",
    }
    manifest_path = output_dir / "regulondb_validation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tju-csv", type=Path, default=Path("data/raw/tju_pan_promoter/reg_and_gen/Datasets/Escherichia_coli/Dataset.csv"))
    parser.add_argument("--regulondb-tsv", type=Path, default=Path("data/processed/m1_regulondb_tss_20260930/promoters_clean.tsv"))
    parser.add_argument("--motif-meme", type=Path, default=Path("results/motif/tjupan_m3_20261003/summary/matrices/escherichia_coli.meme"))
    parser.add_argument("--fimo", type=Path, default=Path("tmp/meme-suite-5.5.9/bin/fimo"))
    parser.add_argument("--work-dir", type=Path, default=Path("tmp/regulondb_external_validation_20261003"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/cross_dataset"))
    args = parser.parse_args()
    manifest = run(
        tju_csv=args.tju_csv,
        regulondb_tsv=args.regulondb_tsv,
        motif_meme=args.motif_meme,
        fimo=args.fimo,
        work_dir=args.work_dir,
        output_dir=args.output_dir,
    )
    print(
        f"regulondb={manifest['regulondb']['unique_81bp_tss61_sequences']}; "
        f"purged={sum(manifest['overlap_filter'][key] for key in ('exact_same_orientation','exact_reverse_complement','ungapped_hamming_identity_at_least_73_of_81'))}; "
        f"independent={manifest['overlap_filter']['independent_sequences_retained']}"
    )


if __name__ == "__main__":
    main()
