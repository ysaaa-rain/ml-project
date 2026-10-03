"""Analyze sequence-level motif co-occurrence and within-sequence spacing."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from itertools import combinations
from pathlib import Path
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES
from experiments.summarize_tjupan_fimo import benjamini_hochberg, fisher_exact_two_sided, _read_hits
from preprocessing.provenance import sha256_file


def _read_fasta(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    identifier: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if identifier is None:
            return
        sequence = "".join(chunks).upper()
        label = identifier.rsplit("|label=", 1)[-1]
        if label not in {"0", "1"} or len(sequence) != 81 or set(sequence) - set("ACGT"):
            raise ValueError(f"invalid labeled 81-bp sequence in {path}: {identifier}")
        if identifier in records:
            raise ValueError(f"duplicate FASTA identifier: {identifier}")
        records[identifier] = {
            "label": label,
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        }

    for line_number, raw in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            identifier, chunks = line[1:], []
        elif identifier is None:
            raise ValueError(f"sequence before header at {path}:{line_number}")
        else:
            chunks.append(line)
    finish()
    return records


def site_relation(first: dict[str, Any], second: dict[str, Any]) -> tuple[str, int]:
    """Return left-to-right order and signed gap; overlap is negative width."""
    if int(first["stop"]) < int(second["start"]):
        return "A_before_B", int(second["start"]) - int(first["stop"]) - 1
    if int(second["stop"]) < int(first["start"]):
        return "B_before_A", int(first["start"]) - int(second["stop"]) - 1
    overlap_bp = min(int(first["stop"]), int(second["stop"])) - max(int(first["start"]), int(second["start"])) + 1
    return "overlapping", -overlap_bp


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table without a schema: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _quartiles(values: list[int]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    if len(values) == 1:
        value = float(values[0])
        return value, value
    lower, upper = statistics.quantiles(values, n=4, method="inclusive")[0], statistics.quantiles(values, n=4, method="inclusive")[2]
    return lower, upper


def analyze(
    *,
    fimo_run_manifest: Path,
    input_manifest: Path,
    meme_summary: Path,
    run_root: Path,
    output_dir: Path,
    q_threshold: float = 0.05,
) -> dict[str, Any]:
    run_manifest = json.loads(fimo_run_manifest.read_text(encoding="utf-8"))
    json.loads(input_manifest.read_text(encoding="utf-8"))
    if run_manifest.get("status") != "complete" or set(run_manifest.get("results", {})) != set(SPECIES):
        raise ValueError("a complete six-species de novo FIMO run is required")
    if run_manifest.get("input_manifest_sha256") != sha256_file(input_manifest):
        raise ValueError("supplied input manifest does not match the de novo FIMO run")
    if run_manifest.get("meme_summary_manifest_sha256") is None:
        raise ValueError("de novo FIMO manifest does not record its MEME summary dependency")
    if not 0 < q_threshold <= 1:
        raise ValueError("q threshold must be in (0, 1]")

    with meme_summary.open(encoding="utf-8", newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    motifs_by_species: dict[str, dict[str, dict[str, str]]] = {slug: {} for slug in SPECIES}
    for row in summary_rows:
        slug = row["species"]
        if slug not in motifs_by_species:
            raise ValueError(f"unexpected species in MEME summary: {slug}")
        motifs_by_species[slug][row["motif_id"]] = row
    if any(len(motifs_by_species[slug]) != 10 for slug in SPECIES):
        raise ValueError("expected exactly 10 de novo motifs per species")

    sequence_class: dict[str, dict[str, str]] = {}
    sequence_sites: dict[str, dict[str, dict[str, set[tuple[int, int, str]]]]] = {}
    source_hashes: dict[str, str] = {}
    for slug in SPECIES:
        run_record = run_manifest["results"][slug]
        fasta_path = Path(run_record["combined_fasta"])
        fimo_path = run_root / slug / "fimo.tsv"
        if not fasta_path.is_file() or sha256_file(fasta_path) != run_record["combined_fasta_sha256"]:
            raise ValueError(f"combined FASTA is missing or changed: {fasta_path}")
        if not fimo_path.is_file() or sha256_file(fimo_path) != run_record["fimo_tsv_sha256"]:
            raise ValueError(f"FIMO output is missing or changed: {fimo_path}")
        source_hashes[f"{slug}/combined_fasta"] = sha256_file(fasta_path)
        source_hashes[f"{slug}/fimo_tsv"] = sha256_file(fimo_path)
        records = _read_fasta(fasta_path)
        seq_class: dict[str, str] = {}
        id_to_hash: dict[str, str] = {}
        for identifier, info in records.items():
            seq_hash = info["sequence_sha256"]
            prior_label = seq_class.setdefault(seq_hash, info["label"])
            if prior_label != info["label"]:
                raise ValueError(f"identical sequence has conflicting labels in {slug}")
            id_to_hash[identifier] = seq_hash
        by_sequence: dict[str, dict[str, set[tuple[int, int, str]]]] = {
            seq_hash: {motif_id: set() for motif_id in motifs_by_species[slug]}
            for seq_hash in seq_class
        }
        for hit in _read_hits(fimo_path, slug):
            if hit["motif_id"] not in motifs_by_species[slug]:
                raise ValueError(f"unexpected motif in FIMO output: {slug} {hit['motif_id']}")
            if hit["q_value"] <= q_threshold:
                seq_hash = id_to_hash[hit["sequence_name"]]
                by_sequence[seq_hash][hit["motif_id"]].add((hit["start"], hit["stop"], hit["strand"]))
        sequence_class[slug] = seq_class
        sequence_sites[slug] = by_sequence

    cooccurrence_rows: list[dict[str, Any]] = []
    spacing_rows: list[dict[str, Any]] = []
    tests: list[tuple[dict[str, Any], float]] = []
    for slug in SPECIES:
        motif_ids = sorted(motifs_by_species[slug], key=lambda value: int(value.rsplit("-", 1)[1]))
        for label, label_name in (("1", "positive"), ("0", "control")):
            seq_hashes = [seq_hash for seq_hash, item_label in sequence_class[slug].items() if item_label == label]
            total = len(seq_hashes)
            for motif_a, motif_b in combinations(motif_ids, 2):
                a_set = {seq_hash for seq_hash in seq_hashes if sequence_sites[slug][seq_hash][motif_a]}
                b_set = {seq_hash for seq_hash in seq_hashes if sequence_sites[slug][seq_hash][motif_b]}
                both = a_set & b_set
                a_only, b_only = a_set - b_set, b_set - a_set
                neither = total - len(a_set | b_set)
                p_value = fisher_exact_two_sided(len(both), len(a_only), len(b_only), neither)
                co_row: dict[str, Any] = {
                    "species": slug,
                    "label_group": label_name,
                    "motif_a": motif_a,
                    "motif_a_consensus": motifs_by_species[slug][motif_a]["consensus"],
                    "motif_b": motif_b,
                    "motif_b_consensus": motifs_by_species[slug][motif_b]["consensus"],
                    "unique_sequences": total,
                    "motif_a_positive_sequences": len(a_set),
                    "motif_b_positive_sequences": len(b_set),
                    "cooccurring_unique_sequences": len(both),
                    "cooccurrence_independence_fisher_p": p_value,
                    "bh_q_value_540_tests": None,
                }
                tests.append((co_row, p_value))

                selected_gaps: list[int] = []
                order_counts = {"A_before_B": 0, "B_before_A": 0, "overlapping": 0}
                for seq_hash in both:
                    sites_a = sorted(sequence_sites[slug][seq_hash][motif_a])
                    sites_b = sorted(sequence_sites[slug][seq_hash][motif_b])
                    candidates: list[tuple[int, int, str]] = []
                    for start_a, stop_a, strand_a in sites_a:
                        for start_b, stop_b, strand_b in sites_b:
                            relation, gap = site_relation(
                                {"start": start_a, "stop": stop_a},
                                {"start": start_b, "stop": stop_b},
                            )
                            separation = 0 if relation == "overlapping" else gap
                            candidates.append((separation, gap, relation))
                    _, nearest_gap, nearest_relation = min(candidates, key=lambda item: (item[0], abs(item[1]), item[2]))
                    selected_gaps.append(nearest_gap)
                    order_counts[nearest_relation] += 1
                q1, q3 = _quartiles(selected_gaps)
                spacing_rows.append({
                    "species": slug,
                    "label_group": label_name,
                    "motif_a": motif_a,
                    "motif_b": motif_b,
                    "cooccurring_unique_sequences": len(both),
                    "sequences_with_selected_overlapping_sites": order_counts["overlapping"],
                    "sequences_motif_a_before_b": order_counts["A_before_B"],
                    "sequences_motif_b_before_a": order_counts["B_before_A"],
                    "nearest_signed_gap_median_bp": statistics.median(selected_gaps) if selected_gaps else None,
                    "nearest_signed_gap_q1_bp": q1,
                    "nearest_signed_gap_q3_bp": q3,
                    "gap_definition": "positive=unmatched bases between intervals; negative=number of shared bases; order stored separately",
                })
    adjusted = benjamini_hochberg([p_value for _, p_value in tests])
    for (row, _), q_value in zip(tests, adjusted, strict=True):
        row["bh_q_value_540_tests"] = q_value

    output_dir.mkdir(parents=True, exist_ok=True)
    cooccurrence_path = output_dir / "motif_pair_cooccurrence.tsv"
    spacing_path = output_dir / "motif_pair_spacing.tsv"
    _write_tsv(cooccurrence_path, cooccurrence_rows := [row for row, _ in tests])
    _write_tsv(spacing_path, spacing_rows)
    manifest = {
        "schema_version": "tjupan-motif-pairs-1.0",
        "project": "PR01-02",
        "status": "complete",
        "purpose": "exploratory sequence-level co-occurrence and nearest-site spacing/order for de novo FIMO hits",
        "dataset": run_manifest.get("dataset"),
        "fimo_version": run_manifest.get("fimo_version"),
        "q_threshold_per_motif": q_threshold,
        "sequence_unit": "unique exact DNA sequence within each species and label; duplicate records are counted once",
        "cooccurrence_test": "two-sided Fisher exact test of motif-presence association within each species and label; BH across 6 species x 2 labels x 45 motif pairs = 540 tests",
        "spacing_rule": "one nearest pair of FIMO sites per co-occurring unique sequence; signed gap is positive for intervening bases and negative for shared bases; order is relative to input-string coordinates",
        "interpretation_limits": [
            "MEME motifs were discovered from the same positive sequences later scanned; all pair results are exploratory and selection-biased",
            "FIMO site q-values are per motif and do not correct for the prior motif-discovery step",
            "control labels are not experimentally verified non-promoters",
            "input strand/orientation is unverified; reported order is only left-to-right in the supplied string",
            "nearest-site spacing is descriptive; no positional randomization test is included in this run",
            "co-occurrence or overlap does not establish separate biological elements or physical interaction",
        ],
        "species": list(SPECIES),
        "motif_pairs_per_species": 45,
        "tests": len(cooccurrence_rows),
        "tests_bh_q_le_0_05": sum(row["bh_q_value_540_tests"] <= 0.05 for row in cooccurrence_rows),
        "input_sha256": {
            "fimo_run_manifest": sha256_file(fimo_run_manifest),
            "meme_summary": sha256_file(meme_summary),
            "input_manifest": sha256_file(input_manifest),
            "per_species_inputs": source_hashes,
        },
        "output_sha256": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in (cooccurrence_path, spacing_path)
        },
    }
    manifest_path = output_dir / "motif_pair_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fimo-run-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo/fimo_run_manifest.json"))
    parser.add_argument("--input-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/input_manifest.json"))
    parser.add_argument("--meme-summary", type=Path, default=Path("results/motif/tjupan_m3_20261003/summary/motif_summary.tsv"))
    parser.add_argument("--run-root", type=Path, default=Path("tmp/fimo_runs/tjupan_promloop_20261003"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/pairs"))
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args()
    manifest = analyze(
        fimo_run_manifest=args.fimo_run_manifest,
        input_manifest=args.input_manifest,
        meme_summary=args.meme_summary,
        run_root=args.run_root,
        output_dir=args.output_dir,
        q_threshold=args.q_threshold,
    )
    print(f"status={manifest['status']}")
    print(f"BH q<=0.05: {manifest['tests_bh_q_le_0_05']}/{manifest['tests']}")


if __name__ == "__main__":
    main()
