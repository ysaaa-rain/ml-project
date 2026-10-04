"""Compare the six TJU Pan de novo motif libraries with one pooled Tomtom run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES
from experiments.summarize_tjupan_fimo import benjamini_hochberg
from preprocessing.provenance import sha256_file


TOMTOM_COLUMNS = (
    "Query_ID", "Target_ID", "Optimal_offset", "p-value", "E-value", "q-value",
    "Overlap", "Query_consensus", "Target_consensus", "Orientation",
)


def _matrix_blocks(path: Path) -> list[str]:
    blocks: list[str] = []
    current: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MOTIF "):
            if current:
                blocks.append("\n".join(current))
            current = [line]
        elif current:
            current.append(line)
    if current:
        blocks.append("\n".join(current))
    if len(blocks) != 10:
        raise ValueError(f"expected 10 matrices in {path}, found {len(blocks)}")
    return blocks


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(
    *,
    matrix_dir: Path,
    fimo_manifest_path: Path,
    tomtom: Path,
    work_dir: Path,
    output_dir: Path,
    q_threshold: float = 0.05,
) -> dict[str, Any]:
    fimo_manifest = json.loads(fimo_manifest_path.read_text(encoding="utf-8"))
    if fimo_manifest.get("fimo_version") != "5.5.9" or set(fimo_manifest.get("results", {})) != set(SPECIES):
        raise ValueError("complete six-species FIMO manifest required")
    pooled_counts = {base: 0 for base in "ACGT"}
    for slug in SPECIES:
        counts = fimo_manifest["results"][slug]["background_base_counts"]
        for base in pooled_counts:
            pooled_counts[base] += int(counts[base])
    total = sum(pooled_counts.values())
    background = {base: pooled_counts[base] / total for base in "ACGT"}
    work_dir.mkdir(parents=True, exist_ok=True)
    pooled_meme = work_dir / "six_species_motifs_pooled_background.meme"
    sections = [
        "MEME version 5",
        "",
        "ALPHABET= ACGT",
        "",
        "strands: + -",
        "",
        "Background letter frequencies",
        " ".join(f"{base} {background[base]:.12g}" for base in "ACGT"),
        "",
    ]
    matrix_hashes: dict[str, str] = {}
    for slug in SPECIES:
        matrix_path = matrix_dir / f"{slug}.meme"
        if not matrix_path.is_file():
            raise FileNotFoundError(matrix_path)
        matrix_hashes[slug] = sha256_file(matrix_path)
        sections.extend(_matrix_blocks(matrix_path))
        sections.append("")
    pooled_meme.write_text("\n".join(sections), encoding="utf-8")

    version_run = subprocess.run([str(tomtom), "-version"], capture_output=True, text=True, check=True, timeout=15)
    version = (version_run.stdout or version_run.stderr).strip().splitlines()[0]
    if version != "5.5.9":
        raise ValueError(f"expected Tomtom 5.5.9, got {version}")
    command = [
        str(tomtom), "-text", "-dist", "ed", "-min-overlap", "6", "-thresh", "1",
        str(pooled_meme), str(pooled_meme),
    ]
    process = subprocess.run(command, capture_output=True, text=True, check=False, timeout=300)
    if process.returncode != 0:
        raise RuntimeError(f"Tomtom failed: {process.stderr[-4000:]}")
    raw_path = work_dir / "tomtom_all_pairs.tsv"
    raw_path.write_text(process.stdout, encoding="utf-8")
    lines = [line for line in process.stdout.splitlines() if line.strip() and not line.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if not rows or tuple(rows[0]) != TOMTOM_COLUMNS:
        raise ValueError("unexpected or empty Tomtom output schema")
    expected_ids = {f"{slug}_MEME-{i}" for slug in SPECIES for i in range(1, 11)}
    query_ids = {row["Query_ID"] for row in rows}
    if query_ids != expected_ids:
        raise ValueError("Tomtom query set does not contain all 60 motifs")
    if any(row["Target_ID"] not in expected_ids for row in rows):
        raise ValueError("Tomtom target set does not contain the expected 60 motifs")

    significant: list[dict[str, Any]] = []
    by_species_pair: dict[tuple[str, str], list[dict[str, Any]]] = {}
    all_cross_species: list[dict[str, Any]] = []
    for row in rows:
        query_species = row["Query_ID"].rsplit("_MEME-", 1)[0]
        target_species = row["Target_ID"].rsplit("_MEME-", 1)[0]
        if query_species == target_species:
            continue
        record = {
            "query_species": query_species,
            "query_motif": row["Query_ID"].rsplit("_", 1)[-1],
            "target_species": target_species,
            "target_motif": row["Target_ID"].rsplit("_", 1)[-1],
            "p_value": float(row["p-value"]),
            "e_value": float(row["E-value"]),
            "q_value_within_query_60_targets": float(row["q-value"]),
            "global_bh_q_3000_cross_species_comparisons": None,
            "overlap_bp": int(row["Overlap"]),
            "orientation": row["Orientation"],
            "query_consensus": row["Query_consensus"],
            "target_consensus": row["Target_consensus"],
        }
        key = (query_species, target_species)
        by_species_pair.setdefault(key, []).append(record)
        all_cross_species.append(record)
    global_q = benjamini_hochberg([record["p_value"] for record in all_cross_species])
    for record, q_value in zip(all_cross_species, global_q, strict=True):
        record["global_bh_q_3000_cross_species_comparisons"] = q_value
        if q_value <= q_threshold:
            significant.append(record)

    species_pair_rows: list[dict[str, Any]] = []
    species_order = list(SPECIES)
    for i, species_a in enumerate(species_order):
        for species_b in species_order[i + 1:]:
            forward = by_species_pair.get((species_a, species_b), [])
            reverse = by_species_pair.get((species_b, species_a), [])
            significant_forward = [r for r in forward if r["q_value_within_query_60_targets"] <= q_threshold]
            significant_reverse = [r for r in reverse if r["q_value_within_query_60_targets"] <= q_threshold]
            global_forward = [r for r in forward if r["global_bh_q_3000_cross_species_comparisons"] <= q_threshold]
            global_reverse = [r for r in reverse if r["global_bh_q_3000_cross_species_comparisons"] <= q_threshold]
            combined = forward + reverse
            best = min(combined, key=lambda r: (r["global_bh_q_3000_cross_species_comparisons"], r["p_value"]))
            species_pair_rows.append({
                "species_a": species_a,
                "species_b": species_b,
                "per_query_q_candidates_a_to_b": len(significant_forward),
                "per_query_q_candidates_b_to_a": len(significant_reverse),
                "global_bh_significant_pairs_a_to_b": len(global_forward),
                "global_bh_significant_pairs_b_to_a": len(global_reverse),
                "best_cross_species_query_motif": f"{best['query_species']}:{best['query_motif']}",
                "best_cross_species_target_motif": f"{best['target_species']}:{best['target_motif']}",
                "best_p_value": best["p_value"],
                "best_e_value": best["e_value"],
                "best_q_value_within_query_60_targets": best["q_value_within_query_60_targets"],
                "best_global_bh_q_3000_cross_species_comparisons": best["global_bh_q_3000_cross_species_comparisons"],
                "best_overlap_bp": best["overlap_bp"],
                "best_orientation": best["orientation"],
            })

    output_dir.mkdir(parents=True, exist_ok=True)
    significant_path = output_dir / "significant_cross_species_motif_matches.tsv"
    sig_columns = [
        "query_species", "query_motif", "target_species", "target_motif", "p_value", "e_value",
        "q_value_within_query_60_targets", "global_bh_q_3000_cross_species_comparisons",
        "overlap_bp", "orientation", "query_consensus", "target_consensus",
    ]
    _write_tsv(significant_path, sig_columns, significant)
    pair_path = output_dir / "species_pair_similarity.tsv"
    _write_tsv(pair_path, list(species_pair_rows[0]), species_pair_rows)
    manifest = {
        "schema_version": "tjupan-cross-species-motif-similarity-1.0",
        "status": "complete",
        "tool": "MEME Suite Tomtom",
        "version": version,
        "query_and_target_motifs": 60,
        "ordered_query_target_comparisons": len(rows),
        "cross_species_comparisons": sum(len(v) for v in by_species_pair.values()),
        "q_threshold": q_threshold,
        "q_value_scope": "Tomtom q-value reported per query motif across the 60 target motifs, including same-species targets",
        "global_test_correction": "Benjamini-Hochberg across all 3,000 directed cross-species query-target comparisons; used for confirmatory significance",
        "parameters": ["-dist ed", "-min-overlap 6", "-thresh 1", "reverse-complement target scoring enabled"],
        "background": "pooled mononucleotide frequencies from all six species label=0 controls",
        "pooled_background_frequencies": background,
        "matrix_sha256": matrix_hashes,
        "fimo_manifest_sha256": sha256_file(fimo_manifest_path),
        "pooled_meme_sha256": sha256_file(pooled_meme),
        "raw_tomtom_sha256": sha256_file(raw_path),
        "per_query_q_candidate_cross_species_pairs": sum(
            record["q_value_within_query_60_targets"] <= q_threshold for record in all_cross_species
        ),
        "globally_significant_cross_species_pairs_bh_q_le_threshold": len(significant),
        "species_pair_count": len(species_pair_rows),
        "biological_interpretation": "PWM similarity identifies sequence-pattern resemblance; it does not establish shared function, genomic position, or common regulatory mechanism.",
        "output_sha256": {p.name: sha256_file(p) for p in (significant_path, pair_path)},
    }
    manifest_path = output_dir / "cross_species_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/summary/matrices"))
    parser.add_argument("--fimo-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo/fimo_run_manifest.json"))
    parser.add_argument("--tomtom", type=Path, default=Path("tmp/meme-suite-5.5.9/bin/tomtom"))
    parser.add_argument("--work-dir", type=Path, default=Path("tmp/tjupan_cross_species_20261003"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/cross_species"))
    args = parser.parse_args()
    manifest = run(
        matrix_dir=args.matrix_dir,
        fimo_manifest_path=args.fimo_manifest,
        tomtom=args.tomtom,
        work_dir=args.work_dir,
        output_dir=args.output_dir,
    )
    print(
        f"comparisons={manifest['ordered_query_target_comparisons']}; "
        f"per_query_candidates={manifest['per_query_q_candidate_cross_species_pairs']}; "
        f"global_significant={manifest['globally_significant_cross_species_pairs_bh_q_le_threshold']}"
    )


if __name__ == "__main__":
    main()
