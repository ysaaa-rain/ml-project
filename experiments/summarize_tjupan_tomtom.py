"""Summarize per-species TJU Pan Tomtom comparisons against the reference library."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES
from experiments.summarize_motif_comparison import validate_reference_consensus


REFERENCE_IDS = {
    "Ecoli_sigma70_minus10",
    "Ecoli_sigma70_minus35",
    "Ecoli_UP_proximal",
    "Ecoli_UP_distal",
}
TOMTOM_COLUMNS = (
    "Query_ID", "Target_ID", "Optimal_offset", "p-value", "E-value", "q-value",
    "Overlap", "Query_consensus", "Target_consensus", "Orientation",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    if not rows or set(rows[0]) != set(TOMTOM_COLUMNS):
        raise ValueError(f"invalid or empty Tomtom TSV: {path}")
    return rows


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(
    tomtom_root: Path,
    matrix_summary_dir: Path,
    reference: Path,
    output_dir: Path,
    species: list[str] | None = None,
) -> dict[str, Any]:
    selected_species = species or list(SPECIES)
    if len(set(selected_species)) != len(selected_species) or any(slug not in SPECIES for slug in selected_species):
        raise ValueError(f"invalid or repeated species list: {selected_species}")
    matrix_manifest_path = matrix_summary_dir / "summary_manifest.json"
    matrix_manifest = json.loads(matrix_manifest_path.read_text(encoding="utf-8"))
    matrix_summary_path = matrix_summary_dir / "motif_summary.tsv"
    with matrix_summary_path.open(encoding="utf-8", newline="") as handle:
        motif_rows = list(csv.DictReader(handle, delimiter="\t"))
    if matrix_manifest.get("included_species") != selected_species:
        raise ValueError("selected species do not match the matrix-summary manifest")
    reference_consensus = validate_reference_consensus(reference)

    all_rows: list[dict[str, Any]] = []
    best_rows: list[dict[str, Any]] = []
    tomtom_hashes: dict[str, str] = {}
    for slug in selected_species:
        matrix_file = matrix_summary_dir / "matrices" / f"{slug}.meme"
        if not matrix_file.is_file():
            raise FileNotFoundError(f"missing exported motif matrix: {matrix_file}")
        group_motifs = [row for row in motif_rows if row["species"] == slug]
        if len(group_motifs) != matrix_manifest["species_motif_counts"][slug]:
            raise ValueError(f"{slug}: motif-summary count does not match its manifest")
        tomtom_path = tomtom_root / slug / "tomtom.tsv"
        rows = _read_tsv(tomtom_path)
        by_query: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            by_query.setdefault(row["Query_ID"], []).append(row)
        expected_queries = {f"{slug}_{row['motif_id']}" for row in group_motifs}
        if set(by_query) != expected_queries:
            raise ValueError(f"{slug}: Tomtom queries do not match exported motifs")
        expected_rows = len(group_motifs) * len(REFERENCE_IDS)
        if len(rows) != expected_rows or any({match["Target_ID"] for match in values} != REFERENCE_IDS for values in by_query.values()):
            raise ValueError(f"{slug}: expected exactly {expected_rows} query-reference pairs")
        tomtom_hashes[f"{slug}/tomtom.tsv"] = sha256_file(tomtom_path)

        motif_by_id = {row["motif_id"]: row for row in group_motifs}
        for query_id in sorted(by_query, key=lambda item: int(item.rsplit("-", 1)[1])):
            meme_id = query_id.rsplit("_", 1)[1]
            motif = motif_by_id[meme_id]
            matches = sorted(by_query[query_id], key=lambda row: (float(row["p-value"]), float(row["q-value"])))
            for match in matches:
                all_rows.append({
                    "species": slug,
                    "query_meme_id": meme_id,
                    "query_consensus": motif["consensus"],
                    "query_meme_evalue": motif["e_value"],
                    "target_id": match["Target_ID"],
                    "tomtom_p_value_ranking_only": match["p-value"],
                    "tomtom_e_value": match["E-value"],
                    "tomtom_q_value_unreliable_small_target_db": match["q-value"],
                    "overlap": match["Overlap"],
                    "alignment_offset": match["Optimal_offset"],
                    "orientation": match["Orientation"],
                    "query_aligned_consensus": match["Query_consensus"],
                    "target_aligned_consensus": match["Target_consensus"],
                })
            best = matches[0]
            best_rows.append({
                "species": slug,
                "query_meme_id": meme_id,
                "query_consensus": motif["consensus"],
                "query_width_bp": motif["width_bp"],
                "query_sites_reported": motif["sites_reported"],
                "query_meme_evalue": motif["e_value"],
                "best_reference_match_descriptive_only": best["Target_ID"],
                "tomtom_p_value_ranking_only": best["p-value"],
                "tomtom_q_value_unreliable_small_target_db": best["q-value"],
                "overlap": best["Overlap"],
                "alignment_offset": best["Optimal_offset"],
                "orientation": best["Orientation"],
            })

    output_dir.mkdir(parents=True, exist_ok=True)
    all_path = output_dir / "all_species_reference_matches.tsv"
    best_path = output_dir / "best_reference_match_per_motif.tsv"
    _write_tsv(all_path, list(all_rows[0]), all_rows)
    _write_tsv(best_path, list(best_rows[0]), best_rows)
    outputs = {path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in (all_path, best_path)}
    manifest = {
        "purpose": "TJU Pan motif matrix vs idealized -10/-35/UP reference similarity; not genomic locus overlap",
        "status": "complete" if selected_species == list(SPECIES) else "partial_validation_subset",
        "expected_species": list(SPECIES),
        "included_species": selected_species,
        "tool": "MEME Suite Tomtom 5.5.9",
        "parameters": ["-dist ed", "-min-overlap 6", "-thresh 1", "target reverse-complement search enabled"],
        "reference_consensus_validation": reference_consensus,
        "reference_motif_ids": sorted(REFERENCE_IDS),
        "matrix_summary_manifest_sha256": sha256_file(matrix_manifest_path),
        "motif_matrix_sha256": {
            slug: sha256_file(matrix_summary_dir / "matrices" / f"{slug}.meme") for slug in selected_species
        },
        "tomtom_output_sha256": tomtom_hashes,
        "pair_count": len(all_rows),
        "statistical_limit": "Only four target motifs and eight strand-aware target comparisons per query; Tomtom p/q values are unreliable and are descriptive ranking fields only.",
        "biological_limit": "The -10/-35 templates are E. coli sigma70 consensus boxes; matrix similarity does not establish genomic position, co-occurrence, activity, or sigma specificity.",
        "output_sha256": outputs,
    }
    manifest_path = output_dir / "comparison_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tomtom-root", type=Path, required=True)
    parser.add_argument("--matrix-summary-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--species", nargs="+", help="默认汇总六物种；调试时可指定完成的部分物种")
    args = parser.parse_args()
    manifest = summarize(args.tomtom_root, args.matrix_summary_dir, args.reference, args.output_dir, args.species)
    print(f"summarized {manifest['pair_count']} Tomtom pairs for {len(manifest['included_species'])} species")


if __name__ == "__main__":
    main()
