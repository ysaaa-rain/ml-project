"""Export sequence-free MEME motif matrices and audit-ready summary tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES
from experiments.render_meme_logos import LETTERS, parse_meme_xml, sha256_file


def _background_frequencies(path: Path) -> dict[str, float]:
    root = ET.parse(path).getroot()
    matrix = root.find(".//background_frequencies/alphabet_array")
    if matrix is None:
        raise ValueError(f"MEME XML lacks background frequencies: {path}")
    frequencies = {node.attrib["letter_id"]: float(node.text or "nan") for node in matrix.findall("./value")}
    if set(frequencies) != set(LETTERS) or any(not math.isfinite(value) or value < 0 for value in frequencies.values()):
        raise ValueError(f"invalid background frequencies in {path}: {frequencies}")
    if not math.isclose(sum(frequencies.values()), 1.0, rel_tol=0, abs_tol=2e-5):
        raise ValueError(f"background frequencies do not sum to one in {path}: {frequencies}")
    return frequencies


def _max_high_confidence_run(columns: list[dict[str, float]], threshold: float = 0.8) -> int:
    best = current = 0
    previous: str | None = None
    for column in columns:
        letter, probability = max(column.items(), key=lambda item: item[1])
        if probability >= threshold and letter == previous:
            current += 1
        elif probability >= threshold:
            current = 1
        else:
            current = 0
        best = max(best, current)
        previous = letter if probability >= threshold else None
    return best


def _mean_information(columns: list[dict[str, float]]) -> float:
    values = []
    for column in columns:
        entropy = -sum(value * math.log2(value) for value in column.values() if value > 0)
        values.append(max(0.0, 2.0 - entropy))
    return sum(values) / len(values)


def _write_matrix(path: Path, motifs: list[dict[str, Any]], background: dict[str, float], species: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\n")
        handle.write("Background letter frequencies\n")
        handle.write(" ".join(f"{letter} {background[letter]:.8g}" for letter in LETTERS) + "\n\n")
        for motif in motifs:
            handle.write(f"MOTIF {species}_{motif['id']} {motif['consensus']}\n")
            handle.write(
                f"letter-probability matrix: alength= 4 w= {motif['width']} "
                f"nsites= {motif['sites']} E= {motif['e_value']:.8g}\n"
            )
            for column in motif["columns"]:
                handle.write(" ".join(f"{column[letter]:.8g}" for letter in LETTERS) + "\n")
            handle.write("\n")


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = [
        "species", "positive_count", "control_count", "motif_id", "consensus", "width_bp",
        "sites_reported", "p_value", "e_value", "mean_information_bits", "mean_gc_probability",
        "max_same_high_probability_run",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def export(
    meme_root: Path,
    run_manifest_path: Path,
    output_dir: Path,
    species: list[str] | None = None,
) -> dict[str, Any]:
    run = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    selected_species = species or list(SPECIES)
    if len(set(selected_species)) != len(selected_species) or any(slug not in SPECIES for slug in selected_species):
        raise ValueError(f"invalid or repeated species list: {selected_species}")
    rows: list[dict[str, Any]] = []
    source_hashes: dict[str, str] = {}
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_dir = output_dir / "matrices"
    matrix_dir.mkdir(parents=True, exist_ok=True)

    for slug in selected_species:
        record = run.get("results", {}).get(slug)
        if not record or record.get("status") != "success" or record.get("return_code") != 0:
            raise ValueError(f"{slug} is not a successful MEME run in {run_manifest_path}")
        xml_path = meme_root / slug / "meme.xml"
        text_path = meme_root / slug / "meme.txt"
        if not xml_path.is_file() or not text_path.is_file():
            raise FileNotFoundError(f"missing MEME XML/text output for {slug}")
        recorded_outputs = record.get("output_files", {})
        for name, path in (("meme.xml", xml_path), ("meme.txt", text_path)):
            expected_hash = recorded_outputs.get(name, {}).get("sha256")
            if expected_hash != sha256_file(path):
                raise ValueError(f"{slug}: {name} hash does not match the successful run manifest")
        motifs = parse_meme_xml(xml_path)
        text_motif_count = sum(line.startswith("MOTIF ") for line in text_path.read_text(encoding="utf-8").splitlines())
        if len(motifs) != text_motif_count:
            raise ValueError(f"{slug}: XML motif count {len(motifs)} != text motif count {text_motif_count}")
        expected_count = int(record["positive_count"])
        control_count = int(record["control_count"])
        if any(motif["sites"] > expected_count for motif in motifs):
            raise ValueError(f"{slug}: a motif reports more sites than positive input sequences")
        background = _background_frequencies(xml_path)
        matrix_path = matrix_dir / f"{slug}.meme"
        _write_matrix(matrix_path, motifs, background, slug)
        source_hashes[f"{slug}/meme.xml"] = sha256_file(xml_path)
        source_hashes[f"{slug}/meme.txt"] = sha256_file(text_path)
        source_hashes[f"{slug}/positive_fasta"] = record["positive_sha256"]
        source_hashes[f"{slug}/control_fasta"] = record["control_sha256"]

        for motif in motifs:
            mean_gc = sum(column["G"] + column["C"] for column in motif["columns"]) / motif["width"]
            rows.append({
                "species": slug,
                "positive_count": expected_count,
                "control_count": control_count,
                "motif_id": motif["id"],
                "consensus": motif["consensus"],
                "width_bp": motif["width"],
                "sites_reported": motif["sites"],
                "p_value": motif["p_value"],
                "e_value": motif["e_value"],
                "mean_information_bits": round(_mean_information(motif["columns"]), 6),
                "mean_gc_probability": round(mean_gc, 6),
                "max_same_high_probability_run": _max_high_confidence_run(motif["columns"]),
            })

    summary_path = output_dir / "motif_summary.tsv"
    _write_tsv(summary_path, rows)
    outputs = [summary_path, *sorted(matrix_dir.glob("*.meme"))]
    manifest = {
        "purpose": "六物种 MEME motif 概要及概率矩阵导出；不含逐条序列或样本 ID",
        "status": "complete" if selected_species == list(SPECIES) else "partial_validation_subset",
        "expected_species": list(SPECIES),
        "included_species": selected_species,
        "dataset": "TJU Pan course cloud drive / reg_and_gen/Datasets",
        "meme_version": run.get("meme_version"),
        "parameters": run.get("parameters"),
        "input_manifest_sha256": run.get("input_manifest_sha256"),
        "species_motif_counts": {
            slug: sum(row["species"] == slug for row in rows) for slug in selected_species
        },
        "e_value_definition": "MEME-reported statistic; a small E-value is a candidate-ranking measure under this search design, not a probability of biological function",
        "sites_limit": "sites_reported is MEME's site count and is not the number of unique promoter sequences",
        "composition_limit": "mean GC probability and high-probability run are descriptive QC fields, not low-complexity or functional significance tests",
        "source_sha256": source_hashes,
        "output_sha256": {
            path.relative_to(output_dir).as_posix(): {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in outputs
        },
    }
    manifest_path = output_dir / "summary_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meme-root", type=Path, required=True)
    parser.add_argument("--run-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--species", nargs="+", help="默认导出六个物种；仅用于验收子集时显式指定")
    args = parser.parse_args()
    result = export(args.meme_root, args.run_manifest, args.output_dir, args.species)
    print(f"exported {sum(result['species_motif_counts'].values())} motifs for {len(result['species_motif_counts'])} species")


if __name__ == "__main__":
    main()
