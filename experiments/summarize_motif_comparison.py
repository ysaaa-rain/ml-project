"""Summarize MEME-to-reference motif comparisons without treating them as loci."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


GROUPS = ("all", "Sigma70", "Sigma24", "Sigma32", "Sigma38", "Sigma28", "Sigma54")
REFERENCE_CONSENSUS = {
    "Ecoli_sigma70_minus10": "TATAAT",
    "Ecoli_sigma70_minus35": "TTGACA",
}
SEQUENCE_COUNTS = {
    "all": 3807,
    "Sigma70": 1942,
    "Sigma24": 513,
    "Sigma32": 300,
    "Sigma38": 176,
    "Sigma28": 140,
    "Sigma54": 96,
}
MOTIF_LINE = re.compile(
    r"^MOTIF\s+(?P<consensus>\S+)\s+(?P<meme_id>MEME-\d+)\s+"
    r"width\s*=\s*(?P<width>\d+)\s+sites\s*=\s*(?P<sites>\d+)"
    r".*?E-value\s*=\s*(?P<evalue>[0-9.eE+-]+)"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_discovered_motifs(path: Path) -> dict[str, dict[str, Any]]:
    motifs: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            match = MOTIF_LINE.match(line.strip())
            if not match:
                continue
            item = match.groupdict()
            consensus = item["consensus"]
            if consensus in motifs:
                raise ValueError(f"重复 consensus，无法唯一映射 Tomtom Query_ID：{path}: {consensus}")
            motifs[consensus] = {
                "meme_id": item["meme_id"],
                "width": int(item["width"]),
                "sites": int(item["sites"]),
                "meme_evalue": float(item["evalue"]),
            }
    if not motifs:
        raise ValueError(f"未能从 MEME 文本解析 motif：{path}")
    return motifs


def read_tomtom(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        content = (line for line in handle if not line.lstrip().startswith("#"))
        rows = [row for row in csv.DictReader(content, delimiter="\t") if row.get("Query_ID")]
    if not rows:
        raise ValueError(f"Tomtom TSV 没有匹配行：{path}")
    return rows


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def motif_number(meme_id: str) -> int:
    return int(meme_id.rsplit("-", 1)[1])


def validate_reference_consensus(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    motifs: dict[str, list[list[float]]] = {}
    index = 0
    while index < len(lines):
        if not lines[index].startswith("MOTIF "):
            index += 1
            continue
        name = lines[index].split()[1]
        index += 1
        while index < len(lines) and not lines[index].startswith("letter-probability matrix:"):
            index += 1
        if index == len(lines):
            raise ValueError(f"reference motif has no matrix: {name}")
        width = int(lines[index].split("w=", 1)[1].split()[0])
        index += 1
        matrix = [[float(value) for value in lines[index + offset].split()] for offset in range(width)]
        motifs[name] = matrix
        index += width

    observed: dict[str, str] = {}
    for name, expected in REFERENCE_CONSENSUS.items():
        if name not in motifs:
            raise ValueError(f"reference database missing required motif: {name}")
        letters = []
        for row in motifs[name]:
            maximum = max(row)
            winners = [letter for letter, value in zip("ACGT", row) if value == maximum]
            letters.append(winners[0] if len(winners) == 1 else "N")
        observed[name] = "".join(letters)
        if observed[name] != expected:
            raise ValueError(f"{name} expected consensus {expected}, got {observed[name]}")
    return observed


def summarize(
    meme_root: Path,
    tomtom_root: Path,
    input_dir: Path,
    reference: Path,
    output_dir: Path,
) -> dict[str, Any]:
    reference_consensus = validate_reference_consensus(reference)
    pair_columns = [
        "source_group", "source_sequence_count", "query_consensus", "query_meme_id",
        "query_width", "query_sites", "query_meme_evalue", "target_id", "optimal_offset",
        "tomtom_p_value", "tomtom_e_value", "tomtom_q_value", "overlap",
        "query_aligned_consensus", "target_aligned_consensus", "orientation",
    ]
    pair_rows: list[dict[str, Any]] = []
    best_rows: list[dict[str, Any]] = []
    input_hashes: dict[str, str] = {"reference_meme": sha256(reference)}
    intermediate_hashes: dict[str, str] = {}
    group_counts: dict[str, dict[str, int]] = {}

    for group in GROUPS:
        meme_file = meme_root / group / "meme.txt"
        tomtom_file = tomtom_root / group / "tomtom.tsv"
        motifs = parse_discovered_motifs(meme_file)
        rows = read_tomtom(tomtom_file)
        by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_query[row["Query_ID"]].append(row)
        if set(by_query) != set(motifs):
            missing = sorted(set(motifs) - set(by_query))
            unexpected = sorted(set(by_query) - set(motifs))
            raise ValueError(f"{group}: MEME/Tomtom motif 对不上；缺少={missing}，多出={unexpected}")

        expected_pairs = len(motifs) * 4
        if len(rows) != expected_pairs:
            raise ValueError(f"{group}: 预期 {expected_pairs} 个 query-target 结果，实际 {len(rows)}")
        group_counts[group] = {"sequences": SEQUENCE_COUNTS[group], "motifs": len(motifs), "matches": len(rows)}
        primary_fasta = input_dir / f"{group}.fasta"
        control_fasta = input_dir / f"{group}_n1.fasta"
        xml_file = meme_root / group / "meme.xml"
        for label, path in ((f"{group}_primary_fasta", primary_fasta), (f"{group}_N1_fasta", control_fasta)):
            if not path.is_file():
                raise FileNotFoundError(f"缺少运行输入，无法写完整 manifest：{path}")
            input_hashes[label] = sha256(path)
        for label, path in (
            (f"{group}_meme_txt", meme_file),
            (f"{group}_meme_xml", xml_file),
            (f"{group}_tomtom_tsv", tomtom_file),
        ):
            if not path.is_file():
                raise FileNotFoundError(f"缺少运行产物，无法写完整 manifest：{path}")
            intermediate_hashes[label] = sha256(path)

        for consensus, matches in by_query.items():
            info = motifs[consensus]
            ordered = sorted(matches, key=lambda row: (float(row["p-value"]), float(row["q-value"])))
            for row in ordered:
                pair_rows.append({
                    "source_group": group,
                    "source_sequence_count": SEQUENCE_COUNTS[group],
                    "query_consensus": consensus,
                    "query_meme_id": info["meme_id"],
                    "query_width": info["width"],
                    "query_sites": info["sites"],
                    "query_meme_evalue": info["meme_evalue"],
                    "target_id": row["Target_ID"],
                    "optimal_offset": row["Optimal_offset"],
                    "tomtom_p_value": row["p-value"],
                    "tomtom_e_value": row["E-value"],
                    "tomtom_q_value": row["q-value"],
                    "overlap": row["Overlap"],
                    "query_aligned_consensus": row["Query_consensus"],
                    "target_aligned_consensus": row["Target_consensus"],
                    "orientation": row["Orientation"],
                })
            best = ordered[0]
            best_rows.append({
                "source_group": group,
                "source_sequence_count": SEQUENCE_COUNTS[group],
                "query_consensus": consensus,
                "query_meme_id": info["meme_id"],
                "query_width": info["width"],
                "query_sites": info["sites"],
                "query_meme_evalue": info["meme_evalue"],
                "best_reference_match": best["Target_ID"],
                "tomtom_p_value_ranking_only": best["p-value"],
                "tomtom_q_value_unreliable_small_target_db": best["q-value"],
                "overlap": best["Overlap"],
                "alignment_offset": best["Optimal_offset"],
                "orientation": best["Orientation"],
                "query_aligned_consensus": best["Query_consensus"],
                "target_aligned_consensus": best["Target_consensus"],
            })

    pair_rows.sort(key=lambda row: (
        GROUPS.index(row["source_group"]), motif_number(row["query_meme_id"]),
        float(row["tomtom_p_value"]),
    ))
    best_rows.sort(key=lambda row: (
        GROUPS.index(row["source_group"]), motif_number(row["query_meme_id"]),
    ))
    all_matches = output_dir / "all_query_reference_matches.tsv"
    best_matches = output_dir / "best_reference_match_per_motif.tsv"
    write_tsv(all_matches, pair_columns, pair_rows)
    write_tsv(best_matches, list(best_rows[0]), best_rows)

    output_hashes = {
        all_matches.name: sha256(all_matches),
        best_matches.name: sha256(best_matches),
    }
    matrices_dir = output_dir / "matrices"
    if matrices_dir.is_dir():
        output_hashes.update({
            f"matrices/{path.name}": sha256(path)
            for path in sorted(matrices_dir.glob("*.meme"))
        })

    manifest: dict[str, Any] = {
        "purpose": "MEME motifs vs idealized literature consensus sequence similarity; not genomic locus overlap",
        "comparison_revision": "2026-10-03: reference matrix consensus audit and full Tomtom rerun",
        "tool": "MEME Suite Tomtom 5.5.9",
        "meme_run": {
            "tool": "MEME Suite MEME 5.5.9",
            "exit_status": "all seven formal jobs exited 0; each XML parsed with xmllint",
            "parameters": ["-dna", "-mod zoops", "-objfun de", "-neg matched N1", "-nmotifs 10", "-minw 6", "-maxw 20", "-revcomp", "-seed 20260930", "-brief 50"],
            "width_search": "MEME default heuristic seed widths; -allw was not used in final jobs",
            "all_group_command": "meme all.fasta -dna -mod zoops -objfun de -neg all_n1.fasta -nmotifs 10 -minw 6 -maxw 20 -revcomp -seed 20260930 -brief 50 -oc all",
            "sigma_group_command": "meme <group>.fasta -dna -mod zoops -objfun de -neg <group>_n1.fasta -nmotifs 10 -minw 6 -maxw 20 -searchsize 100000 -revcomp -seed 20260930 -brief 50 -nostatus -oc <group>",
        },
        "distance": "Euclidean distance (-dist ed)",
        "minimum_overlap": 6,
        "tomtom_threshold": 1.0,
        "orientation": "Tomtom default: target reverse complements are also compared",
        "tomtom_command": "tomtom -dist ed -thresh 1 -min-overlap 6 -oc <output_dir> <group>/meme.txt baselines/known_promoter_elements.meme",
        "reference_motifs": [
            "Ecoli_sigma70_minus10", "Ecoli_sigma70_minus35",
            "Ecoli_UP_proximal", "Ecoli_UP_distal",
        ],
        "validated_reference_consensus": reference_consensus,
        "group_counts": group_counts,
        "total_discovered_motifs": sum(item["motifs"] for item in group_counts.values()),
        "total_query_reference_rows": len(pair_rows),
        "statistical_limit": (
            "The target database contains only four idealized motifs. Tomtom warns that p-values are inaccurate "
            "for such a small target database; reported p/q values are retained for audit but are not used as "
            "confirmatory significance. Best matches are descriptive rankings only."
        ),
        "biological_limit": (
            "These references are E. coli sigma70 -10/-35 consensus and UP subsite consensus templates. "
            "Sequence similarity does not establish promoter-specific presence, TSS-relative position, "
            "functional activity, or sigma-specificity."
        ),
        "input_sha256": input_hashes,
        "intermediate_sha256": intermediate_hashes,
        "output_sha256": output_hashes,
    }
    manifest_path = output_dir / "comparison_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meme-root", type=Path, required=True, help="含 all/Sigma70/.../meme.txt 的目录")
    parser.add_argument("--tomtom-root", type=Path, required=True, help="含 all/Sigma70/.../tomtom.tsv 的目录")
    parser.add_argument("--input-dir", type=Path, required=True, help="含 discovery 与 N1 FASTA 的输入目录")
    parser.add_argument("--reference", type=Path, required=True, help="理想化已知元件 MEME 格式文件")
    parser.add_argument("--output-dir", type=Path, required=True, help="汇总结果目录")
    args = parser.parse_args()
    manifest = summarize(args.meme_root, args.tomtom_root, args.input_dir, args.reference, args.output_dir)
    print(
        f"summarized {manifest['total_discovered_motifs']} motifs and "
        f"{manifest['total_query_reference_rows']} query-reference comparisons "
        f"to {args.output_dir}"
    )


if __name__ == "__main__":
    main()
