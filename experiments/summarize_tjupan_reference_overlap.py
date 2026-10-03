"""Summarize reference-element hits and overlap with de novo FIMO sites."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES
from experiments.summarize_tjupan_fimo import benjamini_hochberg, fisher_exact_two_sided, _read_hits
from preprocessing.provenance import sha256_file


REFERENCES = (
    "Ecoli_sigma70_minus10",
    "Ecoli_sigma70_minus35",
    "Ecoli_UP_proximal",
    "Ecoli_UP_distal",
)


def _read_fasta(path: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    identifier: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if identifier is None:
            return
        sequence = "".join(chunks).upper()
        label = identifier.rsplit("|label=", 1)[-1]
        if label not in {"0", "1"} or len(sequence) != 81 or set(sequence) - set("ACGT"):
            raise ValueError(f"invalid labeled 81-bp sequence in {path}: {identifier}")
        if identifier in result:
            raise ValueError(f"duplicate FASTA identifier: {identifier}")
        result[identifier] = {
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
    return result


def _read_reference_hits(path: Path, allowed_ids: set[str], q_threshold: float) -> list[dict[str, Any]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter="\t")
    required = {"motif_id", "sequence_name", "start", "stop", "strand", "q-value"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise ValueError(f"FIMO output has missing columns: {path}")
    hits: list[dict[str, Any]] = []
    for row in reader:
        start, stop = int(row["start"]), int(row["stop"])
        q_value = float(row["q-value"])
        if row["motif_id"] not in allowed_ids:
            raise ValueError(f"unexpected reference motif in {path}: {row['motif_id']}")
        if start < 1 or stop < start or stop > 81 or row["strand"] not in {"+", "-"}:
            raise ValueError(f"invalid reference hit: {row}")
        if not 0 <= q_value <= 1:
            raise ValueError(f"invalid q value in reference hit: {row}")
        if q_value <= q_threshold:
            hits.append({
                "motif_id": row["motif_id"], "sequence_name": row["sequence_name"],
                "start": start, "stop": stop, "strand": row["strand"], "q_value": q_value,
            })
    return hits


def intervals_overlap(first: dict[str, Any], second: dict[str, Any]) -> bool:
    """Return whether two 1-based inclusive site intervals share a base."""
    return max(int(first["start"]), int(second["start"])) <= min(int(first["stop"]), int(second["stop"]))


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty table without a schema: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(
    *,
    reference_run_manifest: Path,
    denovo_run_manifest: Path,
    input_manifest: Path,
    meme_summary: Path,
    fimo_summary_manifest: Path,
    run_root: Path,
    denovo_run_root: Path,
    output_dir: Path,
    q_threshold: float = 0.05,
) -> dict[str, Any]:
    reference_run = json.loads(reference_run_manifest.read_text(encoding="utf-8"))
    denovo_run = json.loads(denovo_run_manifest.read_text(encoding="utf-8"))
    json.loads(input_manifest.read_text(encoding="utf-8"))
    fimo_summary = json.loads(fimo_summary_manifest.read_text(encoding="utf-8"))
    if reference_run.get("status") != "complete" or denovo_run.get("status") != "complete":
        raise ValueError("both six-species reference and de novo FIMO runs must be complete")
    if fimo_summary.get("status") != "complete":
        raise ValueError("the sequence-free de novo FIMO summary must be complete")
    if reference_run.get("de_novo_fimo_manifest_sha256") != sha256_file(denovo_run_manifest):
        raise ValueError("reference scan was not run against the supplied de novo FIMO manifest")
    if denovo_run.get("input_manifest_sha256") != sha256_file(input_manifest):
        raise ValueError("supplied input manifest does not match the de novo FIMO run")
    if fimo_summary.get("meme_summary_sha256") != sha256_file(meme_summary):
        raise ValueError("supplied MEME summary does not match the FIMO summary")
    if set(reference_run.get("results", {})) != set(SPECIES) or set(denovo_run.get("results", {})) != set(SPECIES):
        raise ValueError("expected all six species in both FIMO run manifests")
    if not 0 < q_threshold <= 1:
        raise ValueError("q threshold must be in (0, 1]")

    with meme_summary.open(encoding="utf-8", newline="") as handle:
        meme_rows = list(csv.DictReader(handle, delimiter="\t"))
    motif_rows: dict[str, dict[str, dict[str, str]]] = {slug: {} for slug in SPECIES}
    for row in meme_rows:
        slug = row["species"]
        if slug not in motif_rows:
            raise ValueError(f"unexpected species in MEME summary: {slug}")
        motif_rows[slug][row["motif_id"]] = row
    if any(len(motif_rows[slug]) != 10 for slug in SPECIES):
        raise ValueError("expected 10 de novo motifs per species")

    known_hits_by_species: dict[str, dict[str, list[dict[str, Any]]]] = {}
    denovo_hits_by_species: dict[str, dict[str, list[dict[str, Any]]]] = {}
    sequences_by_species: dict[str, dict[str, dict[str, str]]] = {}
    input_hashes: dict[str, str] = {}
    for slug in SPECIES:
        known_path = run_root / slug / "fimo.tsv"
        denovo_path = denovo_run_root / slug / "fimo.tsv"
        fasta_path = Path(denovo_run["results"][slug]["combined_fasta"])
        for label, path, manifest_hash in (
            (f"{slug}/reference_fimo", known_path, reference_run["results"][slug]["fimo_tsv_sha256"]),
            (f"{slug}/denovo_fimo", denovo_path, denovo_run["results"][slug]["fimo_tsv_sha256"]),
            (f"{slug}/combined_fasta", fasta_path, denovo_run["results"][slug]["combined_fasta_sha256"]),
        ):
            if not path.is_file() or sha256_file(path) != manifest_hash:
                raise ValueError(f"missing or changed input {label}: {path}")
            input_hashes[label] = sha256_file(path)
        sequences = _read_fasta(fasta_path)
        known_hits = _read_reference_hits(known_path, set(REFERENCES), q_threshold)
        denovo_hits = _read_hits(denovo_path, slug)
        known_hits_by_species[slug] = {name: [] for name in REFERENCES}
        for hit in known_hits:
            if hit["sequence_name"] not in sequences:
                raise ValueError(f"reference FIMO hit references unknown sequence: {hit['sequence_name']}")
            known_hits_by_species[slug][hit["motif_id"]].append(hit)
        denovo_hits_by_species[slug] = {motif_id: [] for motif_id in motif_rows[slug]}
        for hit in denovo_hits:
            if hit["sequence_name"] not in sequences:
                raise ValueError(f"de novo FIMO hit references unknown sequence: {hit['sequence_name']}")
            if hit["q_value"] <= q_threshold:
                denovo_hits_by_species[slug][hit["motif_id"]].append(hit)
        sequences_by_species[slug] = sequences

    association_rows: list[dict[str, Any]] = []
    raw_p_values: list[float] = []
    for slug in SPECIES:
        sequences = sequences_by_species[slug]
        unique_by_label = {
            label: {record["sequence_sha256"] for record in sequences.values() if record["label"] == label}
            for label in ("1", "0")
        }
        shared = unique_by_label["1"] & unique_by_label["0"]
        for reference in REFERENCES:
            hit_ids = {
                label: {
                    sequences[hit["sequence_name"]]["sequence_sha256"]
                    for hit in known_hits_by_species[slug][reference]
                    if sequences[hit["sequence_name"]]["label"] == label
                }
                for label in ("1", "0")
            }
            pos_hits, ctl_hits = hit_ids["1"] - shared, hit_ids["0"] - shared
            pos_total, ctl_total = unique_by_label["1"] - shared, unique_by_label["0"] - shared
            p_value = fisher_exact_two_sided(
                len(pos_hits), len(pos_total) - len(pos_hits),
                len(ctl_hits), len(ctl_total) - len(ctl_hits),
            )
            raw_p_values.append(p_value)
            positive_coverage = len(hit_ids["1"]) / len(unique_by_label["1"])
            control_coverage = len(hit_ids["0"]) / len(unique_by_label["0"])
            association_rows.append({
                "species": slug,
                "reference_motif": reference,
                "q_threshold": q_threshold,
                "positive_reference_site_count": sum(
                    1 for hit in known_hits_by_species[slug][reference]
                    if sequences[hit["sequence_name"]]["label"] == "1"
                ),
                "positive_hit_unique_sequences": len(hit_ids["1"]),
                "positive_unique_sequences": len(unique_by_label["1"]),
                "positive_unique_sequence_coverage": positive_coverage,
                "control_reference_site_count": sum(
                    1 for hit in known_hits_by_species[slug][reference]
                    if sequences[hit["sequence_name"]]["label"] == "0"
                ),
                "control_hit_unique_sequences": len(hit_ids["0"]),
                "control_unique_sequences": len(unique_by_label["0"]),
                "control_unique_sequence_coverage": control_coverage,
                "coverage_difference_positive_minus_control": positive_coverage - control_coverage,
                "cross_label_shared_sequences_excluded_from_fisher": len(shared),
                "fisher_two_sided_p_value": p_value,
                "fisher_bh_q_value_24_tests": None,
            })
    adjusted = benjamini_hochberg(raw_p_values)
    for row, q_value in zip(association_rows, adjusted, strict=True):
        row["fisher_bh_q_value_24_tests"] = q_value

    overlap_rows: list[dict[str, Any]] = []
    for slug in SPECIES:
        sequences = sequences_by_species[slug]
        for motif_id in sorted(denovo_hits_by_species[slug], key=lambda name: int(name.rsplit("-", 1)[1])):
            sites = denovo_hits_by_species[slug][motif_id]
            for reference in REFERENCES:
                reference_by_sequence: dict[str, list[dict[str, Any]]] = {}
                for hit in known_hits_by_species[slug][reference]:
                    reference_by_sequence.setdefault(hit["sequence_name"], []).append(hit)
                overlapped: list[dict[str, Any]] = []
                same_strand: list[dict[str, Any]] = []
                overlap_sequences: set[str] = set()
                for site in sites:
                    candidates = reference_by_sequence.get(site["sequence_name"], [])
                    matches = [ref for ref in candidates if intervals_overlap(site, ref)]
                    if matches:
                        overlapped.append(site)
                        overlap_sequences.add(site["sequence_name"])
                        if any(site["strand"] == ref["strand"] for ref in matches):
                            same_strand.append(site)
                overlap_rows.append({
                    "species": slug,
                    "de_novo_motif_id": motif_id,
                    "de_novo_consensus": motif_rows[slug][motif_id]["consensus"],
                    "de_novo_meme_e_value": motif_rows[slug][motif_id]["e_value"],
                    "reference_motif": reference,
                    "q_threshold_each_motif": q_threshold,
                    "de_novo_significant_site_count": len(sites),
                    "overlapping_de_novo_site_count_any_strand": len(overlapped),
                    "overlapping_de_novo_site_fraction": len(overlapped) / len(sites) if sites else 0.0,
                    "overlapping_de_novo_site_count_same_input_strand": len(same_strand),
                    "overlapping_unique_record_count": len(overlap_sequences),
                })

    output_dir.mkdir(parents=True, exist_ok=True)
    association_path = output_dir / "known_reference_coverage_association.tsv"
    overlap_path = output_dir / "denovo_known_reference_overlap.tsv"
    _write_tsv(association_path, association_rows)
    _write_tsv(overlap_path, overlap_rows)
    output_manifest = {
        "schema_version": "tjupan-known-reference-overlap-summary-1.0",
        "project": "PR01-02",
        "status": "complete",
        "purpose": "aggregate known-reference motif coverage and coordinate overlap with de novo FIMO hits; no record identifiers or sequences are exported",
        "reference_source": "four idealized E. coli sigma70 -10/-35 and UP templates; not species-specific gold standards",
        "fimo_version": reference_run["fimo_version"],
        "q_threshold": q_threshold,
        "q_value_scope": "each reference motif and each de novo motif is thresholded separately by FIMO",
        "overlap_definition": "same input FASTA record and inclusive 1-based intervals share at least one base; strand ignored for any-strand overlap, same input-strand reported separately",
        "association_test": "two-sided Fisher exact test on unique exact sequences, with cross-label shared sequences excluded; Benjamini-Hochberg across 24 species-by-reference comparisons",
        "interpretation_limits": [
            "reference motifs are idealized E. coli templates; matches in other species are exploratory and do not establish orthology or function",
            "de novo motifs were discovered from the same positive sequences later scanned, so overlap is not independent validation",
            "negative controls are source label=0, not experimentally verified non-promoters",
            "input orientation is not independently verified; coordinates remain 1-81 and strand is relative to the input string",
            "coordinate overlap indicates shared bases, not biological interaction or functional equivalence",
        ],
        "species": list(SPECIES),
        "reference_motif_count": len(REFERENCES),
        "de_novo_species_motif_count": sum(len(value) for value in motif_rows.values()),
        "fisher_tests": len(association_rows),
        "fisher_bh_q_le_0_05_count": sum(row["fisher_bh_q_value_24_tests"] <= 0.05 for row in association_rows),
        "source_sha256": {
            "reference_run_manifest": sha256_file(reference_run_manifest),
            "de_novo_run_manifest": sha256_file(denovo_run_manifest),
            "input_manifest": sha256_file(input_manifest),
            "meme_summary": sha256_file(meme_summary),
            "de_novo_fimo_summary_manifest": sha256_file(fimo_summary_manifest),
            "combined_fasta_and_fimo_inputs": input_hashes,
        },
        "output_sha256": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in (association_path, overlap_path)
        },
    }
    manifest_path = output_dir / "reference_overlap_summary_manifest.json"
    manifest_path.write_text(json.dumps(output_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return output_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-run-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/reference_elements/reference_fimo_run_manifest.json"))
    parser.add_argument("--denovo-run-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo/fimo_run_manifest.json"))
    parser.add_argument("--input-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/input_manifest.json"))
    parser.add_argument("--meme-summary", type=Path, default=Path("results/motif/tjupan_m3_20261003/summary/motif_summary.tsv"))
    parser.add_argument("--fimo-summary-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo/fimo_summary_manifest.json"))
    parser.add_argument("--run-root", type=Path, default=Path("tmp/fimo_runs/tjupan_promloop_20261003/reference_elements"))
    parser.add_argument("--denovo-run-root", type=Path, default=Path("tmp/fimo_runs/tjupan_promloop_20261003"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/reference_elements"))
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args()
    manifest = summarize(
        reference_run_manifest=args.reference_run_manifest,
        denovo_run_manifest=args.denovo_run_manifest,
        input_manifest=args.input_manifest,
        meme_summary=args.meme_summary,
        fimo_summary_manifest=args.fimo_summary_manifest,
        run_root=args.run_root,
        denovo_run_root=args.denovo_run_root,
        output_dir=args.output_dir,
        q_threshold=args.q_threshold,
    )
    print(f"status={manifest['status']}")
    print(f"Fisher BH q<=0.05: {manifest['fisher_bh_q_le_0_05_count']}/{manifest['fisher_tests']}")


if __name__ == "__main__":
    main()
