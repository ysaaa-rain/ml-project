"""Run species-wise FIMO scans on TJU Pan promoters and matched controls."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import json
from pathlib import Path
import re
import subprocess
import time
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES
from preprocessing.provenance import sha256_file


DNA = set("ACGT")
EXPECTED_SPECIES = list(SPECIES)


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _read_fasta(path: Path, expected_label: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if header is None:
            return
        sequence = "".join(chunks).upper()
        if not sequence or len(sequence) != 81 or not set(sequence) <= DNA:
            raise ValueError(f"invalid 81-bp A/C/G/T FASTA sequence in {path}: {header}")
        if f"|label={expected_label}" not in header:
            raise ValueError(f"unexpected label in {path}: {header}")
        records.append((header, sequence))

    for line_number, raw_line in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header = line[1:]
            if not header or any(character.isspace() for character in header):
                raise ValueError(f"invalid FASTA header at {path}:{line_number}")
            chunks = []
        elif header is None:
            raise ValueError(f"sequence before FASTA header at {path}:{line_number}")
        else:
            chunks.append(line)
    finish()
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    if len({name for name, _ in records}) != len(records):
        raise ValueError(f"duplicate FASTA identifiers in {path}")
    return records


def _prepare_scan_inputs(
    positive_fasta: Path,
    control_fasta: Path,
    combined_fasta: Path,
    background_file: Path,
) -> dict[str, Any]:
    positives = _read_fasta(positive_fasta, "1")
    controls = _read_fasta(control_fasta, "0")
    identifiers = [name for name, _ in positives + controls]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("positive/control FASTA identifiers overlap")

    combined_fasta.parent.mkdir(parents=True, exist_ok=True)
    with combined_fasta.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in positives + controls:
            handle.write(f">{name}\n{sequence}\n")

    counts = Counter(base for _, sequence in controls for base in sequence)
    total = sum(counts.values())
    if total == 0 or any(counts[base] == 0 for base in "ACGT"):
        raise ValueError(f"control sequences do not provide all four DNA bases: {counts}")
    frequencies = {base: counts[base] / total for base in "ACGT"}
    background_file.parent.mkdir(parents=True, exist_ok=True)
    background_file.write_text(
        "".join(f"{base} {frequencies[base]:.12g}\n" for base in "ACGT"),
        encoding="ascii",
    )
    positive_sequences = {sequence for _, sequence in positives}
    control_sequences = {sequence for _, sequence in controls}
    return {
        "positive_records": len(positives),
        "control_records": len(controls),
        "positive_unique_sequences": len(positive_sequences),
        "control_unique_sequences": len(control_sequences),
        "cross_label_exact_sequence_overlap": len(positive_sequences & control_sequences),
        "background_base_counts": dict(counts),
        "background_frequencies": frequencies,
    }


def _tool_version(fimo: Path) -> str:
    result = subprocess.run([str(fimo), "--version"], capture_output=True, text=True, check=True, timeout=15)
    return (result.stdout or result.stderr).strip().splitlines()[0]


def _motif_pi0_estimates(log_path: Path, species: str) -> dict[str, float]:
    current_motif: str | None = None
    estimates: dict[str, float] = {}
    motif_pattern = re.compile(rf"Using motif \+{re.escape(species)}_(MEME-\d+) of width")
    pi0_pattern = re.compile(r"Estimated pi_0=([0-9.eE+-]+)")
    for line in log_path.read_text(encoding="utf-8").splitlines():
        motif_match = motif_pattern.search(line)
        if motif_match:
            current_motif = motif_match.group(1)
            continue
        pi0_match = pi0_pattern.search(line)
        if pi0_match and current_motif:
            estimates[current_motif] = float(pi0_match.group(1))
    return estimates


def run(
    input_manifest_path: Path,
    summary_manifest_path: Path,
    input_dir: Path,
    matrix_dir: Path,
    fimo: Path,
    scan_input_dir: Path,
    run_root: Path,
    result_dir: Path,
    species_slugs: list[str],
    q_threshold: float = 0.05,
) -> dict[str, Any]:
    prepared = json.loads(input_manifest_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_manifest_path.read_text(encoding="utf-8"))
    if prepared.get("status") != "inputs_validated_and_written":
        raise ValueError("input manifest is not validated")
    if summary.get("status") != "complete" or summary.get("included_species") != EXPECTED_SPECIES:
        raise ValueError("complete six-species MEME summary is required before FIMO")
    if len(set(species_slugs)) != len(species_slugs) or any(slug not in SPECIES for slug in species_slugs):
        raise ValueError(f"invalid or repeated species list: {species_slugs}")
    if not 0 < q_threshold <= 1:
        raise ValueError(f"q-value threshold must be in (0, 1], got {q_threshold}")
    version = _tool_version(fimo)
    if version != "5.5.9":
        raise ValueError(f"expected MEME Suite FIMO 5.5.9, got {version}")

    manifest_path = result_dir / "fimo_run_manifest.json"
    parameters = {
        "qvalue_threshold": q_threshold,
        "threshold_type": "FIMO-estimated q-value; the local 5.5.9 implementation estimates pi0 from a uniformly sampled p-value reservoir",
        "qvalue_scope": "computed separately for each motif across that motif's tested positions, not pooled across the 10 motifs",
        "reverse_complement_search": "enabled by MEME matrix strands + -; FIMO --norc is not set",
        "motif_scope": "all 10 MEME motifs discovered for each species, scanned on that species' positive and control sequences",
        "background_model": "zero-order A/C/G/T frequencies estimated from that species' valid 81-bp negative controls",
        "max_stored_scores": 1000000,
        "coordinates": "1-based positions in the input FASTA; not converted to TSS-relative coordinates",
    }
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        prior_parameters = manifest.get("parameters", {})
        if manifest.get("fimo_version") != version or any(
            prior_parameters.get(key) != parameters.get(key)
            for key in ("qvalue_threshold", "max_stored_scores", "coordinates")
        ):
            raise ValueError("existing FIMO manifest uses a different software version or parameter set")
    else:
        manifest = {
            "schema_version": "tjupan-fimo-runs-1.0",
            "project": "PR01-02",
            "dataset": prepared["source"],
            "input_manifest": str(input_manifest_path),
            "input_manifest_sha256": sha256_file(input_manifest_path),
            "meme_summary_manifest": str(summary_manifest_path),
            "meme_summary_manifest_sha256": sha256_file(summary_manifest_path),
            "fimo_tool": str(fimo),
            "fimo_version": version,
            "parameters": parameters,
            "results": {},
        }
    manifest.update(
        {
            "input_manifest": str(input_manifest_path),
            "input_manifest_sha256": sha256_file(input_manifest_path),
            "meme_summary_manifest": str(summary_manifest_path),
            "meme_summary_manifest_sha256": sha256_file(summary_manifest_path),
            "parameters": parameters,
        }
    )
    result_dir.mkdir(parents=True, exist_ok=True)
    _write_json(manifest_path, manifest)

    for slug in species_slugs:
        species_info = prepared["species"][slug]
        positive = input_dir / f"{slug}.positive.fasta"
        control = input_dir / f"{slug}.control.fasta"
        expected_inputs = (
            (positive, species_info["positive_samples_csv"]["fasta_sha256"]),
            (control, species_info["negative_control"]["fasta_sha256"]),
        )
        for path, expected_hash in expected_inputs:
            if not path.is_file() or sha256_file(path) != expected_hash:
                raise ValueError(f"input FASTA is missing or differs from its manifest: {path}")

        matrix = matrix_dir / f"{slug}.meme"
        matrix_key = f"matrices/{slug}.meme"
        matrix_hash = summary.get("output_sha256", {}).get(matrix_key, {}).get("sha256")
        if not matrix.is_file() or not matrix_hash or sha256_file(matrix) != matrix_hash:
            raise ValueError(f"MEME matrix is missing or differs from its complete summary manifest: {matrix}")

        prior = manifest["results"].get(slug)
        if prior and prior.get("status") == "success":
            fimo_tsv = run_root / slug / "fimo.tsv"
            if not fimo_tsv.is_file() or sha256_file(fimo_tsv) != prior.get("fimo_tsv_sha256"):
                raise ValueError(f"successful FIMO output is missing or changed: {fimo_tsv}")
            log_path = Path(prior["log"])
            if log_path.is_file():
                prior["log_sha256"] = sha256_file(log_path)
                prior["motif_pi0_estimates"] = _motif_pi0_estimates(log_path, slug)
            print(f"SKIP {slug}: successful output already recorded", flush=True)
            continue
        attempt_history = list((prior or {}).get("attempt_history", []))
        if prior:
            attempt_history.append({key: value for key, value in prior.items() if key != "attempt_history"})

        combined = scan_input_dir / f"{slug}.combined.fasta"
        background = scan_input_dir / f"{slug}.control.background"
        counts = _prepare_scan_inputs(positive, control, combined, background)
        output_dir = run_root / slug
        if output_dir.exists() and any(output_dir.iterdir()):
            raise FileExistsError(f"refusing to overwrite nonempty output {output_dir}; use a fresh run root")
        log_path = run_root / "logs" / f"{slug}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        command = [
            str(fimo), "--oc", str(output_dir), "--bfile", str(background),
            "--qv-thresh", "--thresh", str(q_threshold), "--max-stored-scores", "1000000",
            str(matrix), str(combined),
        ]
        started_at = datetime.now().astimezone().isoformat(timespec="seconds")
        started_clock = time.monotonic()
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        elapsed = time.monotonic() - started_clock
        finished_at = datetime.now().astimezone().isoformat(timespec="seconds")
        log_path.write_text(
            f"COMMAND\n{json.dumps(command, ensure_ascii=False)}\n\nSTDOUT\n{process.stdout}\n\nSTDERR\n{process.stderr}",
            encoding="utf-8",
        )
        record: dict[str, Any] = {
            "status": "success" if process.returncode == 0 else "failed",
            "positive_count": counts["positive_records"],
            "control_count": counts["control_records"],
            "positive_unique_sequences": counts["positive_unique_sequences"],
            "control_unique_sequences": counts["control_unique_sequences"],
            "cross_label_exact_sequence_overlap": counts["cross_label_exact_sequence_overlap"],
            "combined_fasta": str(combined),
            "combined_fasta_sha256": sha256_file(combined),
            "background_file": str(background),
            "background_file_sha256": sha256_file(background),
            "background_base_counts": counts["background_base_counts"],
            "background_frequencies": counts["background_frequencies"],
            "matrix_file": str(matrix),
            "matrix_sha256": matrix_hash,
            "command": command,
            "output_directory": str(output_dir),
            "log": str(log_path),
            "started_at": started_at,
            "finished_at": finished_at,
            "elapsed_seconds": round(elapsed, 2),
            "return_code": process.returncode,
        }
        if attempt_history:
            record["attempt_history"] = attempt_history
        fimo_tsv = output_dir / "fimo.tsv"
        if process.returncode == 0:
            required_outputs = (fimo_tsv, output_dir / "fimo.xml", output_dir / "cisml.xml")
            missing = [str(path) for path in required_outputs if not path.is_file()]
            if missing:
                record["status"] = "failed"
                record["missing_outputs"] = missing
            else:
                record["fimo_tsv_sha256"] = sha256_file(fimo_tsv)
                record["fimo_tsv_bytes"] = fimo_tsv.stat().st_size
                record["log_sha256"] = sha256_file(log_path)
                record["motif_pi0_estimates"] = _motif_pi0_estimates(log_path, slug)
                if len(record["motif_pi0_estimates"]) != 10:
                    record["status"] = "failed"
                    record["qvalue_pi0_audit_error"] = (
                        "expected one pi0 estimate per motif; found "
                        f"{len(record['motif_pi0_estimates'])}"
                    )
                record["output_files"] = {
                    path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
                    for path in required_outputs
                }
        manifest["results"][slug] = record
        manifest["status"] = "complete" if all(
            manifest["results"].get(name, {}).get("status") == "success" for name in EXPECTED_SPECIES
        ) else "partial_validation_subset"
        _write_json(manifest_path, manifest)
        print(
            f"{record['status'].upper()} {slug}: rc={process.returncode}, "
            f"elapsed={elapsed:.2f}s, q≤{q_threshold}",
            flush=True,
        )
        if record["status"] != "success":
            raise RuntimeError(f"FIMO failed for {slug}; see {log_path} and {manifest_path}")

    manifest["status"] = "complete" if all(
        manifest["results"].get(name, {}).get("status") == "success" for name in EXPECTED_SPECIES
    ) else "partial_validation_subset"
    _write_json(manifest_path, manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/input_manifest.json"))
    parser.add_argument("--summary-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/summary/summary_manifest.json"))
    parser.add_argument("--input-dir", type=Path, default=Path("tmp/tjupan_motif_inputs_20261003"))
    parser.add_argument("--matrix-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/summary/matrices"))
    parser.add_argument("--fimo", type=Path, default=Path("tmp/meme-suite-5.5.9/bin/fimo"))
    parser.add_argument("--scan-input-dir", type=Path, default=Path("tmp/tjupan_fimo_inputs_20261003"))
    parser.add_argument("--run-root", type=Path, default=Path("tmp/fimo_runs/tjupan_promloop_20261003"))
    parser.add_argument("--result-dir", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo"))
    parser.add_argument("--species", nargs="*", default=EXPECTED_SPECIES)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args()
    manifest = run(
        args.input_manifest, args.summary_manifest, args.input_dir, args.matrix_dir, args.fimo,
        args.scan_input_dir, args.run_root, args.result_dir, args.species, args.q_threshold,
    )
    print(f"FIMO manifest status: {manifest['status']}")


if __name__ == "__main__":
    main()
