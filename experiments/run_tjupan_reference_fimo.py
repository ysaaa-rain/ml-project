"""Scan TJU Pan promoters with the four idealized known-element PWMs."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import subprocess
import time
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES
from preprocessing.provenance import sha256_file


EXPECTED_REFERENCES = {
    "Ecoli_sigma70_minus10",
    "Ecoli_sigma70_minus35",
    "Ecoli_UP_proximal",
    "Ecoli_UP_distal",
}


def _read_reference_names(path: Path) -> set[str]:
    names = {line.split()[1] for line in path.read_text(encoding="utf-8").splitlines() if line.startswith("MOTIF ")}
    if names != EXPECTED_REFERENCES:
        raise ValueError(f"expected exactly four audited reference motifs, got {sorted(names)}")
    return names


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def run(
    *,
    reference_meme: Path,
    fimo_run_manifest: Path,
    fimo: Path,
    run_root: Path,
    output_manifest: Path,
    q_threshold: float = 0.05,
) -> dict[str, Any]:
    source_manifest = json.loads(fimo_run_manifest.read_text(encoding="utf-8"))
    if source_manifest.get("status") != "complete" or set(source_manifest.get("results", {})) != set(SPECIES):
        raise ValueError("a complete six-species de novo FIMO run is required")
    if not 0 < q_threshold <= 1:
        raise ValueError(f"q threshold must be in (0, 1], got {q_threshold}")
    reference_names = sorted(_read_reference_names(reference_meme))
    version_result = subprocess.run([str(fimo), "--version"], capture_output=True, text=True, check=True, timeout=15)
    version = (version_result.stdout or version_result.stderr).strip().splitlines()[0]
    if version != "5.5.9":
        raise ValueError(f"expected MEME Suite FIMO 5.5.9, got {version}")

    parameters = {
        "reference_motifs": reference_names,
        "qvalue_threshold": q_threshold,
        "qvalue_scope": "FIMO computes q-values separately for each reference motif across its tested positions",
        "background_model": "same species-specific zero-order frequencies estimated from valid label=0 controls in the de novo FIMO run",
        "reverse_complement_search": "enabled by strands + - in the reference MEME file",
        "coordinates": "1-based inclusive positions in the source FASTA; not converted to TSS-relative coordinates",
        "max_stored_scores": 1000000,
        "purpose": "reference-site localization and coordinate overlap with de novo FIMO sites; not proof of conserved function",
    }
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    if output_manifest.is_file():
        manifest = json.loads(output_manifest.read_text(encoding="utf-8"))
        if (
            manifest.get("fimo_version") != version
            or manifest.get("parameters") != parameters
            or manifest.get("de_novo_fimo_manifest_sha256") != sha256_file(fimo_run_manifest)
            or manifest.get("reference_meme_sha256") != sha256_file(reference_meme)
        ):
            raise ValueError("existing reference-FIMO manifest uses different inputs or parameters")
    else:
        manifest = {
            "schema_version": "tjupan-known-reference-fimo-1.0",
            "project": "PR01-02",
            "dataset": source_manifest.get("dataset"),
            "status": "running",
            "fimo_version": version,
            "fimo_executable": str(fimo),
            "reference_meme": str(reference_meme),
            "reference_meme_sha256": sha256_file(reference_meme),
            "de_novo_fimo_manifest": str(fimo_run_manifest),
            "de_novo_fimo_manifest_sha256": sha256_file(fimo_run_manifest),
            "parameters": parameters,
            "results": {},
        }
    manifest["status"] = "running"
    _write_json(output_manifest, manifest)

    for slug in SPECIES:
        source = source_manifest["results"][slug]
        fasta = Path(source["combined_fasta"])
        background = Path(source["background_file"])
        if not fasta.is_file() or sha256_file(fasta) != source["combined_fasta_sha256"]:
            raise ValueError(f"combined FASTA is missing or changed: {fasta}")
        if not background.is_file() or sha256_file(background) != source["background_file_sha256"]:
            raise ValueError(f"control background is missing or changed: {background}")

        prior = manifest["results"].get(slug)
        output_dir = run_root / slug
        tsv_path = output_dir / "fimo.tsv"
        if prior and prior.get("status") == "success":
            if not tsv_path.is_file() or sha256_file(tsv_path) != prior.get("fimo_tsv_sha256"):
                raise ValueError(f"recorded successful FIMO output is missing or changed: {tsv_path}")
            print(f"SKIP {slug}: verified successful output", flush=True)
            continue
        attempt_history = list((prior or {}).get("attempt_history", []))
        if prior:
            attempt_history.append({key: value for key, value in prior.items() if key != "attempt_history"})
        if output_dir.exists() and any(output_dir.iterdir()):
            raise FileExistsError(f"refusing to overwrite nonempty directory: {output_dir}")

        log_path = run_root / "logs" / f"{slug}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        command = [
            str(fimo), "--oc", str(output_dir), "--bfile", str(background),
            "--qv-thresh", "--thresh", str(q_threshold),
            "--max-stored-scores", "1000000", str(reference_meme), str(fasta),
        ]
        started_at = datetime.now().astimezone().isoformat(timespec="seconds")
        started = time.monotonic()
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        elapsed = round(time.monotonic() - started, 2)
        log_path.write_text(
            "COMMAND\n" + json.dumps(command, ensure_ascii=False) + "\n\nSTDOUT\n" + process.stdout + "\n\nSTDERR\n" + process.stderr,
            encoding="utf-8",
        )
        required = (tsv_path, output_dir / "fimo.xml", output_dir / "cisml.xml")
        missing = [str(path) for path in required if not path.is_file()] if process.returncode == 0 else []
        status = "success" if process.returncode == 0 and not missing else "failed"
        record: dict[str, Any] = {
            "status": status,
            "combined_fasta_sha256": source["combined_fasta_sha256"],
            "background_sha256": source["background_file_sha256"],
            "command": command,
            "output_directory": str(output_dir),
            "log": str(log_path),
            "started_at": started_at,
            "elapsed_seconds": elapsed,
            "return_code": process.returncode,
            "missing_outputs": missing,
        }
        if attempt_history:
            record["attempt_history"] = attempt_history
        if status == "success":
            record["fimo_tsv_sha256"] = sha256_file(tsv_path)
            record["output_sha256"] = {path.name: sha256_file(path) for path in required}
        manifest["results"][slug] = record
        _write_json(output_manifest, manifest)
        print(f"{status.upper()} {slug}: return_code={process.returncode}, {elapsed}s", flush=True)
        if status != "success":
            raise RuntimeError(f"FIMO reference scan failed for {slug}; see {log_path}")

    manifest["status"] = "complete" if all(
        manifest["results"].get(slug, {}).get("status") == "success" for slug in SPECIES
    ) else "incomplete"
    manifest["completed_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
    _write_json(output_manifest, manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-meme", type=Path, default=Path("baselines/known_promoter_elements.meme"))
    parser.add_argument("--fimo-run-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/fimo/fimo_run_manifest.json"))
    parser.add_argument("--fimo", type=Path, default=Path("tmp/meme-suite-5.5.9/bin/fimo"))
    parser.add_argument("--run-root", type=Path, default=Path("tmp/fimo_runs/tjupan_promloop_20261003/reference_elements"))
    parser.add_argument("--output-manifest", type=Path, default=Path("results/motif/tjupan_m3_20261003/reference_elements/reference_fimo_run_manifest.json"))
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args()
    manifest = run(
        reference_meme=args.reference_meme,
        fimo_run_manifest=args.fimo_run_manifest,
        fimo=args.fimo,
        run_root=args.run_root,
        output_manifest=args.output_manifest,
        q_threshold=args.q_threshold,
    )
    print(f"status={manifest['status']}")


if __name__ == "__main__":
    main()
