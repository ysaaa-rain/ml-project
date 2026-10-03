"""Run reproducible, species-wise MEME searches on TJU Pan positive promoters."""

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


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _tool_version(meme: Path) -> str:
    result = subprocess.run(
        [str(meme), "-version"], capture_output=True, text=True, check=True, timeout=15
    )
    return (result.stdout or result.stderr).strip().splitlines()[0]


def _attempt_history(prior: dict[str, Any] | None) -> list[dict[str, Any]]:
    if prior is None:
        return []
    history = list(prior.get("attempt_history", []))
    history.append(dict(prior))
    return history


def run(
    input_manifest: Path,
    meme: Path,
    input_dir: Path,
    run_root: Path,
    result_dir: Path,
    species_slugs: list[str],
) -> dict[str, Any]:
    prepared = json.loads(input_manifest.read_text(encoding="utf-8"))
    version = _tool_version(meme)
    expected = {
        "-dna": True,
        "-mod": "zoops",
        "-objfun": "de",
        "-nmotifs": 10,
        "-minw": 6,
        "-maxw": 20,
        "-searchsize": 0,
        "-revcomp": True,
        "-seed": 20261003,
        "-brief": 50,
    }
    manifest_path = result_dir / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        legacy_expected = {key: value for key, value in expected.items() if key != "-brief"}
        if manifest.get("meme_version") != version or manifest.get("parameters") not in (expected, legacy_expected):
            raise ValueError(f"existing run manifest uses a different MEME version or parameter set: {manifest_path}")
    else:
        manifest = {
            "schema_version": "tjupan-meme-runs-1.0",
            "project": "PR01-02",
            "dataset": prepared["source"],
            "input_manifest": str(input_manifest),
            "input_manifest_sha256": sha256_file(input_manifest),
            "motif_tool": str(meme),
            "meme_version": version,
            "parameters": expected,
            "search_note": "searchsize=0 disables sequence sampling and uses the full primary dataset; objective-specific holdout behavior still applies; runtime may be long",
            "results": {},
        }
    result_dir.mkdir(parents=True, exist_ok=True)
    manifest["parameters"] = expected
    manifest["search_note"] = "searchsize=0 disables sequence sampling and uses the full primary dataset; objective-specific holdout behavior still applies; runtime may be long"
    manifest["input_manifest"] = str(input_manifest)
    manifest["input_manifest_sha256"] = sha256_file(input_manifest)
    _write_json(manifest_path, manifest)

    for slug in species_slugs:
        if slug not in SPECIES:
            raise ValueError(f"unknown species slug: {slug}")
        species_info = prepared["species"][slug]
        positive = input_dir / f"{slug}.positive.fasta"
        control = input_dir / f"{slug}.control.fasta"
        for path, expected_hash in (
            (positive, species_info["positive_samples_csv"]["fasta_sha256"]),
            (control, species_info["negative_control"]["fasta_sha256"]),
        ):
            if not path.is_file() or sha256_file(path) != expected_hash:
                raise ValueError(f"input FASTA is missing or differs from the prepared manifest: {path}")

        prior = manifest["results"].get(slug)
        if prior and prior.get("status") == "success":
            print(f"SKIP {slug}: successful output already recorded", flush=True)
            continue
        attempt_history = _attempt_history(prior)

        output_dir = run_root / slug
        if output_dir.exists() and any(output_dir.iterdir()):
            raise FileExistsError(
                f"refusing to overwrite nonempty output {output_dir}; choose a fresh run root"
            )
        log_path = run_root / "logs" / f"{slug}.stdout.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        command = [
            str(meme),
            str(positive),
            "-dna",
            "-mod",
            "zoops",
            "-objfun",
            "de",
            "-neg",
            str(control),
            "-nmotifs",
            "10",
            "-minw",
            "6",
            "-maxw",
            "20",
            "-searchsize",
            "0",
            "-revcomp",
            "-seed",
            "20261003",
            "-brief",
            "50",
            "-nostatus",
            "-oc",
            str(output_dir),
        ]
        record: dict[str, Any] = {
            "status": "running",
            "species_directory": SPECIES[slug],
            "positive_count": species_info["positive_samples_csv"]["records"],
            "control_count": species_info["negative_control"]["records"],
            "positive_fasta": str(positive),
            "positive_sha256": sha256_file(positive),
            "control_fasta": str(control),
            "control_sha256": sha256_file(control),
            "command": command,
            "output_directory": str(output_dir),
            "log": str(log_path),
            "started_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        }
        if attempt_history:
            record["attempt_history"] = attempt_history
        manifest["results"][slug] = record
        _write_json(manifest_path, manifest)
        print(f"START {slug}: primary={record['positive_count']} control={record['control_count']}", flush=True)
        started = time.monotonic()
        try:
            with log_path.open("w", encoding="utf-8") as log:
                log.write("COMMAND\t" + json.dumps(command, ensure_ascii=False) + "\n")
                log.flush()
                completed = subprocess.run(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
        except OSError as exc:
            record.update(status="failed", error=str(exc))
            completed_returncode = None
        else:
            completed_returncode = completed.returncode
        record["return_code"] = completed_returncode
        record["elapsed_seconds"] = round(time.monotonic() - started, 2)
        record["finished_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
        if completed_returncode == 0:
            required = (output_dir / "meme.txt", output_dir / "meme.xml")
            missing = [str(path) for path in required if not path.is_file()]
            if missing:
                record.update(status="failed", error=f"missing expected MEME output: {missing}")
            else:
                output_files = sorted(
                    path for path in output_dir.iterdir() if path.is_file() and path != log_path
                )
                record.update(
                    status="success",
                    output_files={
                        path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
                        for path in output_files
                    },
                )
        elif completed_returncode is not None:
            record.update(status="failed", error=f"MEME exited with code {completed_returncode}")
        _write_json(manifest_path, manifest)
        print(
            f"{record['status'].upper()} {slug}: elapsed={record['elapsed_seconds']}s "
            f"return_code={completed_returncode}",
            flush=True,
        )
        if record["status"] != "success":
            raise RuntimeError(f"MEME failed for {slug}; inspect {log_path}")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meme", default="tmp/meme-suite-5.5.9/bin/meme")
    parser.add_argument("--input-dir", default="tmp/tjupan_motif_inputs_20261003")
    parser.add_argument("--input-manifest", default="results/motif/tjupan_m3_20261003/input_manifest.json")
    parser.add_argument("--run-root", default="tmp/meme_runs/tjupan_promloop_20261003")
    parser.add_argument("--result-dir", default="results/motif/tjupan_m3_20261003")
    parser.add_argument("--species", nargs="+", default=list(SPECIES))
    args = parser.parse_args()
    run(
        Path(args.input_manifest),
        Path(args.meme),
        Path(args.input_dir),
        Path(args.run_root),
        Path(args.result_dir),
        args.species,
    )


if __name__ == "__main__":
    main()
