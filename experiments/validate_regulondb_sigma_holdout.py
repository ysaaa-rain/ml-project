"""Train RegulonDB sigma-group motifs on a fixed split and test on held-out pairs."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
import math
from pathlib import Path
import random
import shutil
import subprocess
from typing import Any

from experiments.analyze_regulondb_sigma_fimo import GROUPS, run as run_sigma_fimo
from experiments.render_meme_logos import parse_meme_xml
from preprocessing.provenance import sha256_file


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = ROOT / "tmp/motif_inputs_20260930"
DEFAULT_INPUT_MANIFEST = DEFAULT_INPUT_DIR / "motif_input_manifest.json"
DEFAULT_MEME = ROOT / "tmp/meme-suite-5.5.9/bin/meme"
DEFAULT_FIMO = ROOT / "tmp/meme-suite-5.5.9/bin/fimo"
DEFAULT_WORK_DIR = ROOT / "tmp/regulondb_sigma_holdout_20261004"
DEFAULT_RESULT_DIR = ROOT / "results/motif/regulondb_sigma_holdout_20261004"
SEED = 20261004
HOLDOUT_FRACTION = 0.20


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence: list[str] = []

    def finish() -> None:
        if identifier is None:
            return
        value = "".join(sequence).upper()
        if len(value) != 81 or set(value) - set("ACGT"):
            raise ValueError(f"expected an 81-nt A/C/G/T sequence in {path}: {identifier}")
        records.append((identifier, value))

    for line_number, raw in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            identifier, sequence = line[1:].split()[0], []
        elif identifier is None:
            raise ValueError(f"sequence before FASTA header at {path}:{line_number}")
        else:
            sequence.append(line)
    finish()
    if len({name for name, _ in records}) != len(records):
        raise ValueError(f"duplicate FASTA identifiers in {path}")
    return records


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")


def _write_fimo_matrix(source: Path, target: Path, motif_ids: list[str]) -> None:
    """Normalize MEME text motif headers to the ID layout expected by FIMO analysis."""
    lines = source.read_text(encoding="utf-8").splitlines()
    motif_index = 0
    normalized: list[str] = []
    for line in lines:
        line = line.replace(f"{ROOT.resolve()}/", "").rstrip()
        if line.startswith("MOTIF "):
            if motif_index >= len(motif_ids):
                raise ValueError(f"more text motifs than XML motifs in {source}")
            motif_index += 1
            line = f"MOTIF {motif_index} {motif_ids[motif_index - 1]}"
        normalized.append(line)
    if motif_index != len(motif_ids):
        raise ValueError(f"text/XML motif counts differ in {source}: {motif_index} != {len(motif_ids)}")
    target.write_text("\n".join(normalized) + "\n", encoding="utf-8")


def _split_paired_records(
    primary: list[tuple[str, str]],
    controls: list[tuple[str, str]],
    *,
    fraction: float = HOLDOUT_FRACTION,
    seed: int = SEED,
) -> dict[str, dict[str, list[tuple[str, str]]]]:
    if not 0 < fraction < 1:
        raise ValueError("holdout fraction must be between 0 and 1")
    primary_by_id = dict(primary)
    controls_by_id = {
        identifier.removesuffix("__dinucleotide1"): sequence
        for identifier, sequence in controls
    }
    if len(primary_by_id) != len(primary) or len(controls_by_id) != len(controls):
        raise ValueError("duplicate primary or control identifiers")
    if set(primary_by_id) != set(controls_by_id):
        raise ValueError("primary and N1 control IDs are not paired one-to-one")
    if len(primary) < 2:
        raise ValueError("each sigma group needs at least two paired records")

    ordered_ids = [identifier for identifier, _ in primary]
    holdout_count = min(len(ordered_ids) - 1, max(1, math.ceil(len(ordered_ids) * fraction)))
    holdout_ids = set(random.Random(seed).sample(ordered_ids, holdout_count))
    partitions: dict[str, dict[str, list[tuple[str, str]]]] = {
        "train": {"primary": [], "control": []},
        "holdout": {"primary": [], "control": []},
    }
    for identifier, sequence in primary:
        partition = "holdout" if identifier in holdout_ids else "train"
        partitions[partition]["primary"].append((identifier, sequence))
        control_id = f"{identifier}__dinucleotide1"
        partitions[partition]["control"].append((control_id, controls_by_id[identifier]))
    return partitions


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def _require_project_output(path: Path) -> None:
    try:
        path.resolve().relative_to(ROOT.resolve())
    except ValueError as error:
        raise ValueError(f"experiment outputs must remain inside the project directory: {path}") from error


def run(
    *,
    input_dir: Path = DEFAULT_INPUT_DIR,
    input_manifest_path: Path = DEFAULT_INPUT_MANIFEST,
    meme: Path = DEFAULT_MEME,
    fimo: Path = DEFAULT_FIMO,
    work_dir: Path = DEFAULT_WORK_DIR,
    fimo_work_dir: Path | None = None,
    result_dir: Path = DEFAULT_RESULT_DIR,
    reuse_meme_runs: bool = False,
) -> dict[str, Any]:
    if not meme.is_file() or not fimo.is_file():
        raise FileNotFoundError("MEME and FIMO executables must both exist")
    _require_project_output(work_dir)
    _require_project_output(result_dir)
    if fimo_work_dir is not None:
        _require_project_output(fimo_work_dir)
    if work_dir.exists() and any(work_dir.iterdir()) and not reuse_meme_runs:
        raise FileExistsError(f"work directory is not empty: {work_dir}")
    if result_dir.exists() and any(result_dir.iterdir()):
        raise FileExistsError(f"result directory is not empty: {result_dir}")

    input_manifest = json.loads(input_manifest_path.read_text(encoding="utf-8"))
    if input_manifest.get("source_dataset") != "regulondb" or input_manifest.get("split") != "discovery":
        raise ValueError("expected RegulonDB discovery inputs")
    work_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    training_dir = work_dir / "training_inputs"
    holdout_dir = work_dir / "holdout_inputs"
    meme_run_dir = work_dir / "meme_runs"
    matrix_dir = work_dir / "matrices"
    for directory in (training_dir, holdout_dir, meme_run_dir, matrix_dir):
        directory.mkdir(parents=True, exist_ok=True)

    split_groups: dict[str, Any] = {}
    test_outputs: dict[str, Any] = {}
    train_hashes: dict[str, str] = {}
    holdout_hashes: dict[str, str] = {}
    training_runs: dict[str, Any] = {}
    matrix_hashes: dict[str, str] = {}

    for index, group in enumerate(GROUPS):
        group_seed = SEED + index
        source = input_manifest["outputs"][group]
        primary_path = input_dir / source["primary_fasta"]
        control_path = input_dir / source["n1_fasta"]
        if sha256_file(primary_path) != source["primary_sha256"] or sha256_file(control_path) != source["n1_sha256"]:
            raise ValueError(f"{group} source FASTA hash differs from its manifest")
        primary, controls = _read_fasta(primary_path), _read_fasta(control_path)
        if len(primary) != source["sequence_count"] or len(controls) != source["sequence_count"]:
            raise ValueError(f"{group} counts differ from the input manifest")
        partitions = _split_paired_records(primary, controls, seed=group_seed)

        input_outputs: dict[str, Any] = {}
        for partition_name, partition in (("train", partitions["train"]), ("holdout", partitions["holdout"])):
            directory = training_dir if partition_name == "train" else holdout_dir
            primary_out = directory / f"{group}.fasta"
            control_out = directory / f"{group}_n1.fasta"
            _write_fasta(partition["primary"], primary_out)
            _write_fasta(partition["control"], control_out)
            input_outputs[f"{partition_name}_primary"] = {
                "records": len(partition["primary"]),
                "sha256": sha256_file(primary_out),
            }
            input_outputs[f"{partition_name}_N1"] = {
                "records": len(partition["control"]),
                "sha256": sha256_file(control_out),
            }
        train_hashes[group] = input_outputs["train_primary"]["sha256"]
        train_hashes[f"{group}_N1"] = input_outputs["train_N1"]["sha256"]
        holdout_hashes[group] = input_outputs["holdout_primary"]["sha256"]
        holdout_hashes[f"{group}_N1"] = input_outputs["holdout_N1"]["sha256"]
        split_groups[group] = {
            "total_pairs": len(primary),
            "train_pairs": input_outputs["train_primary"]["records"],
            "holdout_pairs": input_outputs["holdout_primary"]["records"],
            "split_seed": group_seed,
            "train_primary_sha256": input_outputs["train_primary"]["sha256"],
            "train_N1_sha256": input_outputs["train_N1"]["sha256"],
            "holdout_primary_sha256": input_outputs["holdout_primary"]["sha256"],
            "holdout_N1_sha256": input_outputs["holdout_N1"]["sha256"],
        }
        test_outputs[group] = {
            "sequence_count": input_outputs["holdout_primary"]["records"],
            "primary_fasta": f"{group}.fasta",
            "primary_sha256": input_outputs["holdout_primary"]["sha256"],
            "n1_fasta": f"{group}_n1.fasta",
            "n1_sha256": input_outputs["holdout_N1"]["sha256"],
        }

        output_dir = meme_run_dir / group
        command = [
            str(meme), str(training_dir / f"{group}.fasta"),
            "-dna", "-mod", "zoops", "-objfun", "de",
            "-neg", str(training_dir / f"{group}_n1.fasta"),
            "-nmotifs", "10", "-minw", "6", "-maxw", "20",
            "-searchsize", "100000", "-revcomp", "-seed", str(SEED),
            "-brief", "50", "-nostatus", "-oc", str(output_dir),
        ]
        stdout_path, stderr_path = work_dir / f"meme_{group}.stdout.log", work_dir / f"meme_{group}.stderr.log"
        xml_path, text_path = output_dir / "meme.xml", output_dir / "meme.txt"
        if reuse_meme_runs:
            if not xml_path.is_file() or not text_path.is_file():
                raise FileNotFoundError(f"cannot reuse incomplete MEME output for {group}: {output_dir}")
            return_code = 0
        else:
            process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False, timeout=3600)
            stdout_path.write_text(process.stdout, encoding="utf-8")
            stderr_path.write_text(process.stderr, encoding="utf-8")
            return_code = process.returncode
        if return_code != 0 or not xml_path.is_file() or not text_path.is_file():
            raise RuntimeError(f"MEME failed for {group}; inspect {stderr_path}")
        motifs = parse_meme_xml(xml_path)
        if len(motifs) != 10:
            raise RuntimeError(f"MEME returned {len(motifs)} motifs for {group}, expected 10")
        matrix_path = matrix_dir / f"{group}.meme"
        _write_fimo_matrix(text_path, matrix_path, [motif["id"] for motif in motifs])
        matrix_hash = sha256_file(matrix_path)
        matrix_hashes[f"matrices/{group}.meme"] = matrix_hash
        training_runs[group] = {
            "return_code": return_code,
            "command": [_relative(Path(part)) if part.startswith("/") else part for part in command],
            "motif_count": len(motifs),
            "meme_txt_sha256": sha256_file(text_path),
            "meme_xml_sha256": sha256_file(xml_path),
            "matrix_sha256": matrix_hash,
            "stdout_sha256": sha256_file(stdout_path) if stdout_path.is_file() else None,
            "stderr_sha256": sha256_file(stderr_path) if stderr_path.is_file() else None,
        }

    split_manifest = {
        "schema_version": "pr01-02-regulondb-sigma-holdout-inputs-1.0",
        "source_dataset": "regulondb",
        "split": "discovery",
        "split_design": "Within each sigma group, paired promoter/N1 records were assigned together to an 80/20 train/holdout split; holdout size is ceil(0.20*N).",
        "seed": SEED,
        "source_manifest_sha256": sha256_file(input_manifest_path),
        "outputs": test_outputs,
        "groups": split_groups,
    }
    split_manifest_path = holdout_dir / "holdout_input_manifest.json"
    split_manifest_path.write_text(json.dumps(split_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    comparison_manifest_path = work_dir / "holdout_matrix_manifest.json"
    comparison_manifest_path.write_text(json.dumps({"output_sha256": matrix_hashes}, indent=2) + "\n", encoding="utf-8")

    selected_fimo_work_dir = fimo_work_dir or work_dir / ("fimo_analysis_reuse" if reuse_meme_runs else "fimo_analysis")
    run_sigma_fimo(
        input_dir=holdout_dir,
        input_manifest_path=split_manifest_path,
        matrix_dir=matrix_dir,
        meme_run_dir=meme_run_dir,
        comparison_manifest_path=comparison_manifest_path,
        fimo=fimo,
        work_dir=selected_fimo_work_dir,
        result_dir=result_dir,
    )
    fimo_manifest_path = result_dir / "sigma_fimo_manifest.json"
    fimo_manifest = json.loads(fimo_manifest_path.read_text(encoding="utf-8"))
    fimo_manifest["schema_version"] = "pr01-02-regulondb-sigma-holdout-1.0"
    fimo_manifest["purpose"] = "在未参与 MEME 发现的 RegulonDB 六 sigma 组留出序列上，复核 motif 命中率与配对背景共现。"
    fimo_manifest["holdout_design"] = {
        "seed": SEED,
        "seed_scheme": "base seed plus zero-based index in GROUPS",
        "train_fraction_approx": 1 - HOLDOUT_FRACTION,
        "holdout_fraction": HOLDOUT_FRACTION,
        "holdout_size_rule": "ceil(0.20*N) per sigma group",
        "split_unit": "paired promoter and its own N1 control",
        "groups": split_groups,
        "training_meme_runs": training_runs,
        "training_fasta_sha256": train_hashes,
        "holdout_fasta_sha256": holdout_hashes,
        "source_manifest_sha256": sha256_file(input_manifest_path),
        "source_input_manifest_path": _relative(input_manifest_path),
        "holdout_input_manifest_sha256": sha256_file(split_manifest_path),
        "reproduction_command": ".venv/bin/python -m experiments.validate_regulondb_sigma_holdout",
        "meme_version": subprocess.run([str(meme), "-version"], capture_output=True, text=True, check=True).stdout.strip(),
        "meme_parameters": ["-dna", "-mod zoops", "-objfun de", "-nmotifs 10", "-minw 6", "-maxw 20", "-searchsize 100000", "-revcomp", f"-seed {SEED}"],
    }
    fimo_manifest["parameters"]["selection_limit"] = (
        "Each motif was discovered only on its source group's train partition and tested on disjoint holdout records. "
        "This is an internal split of one RegulonDB dataset, not an independent external-dataset replication."
    )
    fimo_manifest["parameters"]["position_limit"] = (
        "RegulonDB inputs are TSS-oriented; site coordinates are 1-81 in input direction. "
        "This auxiliary result does not replace six-species TJU Pan analysis."
    )
    fimo_manifest["training_motif_matrix_sha256"] = matrix_hashes
    public_matrix_dir = result_dir / "matrices"
    public_matrix_dir.mkdir(parents=True, exist_ok=True)
    for group in GROUPS:
        source_path = matrix_dir / f"{group}.meme"
        target_path = public_matrix_dir / source_path.name
        shutil.copyfile(source_path, target_path)
        fimo_manifest["output_sha256"][f"matrices/{target_path.name}"] = {
            "bytes": target_path.stat().st_size,
            "sha256": sha256_file(target_path),
        }
    fimo_manifest["created_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
    fimo_manifest_path.write_text(json.dumps(fimo_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return fimo_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--input-manifest", type=Path, default=DEFAULT_INPUT_MANIFEST)
    parser.add_argument("--meme", type=Path, default=DEFAULT_MEME)
    parser.add_argument("--fimo", type=Path, default=DEFAULT_FIMO)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--fimo-work-dir", type=Path)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--reuse-meme-runs", action="store_true", help="复用工作目录中已通过验收的六组 MEME 原生输出")
    args = parser.parse_args()
    result = run(
        input_dir=args.input_dir,
        input_manifest_path=args.input_manifest,
        meme=args.meme,
        fimo=args.fimo,
        work_dir=args.work_dir,
        fimo_work_dir=args.fimo_work_dir,
        result_dir=args.result_dir,
        reuse_meme_runs=args.reuse_meme_runs,
    )
    print(f"Completed six-group held-out validation: {result['rows']}")


if __name__ == "__main__":
    main()
