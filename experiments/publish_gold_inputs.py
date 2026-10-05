"""Certify and describe the already-frozen collaboration inputs, without rebuilding."""
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run():
    config_path = ROOT/"configs/g0_rebuild_20261004.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    directory = ROOT/config["data_output"]
    summaries = ROOT/config["summary_output"]
    destination = directory/"share_manifest.json"
    if destination.exists():
        raise FileExistsError(destination)
    original_path = directory/"run_manifest.json"
    original = json.loads(original_path.read_text(encoding="utf-8"))
    gate_path = summaries/"independent_gate_complete.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    supplement_path = summaries/"supplement_manifest.json"
    supplement = json.loads(supplement_path.read_text(encoding="utf-8"))
    if gate["status"] != "passed_for_stated_cross_split_criteria" or sha(original_path) != gate["input_manifest_sha256"]:
        raise ValueError("original data gate does not match")
    if sha(config_path) != gate["config_sha256"]:
        raise ValueError("frozen config changed")
    expected = {directory/name: digest for name, digest in original["output_sha256"].items()}
    for name, digest in gate["supplement_output_sha256"].items():
        expected[directory/name if name == "dbtbs_audited_candidates.tsv" else summaries/name] = digest
    for name, digest in supplement["output_sha256"].items():
        expected[ROOT/name] = digest
    for path, digest in expected.items():
        if sha(path) != digest:
            raise ValueError(f"frozen artifact changed: {path.relative_to(ROOT)}")
    files = sorted(p for p in directory.rglob("*") if p.is_file() and p != original_path)
    publication = {
        "publication_date": "2026-10-05", "dataset_run_id": config["run_id"],
        "status": "frozen_input_publication_only; no_new_motif_or_ML_experiment",
        "command": "python -B -m experiments.publish_gold_inputs",
        "script_sha256": sha(Path(__file__)), "config_sha256": sha(config_path),
        "git_base_at_publication": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "original_local_manifest_sha256": sha(original_path),
        "original_run_metadata": {key: value for key, value in original.items() if key != "software"},
        "audit_provenance_sha256": {gate_path.relative_to(ROOT).as_posix(): sha(gate_path), supplement_path.relative_to(ROOT).as_posix(): sha(supplement_path)},
        "files": {path.relative_to(directory).as_posix(): {"bytes": path.stat().st_size, "sha256": sha(path)} for path in files},
        "limits": ["Original machine environment manifest stays local; its byte hash is recorded here.",
                   "DBTBS candidates are not gold-certified; no verified TSS/strand.",
                   "Only RegulonDB internal cross-split criteria passed; external-source separation and motif experiments remain pending."]}
    destination.write_bytes((json.dumps(publication, ensure_ascii=False, indent=2)+"\n").encode("utf-8"))
    print(json.dumps({"verified_frozen_artifacts": len(expected), "published_payload_files": len(files),
                      "payload_bytes": sum(p.stat().st_size for p in files), "share_manifest_sha256": sha(destination)}))

if __name__ == "__main__":
    run()
