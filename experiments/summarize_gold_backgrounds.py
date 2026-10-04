"""Measure matched-background GC differences without altering frozen inputs."""
import hashlib
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def read_fasta(path):
    records = []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            records.append([line[1:], ""])
        elif line:
            records[-1][1] += line
    return records

def run():
    config_path = ROOT/"configs/g0_rebuild_20261004.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    data_dir, result_dir = ROOT/config["data_output"], ROOT/config["summary_output"]
    destination = result_dir/"background_gc_summary.tsv"
    if destination.exists():
        raise FileExistsError(destination)
    gc = lambda sequence: (sequence.count("G")+sequence.count("C"))/len(sequence)
    rows = []
    for upstream in config["windows_upstream"]:
        directory = data_dir/f"window_{upstream}_{config['downstream']}"/"train"
        for path in sorted(directory.glob("Sigma*_other_matched_primary.fasta")):
            sigma = path.name.split("_")[0]
            positive = read_fasta(path)
            control = read_fasta(directory/f"{sigma}_other_matched.fasta")
            if len(positive) != len(control):
                raise ValueError("unpaired controls")
            ids = [item[0].split("__draw")[0] for item in control]
            if len(ids) != len(set(ids)):
                raise ValueError("reused controls")
            differences = [abs(gc(a[1])-gc(b[1])) for a,b in zip(positive, control)]
            rows.append({"upstream": upstream, "sigma": sigma, "matched_pairs": len(control),
                         "primary_mean_gc": sum(gc(s) for _,s in positive)/len(positive) if positive else None,
                         "control_mean_gc": sum(gc(s) for _,s in control)/len(control) if control else None,
                         "mean_absolute_paired_gc_difference": sum(differences)/len(differences) if differences else None,
                         "max_absolute_paired_gc_difference": max(differences) if differences else None})
    pd.DataFrame(rows).to_csv(destination, sep="\t", index=False, lineterminator="\n")
    metadata = pd.read_csv(data_dir/"standard_metadata.tsv", sep="\t")
    index = metadata[["sequence_id", "upstream", "cluster_id", "split", "gold_representative"]].copy()
    index["sequence_sha256"] = metadata.sequence.map(lambda s: hashlib.sha256(s.encode("ascii")).hexdigest())
    index["selection_reason"] = metadata.gold_representative.map({True: "gold_cluster_representative_by_annotation_quality_then_ID", False: "retained_in_audit_only_no_label_transfer"})
    index_path = data_dir/"sequence_hash_and_selection.tsv"
    if index_path.exists():
        raise FileExistsError(index_path)
    index.to_csv(index_path, sep="\t", index=False, lineterminator="\n")
    manifest = {"command": "python -B -m experiments.summarize_gold_backgrounds", "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "input_metadata_sha256": hashlib.sha256((data_dir/"standard_metadata.tsv").read_bytes()).hexdigest(),
                "output_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in [destination, index_path]}}
    (result_dir/"supplement_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")

if __name__ == "__main__":
    run()
