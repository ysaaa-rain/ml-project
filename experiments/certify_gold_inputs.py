"""Independent exhaustive cross-split gate for a frozen G0 data run."""
from __future__ import annotations
import argparse
from datetime import datetime
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def certify(config_path, output_name):
    config = json.loads(config_path.read_text(encoding="utf-8"))
    data_dir = ROOT/config["data_output"]
    result_dir = ROOT/config["summary_output"]
    if Path(output_name).name != output_name:
        raise ValueError("output name must be a filename")
    destination = result_dir/output_name
    if destination.exists():
        raise FileExistsError(destination)
    manifest = json.loads((data_dir/"run_manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["output_sha256"].items():
        if sha(data_dir/name) != expected:
            raise ValueError(f"frozen output changed: {name}")
    metadata = pd.read_csv(data_dir/"standard_metadata.tsv", sep="\t")
    audits = []
    for upstream in config["windows_upstream"]:
        window = metadata[metadata.upstream == upstream]
        train = window[window.split == "train"].sequence.tolist()
        holdout = window[window.split == "holdout"].sequence.tolist()
        length = upstream+config["downstream"]+1
        cutoff = length-math.ceil(config["similarity"]*length)
        # This gate has no seed index: every train x holdout pair is examined.
        hits_edit = hits_shift = 0
        for orientation in ("direct", "reverse_complement"):
            controls = holdout if orientation == "direct" else [s.translate(str.maketrans("ACGT", "TGCA"))[::-1] for s in holdout]
            b = np.asarray([list(s.encode("ascii")) for s in controls], dtype=np.uint8)
            for start in range(0, len(train), 64):
                queries = train[start:start+64]
                distances = process.cdist(queries, controls, scorer=Levenshtein.distance,
                                          score_cutoff=cutoff, dtype=np.uint8, workers=1)
                hits_edit += int(np.count_nonzero(distances <= cutoff))
                a = np.asarray([list(s.encode("ascii")) for s in queries], dtype=np.uint8)
                for shift in range(1, cutoff+1):
                    for left, right in ((a[:, shift:], b[:, :-shift]), (a[:, :-shift], b[:, shift:])):
                        mismatch = np.count_nonzero(left[:, None, :] != right[None, :, :], axis=2)
                        hits_shift += int(np.count_nonzero(mismatch+shift <= cutoff))
        audits.append({"upstream": upstream, "length": length, "train_loci": len(train),
                       "holdout_loci": len(holdout), "pairs_per_orientation": len(train)*len(holdout),
                       "edit_distance_hits": hits_edit, "shifted_overlap_hits": hits_shift})
        print(f"independent gate {upstream}: edit={hits_edit}, shifted={hits_shift}", flush=True)
    if any(a["edit_distance_hits"] or a["shifted_overlap_hits"] for a in audits):
        raise ValueError("cross-split homology gate failed")
    primary = metadata[metadata.upstream == 80]
    labels = []
    for sigma in sorted({s for value in primary.gold_sigma_labels for s in json.loads(value)}):
        members = primary[primary.gold_sigma_labels.map(lambda v: sigma in json.loads(v))]
        representatives = members[members.gold_representative]
        labels.append({"sigma": sigma, "gold_loci_before_deduplication": len(members),
                       "clusters_containing_gold_label": int(members.cluster_id.nunique()),
                       "representatives_with_actual_gold_label": len(representatives),
                       "clusters_without_representative_label": int(members.cluster_id.nunique()-representatives.cluster_id.nunique())})
    pd.DataFrame(labels).to_csv(result_dir/"label_retention.tsv", sep="\t", index=False, lineterminator="\n")
    dbtbs_report = json.loads((ROOT/"data/interim/dbtbs_download_report_20260916.json").read_text(encoding="utf-8"))
    dbtbs_dir = ROOT/"data/raw/dbtbs_v4.1_20260916/pages"
    matches = missing = changed = 0
    for name, expected in dbtbs_report["page_hashes"].items():
        path = dbtbs_dir/name
        if not path.exists():
            missing += 1
        elif sha(path) != expected:
            changed += 1
        else:
            matches += 1
    dbtbs = pd.read_csv(ROOT/"data/interim/dbtbs_metadata_20260916.tsv", sep="\t")
    corrected = dbtbs.rename(columns={"source_tss_coordinate": "source_region_start", "source_tss_end_coordinate": "source_region_end"}).copy()
    corrected["evidence_level"] = "source_reference_present_requires_review"
    corrected["coordinate_status"] = "TSS_and_strand_unverified"
    corrected["original_evidence_level"] = dbtbs.evidence_level
    corrected_path = data_dir/"dbtbs_audited_candidates.tsv"
    if corrected_path.exists():
        raise FileExistsError(corrected_path)
    corrected.to_csv(corrected_path, sep="\t", index=False, lineterminator="\n")
    output = {"status": "passed_for_stated_cross_split_criteria", "created_at": datetime.now().astimezone().isoformat(),
              "command": "python -B -m experiments.certify_gold_inputs "+" ".join(sys.argv[1:]), "cross_split_audit": audits,
              "input_manifest_sha256": sha(data_dir/"run_manifest.json"), "config_sha256": sha(config_path),
              "builder_sha256": sha(ROOT/"experiments/rebuild_gold_inputs.py"), "gate_script_sha256": sha(Path(__file__)),
              "build_command": "python -B -m experiments.rebuild_gold_inputs", "software": {"numpy": np.__version__, "rapidfuzz": __import__("rapidfuzz").__version__},
              "multi_sigma_loci": int(primary.sigma_label_status.eq("multi_sigma").sum()),
              "supplement_output_sha256": {corrected_path.name: sha(corrected_path), "label_retention.tsv": sha(result_dir/"label_retention.tsv")},
              "single_sigma_ml_representatives": int(primary.ml_eligible.sum()),
              "DBTBS": {"downloaded_pages_matching_hash": matches, "missing_pages": missing, "changed_pages": changed,
                        "nonempty_source_evidence": int(dbtbs.source_evidence.fillna("").ne("").sum()),
                        "evidence_gate": "references retained; experimental label assigned by old importer is not independently verified",
                        "coordinate_gate": "not_passed; source_tss_coordinate is a promoter-region boundary, not a verified TSS"},
              "limits": ["No motifs discovered or validated.", "Each cluster has one representative; sigma labels of other members are not transferred to it.",
                         "The criteria do not certify arbitrary shorter local alignments or external-source separation."]}
    destination.write_text(json.dumps(output, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/g0_rebuild_20261004.json")
    parser.add_argument("--output-name", default="independent_gate_complete.json")
    args = parser.parse_args()
    certify(args.config, args.output_name)
