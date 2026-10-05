"""Read-only M1 source-to-published-input acceptance checks."""
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import random
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rc(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def run():
    directory = ROOT/"data/processed/g0_rebuild_20261004"
    config = json.loads((ROOT/"configs/g0_rebuild_20261004.json").read_text(encoding="utf-8"))
    share = json.loads((directory/"share_manifest.json").read_text(encoding="utf-8"))
    for name, entry in share["files"].items():
        assert sha(directory/name) == entry["sha256"], f"shared file changed: {name}"
    for name, expected in share["original_run_metadata"]["source_sha256"].items():
        assert sha(ROOT/name) == expected, f"source changed: {name}"
    for name, expected in share["audit_provenance_sha256"].items():
        assert sha(ROOT/name) == expected, f"gate evidence changed: {name}"
    genome_lines = (ROOT/config["genome_fasta"]).read_text(encoding="utf-8").splitlines()
    genome = "".join(s.strip() for s in genome_lines if not s.startswith(">"))
    expected_loci = defaultdict(list)
    excluded = []
    for number, line in enumerate((ROOT/config["regulondb_gff3"]).read_text(encoding="utf-8").splitlines(), 1):
        if line.startswith("#") or not line:
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] not in {"promoter", "transcription_start_site"}:
            continue
        attrs = {unquote(item.split("=", 1)[0].strip()): unquote(item.split("=", 1)[1].strip()) for item in fields[8].split(";") if "=" in item}
        sequence = attrs.get("Sequence", "")
        if not sequence:
            excluded.append(attrs["name"])
            continue
        marker = [i for i, base in enumerate(sequence) if base in "ACGT"]
        assert len(marker) == 1 and fields[3] == fields[4] and fields[6] in {"+", "-"}
        tss = int(fields[3])
        left = tss-1-marker[0] if fields[6] == "+" else tss-1-(len(sequence)-1-marker[0])
        fragment = genome[left:left+len(sequence)]
        oriented = fragment if fields[6] == "+" else rc(fragment)
        if oriented != sequence.upper():
            excluded.append(attrs["name"])
            continue
        expected_loci[(tss, fields[6])].append({"id": attrs["name"], "sigma": attrs.get("SigmaFactor", ""),
                                               "confidence": attrs.get("Confidence"), "evidence": attrs.get("Evidence", ""), "line": number})
    with (directory/"standard_metadata.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == len(expected_loci)*3
    observed = set()
    for row in rows:
        tss, upstream, downstream = int(row["TSS"]), int(row["upstream"]), int(row["downstream"])
        key = tss, row["strand"]
        evidence = expected_loci[key]
        assert json.loads(row["evidence"]) == evidence
        labels = sorted({r["sigma"] for r in evidence if r["sigma"]})
        gold = sorted({r["sigma"] for r in evidence if r["sigma"] and r["confidence"] in config["accepted_confidence"]})
        assert json.loads(row["sigma_raw"]) == labels and json.loads(row["gold_sigma_labels"]) == gold
        start, end = (tss-upstream, tss+downstream) if key[1] == "+" else (tss-downstream, tss+upstream)
        assert start >= 1 and end <= len(genome)
        actual = genome[start-1:end]
        actual = actual if key[1] == "+" else rc(actual)
        assert row["sequence"] == actual and len(actual) == upstream+downstream+1
        assert int(row["genomic_start"]) == start and int(row["genomic_end"]) == end and int(row["tss_index_1based"]) == upstream+1
        assert row["sequence_orientation"] == "transcription"
        assert (key, upstream) not in observed
        observed.add((key, upstream))
    clusters = defaultdict(list)
    for row in rows:
        if int(row["upstream"]) == 80:
            clusters[row["cluster_id"]].append(row)
    keys = sorted(clusters, key=lambda key: min((int(r["TSS"]), r["strand"]) for r in clusters[key]))
    shuffled = keys.copy()
    random.Random(config["seed"]).shuffle(shuffled)
    heldout = set(shuffled[:math.ceil(len(keys)*config["holdout_fraction"])])
    sequence_splits = {}
    for key, members in clusters.items():
        assert {r["split"] for r in members} == ({"holdout"} if key in heldout else {"train"})
        cluster_labels = sorted({s for r in members for s in json.loads(r["sigma_raw"])})
        gold_members = [r for r in members if json.loads(r["gold_sigma_labels"])]
        selected = min(gold_members, key=lambda r: (-sum(e["confidence"] == "Confirmed" for e in json.loads(r["evidence"])), r["sequence_id"])) if gold_members else None
        assert [r["sequence_id"] for r in members if r["gold_representative"] == "True"] == ([selected["sequence_id"]] if selected else [])
        for row in members:
            assert json.loads(row["cluster_sigma_labels"]) == cluster_labels
            labels, gold = json.loads(row["sigma_raw"]), json.loads(row["gold_sigma_labels"])
            eligible = row is selected and len(labels) == 1 and labels == gold and len(cluster_labels) == 1
            assert (row["ml_eligible"] == "True") == eligible
            canonical = min(row["sequence"], rc(row["sequence"]))
            assert canonical not in sequence_splits or sequence_splits[canonical] == row["split"]
            sequence_splits[canonical] = row["split"]
    primary = [r for r in rows if int(r["upstream"]) == 80]
    published = {"audit_date": "2026-10-05", "status": "RegulonDB_internal_M1_checks_passed; overall_M1_partial",
                 "command": "python -B -m experiments.audit_shared_m1", "script_sha256": sha(Path(__file__)),
                 "share_manifest_sha256": sha(directory/"share_manifest.json"), "checked_windows": len(rows),
                 "verified_loci": len(expected_loci), "clusters": len(clusters), "excluded_source_records": len(excluded),
                 "gold_representatives": sum(r["gold_representative"] == "True" for r in primary),
                 "single_sigma_ml": sum(r["ml_eligible"] == "True" for r in primary),
                 "checks": ["shared payload byte hashes", "original source hashes", "all genomic windows and bounds", "all source evidence and sigma labels", "one deterministic representative per gold cluster", "seeded global split reconstruction", "strict ML eligibility", "exact/RC overlap across split"],
                 "near_homology_gate": "reuse exhaustive 2026-10-04 gate after verifying byte-identical inputs; not rerun here",
                 "remaining_M1_requirements": ["DBTBS record-level evidence and reliable TSS/strand", "external-source homology separation", "explicit data reuse/redistribution terms not verified; do not substitute software license"],
                 "experiment_status": "No new motif discovery, scanning or classifier training"}
    destination = ROOT/"results/data_audit/g0_rebuild_20261004/m1_acceptance_20261005.json"
    if destination.exists():
        raise FileExistsError(destination)
    destination.write_bytes((json.dumps(published, ensure_ascii=False, indent=2)+"\n").encode("utf-8"))
    print(json.dumps(published, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    run()
