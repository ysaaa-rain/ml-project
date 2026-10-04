"""Reconstruct genomic promoter windows and freeze conservative cluster splits.

No motif tool is invoked. Source files and all historical outputs are untouched.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess

import pandas as pd
from rapidfuzz.distance import Levenshtein
from experiments.preflight import collect_environment
from preprocessing.source_adapters import parse_gff3_attributes
from preprocessing.shuffle import shuffled_dinucleotide, dinucleotide_counts

ROOT = Path(__file__).resolve().parents[1]
SIGMAS = {"Sigma70", "Sigma24", "Sigma32", "Sigma38", "Sigma28", "Sigma54"}


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write_json(p, value):
    p.write_bytes((json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))


def write_tsv(p, rows, columns=None):
    pd.DataFrame(rows, columns=columns).to_csv(p, sep="\t", index=False, lineterminator="\n")


def fasta(p, rows):
    p.write_bytes("".join(f">{r['sequence_id']}\n{r['sequence']}\n" for r in rows).encode("ascii"))


def near_sequences(a, b, cutoff):
    """Equal-length, full coverage edit distance or shifted >=90% overlap.

    The overlap score includes unaligned end bases as errors. For shift d,
    errors = d + mismatches in overlap; thus coverage is >=1-cutoff/L.
    """
    if len(a) != len(b):
        raise ValueError("cluster comparisons must use the same genomic window")
    if Levenshtein.distance(a, b, score_cutoff=cutoff) <= cutoff:
        return True
    for shift in range(1, cutoff+1):
        for x, y in ((a[shift:], b[:-shift]), (a[:-shift], b[shift:])):
            if shift + sum(left != right for left, right in zip(x, y)) <= cutoff:
                return True
    return False


def cluster_window(sequences, similarity):
    """Complete exact-seed candidates for bounded edits; no frequency pruning.

    Split a query into k+1 disjoint pieces for k=floor((1-similarity)*L).
    At most k edits/end bases can disrupt those pieces, so at least one
    remains intact. Its first floor(L/(k+1)) bases occur in the other string.
    Index every substring of that size, including the reverse complement.
    Candidate pairs are confirmed explicitly; low complexity is never pruned.
    """
    if not sequences:
        return [], 0
    length = len(sequences[0])
    if any(len(s) != length for s in sequences):
        raise ValueError("mixed window lengths")
    cutoff = length - math.ceil(similarity*length)
    parts = cutoff+1
    seed_length = length//parts
    index = defaultdict(set)
    for i, s in enumerate(sequences):
        for variant in (s, rc(s)):
            for start in range(length-seed_length+1):
                index[variant[start:start+seed_length]].add(i)
    edges = []
    candidate_count = 0
    for i, s in enumerate(sequences):
        candidates = set()
        for part in range(parts):
            start = part*length//parts
            candidates.update(index[s[start:start+seed_length]])
        for j in sorted(candidates):
            if j <= i:
                continue
            candidate_count += 1
            if near_sequences(s, sequences[j], cutoff) or near_sequences(s, rc(sequences[j]), cutoff):
                edges.append((i, j))
    return edges, candidate_count


def run(config_path):
    config = json.loads(config_path.read_text(encoding="utf-8"))
    data_dir, result_dir = ROOT/config["data_output"], ROOT/config["summary_output"]
    for p in (data_dir, result_dir):
        p.resolve().relative_to(ROOT)
        if p.exists() and any(p.iterdir()):
            raise FileExistsError(f"use new output directories: {p}")
    gff, genome_path = ROOT/config["regulondb_gff3"], ROOT/config["genome_fasta"]
    if sha(gff) != config["expected_gff3_sha256"] or sha(genome_path) != config["expected_genome_sha256"]:
        raise ValueError("source bytes differ from the frozen snapshot")
    genome_lines = genome_path.read_text(encoding="utf-8").splitlines()
    headers = [s for s in genome_lines if s.startswith(">")]
    if len(headers) != 1:
        raise ValueError("expected a single-contig genome")
    accession = headers[0][1:].split()[0]
    genome = "".join(s.strip() for s in genome_lines if not s.startswith(">")).upper()
    loci = {}
    excluded = []
    raw_count = 0
    max_upstream = max(config["windows_upstream"])
    downstream = config["downstream"]
    for line_no, line in enumerate(gff.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) != 9 or f[2] not in {"transcription_start_site", "promoter"}:
            continue
        raw_count += 1
        attrs = parse_gff3_attributes(f[8])
        record_id = attrs.get("name", f"line_{line_no}")
        s = attrs.get("Sequence", "")
        marker = [i for i, base in enumerate(s) if base in "ACGT"]
        tss, strand = int(f[3]), f[6]
        reason = None
        if not s:
            reason = "missing_source_sequence"
        elif len(marker) != 1:
            reason = "ambiguous_source_tss_marker"
        elif strand not in {"+", "-"} or f[0] != accession or f[3] != f[4]:
            reason = "invalid_genome_coordinate"
        elif set(s.upper())-set("ACGT"):
            reason = "ambiguous_bases"
        else:
            offset = marker[0]
            start = tss-offset if strand == "+" else tss-(len(s)-1-offset)
            stop = start+len(s)-1
            original = genome[start-1:stop]
            original = original if strand == "+" else rc(original)
            if start < 1 or stop > len(genome) or original != s.upper():
                reason = "source_genome_orientation_disagreement"
        left = tss-max_upstream if strand == "+" else tss-downstream
        right = tss+downstream if strand == "+" else tss+max_upstream
        if reason is None and (left < 1 or right > len(genome)):
            reason = "genomic_window_out_of_bounds"
        if reason:
            excluded.append({"source_record_id": record_id, "line": line_no, "reason": reason})
            continue
        key = (tss, strand)
        if key not in loci:
            loci[key] = {"sequence_id": f"regulondb_{tss}_{'plus' if strand == '+' else 'minus'}",
                         "TSS": tss, "strand": strand, "records": [], "sigma_labels": set(), "gold_sigma_labels": set()}
        locus = loci[key]
        sigma = attrs.get("SigmaFactor", "")
        if sigma and sigma not in SIGMAS:
            raise ValueError(f"unmapped sigma label: {sigma}")
        locus["records"].append({"id": record_id, "sigma": sigma, "confidence": attrs.get("Confidence"),
                                 "evidence": attrs.get("Evidence", ""), "line": line_no})
        if sigma:
            locus["sigma_labels"].add(sigma)
            if attrs.get("Confidence") in config["accepted_confidence"]:
                locus["gold_sigma_labels"].add(sigma)
    ordered = [loci[key] for key in sorted(loci)]
    sequences = {}
    for upstream in config["windows_upstream"]:
        current = []
        for locus in ordered:
            tss, strand = locus["TSS"], locus["strand"]
            start = tss-upstream if strand == "+" else tss-downstream
            stop = tss+downstream if strand == "+" else tss+upstream
            s = genome[start-1:stop]
            s = s if strand == "+" else rc(s)
            if len(s) != upstream+downstream+1 or set(s)-set("ACGT"):
                raise ValueError("genomic window failed length/alphabet validation")
            current.append(s)
        sequences[upstream] = current
    parents = list(range(len(ordered)))
    def find(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i
    def union(i, j):
        a, b = find(i), find(j)
        parents[max(a,b)] = min(a,b)
    cluster_stats = []
    for upstream, current in sequences.items():
        edges, candidates = cluster_window(current, config["similarity"])
        for i,j in edges:
            union(i,j)
        cluster_stats.append({"upstream": upstream, "length": len(current[0]), "candidate_pairs": candidates,
                              "confirmed_pairs": len(edges), "cutoff": len(current[0])-math.ceil(config["similarity"]*len(current[0]))})
        print(f"clustered {upstream}: {len(edges)} confirmed pairs", flush=True)
    clusters = defaultdict(list)
    for i, locus in enumerate(ordered):
        clusters[find(i)].append(i)
    # Freeze one split per global cluster. No outcome-driven resampling.
    cluster_keys = sorted(clusters)
    shuffled = list(cluster_keys)
    random.Random(config["seed"]).shuffle(shuffled)
    holdout_keys = set(shuffled[:math.ceil(len(shuffled)*config["holdout_fraction"])])
    data_dir.mkdir(parents=True)
    result_dir.mkdir(parents=True)
    rows, members, representatives = [], [], []
    for cluster_key, indices in sorted(clusters.items()):
        cluster_id = ordered[cluster_key]["sequence_id"]
        split = "holdout" if cluster_key in holdout_keys else "train"
        gold_indices = [i for i in indices if ordered[i]["gold_sigma_labels"]]
        # Representative is selected using annotation quality/ID only, never motif output.
        representative = min(gold_indices, key=lambda i: (-sum(r["confidence"] == "Confirmed" for r in ordered[i]["records"]), ordered[i]["sequence_id"])) if gold_indices else None
        cluster_labels = sorted(set().union(*(ordered[i]["sigma_labels"] for i in indices)))
        for i in indices:
            locus = ordered[i]
            labels, gold_labels = sorted(locus["sigma_labels"]), sorted(locus["gold_sigma_labels"])
            for upstream in config["windows_upstream"]:
                row = {"sequence_id": locus["sequence_id"], "species": "Escherichia coli", "strain": "K-12 MG1655",
                       "source": "RegulonDB", "source_version": config["source_snapshot"], "genome_accession": accession,
                       "TSS": locus["TSS"], "strand": locus["strand"], "sequence_orientation": "transcription",
                       "coordinate_status": "source_sequence_and_genome_verified", "upstream": upstream, "downstream": downstream,
                       "genomic_start": locus["TSS"]-upstream if locus["strand"] == "+" else locus["TSS"]-downstream,
                       "genomic_end": locus["TSS"]+downstream if locus["strand"] == "+" else locus["TSS"]+upstream,
                       "tss_index_1based": upstream+1, "sequence": sequences[upstream][i],
                       "sigma_raw": json.dumps(labels), "gold_sigma_labels": json.dumps(gold_labels),
                       "sigma_family": json.dumps({label: "sigma54" if label == "Sigma54" else "sigma70" for label in labels}),
                       "sigma_label_status": "multi_sigma" if len(labels)>1 else "single_sigma" if labels else "unknown",
                       "evidence": json.dumps(locus["records"]), "cluster_id": cluster_id, "split": split,
                       "cluster_sigma_labels": json.dumps(cluster_labels), "gold_representative": i == representative,
                       "ml_eligible": i == representative and len(labels) == 1 and labels == gold_labels and len(cluster_labels) == 1}
                rows.append(row)
            members.append({"sequence_id": locus["sequence_id"], "cluster_id": cluster_id, "split": split,
                            "gold_representative": i == representative})
        if representative is not None:
            representatives.append(representative)
    write_tsv(data_dir/"standard_metadata.tsv", rows)
    write_tsv(data_dir/"cluster_members.tsv", members)
    write_tsv(data_dir/"excluded_source_records.tsv", excluded)
    pd.DataFrame(rows).query("ml_eligible and upstream == 80").to_csv(data_dir/"single_sigma_ml.tsv",sep="\t",index=False,lineterminator="\n")
    summaries = []
    background_counts = []
    for upstream in config["windows_upstream"]:
        window_rows = [r for r in rows if r["upstream"] == upstream and r["gold_representative"]]
        for split in ("train", "holdout"):
            partition = [r for r in window_rows if r["split"] == split]
            target_dir = data_dir/f"window_{upstream}_{downstream}"/split
            target_dir.mkdir(parents=True)
            fasta(target_dir/"all_gold.fasta", partition)
            for sigma in sorted(SIGMAS):
                # Cluster member uncertainty excludes a target from other-sigma controls.
                positive = [r for r in partition if sigma in json.loads(r["gold_sigma_labels"])]
                other = [r for r in partition if sigma not in json.loads(r["cluster_sigma_labels"])]
                fasta(target_dir/f"{sigma}.fasta", positive)
                summaries.append({"upstream": upstream, "split": split, "sigma": sigma, "unique_cluster_representatives": len(positive),
                                  "method": "STREME" if len(positive)>=50 else "MEME" if len(positive)>=10 else "not_formal",
                                  "ml_single_sigma": sum(r["ml_eligible"] for r in positive)})
                if split != "train":
                    continue
                rng = random.Random(config["seed"]+upstream+sum(map(ord,sigma)))
                controls = []
                for r in positive:
                    shuffled_seq = r["sequence"]
                    for _ in range(100):
                        shuffled_seq = shuffled_dinucleotide(r["sequence"], rng)
                        if shuffled_seq != r["sequence"]:
                            break
                    if shuffled_seq == r["sequence"]:
                        raise ValueError("composition-preserving background cannot be changed")
                    if dinucleotide_counts(shuffled_seq) != dinucleotide_counts(r["sequence"]):
                        raise ValueError("background composition mismatch")
                    controls.append({"sequence_id": r["sequence_id"]+"__N1", "sequence": shuffled_seq})
                fasta(target_dir/f"{sigma}_N1.fasta", controls)
                # One-to-one matching without replacement in frozen 5% GC strata.
                bin_width = config["gc_bin_width"]
                gc_bin = lambda r: int((r["sequence"].count("G")+r["sequence"].count("C"))/len(r["sequence"])/bin_width)
                other_bins = defaultdict(list)
                for r in other:
                    other_bins[gc_bin(r)].append(r)
                matched_positive, matched_other = [], []
                for r in positive:
                    pool = other_bins.get(gc_bin(r), [])
                    if not pool:
                        continue
                    chosen = pool.pop(rng.randrange(len(pool)))
                    matched_positive.append(r)
                    matched_other.append(chosen)
                fasta(target_dir/f"{sigma}_other_matched_primary.fasta", matched_positive)
                # Draw IDs remain distinct and traceable; no control is reused.
                fasta(target_dir/f"{sigma}_other_matched.fasta", [{"sequence_id": r["sequence_id"]+f"__draw{i}", "sequence": r["sequence"]} for i,r in enumerate(matched_other)])
                background_counts.append({"upstream": upstream, "sigma": sigma, "primary": len(positive), "N1": len(controls),
                                          "other_matched": len(matched_other), "other_unique": len({r["sequence_id"] for r in matched_other}),
                                          "unmatched_primary": len(positive)-len(matched_positive), "gc_bin_width": bin_width,
                                          "sampling": "without_replacement; one control per matched primary",
                                          "other_background_method": "STREME" if len(matched_other)>=50 else "MEME" if len(matched_other)>=10 else "not_formal"})
    write_tsv(result_dir/"sigma_split_summary.tsv", summaries)
    write_tsv(result_dir/"background_summary.tsv", background_counts)
    write_tsv(result_dir/"exclusion_summary.tsv", [{"reason": k, "records": n} for k,n in sorted(Counter(r["reason"] for r in excluded).items())])
    write_tsv(result_dir/"cluster_summary.tsv", [{"split": split, "clusters": sum((k in holdout_keys)==(split=="holdout") for k in clusters),
                                                "gold_representatives": sum((find(i) in holdout_keys)==(split=="holdout") for i in representatives)} for split in ["train","holdout"]])
    assert pd.DataFrame(members).groupby("cluster_id").split.nunique().max() == 1
    dbtbs = pd.read_csv(ROOT/"data/interim/dbtbs_metadata_20260916.tsv",sep="\t")
    dbtbs_status = {"records": len(dbtbs), "missing_tss": int(dbtbs.tss_position.isna().sum()), "missing_strand": int(dbtbs.strand.isna().sum()),
                    "coordinate_gate": "not_passed", "permitted": "nonpositional candidate discovery after evidence audit; no strict positional validation"}
    manifest = {"run_id": config["run_id"], "created_at": datetime.now().astimezone().isoformat(), "status": "RegulonDB genomic inputs rebuilt; no motif experiment",
                "source_records": raw_count, "excluded_records": len(excluded), "verified_loci": len(ordered), "clusters": len(clusters),
                "gold_representatives": len(representatives), "split_seed": config["seed"], "holdout_fraction_requested": config["holdout_fraction"],
                "split_unit": "global homology cluster across all three windows, including weak/unknown records", "cluster_stats": cluster_stats,
                "cross_split_cluster_overlap": 0, "DBTBS": dbtbs_status, "config": config,
                "source_sha256": {config["regulondb_gff3"]:sha(gff),config["genome_fasta"]:sha(genome_path)},
                "git_head": subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                "git_dirty": bool(subprocess.check_output(["git","status","--porcelain"],cwd=ROOT)),
                "software": {"rapidfuzz": __import__("rapidfuzz").__version__, "python_environment": collect_environment()},
                "limits": ["No motif discovery/scan/classifier was run.", "Gold confidence supports source annotation, not wet-lab confirmation of each new motif.",
                           "Clusters use stated global-edit/shifted-overlap criteria, not arbitrary local alignments.", "Other-sigma controls use frozen GC strata without replacement; unmatched primaries are reported.",
                           "External datasets have not yet been certified disjoint from discovery; DBTBS coordinates remain missing."],
                "output_sha256": {p.relative_to(data_dir).as_posix():sha(p) for p in data_dir.rglob("*") if p.is_file()}}
    write_json(data_dir/"run_manifest.json", manifest)
    public = {k:v for k,v in manifest.items() if k not in {"software","output_sha256"}}
    public["data_manifest_sha256"] = sha(data_dir/"run_manifest.json")
    public["summary_sha256"] = {p.name:sha(p) for p in result_dir.glob("*.tsv")}
    write_json(result_dir/"rebuild_summary.json", public)
    print(json.dumps({k:manifest[k] for k in ["source_records","excluded_records","verified_loci","clusters","gold_representatives","cross_split_cluster_overlap"]},indent=2))
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=ROOT/"configs/g0_rebuild_20261004.json")
    run(parser.parse_args().config)
