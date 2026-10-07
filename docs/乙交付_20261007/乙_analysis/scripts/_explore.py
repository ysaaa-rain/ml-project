# -*- coding: utf-8 -*-
"""乙 - 数据勘察：确认哪些列能给哪些分析提供依据（探索用，非交付脚本）"""
import collections
import csv
import io
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
DATA = os.path.join(BASE, "data", "data", "processed", "pr01_02_data_v2")


def rows(path):
    with io.open(path, encoding="utf-8", errors="replace", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            yield r


def main():
    sm = list(rows(os.path.join(DATA, "sequence_master.tsv")))
    print("sequence_master 行数 = %d" % len(sm))

    srcs = collections.Counter(r["source"] for r in sm)
    print("\n=== 按 source ===")
    for k, v in sorted(srcs.items(), key=lambda x: -x[1]):
        print("  %-28s %6d" % (k, v))

    print("\n=== source × species ===")
    for k, v in sorted(collections.Counter((r["source"], r["species"]) for r in sm).items()):
        print("  %-28s %-22s %6d" % (k[0], k[1], v))

    print("\n=== source × split ===")
    for k, v in sorted(collections.Counter((r["source"], r["split"]) for r in sm).items()):
        print("  %-28s %-22s %6d" % (k[0], k[1], v))

    print("\n=== source × tier ===")
    for k, v in sorted(collections.Counter((r["source"], r["tier"]) for r in sm).items()):
        print("  %-28s %-10s %6d" % (k[0], k[1], v))

    print("\n=== source × label ===")
    for k, v in sorted(collections.Counter((r["source"], r["label"]) for r in sm).items()):
        print("  %-28s %-4s %6d" % (k[0], k[1], v))

    print("\n=== 坐标证据可用性（按 source）===")
    by = collections.defaultdict(lambda: {"n": 0, "tss_nonempty": 0, "strand_nonempty": 0,
                                          "assembly_nonempty": 0, "ev": collections.Counter()})
    for r in sm:
        d = by[r["source"]]
        d["n"] += 1
        gtv = (r["genomic_tss_values"] or "").strip()
        gs = (r["genomic_strands"] or "").strip()
        asm = (r["assembly"] or "").strip()
        if gtv and gtv not in ("[]", '[""]', ""):
            d["tss_nonempty"] += 1
        if gs and gs not in ("[]", '[""]', ""):
            d["strand_nonempty"] += 1
        if asm:
            d["assembly_nonempty"] += 1
        d["ev"][r["coordinate_evidence"]] += 1
    for s in sorted(by):
        d = by[s]
        print("  %-28s n=%6d  TSS非空=%6d (%.1f%%)  链非空=%6d (%.1f%%)  assembly非空=%6d"
              % (s, d["n"], d["tss_nonempty"], 100.0 * d["tss_nonempty"] / d["n"],
                 d["strand_nonempty"], 100.0 * d["strand_nonempty"] / d["n"], d["assembly_nonempty"]))
        for e, c in d["ev"].most_common():
            print("        evidence: %-52s %6d" % (e, c))

    print("\n=== length / local_tss_index / window_offsets 取值 ===")
    print("  length:", collections.Counter(r["length"] for r in sm).most_common(5))
    print("  local_tss_index:", collections.Counter(r["local_tss_index"] for r in sm).most_common(5))
    print("  window_offsets:", collections.Counter(r["window_offsets"] for r in sm).most_common(5))

    print("\n=== 序列碱基组成（非法字符检查）===")
    bad = collections.Counter()
    lens = collections.Counter()
    for r in sm:
        s = r["sequence"].upper()
        lens[len(s)] += 1
        for ch in set(s):
            if ch not in "ACGT":
                bad[ch] += 1
    print("  长度分布:", lens.most_common(5))
    print("  非法字符:", bad.most_common(10) or "无")

    print("\n=== sigma_associations 多标签情况 ===")
    sa = list(rows(os.path.join(DATA, "sigma_associations.tsv")))
    print("  行数 = %d" % len(sa))
    per_seq = collections.defaultdict(set)
    for r in sa:
        per_seq[r["sequence_id"]].add(r["sigma_id"])
    dist = collections.Counter(len(v) for v in per_seq.values())
    print("  每条序列的 sigma 数分布:", sorted(dist.items()))
    print("  sigma_raw 取值:", collections.Counter(r["sigma_raw"] for r in sa).most_common(30))
    print("  sigma_association_evidence:",
          collections.Counter(r["sigma_association_evidence"] for r in sa).most_common(5))
    multi = sum(1 for v in per_seq.values() if len(v) > 1)
    print("  多标签序列数 = %d，单标签序列数 = %d" % (multi, len(per_seq) - multi))


if __name__ == "__main__":
    main()
