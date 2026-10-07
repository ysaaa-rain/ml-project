# -*- coding: utf-8 -*-
"""
乙 - 步骤 8：按 σ 组的位置与间距明细（补 M1/M2 报告框架 表 4-6 / 4-7 / 4-8）

为什么单独一个脚本：
  脚本 04 的产出（T8~T15）已经通过 174 项数字核对，不再改动。
  本脚本只**新增** T19 / T20 与 F10 / F11，不覆盖任何既有产物。

对应报告框架的哪几张表：
  表 4-6 -10 box 位置分布（仅 RegulonDB）→ T19 中 element = -10 box
  表 4-7 -35 box 位置分布（仅 RegulonDB）→ T19 中 element = -35 box
  表 4-8 -10/-35 间距分布（仅 RegulonDB）→ T20
  图 4-1 序列长度分布                   → F10（长度是构建时固定值，图上会注明）
  图 4-4 -10/-35 间距分布               → F11

口径（与脚本 04 完全一致，不新立标准）：
  * 元件串：-35 = TTGACA，-10 = TATAAT；期望窗口 -35 `[-37,-32]`、-10 `[-12,-7]`
  * 主口径精确匹配（0 错配）；≤1 错配仅作敏感性
  * 只扫正向；只用 discovery；背景仍用配对的 dinucleotide_null（本脚本只在需要时用）
  * **每条序列每个元件只取一个代表命中**：先取落在期望窗口内的，若没有则取全体，
    再按「离期望窗口中心最近」挑选（-35 中心 -34.5，-10 中心 -9.5；并列取更靠上游者）。
    这条规则是为了让"位置"的均值/中位数不会被同一条序列的多个偶然命中拉偏。
  * 间距 = 两个代表命中之间夹的碱基数 = s10 - s35 - 6（教科书 ≈17）
  * σ 分组仍是 core 层单标签、discovery 划分

产出：乙_产出/tables/T19、T20；figures/F10、F11；out/positions_by_sigma.json
"""
import collections
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 130

# (显示名, 元件串, 期望窗口起, 期望窗口止, 期望窗口中心)
ELEMENTS = [("-35 box (TTGACA)", "TTGACA", -37, -32, -34.5),
            ("-10 box (TATAAT)", "TATAAT", -12, -7, -9.5)]

SPACER_LO, SPACER_HI = -5, 40


def mismatch_variants(motif, k):
    out, cur = {motif}, {motif}
    for _ in range(k):
        nxt = set()
        for s in cur:
            for i in range(len(s)):
                for b in "ACGT":
                    if b != s[i]:
                        nxt.add(s[:i] + b + s[i + 1:])
        cur = nxt
        out |= cur
    return out


def hit_offsets(seq, motifs):
    offs = set()
    for m in motifs:
        for mt in re.finditer("(?=(%s))" % m, seq):
            offs.add(C.offset_of_index(mt.start() + 1))
    return sorted(offs)


def representative(offsets, lo, hi, center):
    if not offsets:
        return None
    inwin = [o for o in offsets if lo <= o <= hi]
    pool = inwin if inwin else offsets
    return min(pool, key=lambda o: (abs(o - center), o))


def desc(vals):
    a = np.asarray([v for v in vals if v is not None], dtype=float)
    if a.size == 0:
        return {"n": 0, "mean": None, "sd": None, "median": None, "q1": None,
                "q3": None, "iqr": None, "min": None, "max": None, "range": None}
    return {"n": int(a.size), "mean": float(a.mean()),
            "sd": float(a.std(ddof=1)) if a.size > 1 else 0.0,
            "median": float(np.median(a)),
            "q1": float(np.percentile(a, 25)), "q3": float(np.percentile(a, 75)),
            "iqr": float(np.percentile(a, 75) - np.percentile(a, 25)),
            "min": float(a.min()), "max": float(a.max()),
            "range": float(a.max() - a.min())}


def main():
    base = os.path.join(C.DATA, "regulondb_ecoli", "core", "sigma")
    sigs = sorted(d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d)))

    rows19, rows20, result = [], [], {"elements": {}, "per_sigma": {}}

    for sig in sigs:
        pp = os.path.join(base, sig, "discovery.positive.fasta")
        if not os.path.exists(pp):
            continue
        seqs = [s for _h, s in C.read_fasta(pp)]
        rec = {"n_sequences": len(seqs), "elements": {}}

        # ---- 位置统计（两种匹配口径）----
        for elem, motif, lo, hi, center in ELEMENTS:
            for tag, k in [("exact", 0), ("mis1", 1)]:
                pats = mismatch_variants(motif, k)
                reps, n_all_matches = [], 0
                for s in seqs:
                    offs = hit_offsets(s, pats)
                    n_all_matches += len(offs)
                    r = representative(offs, lo, hi, center)
                    if r is not None:
                        reps.append(r)
                d = desc(reps)
                rec["elements"]["%s|%s" % (elem, tag)] = dict(d, n_matches=n_all_matches,
                                                              n_sequences=len(seqs))
                rows19.append([sig, elem, tag, len(seqs), d["n"], n_all_matches,
                               C.fmt(d["mean"], 2), C.fmt(d["median"], 1),
                               C.fmt(d["q1"], 1), C.fmt(d["q3"], 1),
                               C.fmt(d["sd"], 2), C.fmt(d["iqr"], 1),
                               C.fmt(d["min"], 0), C.fmt(d["max"], 0),
                               C.fmt(d["range"], 0),
                               C.fmt(100.0 * d["n"] / len(seqs), 2) if seqs else "NA"])

        # ---- 间距统计（两种匹配口径）----
        for tag, k in [("exact", 0), ("mis1", 1)]:
            p35 = mismatch_variants("TTGACA", k)
            p10 = mismatch_variants("TATAAT", k)
            spacers, both, only35, only10, neither = [], 0, 0, 0, 0
            for s in seqs:
                r35 = representative(hit_offsets(s, p35), -37, -32, -34.5)
                r10 = representative(hit_offsets(s, p10), -12, -7, -9.5)
                if r35 is not None and r10 is not None:
                    both += 1
                    sp = r10 - r35 - 6
                    if SPACER_LO <= sp <= SPACER_HI:
                        spacers.append(sp)
                elif r35 is not None:
                    only35 += 1
                elif r10 is not None:
                    only10 += 1
                else:
                    neither += 1
            d = desc(spacers)
            typical = sum(1 for x in spacers if 16 <= x <= 18)      # 教科书 17±1
            wide = sum(1 for x in spacers if 15 <= x <= 19)         # 教科书 17±2
            rec.setdefault("spacing", {})[tag] = dict(
                d, both=both, only_35=only35, only_10=only10, neither=neither,
                typical_16_18=typical, typical_15_19=wide,
                atypical_frac=(1.0 - typical / d["n"]) if d["n"] else None,
                atypical_frac_wide=(1.0 - wide / d["n"]) if d["n"] else None)
            rows20.append([sig, tag, len(seqs), both, only35, only10, neither, d["n"],
                           C.fmt(d["mean"], 2), C.fmt(d["median"], 1),
                           C.fmt(d["q1"], 1), C.fmt(d["q3"], 1), C.fmt(d["sd"], 2),
                           C.fmt(d["min"], 0), C.fmt(d["max"], 0),
                           typical, wide,
                           C.fmt(100.0 * (1.0 - typical / d["n"]), 1) if d["n"] else "NA",
                           C.fmt(100.0 * (1.0 - wide / d["n"]), 1) if d["n"] else "NA"])
        result["per_sigma"][sig] = rec

    C.write_tsv(os.path.join(C.TABLES, "T19_位置统计_按sigma组发现集.tsv"),
                ["sigma", "element", "match_rule", "n_sequences", "n_with_hit",
                 "n_matches_total", "mean_offset", "median_offset", "q1", "q3",
                 "sd", "iqr", "min", "max", "range", "hit_rate_pct"], rows19)
    C.write_tsv(os.path.join(C.TABLES, "T20_间距统计_按sigma组发现集.tsv"),
                ["sigma", "match_rule", "n_sequences", "both_elements", "only_35",
                 "only_10", "neither", "n_spacers", "mean_spacer", "median_spacer",
                 "q1", "q3", "sd", "min", "max",
                 "n_in_16_18", "n_in_15_19", "atypical_pct_17p1", "atypical_pct_17p2"],
                rows20)

    # ================= F10 序列长度分布 =================
    sm = C.load_sequence_master()
    lens = collections.Counter(len(r["sequence"]) for r in sm)
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    ax.bar(list(lens.keys()), list(lens.values()), width=0.8, color="#2E75B6")
    ax.set_xlim(60, 100)
    ax.set_xlabel("sequence length (bp)")
    ax.set_ylabel("sequences")
    ax.set_title("Sequence length distribution (all %d records)\n"
                 "single value 81 bp - length is fixed by construction, not observed"
                 % len(sm), fontsize=10)
    ax.annotate("all %d records = 81 bp\n(window [-60,+20] around TSS)" % len(sm),
                xy=(81, lens.get(81, 0)), xytext=(86, len(sm) * 0.65),
                arrowprops=dict(arrowstyle="->", color="crimson"), color="crimson",
                fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    f10 = os.path.join(C.FIGURES, "F10_序列长度分布_常量为81bp.png")
    fig.savefig(f10); plt.close(fig)

    # ================= F11 间距分布 =================
    # 图里用**全部** σ 组，与 T20 的分组合计一致（避免图与表对不上）
    ref_sigs = [s for s in sigs if s in result["per_sigma"]]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    allsp = {}
    for tag, label, ax in [("exact", "exact match (0 mismatch)", axes[0]),
                           ("mis1", "<=1 mismatch (sensitivity)", axes[1])]:
        pool = []
        p35 = mismatch_variants("TTGACA", 0 if tag == "exact" else 1)
        p10 = mismatch_variants("TATAAT", 0 if tag == "exact" else 1)
        for sig in ref_sigs:
            pp = os.path.join(base, sig, "discovery.positive.fasta")
            for _h, s in C.read_fasta(pp):
                r35 = representative(hit_offsets(s, p35), -37, -32, -34.5)
                r10 = representative(hit_offsets(s, p10), -12, -7, -9.5)
                if r35 is not None and r10 is not None:
                    sp = r10 - r35 - 6
                    if SPACER_LO <= sp <= SPACER_HI:
                        pool.append(sp)
        allsp[tag] = pool
        ax.hist(pool, bins=range(SPACER_LO, SPACER_HI + 2), color="#548235",
                edgecolor="white", linewidth=0.4)
        ax.axvline(17, color="crimson", linestyle="--", linewidth=1.4,
                   label="textbook 17 bp")
        ax.set_title("%s\nn = %d spacers" % (label, len(pool)), fontsize=10)
        ax.set_xlabel("spacer between -35 and -10 hexamers (bp)")
        ax.set_ylabel("sequences")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("RegulonDB E. coli core, discovery: -35/-10 spacing "
                 "(one representative hit per element per sequence)")
    fig.tight_layout()
    f11 = os.path.join(C.FIGURES, "F11_间距分布_RegulonDBcore.png")
    fig.savefig(f11); plt.close(fig)

    result["length_distribution"] = {str(k): v for k, v in sorted(lens.items())}
    result["spacing_overall"] = {k: {"n": len(v), "min": min(v) if v else None,
                                     "max": max(v) if v else None,
                                     "median": float(np.median(v)) if v else None}
                                 for k, v in allsp.items()}
    result["figures"] = [os.path.basename(f10), os.path.basename(f11)]
    with open(os.path.join(C.OUT, "positions_by_sigma.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)

    # ================= 控制台 =================
    print("=== 长度分布（图 4-1）===")
    print("  取值: %s（%d 条全部相同）" % (dict(lens), sum(lens.values())))
    print("  结论：长度是构建时固定的窗口，没有分布。图上已注明。")
    print("\n=== 位置统计（精确口径，代表性命中；图例：偏移 bp）===")
    print("  %-9s %-18s %5s %5s %8s %8s %8s %8s" %
          ("sigma", "element", "n", "hit", "mean", "median", "IQR", "range"))
    for r in rows19:
        if r[2] != "exact":
            continue
        print("  %-9s %-18s %5d %5d %8s %8s %8s %8s"
              % (r[0], r[1], int(r[3]), int(r[4]), r[6], r[7], r[11], r[14]))
    print("\n=== 间距统计（两六聚体之间夹的碱基数，教科书≈17）===")
    print("  %-9s %-6s %5s %5s %8s %8s %8s %10s %10s" %
          ("sigma", "rule", "n_seq", "both", "mean", "median", "IQR",
           "非典型(17±1)", "非典型(17±2)"))
    for r in rows20:
        iqr = "%s-%s" % (r[10], r[11]) if r[7] != "0" else "NA"
        print("  %-9s %-6s %5d %5d %8s %8s %8s %10s %10s"
              % (r[0], r[1], int(r[2]), int(r[3]), r[8], r[9], iqr, r[17], r[18]))
    print("\n  全体（%d 个 RegulonDB σ 组合计）：精确 = %d 对，≤1 错配 = %d 对"
          % (len(ref_sigs), len(allsp["exact"]), len(allsp["mis1"])))
    print("\n图：%s / %s" % (os.path.basename(f10), os.path.basename(f11)))


if __name__ == "__main__":
    main()
