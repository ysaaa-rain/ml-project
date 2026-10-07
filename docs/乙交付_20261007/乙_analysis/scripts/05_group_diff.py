# -*- coding: utf-8 -*-
"""
乙 - 步骤 5：σ 组间差异比较（必须报"差多少"，不能只报显著性）

口径：
  * 只用 discovery 划分（holdout 一律不碰）
  * 只用单标签 σ 序列（core 层），与 sigma_group_summary.tsv 口径一致
  * 只在**同一来源内部**比较。跨来源（RegulonDB vs DBTBS）不比：
    两者物种不同、σ 命名体系不同，直接比会把物种效应当成 σ 效应
    （甲_σ组可用性清单 口径说明 2）
  * 参照组：RegulonDB 内用 Sigma70，DBTBS 内用 SigA（各自样本量最大）
  * 可用性分级沿用甲_σ组可用性清单：
        可做正式结论 → 报完整统计量
        只能探索     → 报统计量，但全部标"探索性"
        只能描述     → 只报样本量与描述统计（n/均值/中位数），**不报 p 值、不报效应量**
  * 效应量（回答"差多少"）：
        GC 含量（连续）  → 中位数差（百分点）+ bootstrap 95%CI + 秩双列相关 r
        元件命中率（比例）→ 风险差（百分点）+ 95%CI + 比值比
  * 多重比较用 Benjamini–Hochberg 校正，同时给出原始 p 与 q

产出：乙_产出/tables/T16~T18、figures/F9、out/group_diff.json
"""
import collections
import json
import os
import re
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 130

RNG = np.random.default_rng(20261005)
N_BOOT = 5000

ELEMENTS = [("-35 box (TTGACA)", "TTGACA"), ("-10 box (TATAAT)", "TATAAT")]
REFERENCE = {"regulondb_ecoli": "Sigma70", "dbtbs_bsub": "SigA"}
INFERENTIAL = {"confirmation_candidate": True, "exploratory_only": True,
               "descriptive_only": False}


def bh(pvals):
    """Benjamini-Hochberg q 值。"""
    p = np.asarray(pvals, float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    q = np.empty(n, float)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = n - rank
        val = min(prev, p[idx] * n / i)
        q[idx] = val
        prev = val
    return q


def rank_biserial(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.size == 0 or y.size == 0:
        return float("nan")
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return float(2.0 * u / (x.size * y.size) - 1.0)


def median_diff_boot(x, y, n_boot=N_BOOT):
    """中位数差（x − y）与 bootstrap 95%CI。"""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.size < 2 or y.size < 2:
        return float("nan"), float("nan"), float("nan")
    obs = float(np.median(x) - np.median(y))
    bx = RNG.choice(x, size=(n_boot, x.size), replace=True)
    by = RNG.choice(y, size=(n_boot, y.size), replace=True)
    diffs = np.median(bx, axis=1) - np.median(by, axis=1)
    return obs, float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def risk_diff_ci(a_hit, a_n, b_hit, b_n):
    ra = a_hit / a_n if a_n else float("nan")
    rb = b_hit / b_n if b_n else float("nan")
    rd = ra - rb
    se = float("nan")
    if a_n and b_n:
        se = np.sqrt(ra * (1 - ra) / a_n + rb * (1 - rb) / b_n)
    return rd, rd - 1.96 * se, rd + 1.96 * se


def hits(recs, motif):
    pat = re.compile("(?=(%s))" % motif)
    return [1 if pat.search(r["sequence"]) else 0 for r in recs]


def main():
    sm = C.load_sequence_master()
    sig_idx = C.build_sigma_index(C.load_sigma_associations())
    rec_by_id = {r["sequence_id"]: r for r in sm}
    summary = {(r["source"], r["sigma"]): r for r in C.load_sigma_group_summary()}

    # ---------- 组装分组 ----------
    groups = collections.defaultdict(list)
    for sid, v in sig_idx.items():
        rec = rec_by_id.get(sid)
        if rec is None or rec["source"].startswith("tjupan_"):
            continue
        if rec["tier"] != "core" or rec["split"] != "discovery":
            continue
        lab = C.single_sigma_label(v["sigma_raw"])
        if lab is None:
            continue
        groups[(rec["source"], lab)].append(rec)

    elig = lambda src, sig: summary.get((src, sig), {}).get("eligibility", "")
    elig_cn = lambda src, sig: C.ELIGIBILITY_CN.get(elig(src, sig), "")

    # ---------- T18 各组描述统计（全部 19 组） ----------
    rows18 = []
    for (src, sig) in sorted(groups, key=lambda k: (k[0], -len(groups[k]))):
        recs = groups[(src, sig)]
        gc = [C.gc_percent(r["sequence"]) for r in recs]
        row = [src, C.SOURCE_CN.get(src, src), sig, len(recs),
               C.fmt(np.mean(gc), 2), C.fmt(np.median(gc), 2),
               C.fmt(np.percentile(gc, 25), 2), C.fmt(np.percentile(gc, 75), 2)]
        for name, motif in ELEMENTS:
            h = hits(recs, motif)
            row.append("%d/%d" % (sum(h), len(h)))
            row.append(C.fmt(100.0 * sum(h) / len(h), 2))
        row += [elig(src, sig), elig_cn(src, sig)]
        rows18.append(row)
    C.write_tsv(os.path.join(C.TABLES, "T18_sigma组描述统计_发现集.tsv"),
                ["source", "source_cn", "sigma", "n_discovery", "gc_mean", "gc_median",
                 "gc_q1", "gc_q3",
                 "minus35_hit_seqs", "minus35_hit_pct", "minus10_hit_seqs",
                 "minus10_hit_pct", "eligibility", "eligibility_cn"], rows18)

    # ---------- T16 GC 组间比较（组 vs 同来源参照组） ----------
    rows16, gc_tests = [], []
    for src in ("regulondb_ecoli", "dbtbs_bsub"):
        ref = REFERENCE[src]
        if (src, ref) not in groups:
            continue
        rgc = [C.gc_percent(r["sequence"]) for r in groups[(src, ref)]]
        for (s2, sig) in sorted(groups):
            if s2 != src or sig == ref:
                continue
            gc = [C.gc_percent(r["sequence"]) for r in groups[(src, sig)]]
            e = elig(src, sig)
            do_inf = INFERENTIAL.get(e, False)
            if do_inf:
                obs, lo, hi = median_diff_boot(gc, rgc)
                u, p = stats.mannwhitneyu(gc, rgc, alternative="two-sided")
                r = rank_biserial(gc, rgc)
                gc_tests.append((src, sig, p))
            else:
                obs = float(np.median(gc) - np.median(rgc))
                lo = hi = p = r = float("nan")
            rows16.append([src, C.SOURCE_CN.get(src, src), sig, len(gc),
                           C.fmt(np.median(gc), 2), ref, len(rgc),
                           C.fmt(np.median(rgc), 2),
                           C.fmt(obs, 2), C.fmt(lo, 2), C.fmt(hi, 2),
                           C.fmt(r, 3) if do_inf else "—",
                           C.fmt_p(p) if do_inf else "—", "", e, elig_cn(src, sig),
                           "是" if do_inf else "否"])
    # BH 校正（按族：GC 比较）
    if gc_tests:
        qs = bh([t[2] for t in gc_tests])
        qmap = {(t[0], t[1]): qs[i] for i, t in enumerate(gc_tests)}
        for row in rows16:
            key = (row[0], row[2])
            if key in qmap:
                row[13] = C.fmt_p(qmap[key])
    C.write_tsv(os.path.join(C.TABLES, "T16_GC组间比较_vs参照组.tsv"),
                ["source", "source_cn", "sigma", "n", "gc_median",
                 "reference_sigma", "reference_n", "reference_gc_median",
                 "median_diff_pp", "ci95_lo_pp", "ci95_hi_pp",
                 "rank_biserial_r", "mannwhitney_p", "bh_q", "eligibility",
                 "eligibility_cn", "inferential_reported"], rows16)

    # ---------- T17 元件命中率组间比较 ----------
    rows17, hit_tests = [], []
    for src in ("regulondb_ecoli", "dbtbs_bsub"):
        ref = REFERENCE[src]
        if (src, ref) not in groups:
            continue
        for name, motif in ELEMENTS:
            rh = hits(groups[(src, ref)], motif)
            for (s2, sig) in sorted(groups):
                if s2 != src or sig == ref:
                    continue
                h = hits(groups[(src, sig)], motif)
                e = elig(src, sig)
                do_inf = INFERENTIAL.get(e, False)
                ah, an = sum(h), len(h)
                bh_, bn = sum(rh), len(rh)
                if do_inf:
                    rd, lo, hi = risk_diff_ci(ah, an, bh_, bn)
                    odds, p = stats.fisher_exact([[ah, an - ah], [bh_, bn - bh_]])
                    hit_tests.append((src, sig, name, p))
                else:
                    rd = (ah / an if an else float("nan")) - (bh_ / bn if bn else float("nan"))
                    lo = hi = p = odds = float("nan")
                rows17.append([src, C.SOURCE_CN.get(src, src), sig, name, an, ah,
                               C.fmt(100.0 * ah / an, 2) if an else "NA",
                               ref, bn, bh_,
                               C.fmt(100.0 * bh_ / bn, 2) if bn else "NA",
                               C.fmt(100 * rd, 2), C.fmt(100 * lo, 2), C.fmt(100 * hi, 2),
                               C.fmt(odds, 3) if do_inf else "—",
                               C.fmt_p(p) if do_inf else "—", "",
                               e, elig_cn(src, sig), "是" if do_inf else "否"])
    if hit_tests:
        qs = bh([t[3] for t in hit_tests])
        qmap = {(t[0], t[1], t[2]): qs[i] for i, t in enumerate(hit_tests)}
        for row in rows17:
            key = (row[0], row[2], row[3])
            if key in qmap:
                row[16] = C.fmt_p(qmap[key])
    C.write_tsv(os.path.join(C.TABLES, "T17_元件命中率组间比较_vs参照组.tsv"),
                ["source", "source_cn", "sigma", "element", "n", "hit_n", "hit_pct",
                 "reference_sigma", "reference_n", "reference_hit_n", "reference_hit_pct",
                 "risk_diff_pp", "ci95_lo_pp", "ci95_hi_pp", "odds_ratio",
                 "fisher_p", "bh_q", "eligibility", "eligibility_cn",
                 "inferential_reported"], rows17)

    # ---------- F9 效应量图（GC 中位数差 + 命中率风险差） ----------
    plot_rows = [r for r in rows16 if r[16] == "是"]
    src_abbr = {"regulondb_ecoli": "RDB", "dbtbs_bsub": "DBTBS"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    if plot_rows:
        y = np.arange(len(plot_rows))
        obs = [float(r[8]) for r in plot_rows]
        lo = [float(r[9]) for r in plot_rows]
        hi = [float(r[10]) for r in plot_rows]
        axes[0].errorbar(obs, y, xerr=[np.array(obs) - np.array(lo),
                                       np.array(hi) - np.array(obs)],
                         fmt="o", color="#2E75B6", capsize=3)
        axes[0].axvline(0, color="crimson", linestyle="--", linewidth=1.1)
        axes[0].set_yticks(y)
        axes[0].set_yticklabels(["%s %s (%s)\nn=%s vs %s %s"
                                 % (src_abbr.get(r[0], r[0]), r[2], r[15], r[3],
                                    src_abbr.get(r[0], r[0]), r[5]) for r in plot_rows],
                                fontsize=8)
        axes[0].set_xlabel("median GC difference vs reference (percentage points)")
        axes[0].set_title("GC content: effect size with 95% CI")
        axes[0].grid(axis="x", alpha=0.3)
    p17 = [r for r in rows17 if r[19] == "是" and r[3].startswith("-10")]
    if p17:
        y = np.arange(len(p17))
        obs = [float(r[11]) for r in p17]
        lo = [float(r[12]) for r in p17]
        hi = [float(r[13]) for r in p17]
        axes[1].errorbar(obs, y, xerr=[np.array(obs) - np.array(lo),
                                       np.array(hi) - np.array(obs)],
                         fmt="s", color="#548235", capsize=3)
        axes[1].axvline(0, color="crimson", linestyle="--", linewidth=1.1)
        axes[1].set_yticks(y)
        axes[1].set_yticklabels(["%s %s (%s)\nn=%s vs %s %s"
                                 % (src_abbr.get(r[0], r[0]), r[2], r[18], r[4],
                                    src_abbr.get(r[0], r[0]), r[7]) for r in p17],
                                fontsize=8)
        axes[1].set_xlabel("-10 box hit-rate risk difference vs reference (pp)")
        axes[1].set_title("-10 box exact-match hit rate: effect size with 95% CI")
        axes[1].grid(axis="x", alpha=0.3)
    fig.suptitle("Sigma-group differences vs within-source reference "
                 "(discovery split, single-sigma core records)")
    fig.tight_layout()
    f9 = os.path.join(C.FIGURES, "F9_组间效应量.png")
    fig.savefig(f9); plt.close(fig)

    # ---------- 汇总 ----------
    result = {
        "reference": REFERENCE,
        "groups_n": {("%s/%s" % k): len(v) for k, v in groups.items()},
        "gc_comparisons": rows16,
        "hit_comparisons": rows17,
        "note": "跨来源不作比较：RegulonDB(大肠杆菌) 与 DBTBS(枯草芽孢杆菌) 物种与 σ 命名体系不同",
    }
    with open(os.path.join(C.OUT, "group_diff.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)

    # ---------- 控制台 ----------
    print("=== 分组（core 单标签 discovery）===")
    for (src, sig) in sorted(groups, key=lambda k: (k[0], -len(groups[k]))):
        print("  %-16s %-9s n=%4d  [%s]" % (src, sig, len(groups[(src, sig)]),
                                            elig_cn(src, sig)))

    print("\n=== GC 组间比较（vs 同来源参照组；单位：百分点）===")
    for r in rows16:
        if r[16] == "是":
            print("  %-16s %-9s n=%4d 中位数 %6s%%  vs %-8s %6s%%  差 %+7s [%s, %s]  "
                  "r=%6s  p=%-9s q=%-9s  [%s]"
                  % (r[0], r[2], int(r[3]), r[4], r[5], r[7], r[8], r[9], r[10],
                     r[11], r[12], r[13], r[15]))
        else:
            print("  %-16s %-9s n=%4d 中位数 %6s%%  （%s：只描述，不做统计）"
                  % (r[0], r[2], int(r[3]), r[4], r[15]))

    print("\n=== 元件命中率组间比较（精确口径，vs 同来源参照组；单位：百分点）===")
    for r in rows17:
        if r[19] == "是":
            print("  %-16s %-9s %-18s n=%4d %6s%% vs %6s%%  风险差 %+7s [%s, %s]  "
                  "OR=%-7s p=%-9s q=%-9s"
                  % (r[0], r[2], r[3], int(r[4]), r[6], r[10], r[11], r[12], r[13],
                     r[14], r[15], r[16]))
        else:
            print("  %-16s %-9s %-18s n=%4d %6s%% vs %6s%%  （%s：只描述，不做统计）"
                  % (r[0], r[2], r[3], int(r[4]), r[6], r[10], r[18]))
    print("\n图：%s" % os.path.basename(f9))


if __name__ == "__main__":
    main()
