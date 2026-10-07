# -*- coding: utf-8 -*-
"""
乙 - 步骤 4：已知元件（-10 / -35 box）位置分布 + 背景校准

口径（**预先固定**，不因结果调整；全部写进报告）：

  1. 元件串：教科书固定串 -35 box = TTGACA，-10 box = TATAAT。
     这是"已知答案"的比对，不是本实验发现的 motif（M2 不找暗号）。
  2. 位置约定：沿用本项目既有约定（见 repo
     results/motif/tjupan_m3_20261003/reference_elements/ecoli_sigma70_expected_window_tests.tsv）：
         -35 box 期望窗口 = 序列内 1-based 第 24–29 位 → TSS 相对偏移 [-37, -32]
         -10 box 期望窗口 = 序列内 1-based 第 49–54 位 → TSS 相对偏移 [-12, -7]
     与本地 81 bp、TSS 第 61 位的窗口定义完全一致（偏移 = 下标 - 61）。
  3. 匹配：主口径为**精确匹配（0 错配）**；另做 ≤1 错配的敏感性分析。
  4. 链向：**只扫正向**。本版序列已统一为转录方向（甲_数据字典 · 字段口径说明），
     再扫反向互补会把错误方向的命中当证据。
  5. 位置：命中起始位换算为 TSS 相对偏移；TSS 本身为 0，上游为负。
  6. 统计划分：只用 discovery；**holdout 一律不碰**。
  7. 背景：直接用数据自带的 *.dinucleotide_null.fasta —— 与正样本逐文件配对、同划分、
     同长度、二核苷酸组成保持（GC 自检差值 0.000 pp）。不自己另造背景。
  8. 命中率对所有序列乘 100；"超出背景"= 正样本命中率 − 打乱背景命中率（百分点）。

主结果 = RegulonDB core（按问答记录，位置分析只认 RegulonDB）。
备查附录 = DBTBS core 与 TJU（按仓库 2026-10-05 规范本可做，但与问答记录冲突，
           故只作为"不得用于结论"的附录）。

产出：乙_产出/tables/T8~T15、figures/F4~F8、out/positions.json
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

BASES = ["A", "C", "G", "T"]

# (显示名, 元件串, 期望窗口起始偏移, 期望窗口结束偏移)
ELEMENTS = [("-35 box (TTGACA)", "TTGACA", -37, -32),
            ("-10 box (TATAAT)", "TATAAT", -12, -7)]

SPACER_MIN, SPACER_MAX = -5, 40     # 只统计合理范围内的间距


# ---------------------------------------------------------------- 扫描
def mismatch_variants(motif, k):
    """返回与 motif 汉明距离 ≤ k 的全部串。"""
    out = {motif}
    cur = {motif}
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


def match_offsets(seq, motifs):
    """返回所有命中的 TSS 相对偏移（可重叠，去重）。"""
    offs = set()
    for m in motifs:
        for mt in re.finditer("(?=(%s))" % m, seq):
            offs.add(C.offset_of_index(mt.start() + 1))
    return sorted(offs)


def rate_and_calibration(pos, nul, motifs):
    """
    算一个元件在正样本 vs 打乱背景里的命中情况。
    返回 (hit-bearing 计数, 校准量, 每条序列的命中偏移列表)。
    """
    p_hits = [match_offsets(s, motifs) for s in pos]
    n_hits = [match_offsets(s, motifs) for s in nul]
    p_any = sum(1 for h in p_hits if h)
    n_any = sum(1 for h in n_hits if h)
    np_, nn = len(pos), len(nul)
    rp = p_any / np_ if np_ else float("nan")
    rn = n_any / nn if nn else float("nan")
    if np_ and nn:
        odds, pval = stats.fisher_exact([[p_any, np_ - p_any], [n_any, nn - n_any]])
        rd = rp - rn
        se = np.sqrt(rp * (1 - rp) / np_ + rn * (1 - rn) / nn)
        lo, hi = rd - 1.96 * se, rd + 1.96 * se
    else:
        odds = pval = rd = lo = hi = float("nan")
    return {
        "n_positive": np_, "n_null": nn,
        "positive_any": p_any, "null_any": n_any,
        "rate_positive": rp, "rate_null": rn,
        "excess": rd, "excess_lo": lo, "excess_hi": hi,
        "fold": (rp / rn) if rn else float("inf"),
        "fisher_p": pval, "odds_ratio": odds,
        "hits_per_seq_positive": (sum(len(h) for h in p_hits) / np_) if np_ else float("nan"),
        "hits_per_seq_null": (sum(len(h) for h in n_hits) / nn) if nn else float("nan"),
        "per_seq_hits_positive": p_hits,
        "offsets_positive": [o for h in p_hits for o in h],
        "offsets_null": [o for h in n_hits for o in h],
    }


def offset_summary(offsets, win_lo, win_hi):
    if not offsets:
        return {"n_matches": 0, "mode_offset": None, "mode_count": 0,
                "upstream_frac": float("nan"), "window_frac": float("nan"),
                "median_upstream": float("nan")}
    cnt = collections.Counter(offsets)
    mode_off, mode_n = cnt.most_common(1)[0]
    up = [o for o in offsets if o <= -1]
    in_win = [o for o in offsets if win_lo <= o <= win_hi]
    return {"n_matches": len(offsets), "mode_offset": mode_off, "mode_count": mode_n,
            "upstream_frac": len(up) / len(offsets),
            "window_frac": len(in_win) / len(offsets),
            "median_upstream": float(np.median(up)) if up else float("nan")}


def seqs_with_hit_in_window(per_seq_hits, win_lo, win_hi):
    return sum(1 for h in per_seq_hits if any(win_lo <= o <= win_hi for o in h))


# ---------------------------------------------------------------- 位置碱基组成
def position_information(pos_seqs, nul_seqs):
    """
    逐位置碱基组成与信息量。
      info_vs_uniform = 2 - H(位置碱基分布)
      info_vs_null    = sum_b p_b * log2(p_b / q_b)，q 为配对打乱背景在同位置的频率
    （后者控制了整体碱基组成，是判断"这个位置真的有偏好"的关键。）
    """
    L = len(pos_seqs[0])
    rows = []
    for idx in range(L):
        p = collections.Counter(s[idx] for s in pos_seqs)
        n = len(pos_seqs)
        pf = [p.get(b, 0) / n for b in BASES]
        h = -sum(x * np.log2(x) for x in pf if x > 0)
        q = collections.Counter(s[idx] for s in nul_seqs) if nul_seqs else None
        if q:
            m = len(nul_seqs)
            qf = [q.get(b, 0) / m for b in BASES]
        else:
            qf = [0.25] * 4
        info_null = 0.0
        for pb, qb in zip(pf, qf):
            if pb > 0 and qb > 0:
                info_null += pb * np.log2(pb / qb)
        rows.append({
            "offset": C.offset_of_index(idx + 1),
            "a": pf[0], "c": pf[1], "g": pf[2], "t": pf[3],
            "shannon_entropy_bits": h,
            "info_vs_uniform": 2.0 - h,
            "info_vs_null_bits": info_null,
        })
    return rows


def main():
    result = {"elements": [{"name": n, "motif": m, "window_offsets": [lo, hi]}
                           for n, m, lo, hi in ELEMENTS]}

    # ============ 主结果：RegulonDB core 全体（discovery） ============
    base = os.path.join(C.DATA, "regulondb_ecoli", "core")
    pos_all = [s for _h, s in C.read_fasta(os.path.join(base, "discovery.positive.fasta"))]
    nul_all = [s for _h, s in C.read_fasta(
        os.path.join(base, "discovery.dinucleotide_null.fasta"))]
    n_tot = len(pos_all)

    main_res = {}
    for name, motif, lo, hi in ELEMENTS:
        for tag, k in [("exact", 0), ("mis1", 1)]:
            motifs = mismatch_variants(motif, k)
            d = rate_and_calibration(pos_all, nul_all, motifs)
            d["n_patterns"] = len(motifs)
            d["summary"] = offset_summary(d["offsets_positive"], lo, hi)
            d["seqs_hit_in_window"] = seqs_with_hit_in_window(d["per_seq_hits_positive"], lo, hi)
            main_res["%s|%s" % (name, tag)] = d

    rows = []
    for name, motif, lo, hi in ELEMENTS:
        for tag, tag_cn in [("exact", "精确匹配(0错配)"), ("mis1", "≤1错配(敏感性)")]:
            d = main_res["%s|%s" % (name, tag)]
            s = d["summary"]
            rows.append([name, motif, "%d..%d" % (lo, hi), tag_cn, d["n_patterns"],
                         d["n_positive"], d["positive_any"],
                         C.fmt(100 * d["rate_positive"], 2),
                         d["n_null"], d["null_any"], C.fmt(100 * d["rate_null"], 2),
                         C.fmt(100 * d["excess"], 2),
                         "[%s, %s]" % (C.fmt(100 * d["excess_lo"], 2),
                                       C.fmt(100 * d["excess_hi"], 2)),
                         C.fmt(d["fold"], 2), C.fmt(d["odds_ratio"], 3),
                         C.fmt_p(d["fisher_p"]),
                         C.fmt(d["hits_per_seq_positive"], 3),
                         C.fmt(d["hits_per_seq_null"], 3),
                         s["mode_offset"], s["mode_count"],
                         C.fmt(100 * s["upstream_frac"], 1),
                         C.fmt(100 * s["window_frac"], 1),
                         d["seqs_hit_in_window"]])
    C.write_tsv(os.path.join(C.TABLES, "T8_已知元件命中率与背景校准_RegulonDBcore发现集.tsv"),
                ["element", "motif", "expected_window_offsets", "match_rule", "n_patterns",
                 "n_positive", "positive_with_hit", "positive_rate_pct",
                 "n_null", "null_with_hit", "null_rate_pct",
                 "excess_pp", "excess_95CI_pp", "fold", "odds_ratio", "fisher_p",
                 "hits_per_seq_positive", "hits_per_seq_null",
                 "mode_offset", "mode_count", "upstream_pct", "in_expected_window_pct",
                 "n_seqs_with_hit_in_expected_window"], rows)

    # 偏移分布（精确口径，全体）
    rows = []
    for name, motif, lo, hi in ELEMENTS:
        d = main_res["%s|exact" % name]
        cnt = collections.Counter(d["offsets_positive"])
        for o in sorted(cnt):
            rows.append([name, o, cnt[o], "IN_WINDOW" if lo <= o <= hi else ""])
    C.write_tsv(os.path.join(C.TABLES, "T9_命中位置分布_RegulonDBcore发现集.tsv"),
                ["element", "tss_relative_offset", "n_matches", "in_expected_window"], rows)

    # ============ 逐位置碱基组成与信息量（RegulonDB core 全体） ============
    info_rows = position_information(pos_all, nul_all)
    C.write_tsv(os.path.join(C.TABLES, "T10_逐位置碱基组成与信息量_RegulonDBcore发现集.tsv"),
                ["tss_relative_offset", "a_frac", "c_frac", "g_frac", "t_frac",
                 "shannon_entropy_bits", "info_vs_uniform_bits", "info_vs_null_bits"],
                [[r["offset"], C.fmt(r["a"], 4), C.fmt(r["c"], 4), C.fmt(r["g"], 4),
                  C.fmt(r["t"], 4), C.fmt(r["shannon_entropy_bits"], 4),
                  C.fmt(r["info_vs_uniform"], 4), C.fmt(r["info_vs_null_bits"], 4)]
                 for r in info_rows])
    info_by_off = {r["offset"]: r for r in info_rows}

    # ============ 期望窗口内逐位置最优碱基 vs 共识串 ============
    # 精确串命中率低，并不等于"那个位置没有碱基偏好"。
    # 这里把窗口内每个位置实际最常出现的碱基列出来，与共识串逐个位置比对。
    rows, consensus_check = [], {}
    for name, motif, lo, hi in ELEMENTS:
        obs, hits = [], 0
        for o in range(lo, hi + 1):
            r = info_by_off[o]
            freqs = {"A": r["a"], "C": r["c"], "G": r["g"], "T": r["t"]}
            best = max(freqs, key=lambda b: freqs[b])
            expect = motif[o - lo]
            ok = (best == expect)
            hits += 1 if ok else 0
            obs.append(best)
            rows.append([name, o, expect, best, C.fmt(100 * freqs[best], 2),
                         C.fmt(100 * freqs[expect], 2), "一致" if ok else "不一致",
                         C.fmt(r["info_vs_uniform"], 4), C.fmt(r["info_vs_null_bits"], 4)])
        consensus_check[name] = {"consensus": motif, "observed_argmax": "".join(obs),
                                 "positions": hi - lo + 1, "match_positions": hits}
        rows.append([name, "小计", motif, "".join(obs), "", "", "%d/%d" % (hits, hi - lo + 1),
                     "", ""])
    C.write_tsv(os.path.join(C.TABLES, "T10b_期望窗口逐位置最优碱基_vs_共识串.tsv"),
                ["element", "tss_relative_offset", "consensus_base", "observed_argmax_base",
                 "argmax_pct", "consensus_base_pct", "match", "info_vs_uniform_bits",
                 "info_vs_null_bits"], rows)

    # ============ σ 分组：命中率 + 校准 + 位置分布 ============
    rows, sigma_res, sigma_info = [], {}, {}
    sdir = os.path.join(base, "sigma")
    for sig in sorted(os.listdir(sdir)):
        gdir = os.path.join(sdir, sig)
        if not os.path.isdir(gdir):
            continue
        pp = os.path.join(gdir, "discovery.positive.fasta")
        npp = os.path.join(gdir, "discovery.dinucleotide_null.fasta")
        if not os.path.exists(pp):
            continue
        p = [s for _h, s in C.read_fasta(pp)]
        n = [s for _h, s in C.read_fasta(npp)] if os.path.exists(npp) else []
        rec = {"n_positive": len(p)}
        if p:
            sigma_info[sig] = position_information(p, n)
        for name, motif, lo, hi in ELEMENTS:
            d = rate_and_calibration(p, n, mismatch_variants(motif, 0))
            sm = offset_summary(d["offsets_positive"], lo, hi)
            rec[name] = {"rate_positive": d["rate_positive"], "rate_null": d["rate_null"],
                         "excess": d["excess"], "fold": d["fold"], "fisher_p": d["fisher_p"],
                         "odds_ratio": d["odds_ratio"], "positive_any": d["positive_any"],
                         "null_any": d["null_any"], "n_null": d["n_null"],
                         "n_matches": sm["n_matches"], "mode_offset": sm["mode_offset"],
                         "upstream_frac": sm["upstream_frac"],
                         "window_frac": sm["window_frac"],
                         "median_upstream": sm["median_upstream"],
                         "seqs_hit_in_window": seqs_with_hit_in_window(
                             d["per_seq_hits_positive"], lo, hi),
                         "offsets_positive": d["offsets_positive"]}
            rows.append([sig, name, "%d..%d" % (lo, hi), d["n_positive"], d["positive_any"],
                         C.fmt(100 * d["rate_positive"], 2), d["n_null"], d["null_any"],
                         C.fmt(100 * d["rate_null"], 2), C.fmt(100 * d["excess"], 2),
                         C.fmt(d["fold"], 2), C.fmt(d["odds_ratio"], 3),
                         C.fmt_p(d["fisher_p"]), sm["mode_offset"], sm["mode_count"],
                         C.fmt(100 * sm["upstream_frac"], 1),
                         C.fmt(100 * sm["window_frac"], 1), rec[name]["seqs_hit_in_window"]])
        sigma_res[sig] = rec
    C.write_tsv(os.path.join(C.TABLES, "T11_已知元件命中率_按sigma组发现集.tsv"),
                ["sigma", "element", "expected_window_offsets", "n_positive",
                 "positive_with_hit", "positive_rate_pct", "n_null", "null_with_hit",
                 "null_rate_pct", "excess_pp", "fold", "odds_ratio", "fisher_p",
                 "mode_offset", "mode_count", "upstream_pct", "in_expected_window_pct",
                 "n_seqs_with_hit_in_expected_window"], rows)

    rows = []
    for sig, rec in sorted(sigma_res.items()):
        for name, motif, lo, hi in ELEMENTS:
            cnt = collections.Counter(rec[name]["offsets_positive"])
            for o in sorted(cnt):
                rows.append([sig, name, o, cnt[o], "IN_WINDOW" if lo <= o <= hi else ""])
    C.write_tsv(os.path.join(C.TABLES, "T12_命中位置分布_按sigma组.tsv"),
                ["sigma", "element", "tss_relative_offset", "n_matches", "in_expected_window"],
                rows)

    # 每 σ 组在期望窗口内的信息量
    rows = []
    for sig in sorted(sigma_info):
        for r in sigma_info[sig]:
            rows.append([sig, r["offset"], C.fmt(r["info_vs_uniform"], 4),
                         C.fmt(r["info_vs_null_bits"], 4), C.fmt(r["shannon_entropy_bits"], 4)])
    C.write_tsv(os.path.join(C.TABLES, "T13_逐位置信息量_按sigma组.tsv"),
                ["sigma", "tss_relative_offset", "info_vs_uniform_bits",
                 "info_vs_null_bits", "shannon_entropy_bits"], rows)

    # ============ 合并口径：两个元件≥1 个命中（回答"到底多少启动子带教科书元件"） ============
    rows, combined = [], {}
    for tag, k, tag_cn in [("exact", 0, "精确匹配(0错配)"), ("mis1", 1, "≤1错配(敏感性)")]:
        pats = [mismatch_variants(mt, k) for _n, mt, _lo, _hi in ELEMENTS]
        p_hit = [bool(match_offsets(s, pats[0])) or bool(match_offsets(s, pats[1]))
                 for s in pos_all]
        n_hit = [bool(match_offsets(s, pats[0])) or bool(match_offsets(s, pats[1]))
                 for s in nul_all]
        pa, na = sum(p_hit), sum(n_hit)
        rp = pa / len(pos_all)
        rn = na / len(nul_all)
        odds, pval = stats.fisher_exact([[pa, len(pos_all) - pa], [na, len(nul_all) - na]])
        rd = rp - rn
        se = np.sqrt(rp * (1 - rp) / len(pos_all) + rn * (1 - rn) / len(nul_all))
        combined[tag] = {"positive_any": pa, "n_positive": len(pos_all),
                         "null_any": na, "n_null": len(nul_all),
                         "rate_positive": rp, "rate_null": rn, "excess": rd,
                         "excess_lo": rd - 1.96 * se, "excess_hi": rd + 1.96 * se,
                         "fold": (rp / rn) if rn else float("inf"),
                         "fisher_p": pval, "odds_ratio": odds}
        rows.append(["≥1 个元件命中（-35 或 -10）", tag_cn, len(pos_all), pa,
                     C.fmt(100 * rp, 2), len(nul_all), na, C.fmt(100 * rn, 2),
                     C.fmt(100 * rd, 2),
                     "[%s, %s]" % (C.fmt(100 * (rd - 1.96 * se), 2),
                                   C.fmt(100 * (rd + 1.96 * se), 2)),
                     C.fmt((rp / rn) if rn else float("inf"), 2), C.fmt_p(pval)])
    C.write_tsv(os.path.join(C.TABLES, "T8b_合并口径_至少一个元件_RegulonDBcore发现集.tsv"),
                ["scope", "match_rule", "n_positive", "positive_with_hit",
                 "positive_rate_pct", "n_null", "null_with_hit", "null_rate_pct",
                 "excess_pp", "excess_95CI_pp", "fold", "fisher_p"], rows)

    # ============ 共现与间距（两种匹配口径都算；正式间距分析属实验4/M4，此处仅描述） ============
    rows, cooc, spacers_by_rule = [], {}, {}
    for tag, k, tag_cn in [("exact", 0, "精确匹配(0错配)"), ("mis1", 1, "≤1错配(敏感性)")]:
        pats35 = mismatch_variants("-35 box" and "TTGACA", k)
        pats10 = mismatch_variants("TATAAT", k)
        both = only35 = only10 = neither = 0
        spacers = []
        for s in pos_all:
            h35 = match_offsets(s, pats35)
            h10 = match_offsets(s, pats10)
            if h35 and h10:
                both += 1
                r35 = min(h35, key=lambda o: (abs(o + 34.5), o))
                r10 = min(h10, key=lambda o: (abs(o + 9.5), o))
                sp = r10 - r35 - 6
                if SPACER_MIN <= sp <= SPACER_MAX:
                    spacers.append(sp)
            elif h35:
                only35 += 1
            elif h10:
                only10 += 1
            else:
                neither += 1
        spacers_by_rule[tag] = spacers
        cooc[tag] = {"n": len(pos_all), "both": both, "only_35": only35,
                     "only_10": only10, "neither": neither}
        rows.append(["共现计数", tag_cn, len(pos_all), both, only35, only10, neither,
                     C.fmt(100.0 * both / len(pos_all), 2),
                     C.fmt(100.0 * only35 / len(pos_all), 2),
                     C.fmt(100.0 * only10 / len(pos_all), 2),
                     C.fmt(100.0 * neither / len(pos_all), 2)])
        if spacers:
            a = np.asarray(spacers, float)
            rows.append(["间距（两六聚体之间的碱基数，教科书≈17）", tag_cn, int(a.size),
                         C.fmt(a.mean(), 2), C.fmt(np.median(a), 1),
                         C.fmt(np.percentile(a, 25), 1), C.fmt(np.percentile(a, 75), 1),
                         int(a.min()), int(a.max()),
                         int(collections.Counter(spacers).most_common(1)[0][0]), ""])
    C.write_tsv(os.path.join(C.TABLES, "T14_两元件共现与间距_RegulonDBcore发现集.tsv"),
                ["scope", "match_rule", "n_or_pairs", "both", "only_35", "only_10",
                 "neither_or_none", "both_pct_or_mean", "only_35_pct_or_q1",
                 "only_10_pct_or_q3", "neither_pct_or_min", "extra"], rows)
    sp_stat = {}
    if spacers_by_rule.get("mis1"):
        a = np.asarray(spacers_by_rule["mis1"], float)
        sp_stat = {"rule": "mis1", "n_pairs": int(a.size), "mean": float(a.mean()),
                   "sd": float(a.std(ddof=1)) if a.size > 1 else 0.0,
                   "median": float(np.median(a)), "q1": float(np.percentile(a, 25)),
                   "q3": float(np.percentile(a, 75)), "min": float(a.min()),
                   "max": float(a.max()),
                   "mode": int(collections.Counter(spacers_by_rule["mis1"]).most_common(1)[0][0])}
    if spacers_by_rule.get("exact"):
        a = np.asarray(spacers_by_rule["exact"], float)
        sp_stat["exact_n_pairs"] = int(a.size)
        sp_stat["exact_median"] = float(np.median(a))

    # ============ 备查附录：DBTBS core 与 TJU（不得用于结论） ============
    appendix = {}
    targets = [("dbtbs_bsub", "core")] + [(s, "main") for s in
                                          sorted({r["source"] for r in C.load_sequence_master()
                                                  if r["source"].startswith("tjupan_")})]
    for src, layer in targets:
        b = os.path.join(C.DATA, src, layer)
        p = [s for _h, s in C.read_fasta(os.path.join(b, "discovery.positive.fasta"))]
        n = [s for _h, s in C.read_fasta(os.path.join(b, "discovery.dinucleotide_null.fasta"))]
        appendix[src] = {"n": len(p), "elements": {}}
        for name, motif, lo, hi in ELEMENTS:
            d = rate_and_calibration(p, n, mismatch_variants(motif, 0))
            sm = offset_summary(d["offsets_positive"], lo, hi)
            appendix[src]["elements"][name] = {
                "rate_positive": d["rate_positive"], "rate_null": d["rate_null"],
                "excess": d["excess"], "fold": d["fold"], "fisher_p": d["fisher_p"],
                "positive_any": d["positive_any"], "null_any": d["null_any"],
                "n_null": d["n_null"], "mode_offset": sm["mode_offset"],
                "window_frac": sm["window_frac"],
                "seqs_hit_in_window": seqs_with_hit_in_window(
                    d["per_seq_hits_positive"], lo, hi)}
    rows = []
    for src in sorted(appendix):
        a = appendix[src]
        for name, motif, lo, hi in ELEMENTS:
            d = a["elements"][name]
            rows.append([src, C.SOURCE_CN.get(src, src), a["n"], name, "%d..%d" % (lo, hi),
                         C.fmt(100 * d["rate_positive"], 2),
                         C.fmt(100 * d["rate_null"], 2), C.fmt(100 * d["excess"], 2),
                         C.fmt(d["fold"], 2), C.fmt_p(d["fisher_p"]),
                         d["mode_offset"], C.fmt(100 * d["window_frac"], 1),
                         d["seqs_hit_in_window"]])
    C.write_tsv(os.path.join(C.TABLES, "T15_位置附录_DBTBS与TJU_不得用于结论.tsv"),
                ["source", "source_cn", "n_positive", "element", "expected_window_offsets",
                 "positive_rate_pct", "null_rate_pct", "excess_pp", "fold", "fisher_p",
                 "mode_offset", "in_expected_window_pct",
                 "n_seqs_with_hit_in_expected_window"], rows)

    # ============ 图 ============
    # F4 全体偏移分布（精确口径）+ 期望窗口
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    bins = np.arange(-61.5, 21.5, 1)
    for ax, (name, motif, lo, hi) in zip(axes, ELEMENTS):
        offs = main_res["%s|exact" % name]["offsets_positive"]
        ax.hist(offs, bins=bins, color="#5B9BD5", edgecolor="white", linewidth=0.4)
        ax.axvspan(lo - 0.5, hi + 0.5, color="crimson", alpha=0.16,
                   label="expected window %d..%d" % (lo, hi))
        ax.axvline(0, color="black", linewidth=1.0)
        ax.set_title("%s  (%d matches / %d seqs)" % (name, len(offs), n_tot))
        ax.set_xlabel("match start, TSS-relative offset (bp)")
        ax.set_ylabel("matches")
        ax.legend(fontsize=8)
    fig.suptitle("RegulonDB E. coli core, discovery: exact-match positions "
                 "(TTGACA / TATAAT)")
    fig.tight_layout()
    f4 = os.path.join(C.FIGURES, "F4_命中位置分布_RegulonDBcore.png")
    fig.savefig(f4); plt.close(fig)

    # F5 命中率 vs 背景（按 σ 组，精确口径）
    sigs = sorted(sigma_res, key=lambda s: -sigma_res[s]["n_positive"])
    x = np.arange(len(sigs))
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for ax, (name, motif, lo, hi) in zip(axes, ELEMENTS):
        rp = [100 * sigma_res[s][name]["rate_positive"] for s in sigs]
        rn = [100 * sigma_res[s][name]["rate_null"] for s in sigs]
        ax.bar(x - 0.2, rp, 0.4, label="promoters", color="#2E75B6")
        ax.bar(x + 0.2, rn, 0.4, label="dinucleotide-null background", color="#C9C9C9")
        ax.set_ylabel("% seqs with >=1 hit")
        ax.set_title("%s   exact match" % name)
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(["%s\nn=%d" % (s, sigma_res[s]["n_positive"]) for s in sigs],
                            rotation=30, ha="right", fontsize=8)
    fig.suptitle("Known-element exact-match hit rate vs dinucleotide-null background "
                 "(discovery)")
    fig.tight_layout()
    f5 = os.path.join(C.FIGURES, "F5_命中率vs背景_按sigma组.png")
    fig.savefig(f5); plt.close(fig)

    # F6 逐位置信息量（vs 背景），标注期望窗口
    offs = [r["offset"] for r in info_rows]
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.4), sharex=True)
    axes[0].plot(offs, [r["info_vs_uniform"] for r in info_rows], color="#2E75B6", linewidth=1.2)
    axes[0].set_ylabel("bits vs uniform")
    axes[0].set_title("RegulonDB E. coli core, discovery (n=%d): per-position information"
                      % n_tot)
    axes[1].plot(offs, [r["info_vs_null_bits"] for r in info_rows], color="#548235",
                 linewidth=1.2)
    axes[1].set_ylabel("bits vs null")
    axes[1].set_xlabel("TSS-relative offset (bp)")
    for ax in axes:
        for name, motif, lo, hi in ELEMENTS:
            ax.axvspan(lo - 0.5, hi + 0.5, color="crimson", alpha=0.13)
        ax.axvline(0, color="black", linewidth=0.9)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    f6 = os.path.join(C.FIGURES, "F6_逐位置信息量_RegulonDBcore.png")
    fig.savefig(f6); plt.close(fig)

    # F7 主要 σ 组的逐位置信息量（vs 背景）
    major = [s for s in ["Sigma70", "Sigma38", "Sigma24", "Sigma32", "Sigma54", "Sigma28"]
             if s in sigma_info]
    fig, axes = plt.subplots(1, len(major), figsize=(3.0 * len(major), 3.6), sharey=True)
    for ax, sig in zip(np.atleast_1d(axes), major):
        rr = sigma_info[sig]
        ax.plot([r["offset"] for r in rr], [r["info_vs_null_bits"] for r in rr],
                color="#2E75B6", linewidth=1.1)
        for name, motif, lo, hi in ELEMENTS:
            ax.axvspan(lo - 0.5, hi + 0.5, color="crimson", alpha=0.13)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title("%s\nn=%d" % (sig, sigma_res[sig]["n_positive"]), fontsize=9)
        ax.set_xlabel("offset (bp)", fontsize=8)
        ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel("bits vs null", fontsize=8)
    fig.suptitle("Per-position information vs paired dinucleotide-null, by sigma group "
                 "(RegulonDB core, discovery)")
    fig.tight_layout()
    f7 = os.path.join(C.FIGURES, "F7_逐位置信息量_按sigma组.png")
    fig.savefig(f7); plt.close(fig)

    # F8 附录
    app_srcs = sorted(appendix)
    x = np.arange(len(app_srcs))
    fig, ax = plt.subplots(figsize=(11, 5))
    w = 0.2
    for k, (name, motif, lo, hi) in enumerate(ELEMENTS):
        rp = [100 * appendix[s]["elements"][name]["rate_positive"] for s in app_srcs]
        rn = [100 * appendix[s]["elements"][name]["rate_null"] for s in app_srcs]
        ax.bar(x + (2 * k - 1.5) * w, rp, w, color=["#2E75B6", "#548235"][k],
               label="%s promoters" % name.split(" ")[0])
        ax.bar(x + (2 * k - 0.5) * w, rn, w, color=["#BDD7EE", "#C6E0B4"][k],
               label="%s null" % name.split(" ")[0])
    ax.set_xticks(x)
    ax.set_xticklabels(["%s\nn=%d" % (C.SOURCE_CN.get(s, s), appendix[s]["n"])
                        for s in app_srcs], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("% seqs with >=1 exact hit")
    ax.set_title("APPENDIX ONLY - not for conclusions (Q&A rule vs repo spec conflict)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    f8 = os.path.join(C.FIGURES, "F8_附录_位置_DBTBS与TJU_不得用于结论.png")
    fig.savefig(f8); plt.close(fig)

    # ============ 汇总写出 ============
    result["main"] = {k: {kk: vv for kk, vv in v.items()
                          if kk not in ("offsets_positive", "offsets_null",
                                        "per_seq_hits_positive")}
                      for k, v in main_res.items()}
    result["sigma_groups"] = {s: {k: {kk: vv for kk, vv in v.items()
                                      if kk != "offsets_positive"}
                                  for k, v in r.items() if isinstance(v, dict)}
                              for s, r in sigma_res.items()}
    result["sigma_groups_n"] = {s: r["n_positive"] for s, r in sigma_res.items()}
    result["cooccurrence"] = cooc
    result["combined_either_element"] = combined
    result["spacing"] = sp_stat
    result["consensus_position_check"] = consensus_check
    result["position_information_main"] = {
        "max_info_vs_uniform": max(r["info_vs_uniform"] for r in info_rows),
        "max_info_vs_uniform_offset": max(info_rows, key=lambda r: r["info_vs_uniform"])["offset"],
        "max_info_vs_null": max(r["info_vs_null_bits"] for r in info_rows),
        "max_info_vs_null_offset": max(info_rows, key=lambda r: r["info_vs_null_bits"])["offset"],
    }
    result["appendix"] = appendix
    with open(os.path.join(C.OUT, "positions.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)

    # ============ 控制台 ============
    print("=== 主结果：RegulonDB core 全体（discovery，n=%d）===" % n_tot)
    for name, motif, lo, hi in ELEMENTS:
        for tag, tag_cn in [("exact", "精确匹配"), ("mis1", "≤1错配")]:
            d = main_res["%s|%s" % (name, tag)]
            s = d["summary"]
            print("  %s [%s] 期望窗口偏移 %d..%d，%d 个模式串" % (name, tag_cn, lo, hi,
                                                                d["n_patterns"]))
            print("    正样本 %d/%d = %.2f%%  |  打乱背景 %d/%d = %.2f%%  |  超出 %+.2f pp"
                  " (95%%CI [%+.2f, %+.2f])  富集 %.2f 倍  Fisher p=%s"
                  % (d["positive_any"], d["n_positive"], 100 * d["rate_positive"],
                     d["null_any"], d["n_null"], 100 * d["rate_null"],
                     100 * d["excess"], 100 * d["excess_lo"], 100 * d["excess_hi"],
                     d["fold"], C.fmt_p(d["fisher_p"])))
            print("    命中数 %d，众数偏移 %s，上游 %.1f%%，落在期望窗口 %.1f%%，"
                  "每序列命中 %.3f"
                  % (s["n_matches"], s["mode_offset"], 100 * s["upstream_frac"],
                     100 * s["window_frac"], d["hits_per_seq_positive"]))
        print("    期望窗口内至少一个命中的序列数（精确口径）: %d / %d = %.2f%%"
              % (main_res["%s|exact" % name]["seqs_hit_in_window"], n_tot,
                 100.0 * main_res["%s|exact" % name]["seqs_hit_in_window"] / n_tot))

    print("\n=== 逐位置信息量（全体）===")
    print("  最大 vs uniform: %.3f bits @ offset %d"
          % (result["position_information_main"]["max_info_vs_uniform"],
             result["position_information_main"]["max_info_vs_uniform_offset"]))
    print("  最大 vs 背景   : %.3f bits @ offset %d"
          % (result["position_information_main"]["max_info_vs_null"],
             result["position_information_main"]["max_info_vs_null_offset"]))
    print("  期望窗口内 vs 背景信息量：")
    for name, motif, lo, hi in ELEMENTS:
        vals = [info_by_off[o]["info_vs_null_bits"] for o in range(lo, hi + 1)]
        print("    %-18s 偏移 %d..%d  平均 %.3f bits，最大 %.3f bits"
              % (name, lo, hi, float(np.mean(vals)), float(np.max(vals))))
    print("  期望窗口内逐位置最优碱基 vs 共识串：")
    for name, motif, lo, hi in ELEMENTS:
        c = consensus_check[name]
        print("    %-18s 共识 %s | 实测 %s | 一致 %d/%d 位"
              % (name, c["consensus"], c["observed_argmax"],
                 c["match_positions"], c["positions"]))

    print("\n=== 合并口径：至少一个元件命中 ===")
    for tag, tag_cn in [("exact", "精确匹配"), ("mis1", "≤1错配")]:
        c = combined[tag]
        print("  [%s] 正样本 %d/%d = %.2f%%  |  打乱背景 %d/%d = %.2f%%  |  超出 %+.2f pp "
              "(95%%CI [%+.2f, %+.2f])  富集 %.2f 倍  Fisher p=%s"
              % (tag_cn, c["positive_any"], c["n_positive"], 100 * c["rate_positive"],
                 c["null_any"], c["n_null"], 100 * c["rate_null"], 100 * c["excess"],
                 100 * c["excess_lo"], 100 * c["excess_hi"], c["fold"],
                 C.fmt_p(c["fisher_p"])))

    print("\n=== 共现与间距 ===")
    for tag, tag_cn in [("exact", "精确匹配"), ("mis1", "≤1错配")]:
        c = cooc[tag]
        print("  [%s] 两个元件都有 %d/%d = %.2f%%；只有 -35 %d；只有 -10 %d；都没有 %d"
              % (tag_cn, c["both"], c["n"], 100.0 * c["both"] / c["n"],
                 c["only_35"], c["only_10"], c["neither"]))
        sp = spacers_by_rule.get(tag) or []
        if len(sp) >= 2:
            a = np.asarray(sp, float)
            print("        间距（两六聚体之间的碱基数，教科书≈17）：n=%d 均值=%.2f 中位数=%.1f "
                  "IQR=[%.1f, %.1f] 范围=[%.0f, %.0f] 众数=%d"
                  % (a.size, a.mean(), np.median(a), np.percentile(a, 25),
                     np.percentile(a, 75), a.min(), a.max(),
                     collections.Counter(sp).most_common(1)[0][0]))
        elif len(sp) == 1:
            print("        间距：只有 %d 对，数值 %d bp，样本量不足以统计" % (len(sp), sp[0]))
        else:
            print("        间距：满足条件的配对数为 0，无法统计")

    print("\n=== σ 组命中率（精确口径，discovery，正样本 vs 打乱背景）===")
    print("  %-9s %5s | %-14s %-14s" % ("sigma", "n", "-35 正/背景", "-10 正/背景"))
    for s in sigs:
        r = sigma_res[s]
        n35, n10 = "-35 box (TTGACA)", "-10 box (TATAAT)"
        print("  %-9s %5d | %6.1f%% / %5.1f%%  %6.1f%% / %5.1f%%"
              % (s, r["n_positive"], 100 * r[n35]["rate_positive"],
                 100 * r[n35]["rate_null"], 100 * r[n10]["rate_positive"],
                 100 * r[n10]["rate_null"]))

    print("\n=== 附录（不得用于结论）===")
    for s in app_srcs:
        a = appendix[s]
        n35, n10 = "-35 box (TTGACA)", "-10 box (TATAAT)"
        print("  %-28s n=%5d  -35 %5.1f%%/%5.1f%%  -10 %5.1f%%/%5.1f%%"
              % (s, a["n"], 100 * a["elements"][n35]["rate_positive"],
                 100 * a["elements"][n35]["rate_null"],
                 100 * a["elements"][n10]["rate_positive"],
                 100 * a["elements"][n10]["rate_null"]))
    print("\n图：%s" % " / ".join(os.path.basename(f) for f in (f4, f5, f6, f7, f8)))


if __name__ == "__main__":
    main()
