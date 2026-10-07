# -*- coding: utf-8 -*-
"""
乙 - 步骤 3：质量核查 + GC / 长度 / 缺失

干什么：查序列长度、GC 含量、缺失和重复，按来源/物种/σ 组分组，画图。
为什么：GC 是 M2 唯一能跨所有来源比较的量；重复和缺失不写清楚，
        后面所有"组间差异"都可能是重复样本或缺失值造成的假象。
怎么做：
  A. 质检：长度分布、非法字符、关键列缺失、完全重复、反向互补重复、跨划分重复；
  B. GC：按来源×物种、按 σ 组（core 单标签）、按 TJU 正/天然对照分别统计；
  C. GC 做 dinucleotide_null 对照（打乱序列的 GC 应接近原序列，这是打乱是否合格的自检）。
做到什么程度：每个分组都有 n / 均值 / 标准差 / 中位数 / 四分位 / 极值；
              长度无分布可画时必须说明原因，不画假图。
产出：乙_产出/tables/T4~T7、乙_产出/figures/F1~F3、out/qc_gc.json
"""
import collections
import json
import os
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


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def stat_block(vals):
    a = np.asarray(vals, dtype=float)
    if a.size == 0:
        return dict(n=0, mean=float("nan"), sd=float("nan"), median=float("nan"),
                    q1=float("nan"), q3=float("nan"), min=float("nan"), max=float("nan"))
    return dict(n=int(a.size), mean=float(a.mean()),
                sd=float(a.std(ddof=1)) if a.size > 1 else 0.0,
                median=float(np.median(a)), q1=float(np.percentile(a, 25)),
                q3=float(np.percentile(a, 75)),
                min=float(a.min()), max=float(a.max()))


def rank_biserial(x, y):
    """Mann-Whitney U 的效应量（秩双列相关），范围 [-1,1]。"""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.size == 0 or y.size == 0:
        return float("nan")
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return float(2.0 * u / (x.size * y.size) - 1.0)


def main():
    sm = C.load_sequence_master()
    sig_idx = C.build_sigma_index(C.load_sigma_associations())
    rec_by_id = {r["sequence_id"]: r for r in sm}

    # ================= A. 质检 =================
    qc = {}
    qc["n_records"] = len(sm)

    # A1 长度
    lens = collections.Counter(len(r["sequence"]) for r in sm)
    qc["length_distribution"] = dict(lens)
    # A2 非法字符
    bad_chars = collections.Counter()
    for r in sm:
        for ch in set(r["sequence"].upper()):
            if ch not in "ACGT":
                bad_chars[ch] += 1
    qc["illegal_char_records"] = dict(bad_chars)
    # A3 关键列缺失
    key_cols = ["sequence_id", "species", "source", "sequence", "length",
                "local_tss_index", "split", "tier", "label", "orientation",
                "coordinate_evidence"]
    missing = {c: sum(1 for r in sm if not (r.get(c) or "").strip()) for c in key_cols}
    qc["missing_key_columns"] = missing
    # A4 local_tss_index / window_offsets 一致性
    qc["local_tss_index_values"] = dict(collections.Counter(r["local_tss_index"] for r in sm))
    qc["window_offsets_values"] = dict(collections.Counter(r["window_offsets"] for r in sm))
    # A5 完全重复（同一划分内 / 跨划分）
    def dup_report(seq_field):
        within_all, cross = 0, 0
        bucket = collections.defaultdict(set)
        for r in sm:
            bucket[r[seq_field]].add((r["source"], r["split"]))
        for seq, tags in bucket.items():
            if len(tags) > 1:
                cross += 1
        # 同一 (source, split) 内的重复条数
        b2 = collections.Counter((r["source"], r["split"], r[seq_field]) for r in sm)
        within_all = sum(v - 1 for v in b2.values() if v > 1)
        return within_all, cross

    exact_within, exact_cross = dup_report("sequence")
    # 反向互补重复：对每条序列查它的反向互补是否也存在于数据中
    seq_set = {r["sequence"] for r in sm}
    rc_pairs = 0
    rc_palindromes = 0
    for s in seq_set:
        r = revcomp(s)
        if r == s:
            rc_palindromes += 1
        elif r in seq_set:
            rc_pairs += 1
    rc_pairs //= 2  # 每对会被数两次
    qc["exact_duplicate_extra_records_within_source_split"] = exact_within
    qc["exact_duplicate_sequences_across_source_or_split"] = exact_cross
    qc["rc_duplicate_pairs"] = rc_pairs
    qc["rc_palindromic_sequences"] = rc_palindromes

    # A5a' RegulonDB core 与 extended 是否按「序列内容」重叠
    #      （sequence_id 带层级前缀因此不相等，必须用序列本体判断）
    core_seqs = {r["sequence"] for r in sm if r["source"] == "regulondb_ecoli"
                 and r["tier"] == "core"}
    ext_seqs = {r["sequence"] for r in sm if r["source"] == "regulondb_ecoli"
                and r["tier"] == "extended"}
    qc["regulondb_core_vs_extended"] = {
        "core_unique_sequences": len(core_seqs),
        "extended_unique_sequences": len(ext_seqs),
        "shared_sequences": len(core_seqs & ext_seqs),
        "union": len(core_seqs | ext_seqs),
    }

    # A5b 把"跨来源或跨划分重复"拆开看：是同一来源不同划分（真泄漏嫌疑），
    #     还是同一物种不同来源（这是设计允许的同物种跨数据集重叠）
    seq_tags = collections.defaultdict(set)
    for r in sm:
        seq_tags[r["sequence"]].add((r["source"], r["split"]))
    cross_split_same_source = collections.Counter()
    cross_source_same_split = collections.Counter()
    both = 0
    for seq, tags in seq_tags.items():
        if len(tags) < 2:
            continue
        srcs = {t[0] for t in tags}
        splits = {t[1] for t in tags}
        if len(srcs) == 1:
            cross_split_same_source[(list(srcs)[0], tuple(sorted(splits)))] += 1
        elif len(splits) == 1:
            cross_source_same_split[tuple(sorted(srcs))] += 1
        else:
            both += 1
    qc["dup_same_source_cross_split"] = {"%s | %s" % (k[0], "+".join(k[1])): v
                                         for k, v in cross_split_same_source.items()}
    qc["dup_cross_source_same_split"] = {"+".join(k): v
                                         for k, v in cross_source_same_split.items()}
    qc["dup_both_source_and_split"] = both

    # A5c 与甲的泄漏表交叉核对：跨划分的同源对应该为 0
    leak_path = os.path.join(C.DATA, "leakage_links.tsv")
    split_of = {r["sequence_id"]: r["split"] for r in sm}
    leak_cross = 0
    leak_n = 0
    if os.path.exists(leak_path):
        for row in C.read_tsv(leak_path):
            leak_n += 1
            sa, sb = split_of.get(row["a"]), split_of.get(row["b"])
            if sa and sb and sa != sb:
                leak_cross += 1
    qc["leakage_links_rows"] = leak_n
    qc["leakage_links_cross_split"] = leak_cross

    # A5d 找出那 1 条同(来源,划分)内的多余重复，指名道姓
    b2 = collections.Counter((r["source"], r["split"], r["sequence"]) for r in sm)
    dup_detail = []
    for (src, split, seq), n in b2.items():
        if n > 1:
            ids = [r["sequence_id"] for r in sm
                   if r["source"] == src and r["split"] == split and r["sequence"] == seq]
            dup_detail.append({"source": src, "split": split, "copies": n,
                               "sequence_ids": ids})
    qc["duplicate_detail"] = dup_detail

    # ================= B. GC 统计 =================
    gc_all = {}
    for r in sm:
        gc_all[r["sequence_id"]] = C.gc_percent(r["sequence"])

    # B1 来源 × 物种 × 标签
    # 注意：TJU 同时有正样本和天然对照，两者 GC 有系统差异，
    #       混在一起算会得到一个没有意义的"物种 GC"，必须按 label 拆开。
    LABEL_CN = {"1": "正样本", "0": "天然对照"}
    groups = collections.defaultdict(list)
    for r in sm:
        groups[(r["source"], r["species"], r["label"])].append(gc_all[r["sequence_id"]])
    rows = []
    for (src, sp, lab), vals in sorted(groups.items()):
        s = stat_block(vals)
        rows.append([src, C.SOURCE_CN.get(src, src), sp, lab, LABEL_CN.get(lab, lab),
                     s["n"], C.fmt(s["mean"]), C.fmt(s["sd"]), C.fmt(s["median"]),
                     C.fmt(s["q1"]), C.fmt(s["q3"]), C.fmt(s["min"]), C.fmt(s["max"])])
    C.write_tsv(os.path.join(C.TABLES, "T4_GC_按来源x物种x标签.tsv"),
                ["source", "source_cn", "species", "label", "label_cn", "n", "mean", "sd",
                 "median", "q1", "q3", "min", "max"], rows)

    # B2 σ 组（core 单标签），只用 discovery 做统计
    sigma_gc = collections.defaultdict(list)
    for sid, v in sig_idx.items():
        rec = rec_by_id.get(sid)
        if rec is None or rec["source"].startswith("tjupan_"):
            continue
        if rec["tier"] != "core":
            continue
        lab = C.single_sigma_label(v["sigma_raw"])
        if lab is None:
            continue
        sigma_gc[(rec["source"], lab, rec["split"])].append(gc_all[sid])

    summary = {(r["source"], r["sigma"]): r for r in C.load_sigma_group_summary()}
    rows = []
    for (src, sig) in sorted({(s, g) for (s, g, _sp) in sigma_gc}):
        vals = sigma_gc[(src, sig, "discovery")]
        s = stat_block(vals)
        e = summary.get((src, sig), {})
        rows.append([src, C.SOURCE_CN.get(src, src), sig, s["n"], C.fmt(s["mean"]),
                     C.fmt(s["sd"]), C.fmt(s["median"]), C.fmt(s["q1"]), C.fmt(s["q3"]),
                     e.get("eligibility", ""),
                     C.ELIGIBILITY_CN.get(e.get("eligibility", ""), "")])
    C.write_tsv(os.path.join(C.TABLES, "T5_GC_按sigma组_发现集.tsv"),
                ["source", "source_cn", "sigma", "n_discovery", "mean", "sd", "median",
                 "q1", "q3", "eligibility", "eligibility_cn"], rows)

    # B3 TJU 正样本 vs 天然对照
    tju_cmp = collections.defaultdict(lambda: {"pos": [], "ctl": []})
    for r in sm:
        if not r["source"].startswith("tjupan_"):
            continue
        key = r["species"]
        tju_cmp[key]["pos" if r["label"] == "1" else "ctl"].append(gc_all[r["sequence_id"]])
    rows, tju_stats = [], {}
    for sp in sorted(tju_cmp):
        pos, ctl = tju_cmp[sp]["pos"], tju_cmp[sp]["ctl"]
        sps, scs = stat_block(pos), stat_block(ctl)
        u, p = stats.mannwhitneyu(pos, ctl, alternative="two-sided")
        rb = rank_biserial(pos, ctl)
        diff = sps["median"] - scs["median"]
        tju_stats[sp] = {"n_pos": sps["n"], "n_ctl": scs["n"],
                         "median_pos": sps["median"], "median_ctl": scs["median"],
                         "median_diff": diff, "p": float(p), "rank_biserial": rb}
        rows.append([sp, sps["n"], C.fmt(sps["mean"]), C.fmt(sps["median"]),
                     scs["n"], C.fmt(scs["mean"]), C.fmt(scs["median"]),
                     C.fmt(diff), C.fmt(p, 6), C.fmt(rb)])
    C.write_tsv(os.path.join(C.TABLES, "T6_GC_TJU正样本vs天然对照.tsv"),
                ["species", "n_positive", "mean_positive", "median_positive",
                 "n_natural_control", "mean_natural_control", "median_natural_control",
                 "median_diff_pp", "mannwhitney_p", "rank_biserial"], rows)

    # B4 打乱对照自检：GC 应与原序列接近
    null_check = []
    for src in sorted({r["source"] for r in sm}):
        for split in C.SPLITS:
            pos_p = os.path.join(C.DATA, src, "main" if src.startswith("tjupan_") else "core",
                                 "%s.positive.fasta" % split)
            null_p = os.path.join(C.DATA, src, "main" if src.startswith("tjupan_") else "core",
                                  "%s.dinucleotide_null.fasta" % split)
            if not (os.path.exists(pos_p) and os.path.exists(null_p)):
                continue
            gp = [C.gc_percent(s) for _h, s in C.read_fasta(pos_p)]
            gn = [C.gc_percent(s) for _h, s in C.read_fasta(null_p)]
            null_check.append([src, split, len(gp), C.fmt(np.mean(gp)),
                               len(gn), C.fmt(np.mean(gn)),
                               C.fmt(np.mean(gn) - np.mean(gp))])
    C.write_tsv(os.path.join(C.TABLES, "T7_打乱对照GC自检.tsv"),
                ["source", "split", "n_positive", "mean_gc_positive",
                 "n_null", "mean_gc_null", "mean_diff_pp"], null_check)

    # ================= C. 图 =================
    # F1 来源×物种×标签 GC 箱线图（TJU 的正样本与天然对照必须分开画）
    order = sorted(groups, key=lambda k: (k[2], np.median(groups[k])))
    labels = ["%s %s\n%s\nn=%d" % (C.SOURCE_CN.get(k[0], k[0]).split(" ")[0], k[1],
                                   LABEL_CN.get(k[2], k[2]), len(groups[k])) for k in order]
    fig, ax = plt.subplots(figsize=(13.5, 6))
    bp = ax.boxplot([groups[k] for k in order], labels=labels, showfliers=False,
                    patch_artist=True, widths=0.6)
    for patch, k in zip(bp["boxes"], order):
        patch.set_facecolor("#F5B7B1" if k[2] == "0" else "#7FB3D5")
    ax.set_ylabel("GC content (%) of the 81 bp window")
    ax.set_title("GC content by source x species x label (n = %d sequences)" % len(sm))
    ax.tick_params(axis="x", labelsize=6.5)
    ax.grid(axis="y", alpha=0.3)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="#7FB3D5", label="positive (label=1)"),
                       Patch(facecolor="#F5B7B1", label="natural control (label=0)")],
              loc="upper left", fontsize=8)
    fig.tight_layout()
    f1 = os.path.join(C.FIGURES, "F1_GC_按来源x物种x标签.png")
    fig.savefig(f1); plt.close(fig)

    # F2 σ 组 GC（RegulonDB core 发现集）
    rd_groups = sorted({(s, g) for (s, g, _sp) in sigma_gc if s == "regulondb_ecoli"},
                       key=lambda k: -len(sigma_gc[(k[0], k[1], "discovery")]))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    data = [sigma_gc[(s, g, "discovery")] for (s, g) in rd_groups]
    ax.boxplot(data, labels=["%s\nn=%d" % (g, len(sigma_gc[(s, g, "discovery")]))
                             for (s, g) in rd_groups], showfliers=False)
    ax.set_ylabel("GC content (%)")
    ax.set_title("RegulonDB E. coli core, single-sigma, discovery split")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    f2 = os.path.join(C.FIGURES, "F2_GC_按RegulonDB_sigma组.png")
    fig.savefig(f2); plt.close(fig)

    # F3 TJU 正样本 vs 天然对照
    sps_order = sorted(tju_cmp)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    positions = np.arange(len(sps_order))
    bp = ax.boxplot([tju_cmp[s]["pos"] for s in sps_order], positions=positions - 0.18,
                    widths=0.32, showfliers=False, patch_artist=True,
                    boxprops=dict(facecolor="#7FB3D5"))
    bc = ax.boxplot([tju_cmp[s]["ctl"] for s in sps_order], positions=positions + 0.18,
                    widths=0.32, showfliers=False, patch_artist=True,
                    boxprops=dict(facecolor="#F5B7B1"))
    ax.set_xticks(positions); ax.set_xticklabels(sps_order, rotation=20, ha="right")
    ax.set_ylabel("GC content (%)")
    ax.set_title("TJU six species: positive promoters vs natural controls (label=0)")
    ax.legend([bp["boxes"][0], bc["boxes"][0]], ["positive (label=1)", "natural control (label=0)"],
              loc="upper right")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    f3 = os.path.join(C.FIGURES, "F3_GC_TJU正样本vs天然对照.png")
    fig.savefig(f3); plt.close(fig)

    result = {"qc": qc, "tju_gc_comparison": tju_stats,
              "figures": [os.path.basename(f1), os.path.basename(f2), os.path.basename(f3)],
              "n_sigma_groups_core_single": len(rows)}
    with open(os.path.join(C.OUT, "qc_gc.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)

    # ================= 控制台 =================
    print("=== A 质检 ===")
    print("记录数              : %d" % qc["n_records"])
    print("长度分布            : %s" % qc["length_distribution"])
    print("非法字符记录        : %s" % (qc["illegal_char_records"] or "无"))
    print("关键列缺失          : %s" % {k: v for k, v in missing.items() if v})
    print("local_tss_index     : %s" % qc["local_tss_index_values"])
    print("window_offsets      : %s" % qc["window_offsets_values"])
    print("同(来源,划分)内完全重复的多余记录 : %d" % exact_within)
    for d in dup_detail:
        print("      -> %s %s 出现 %d 次：%s" % (d["source"], d["split"], d["copies"],
                                                ", ".join(d["sequence_ids"])))
    print("跨来源或跨划分的重复序列          : %d" % exact_cross)
    print("     其中 同来源跨划分（泄漏嫌疑）: %s" % (qc["dup_same_source_cross_split"] or "无"))
    print("     其中 跨来源同划分（设计允许）: %s" % (qc["dup_cross_source_same_split"] or "无"))
    print("     两者都跨                    : %d" % both)
    print("甲的泄漏表 %d 对，其中跨划分 %d 对（清单登记 cross_split_leakage_groups=0）"
          % (leak_n, leak_cross))
    ce = qc["regulondb_core_vs_extended"]
    print("RegulonDB core/extended 按序列内容重叠 : core %d + extended %d，共享 %d，并集 %d"
          % (ce["core_unique_sequences"], ce["extended_unique_sequences"],
             ce["shared_sequences"], ce["union"]))
    print("反向互补重复对                    : %d（回文序列 %d 条）"
          % (rc_pairs, rc_palindromes))

    print("\n=== B1 GC 按来源 x 物种 x 标签 ===")
    for r in sorted(_read(os.path.join(C.TABLES, "T4_GC_按来源x物种x标签.tsv")),
                    key=lambda x: (x[3], float(x[8]))):
        print("  %-28s %-20s %-6s n=%5d 中位数=%6s%% 均值=%6s%%"
              % (r[1], r[2], r[4], int(r[5]), r[8], r[6]))

    print("\n=== B2 GC 按 σ 组（RegulonDB core 发现集）===")
    for r in sorted(_read(os.path.join(C.TABLES, "T5_GC_按sigma组_发现集.tsv")),
                    key=lambda x: -int(x[3])):
        print("  %-9s n=%4d 中位数=%s%% 均值=%s%%  [%s]"
              % (r[2], int(r[3]), r[6], r[4], r[10]))

    print("\n=== B3 TJU 正样本 vs 天然对照 GC ===")
    for r in _read(os.path.join(C.TABLES, "T6_GC_TJU正样本vs天然对照.tsv")):
        print("  %-18s 正 n=%5s 中位数=%s%% | 对照 n=%5s 中位数=%s%% | 差=%s pp  p=%.2e  效应量=%s"
              % (r[0], r[1], r[3], r[4], r[6], r[7], float(r[8]), r[9]))

    print("\n=== B4 打乱对照 GC 自检（前 8 行）===")
    for r in list(_read(os.path.join(C.TABLES, "T7_打乱对照GC自检.tsv")))[:8]:
        print("  %-28s %-12s 正=%s 打乱=%s 差=%s pp" % (r[0], r[1], r[3], r[5], r[6]))
    print("\n图已写出：%s / %s / %s" % (os.path.basename(f1), os.path.basename(f2),
                                        os.path.basename(f3)))


def _read(path):
    with open(path, "r", encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield line.rstrip("\n").split("\t")


if __name__ == "__main__":
    main()
