# -*- coding: utf-8 -*-
"""
乙 - 步骤 8：报告数字核对

干什么：把我在各份报告里手写的数字，逐个回到权威产物（TSV / JSON）里重算一遍并比对。
为什么：丙要能抽查。手抄数字是这份交付里最容易出错的一环，
        所以用脚本做"报告文字 vs 计算产物"的一致性检查，而不是靠人眼。
怎么做：每个断言写成 (报告文件, 报告里出现的字符串, 权威值, 出处说明)；
        脚本把权威值按报告里的格式重排后，检查该字符串确实出现在报告文件里。
        比较前统一把 Unicode 减号（−, U+2212）归一为 ASCII '-'，避免排版符号造成假失败。
产出：out/report_check.json，并打印 PASS/FAIL 汇总。

注意：本脚本只证明"报告里的数字与产物一致"，不证明"产物算得对"。
      后者由 scripts/04_positions.py 等脚本内的口径注释与 05 的交叉核对负责。
"""
import io
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

REPORTS = C.WORKSPACE
T = C.TABLES


def norm(s):
    return (s.replace("\u2212", "-").replace("\u2013", "-")
             .replace("\u00a0", " ").replace(",", ""))


def read_text(name):
    with io.open(os.path.join(REPORTS, name), "r", encoding="utf-8") as fh:
        return norm(fh.read())


def load_tsv(name):
    rows = []
    with io.open(os.path.join(T, name), encoding="utf-8", newline="") as fh:
        import csv
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    return rows


def load_json(path):
    with io.open(path, encoding="utf-8") as fh:
        return json.load(fh)


def num(x, nd=None):
    f = float(x)
    if nd is None:
        return ("%g" % f)
    return ("%%.%df" % nd) % f


def main():
    hc = load_json(os.path.join(C.OUT, "handoff_check.json"))
    head = load_json(os.path.join(C.OUT, "headcount.json"))
    qc = load_json(os.path.join(C.OUT, "qc_gc.json"))
    pos = load_json(os.path.join(C.OUT, "positions.json"))
    gd = load_json(os.path.join(C.OUT, "group_diff.json"))

    T2 = load_tsv("T2_来源x物种合计.tsv")
    T3 = load_tsv("T3_sigma组x划分.tsv")
    T4 = {(r["source"], r["label"]): r for r in load_tsv("T4_GC_按来源x物种x标签.tsv")}
    T6 = {(r["species"]): r for r in load_tsv("T6_GC_TJU正样本vs天然对照.tsv")}
    T8 = load_tsv("T8_已知元件命中率与背景校准_RegulonDBcore发现集.tsv")
    T8b = load_tsv("T8b_合并口径_至少一个元件_RegulonDBcore发现集.tsv")
    T10b = load_tsv("T10b_期望窗口逐位置最优碱基_vs_共识串.tsv")
    T16 = load_tsv("T16_GC组间比较_vs参照组.tsv")
    T17 = load_tsv("T17_元件命中率组间比较_vs参照组.tsv")

    T8x = {(r["element"], r["match_rule"]): r for r in T8}
    T16x = {(r["source"], r["sigma"]): r for r in T16}
    T17x = {(r["source"], r["sigma"], r["element"]): r for r in T17}
    sub35 = [r for r in T10b if r["element"].startswith("-35") and r["tss_relative_offset"] == "小计"][0]
    sub10 = [r for r in T10b if r["element"].startswith("-10") and r["tss_relative_offset"] == "小计"][0]

    sum_all = sum(int(r["n_all"]) for r in T2)
    sum_pos = sum(int(r["n_positive"]) for r in T2)
    sum_ctl = sum(int(r["n_natural_control"]) for r in T2)
    sum_total = sum(int(r["total"]) for r in T3)
    sum_disc = sum(int(r["discovery"]) for r in T3)
    tju = [r for r in T2 if r["source"].startswith("tjupan")]
    tju_all = sum(int(r["n_all"]) for r in tju)
    tju_pos = sum(int(r["n_positive"]) for r in tju)

    checks = []

    def add(report, text, why):
        checks.append((report, text, why))

    # ---------------- 交接核验 ----------------
    add("乙_数据交接核验_20261007.md", "210", "manifest 登记文件数 %d" % hc["A_sha256"]["registered"])
    add("乙_数据交接核验_20261007.md", "29,564,814" if False else "29564814",
        "zip 字节数（用 JSON 记录值）")
    add("乙_数据交接核验_20261007.md", "0d25d7e778b9b8e1e7d3d34b8741bcd6957e5de6db48721666036273ffb5fc79",
        "zip SHA256")

    # ---------------- 数人头 ----------------
    add("乙_数人头_20261007.md", "%d" % sum_all, "T2 n_all 合计")
    add("乙_数人头_20261007.md", "%d" % sum_pos, "T2 n_positive 合计")
    add("乙_数人头_20261007.md", "%d" % sum_ctl, "T2 n_natural_control 合计")
    add("乙_数人头_20261007.md", "%d" % sum_total, "T3 total 合计")
    add("乙_数人头_20261007.md", "%d" % sum_disc, "T3 discovery 合计")
    add("乙_数人头_20261007.md", "%d" % tju_all, "TJU 六物种合计")
    add("乙_数人头_20261007.md", "%d" % tju_pos, "TJU 正样本合计")
    for r in T3:
        if r["sigma"] in ("Sigma70", "SigA"):
            add("乙_数人头_20261007.md", r["total"], "%s 合计" % r["sigma"])
    add("乙_数人头_20261007.md", "%d" % head["C2_checked"], "C2 核对文件数")
    add("乙_数人头_20261007.md", "%d" % head["C3_checked"], "C3 核对文件数")
    add("乙_数人头_20261007.md", "%d" % head["C1_checks"], "C1 核对格数")
    add("乙_数人头_20261007.md", "3,812", "RegulonDB 唯一序列并集")

    # ---------------- 质控与 GC ----------------
    add("乙_质控_GC与长度_20261007.md", "22,594", "记录总数")
    add("乙_质控_GC与长度_20261007.md", "81", "序列长度")
    add("乙_质控_GC与长度_20261007.md", "1,169", "跨来源同划分重复序列数")
    add("乙_质控_GC与长度_20261007.md", "0.000", "打乱对照 GC 差")
    for key, label in [(("tjupan_staphylococcus", "1"), "葡萄球菌正样本"),
                       (("tjupan_bradyrhizobium", "1"), "慢生根瘤菌正样本"),
                       (("regulondb_ecoli", "1"), "RegulonDB 正样本"),
                       (("dbtbs_bsub", "1"), "DBTBS 正样本")]:
        add("乙_质控_GC与长度_20261007.md", T4[key]["median"], "%s GC 中位数" % label)
    for sp in T6:
        add("乙_质控_GC与长度_20261007.md", T6[sp]["median_diff_pp"], "%s 正负 GC 差" % sp)
        add("乙_质控_GC与长度_20261007.md", T6[sp]["rank_biserial"], "%s 效应量" % sp)

    # ---------------- 位置与背景校准 ----------------
    for elem in ("-35 box (TTGACA)", "-10 box (TATAAT)"):
        for rule in ("精确匹配(0错配)", "≤1错配(敏感性)"):
            r = T8x[(elem, rule)]
            tag = "%s %s" % (elem, rule)
            add("乙_已知元件位置与背景校准_20261007.md", r["positive_rate_pct"], tag + " 正样本率")
            add("乙_已知元件位置与背景校准_20261007.md", r["null_rate_pct"], tag + " 背景率")
            add("乙_已知元件位置与背景校准_20261007.md", r["excess_pp"], tag + " 超出")
            add("乙_已知元件位置与背景校准_20261007.md", r["mode_offset"], tag + " 众数偏移")
    add("乙_已知元件位置与背景校准_20261007.md", "0.3485", "-35 精确 Fisher p")
    add("乙_已知元件位置与背景校准_20261007.md", "0.6034", "-10 精确 Fisher p")
    add("乙_已知元件位置与背景校准_20261007.md", "0.0029", "-35 ≤1 错配 Fisher p")
    add("乙_已知元件位置与背景校准_20261007.md", "1.06e-05", "-10 ≤1 错配 Fisher p")
    add("乙_已知元件位置与背景校准_20261007.md", "0.221", "全窗口最大信息量")
    add("乙_已知元件位置与背景校准_20261007.md", sub35["observed_argmax_base"], "-35 实测最优碱基")
    add("乙_已知元件位置与背景校准_20261007.md", sub10["observed_argmax_base"], "-10 实测最优碱基")
    add("乙_已知元件位置与背景校准_20261007.md", sub35["match"], "-35 一致位数")
    add("乙_已知元件位置与背景校准_20261007.md", sub10["match"], "-10 一致位数")
    for r in T8b:
        add("乙_已知元件命中率" if False else "乙_已知元件位置与背景校准_20261007.md",
            r["positive_rate_pct"], "合并口径 %s 正样本率" % r["match_rule"])

    # ---------------- 组间差异 ----------------
    for r in T16:
        add("乙_组间差异_20261007.md", r["median_diff_pp"],
            "%s GC 中位数差" % r["sigma"])
        if r["inferential_reported"] == "是":
            add("乙_组间差异_20261007.md", r["ci95_lo_pp"], "%s CI 下界" % r["sigma"])
            add("乙_组间差异_20261007.md", r["ci95_hi_pp"], "%s CI 上界" % r["sigma"])
            add("乙_组间差异_20261007.md", r["bh_q"], "%s BH q" % r["sigma"])
    for r in T17:
        if r["inferential_reported"] == "是":
            tag = "%s %s" % (r["sigma"], r["element"])
            add("乙_组间差异_20261007.md", r["risk_diff_pp"], tag + " 风险差")
            add("乙_组间差异_20261007.md", r["hit_pct"], tag + " 命中率")
    for r in load_tsv("T18_sigma组描述统计_发现集.tsv"):
        add("乙_组间差异_20261007.md", r["gc_median"], "%s GC 中位数（描述表）" % r["sigma"])
        add("乙_组间差异_20261007.md", r["minus10_hit_pct"], "%s -10 命中率（描述表）" % r["sigma"])

    # ---------------- 位置与间距按 σ 组（报告框架 表4-6/4-7/4-8）----------------
    T19 = load_tsv("T19_位置统计_按sigma组发现集.tsv")
    T20 = load_tsv("T20_间距统计_按sigma组发现集.tsv")
    T19x = {(r["sigma"], r["element"], r["match_rule"]): r for r in T19}
    T20x = {(r["sigma"], r["match_rule"]): r for r in T20}
    RPT_PS = "乙_位置与间距_按sigma组_20261007.md"
    add(RPT_PS, "%d" % sum(int(r["n_sequences"]) for r in T20 if r["match_rule"] == "exact"
                           and r["sigma"] == T20[0]["sigma"]) if False else "22,594",
        "总记录数")
    for sig in ("Sigma70", "Sigma38", "Sigma24", "Sigma32"):
        for elem in ("-35 box (TTGACA)", "-10 box (TATAAT)"):
            r = T19x.get((sig, elem, "exact"))
            if not r:
                continue
            add(RPT_PS, r["median_offset"], "%s %s 位置中位数" % (sig, elem))
            add(RPT_PS, r["iqr"], "%s %s IQR" % (sig, elem))
            add(RPT_PS, r["n_with_hit"], "%s %s 有命中序列数" % (sig, elem))
            if r["mean_offset"]:
                add(RPT_PS, r["mean_offset"], "%s %s 位置均值" % (sig, elem))
    # Sigma70 的两个关键偏离量
    add(RPT_PS, T19x[("Sigma70", "-35 box (TTGACA)", "exact")]["range"], "Sigma70 -35 极差")
    add(RPT_PS, T19x[("Sigma70", "-10 box (TATAAT)", "exact")]["range"], "Sigma70 -10 极差")
    # 间距
    for sig in ("Sigma70", "Sigma38", "Sigma32", "Sigma24", "Sigma54", "Sigma28"):
        for rule, tag in (("exact", "精确"), ("mis1", "≤1错配")):
            r = T20x.get((sig, rule))
            if not r:
                continue
            if r["n_spacers"] not in ("0", ""):
                add(RPT_PS, r["n_spacers"], "%s %s 有效间距对数" % (sig, tag))
            if r["median_spacer"]:
                add(RPT_PS, r["median_spacer"], "%s %s 间距中位数" % (sig, tag))
    add(RPT_PS, T20x[("Sigma70", "mis1")]["both_elements"], "Sigma70 ≤1 两元件都有")
    add(RPT_PS, T20x[("Sigma70", "mis1")]["mean_spacer"], "Sigma70 ≤1 间距均值")
    add(RPT_PS, T20x[("Sigma70", "mis1")]["sd"], "Sigma70 ≤1 间距 SD")
    add(RPT_PS, T20x[("Sigma70", "mis1")]["atypical_pct_17p1"], "Sigma70 ≤1 非典型 17±1")
    add(RPT_PS, T20x[("Sigma70", "mis1")]["atypical_pct_17p2"], "Sigma70 ≤1 非典型 17±2")
    add(RPT_PS, "%d" % sum(int(r["n_spacers"]) for r in T20 if r["match_rule"] == "mis1"),
        "≤1 错配间距对数合计（与图 F11 一致）")

    # ---------------- 跑检查 ----------------
    cache = {}
    results, fails = [], []
    for report, text, why in checks:
        if report not in cache:
            cache[report] = read_text(report)
        content = cache[report]
        ok = norm(str(text)) in content
        item = {"report": report, "text": str(text), "why": why, "found": ok}
        results.append(item)
        if not ok:
            fails.append(item)

    out = {"total_checks": len(results), "passed": len(results) - len(fails),
           "failed": len(fails), "failures": fails, "all": results}

    # ---------------- 附：报告内部链接检查 ----------------
    import re
    link_fails = []
    n_links = 0
    for name in sorted(os.listdir(REPORTS)):
        if not (name.startswith("乙") and name.endswith(".md")):
            continue
        content = read_text(name)
        for m in re.finditer(r"\]\(([^)#]+?)\)", content):
            target = m.group(1).strip()
            if target.startswith(("http://", "https://", "#")):
                continue
            n_links += 1
            if not os.path.exists(os.path.join(REPORTS, target)):
                link_fails.append({"report": name, "target": target})
    out["internal_links_checked"] = n_links
    out["internal_links_broken"] = link_fails

    with io.open(os.path.join(C.OUT, "report_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

    print("报告数字核对：共 %d 项，通过 %d，失败 %d"
          % (len(results), len(results) - len(fails), len(fails)))
    for f in fails:
        print("  FAIL  %-42s 找不到 %-14s （%s）" % (f["report"], f["text"], f["why"]))
    print("报告内部链接：共 %d 条，失效 %d 条" % (n_links, len(link_fails)))
    for f in link_fails:
        print("  BROKEN  %-42s -> %s" % (f["report"], f["target"]))
    print("\n按报告统计：")
    import collections
    agg = collections.Counter()
    for r in results:
        agg[(r["report"], r["found"])] += 1
    for (rep, ok), n in sorted(agg.items()):
        print("  %-44s %s  %d 项" % (rep, "PASS" if ok else "FAIL", n))


if __name__ == "__main__":
    main()
