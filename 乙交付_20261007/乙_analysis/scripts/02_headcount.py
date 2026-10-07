# -*- coding: utf-8 -*-
"""
乙 - 步骤 2：数人头

干什么：数每个物种、每个 σ 组、每个来源在三个划分里各有多少条。
为什么：样本量决定后面哪些能做统计、哪些只能描述；这一步错，后面全错。
怎么做：
  A. 用 sequence_master.tsv 数「来源 × 物种 × 划分 × 标签」；
  B. 用 sigma_associations.tsv 判定单标签 σ，数「σ 组 × 划分」；
  C. 与 sigma_group_summary.tsv、manifest.json counts、以及各 σ 子目录 FASTA 三方交叉核对。
做到什么程度：每个数字都能追到文件；三方核对必须零差异，有差异就报出来不掩盖。
产出：乙_产出/tables/ 下 5 张表 + out/headcount.json
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402


def main():
    sm = C.load_sequence_master()
    assoc = C.load_sigma_associations()
    summary = C.load_sigma_group_summary()
    manifest = C.load_manifest()
    sig_idx = C.build_sigma_index(assoc)

    # ---------- A. 来源 × 物种 × 划分 × 标签 ----------
    cnt = collections.Counter()
    for r in sm:
        cnt[(r["source"], r["species"], r["split"], r["label"])] += 1

    rows_a = []
    for (src, sp, split, label), n in sorted(cnt.items()):
        rows_a.append([src, C.SOURCE_CN.get(src, src), sp, split, C.SPLIT_CN.get(split, split),
                       label, n])
    C.write_tsv(os.path.join(C.TABLES, "T1_来源x物种x划分x标签.tsv"),
                ["source", "source_cn", "species", "split", "split_cn", "label", "n"], rows_a)

    # ---------- 物种级合计（含正负样本） ----------
    sp_cnt = collections.Counter()
    sp_pos = collections.Counter()
    for r in sm:
        sp_cnt[(r["source"], r["species"])] += 1
        if r["label"] == "1":
            sp_pos[(r["source"], r["species"])] += 1
    rows_species = [[src, C.SOURCE_CN.get(src, src), sp, sp_cnt[(src, sp)], sp_pos[(src, sp)],
                     sp_cnt[(src, sp)] - sp_pos[(src, sp)]]
                    for (src, sp) in sorted(sp_cnt)]
    C.write_tsv(os.path.join(C.TABLES, "T2_来源x物种合计.tsv"),
                ["source", "source_cn", "species", "n_all", "n_positive", "n_natural_control"],
                rows_species)

    # ---------- B. σ 组 × 划分（单标签，core） ----------
    # 单标签 σ 序列：一条序列只挂一个非 unknown 的 σ
    per_seq_sigma = {sid: v["sigma_raw"] for sid, v in sig_idx.items()}
    rec_by_id = {r["sequence_id"]: r for r in sm}

    sigma_split = collections.Counter()      # (source, sigma, split) -> n
    sigma_skip_multi, sigma_skip_unknown, sigma_skip_tju = 0, 0, 0
    for sid, sset in per_seq_sigma.items():
        rec = rec_by_id.get(sid)
        if rec is None:
            continue
        if rec["source"].startswith("tjupan_"):
            sigma_skip_tju += 1
            continue
        named = {s for s in sset if s and s.lower() != "unknown"}
        if not named:
            sigma_skip_unknown += 1
            continue
        if len(named) > 1:
            sigma_skip_multi += 1
            continue
        label = sorted(named)[0]
        # 口径：只用 core 层（extended 是同一批记录的低证据等级版本，混算会重复计数）
        if rec["tier"] != "core":
            continue
        sigma_split[(rec["source"], label, rec["split"])] += 1

    rows_b = []
    elig = {(r["source"], r["sigma"]): r for r in summary}
    for (src, sig) in sorted({(s, g) for (s, g, _sp) in sigma_split}):
        tot = sum(sigma_split[(src, sig, sp)] for sp in C.SPLITS)
        e = elig.get((src, sig), {})
        rows_b.append([src, C.SOURCE_CN.get(src, src), sig, tot,
                       sigma_split[(src, sig, "discovery")],
                       sigma_split[(src, sig, "development")],
                       sigma_split[(src, sig, "holdout")],
                       e.get("single_sigma_sequences", ""),
                       e.get("eligibility", ""),
                       C.ELIGIBILITY_CN.get(e.get("eligibility", ""), "")])
    C.write_tsv(os.path.join(C.TABLES, "T3_sigma组x划分.tsv"),
                ["source", "source_cn", "sigma", "total", "discovery", "development",
                 "holdout", "summary_total", "eligibility", "eligibility_cn"], rows_b)

    # ---------- C1. 与 sigma_group_summary.tsv 逐格核对 ----------
    diffs = []
    c1_checks = 0
    for e in summary:
        src, sig = e["source"], e["sigma"]
        for sp, key in [("discovery", "discovery"), ("development", "development"),
                        ("holdout", "holdout")]:
            mine = sigma_split[(src, sig, sp)]
            theirs = int(e[key])
            c1_checks += 1
            if mine != theirs:
                diffs.append({"source": src, "sigma": sig, "split": sp,
                              "derived": mine, "sigma_group_summary": theirs})
        tot_mine = sum(sigma_split[(src, sig, sp)] for sp in C.SPLITS)
        c1_checks += 1
        if tot_mine != int(e["single_sigma_sequences"]):
            diffs.append({"source": src, "sigma": sig, "split": "total",
                          "derived": tot_mine,
                          "sigma_group_summary": int(e["single_sigma_sequences"])})

    # ---------- C2. 与各 σ 子目录 FASTA 条数核对 ----------
    fasta_check = []
    c2_checked = 0
    for src in ("regulondb_ecoli", "dbtbs_bsub"):
        sdir = os.path.join(C.DATA, src, "core", "sigma")
        if not os.path.isdir(sdir):
            continue
        for sig in sorted(os.listdir(sdir)):
            gdir = os.path.join(sdir, sig)
            if not os.path.isdir(gdir):
                continue
            for split in C.SPLITS:
                p = os.path.join(gdir, "%s.positive.fasta" % split)
                if not os.path.exists(p):
                    continue
                c2_checked += 1
                n_file = sum(1 for _ in C.read_fasta(p))
                n_derived = sigma_split[(src, sig, split)]
                if n_file != n_derived:
                    fasta_check.append({"source": src, "sigma": sig, "split": split,
                                        "fasta": n_file, "derived": n_derived})

    # ---------- C3. 与 manifest counts 核对（正样本 FASTA vs label=1 记录数） ----------
    # 注意：external_validation 在 sequence_master 的 split 列里仍记为 holdout
    # （它是 holdout 中与其他来源不重叠的子集，靠另一条线索标记），
    # 所以对它不做 split 直配，而是单独用序列内容验证"是 holdout 的子集"。
    manifest_check = []
    external_validation_check = []
    c3_checked = 0
    for srcdir, files in manifest["counts"].items():
        src = srcdir.split("/")[1 - 1]
        tier = srcdir.split("/")[1]
        for fn, n in files.items():
            if not fn.endswith(".positive.fasta"):
                continue
            split = fn.split(".")[0]
            if split == "external_validation":
                # 用序列内容判断：external_validation 的序列应全部出现在 holdout 里
                ev_path = os.path.join(C.DATA, srcdir.replace("/", os.sep), fn)
                ho_path = os.path.join(C.DATA, srcdir.replace("/", os.sep),
                                       "holdout.positive.fasta")
                ev_seqs = [s for _h, s in C.read_fasta(ev_path)] if os.path.exists(ev_path) else []
                ho_seqs = set(s for _h, s in C.read_fasta(ho_path)) if os.path.exists(ho_path) else set()
                subset = all(s in ho_seqs for s in ev_seqs)
                external_validation_check.append({
                    "srcdir": srcdir, "file": fn, "manifest_count": n,
                    "n_records": len(ev_seqs), "holdout_n": len(ho_seqs),
                    "is_subset_of_holdout": subset,
                    "count_matches_records": n == len(ev_seqs),
                })
                continue
            mine = sum(1 for r in sm
                       if r["source"] == src and r["tier"] == tier
                       and r["split"] == split and r["label"] == "1")
            c3_checked += 1
            if mine != n:
                manifest_check.append({"srcdir": srcdir, "file": fn,
                                       "manifest": n, "from_sequence_master": mine})

    # ---------- core / extended 关系检查 ----------
    core_ids = {r["sequence_id"] for r in sm if r["source"] == "regulondb_ecoli"
                and r["tier"] == "core"}
    ext_ids = {r["sequence_id"] for r in sm if r["source"] == "regulondb_ecoli"
               and r["tier"] == "extended"}
    overlap = len(core_ids & ext_ids)

    result = {
        "sequence_master_rows": len(sm),
        "sigma_association_rows": len(assoc),
        "sigma_group_summary_rows": len(summary),
        "A_rows": len(rows_a),
        "B_sigma_groups": len(rows_b),
        "C1_summary_vs_derived_diffs": diffs,
        "C1_checks": c1_checks,
        "C2_fasta_vs_derived_diffs": fasta_check,
        "C2_checked": c2_checked,
        "C3_manifest_vs_sequence_master_diffs": manifest_check,
        "C3_checked": c3_checked,
        "C3b_external_validation": external_validation_check,
        "sigma_total": sum(sigma_split.values()),
        "sigma_total_by_split": {sp: sum(v for (s, g, sp2), v in sigma_split.items()
                                        if sp2 == sp) for sp in C.SPLITS},
        "sigma_skipped": {"tju_no_sigma": sigma_skip_tju,
                          "unknown": sigma_skip_unknown,
                          "multi_sigma": sigma_skip_multi},
        "regulondb_core_extended_overlap": overlap,
        "regulondb_core_n": len(core_ids),
        "regulondb_extended_n": len(ext_ids),
    }
    with open(os.path.join(C.OUT, "headcount.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)

    # ---------- 控制台 ----------
    print("sequence_master 行数        : %d" % len(sm))
    print("sigma_associations 行数     : %d" % len(assoc))
    print("σ 组数（推导，core 单标签） : %d" % len(rows_b))
    print("被排除在 σ 分组之外         : TJU 无标签 %d / unknown %d / 多标签 %d"
          % (sigma_skip_tju, sigma_skip_unknown, sigma_skip_multi))
    print("RegulonDB core∩extended 重叠 : %d（core %d，extended %d）"
          % (overlap, len(core_ids), len(ext_ids)))
    print("C1 sigma_group_summary 核对 : %d 处不一致（核对 %d 格）" % (len(diffs), c1_checks))
    for d in diffs[:12]:
        print("   %s %s %s 推导=%s 清单=%s"
              % (d["source"], d["sigma"], d["split"], d["derived"], d["sigma_group_summary"]))
    print("C2 σ 子目录 FASTA 核对      : %d 处不一致（核对 %d 个文件）"
          % (len(fasta_check), c2_checked))
    for d in fasta_check[:12]:
        print("   %s %s %s FASTA=%s 推导=%s"
              % (d["source"], d["sigma"], d["split"], d["fasta"], d["derived"]))
    print("C3 manifest counts 核对     : %d 处不一致（核对 %d 个文件）"
          % (len(manifest_check), c3_checked))
    for d in manifest_check[:12]:
        print("   %s 清单=%s 推导=%s" % (d["file"], d["manifest"], d["from_sequence_master"]))
    print("C3b external_validation     : %d 个文件" % len(external_validation_check))
    for d in external_validation_check:
        print("   %-24s 清单=%d 实测=%d 是holdout子集=%s"
              % (d["srcdir"], d["manifest_count"], d["n_records"], d["is_subset_of_holdout"]))

    print("\n=== 物种合计（按来源）===")
    for r in rows_species:
        print("  %-30s %-22s 全部=%5d 正=%5d 对照=%5d" % (r[1], r[2], r[3], r[4], r[5]))

    print("\n=== σ 组 × 划分 ===")
    for r in sorted(rows_b, key=lambda x: (x[0], -x[3])):
        print("  %-16s %-9s 合计=%4d 发现=%4d 开发=%3d 留出=%4d  [%s]"
              % (r[0], r[2], r[3], r[4], r[5], r[6], r[9]))


if __name__ == "__main__":
    main()
