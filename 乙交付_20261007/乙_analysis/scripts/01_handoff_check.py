# -*- coding: utf-8 -*-
"""
乙 - 步骤 1：接手核验（甲→乙 交接闸门）

目的：确认乙手上这份数据，与 甲_验收记录_20261006.md 验收过的是同一版。
方法：
  A. 复用 甲 第 1 步的方法——按 manifest.json 的 output_sha256 逐文件校验 SHA256；
  B. 独立复算条数：把 manifest.counts 展平成文件→条数，逐个 FASTA 数 '>' 行；
  C. 独立复算打乱对照的配对不变式：dinucleotide_null 条数应等于同前缀 positive 条数。
乙不直接采信甲的数字，全部自己再数一遍。
产出：out/handoff_check.json
"""
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
DATA = os.path.join(BASE, "data", "data", "processed", "pr01_02_data_v2")
OUT = os.path.join(BASE, "out")
os.makedirs(OUT, exist_ok=True)


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def count_fasta_records(path):
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def main():
    with open(os.path.join(DATA, "manifest.json"), "r", encoding="utf-8") as fh:
        manifest = json.load(fh)

    # --- A. SHA256 ---
    registered = manifest["output_sha256"]
    ok, bad, missing = 0, [], []
    for rel, expected in registered.items():
        p = os.path.join(DATA, rel.replace("/", os.sep))
        if not os.path.exists(p):
            missing.append(rel)
            continue
        actual = sha256(p)
        if actual.lower() == expected.lower():
            ok += 1
        else:
            bad.append({"file": rel, "expected": expected, "actual": actual})

    on_disk = []
    for dirpath, _d, filenames in os.walk(DATA):
        for fn in filenames:
            on_disk.append(os.path.relpath(os.path.join(dirpath, fn), DATA).replace(os.sep, "/"))

    # --- B. counts 展平后逐文件复算 ---
    flat = {}
    for srcdir, files in manifest["counts"].items():
        for fn, n in files.items():
            flat["%s/%s" % (srcdir, fn)] = n

    count_mismatch = []
    for rel, expected in sorted(flat.items()):
        p = os.path.join(DATA, rel.replace("/", os.sep))
        if not os.path.exists(p):
            count_mismatch.append({"file": rel, "expected": expected, "actual": None,
                                   "note": "文件不存在"})
            continue
        actual = count_fasta_records(p)
        if actual != expected:
            count_mismatch.append({"file": rel, "expected": expected, "actual": actual})

    # --- C. 打乱对照配对不变式 ---
    null_files = [r for r in on_disk if r.endswith(".dinucleotide_null.fasta")]
    null_check = []
    for rel in sorted(null_files):
        pos_rel = rel.replace(".dinucleotide_null.fasta", ".positive.fasta")
        p_null = os.path.join(DATA, rel.replace("/", os.sep))
        p_pos = os.path.join(DATA, pos_rel.replace("/", os.sep))
        n_null = count_fasta_records(p_null)
        n_pos = count_fasta_records(p_pos) if os.path.exists(p_pos) else None
        null_check.append({"null_file": rel, "null_records": n_null,
                           "positive_file": pos_rel, "positive_records": n_pos,
                           "paired": n_pos is not None and n_null == n_pos})

    # --- D. 汇总表行数（不含表头） ---
    table_rows = {}
    for name, field in [("sequence_master.tsv", "statistical_units"),
                        ("source_record_mapping.tsv", "source_records_retained"),
                        ("sigma_associations.tsv", "sigma_associations"),
                        ("leakage_links.tsv", "leakage_graph_links"),
                        ("tju_sigma_annotation_candidates.tsv", "annotation_candidate_links")]:
        p = os.path.join(DATA, name)
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            n = sum(1 for _ in fh) - 1
        table_rows[name] = {"field": field, "manifest": manifest.get(field), "actual": n,
                            "match": manifest.get(field) == n}
    # exclusions 行数
    p = os.path.join(DATA, "exclusions.tsv")
    with open(p, "r", encoding="utf-8", errors="replace") as fh:
        n_excl = sum(1 for _ in fh) - 1
    table_rows["exclusions.tsv"] = {"field": "excluded_reasons(sum)",
                                    "manifest": sum(manifest["excluded_reasons"].values()),
                                    "actual": n_excl,
                                    "match": sum(manifest["excluded_reasons"].values()) == n_excl}

    report = {
        "data_dir": DATA,
        "schema": manifest["schema"],
        "status": manifest["status"],
        "seed": manifest["seed"],
        "tss_index_1based": manifest["tss_index_1based"],
        "primary": manifest["primary"],
        "auxiliary": manifest["auxiliary"],
        "A_sha256": {
            "registered": len(registered), "ok": ok,
            "bad_count": len(bad), "bad": bad,
            "missing_count": len(missing), "missing": missing,
            "files_on_disk": len(on_disk),
            "disk_minus_registered": len(on_disk) - len(registered),
        },
        "B_counts": {
            "leaf_items": len(flat), "mismatch_count": len(count_mismatch),
            "mismatch": count_mismatch,
        },
        "C_null_pairing": {
            "null_files": len(null_check),
            "paired_ok": sum(1 for x in null_check if x["paired"]),
            "unpaired": [x for x in null_check if not x["paired"]],
        },
        "D_tables": table_rows,
        "known_limits": manifest["known_limits"],
        "excluded_reasons": manifest["excluded_reasons"],
        "sigma_group_count": len(manifest["sigma_groups"]),
    }

    with open(os.path.join(OUT, "handoff_check.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=2)

    print("schema/status      : %s / %s" % (report["schema"], report["status"]))
    print("seed               : %s" % report["seed"])
    print("tss_index_1based   : %s" % report["tss_index_1based"])
    print("A. SHA256 登记 %d，一致 %d，不匹配 %d，缺失 %d；磁盘 %d 个文件（多 %d）"
          % (len(registered), ok, len(bad), len(missing), len(on_disk),
             len(on_disk) - len(registered)))
    print("B. counts 展平 %d 项，不一致 %d 项" % (len(flat), len(count_mismatch)))
    for m in count_mismatch[:10]:
        print("   %s 登记=%s 实测=%s" % (m["file"], m["expected"], m["actual"]))
    print("C. dinucleotide_null %d 个，配对成立 %d 个"
          % (len(null_check), report["C_null_pairing"]["paired_ok"]))
    print("D. 汇总表行数：")
    for name, info in table_rows.items():
        print("   %-38s 登记=%-7s 实测=%-7s %s"
              % (name, info["manifest"], info["actual"], "OK" if info["match"] else "不一致"))
    print("σ 组数（manifest.sigma_groups）: %d" % len(manifest["sigma_groups"]))


if __name__ == "__main__":
    main()
