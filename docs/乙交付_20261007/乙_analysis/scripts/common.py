# -*- coding: utf-8 -*-
"""
乙 - 公共模块：统一的数据读取与口径

口径（固定，不许在分析脚本里临时改）：
  * 窗口       : 81 bp，相对 TSS 的 [-60, +20]
  * TSS 位置   : 序列内 1-based 第 61 位；position index i 的 TSS 相对偏移 = i - 61
  * 正/负样本  : label=1 正样本；label=0 天然对照（只有 TJU 有）
  * 打乱对照   : *.dinucleotide_null.fasta（二核苷酸保持打乱）
  * σ 分组     : 只取"单标签"序列（一条序列只挂一个非 unknown 的 σ），多标签不硬选
  * 统计用划分 : holdout 一律不做统计；σ 组间比较只用 discovery
"""
import collections
import csv
import io
import os

# ---------------------------------------------------------------- 路径
HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)
WORKSPACE = os.path.dirname(ANALYSIS)
DATA = os.path.join(ANALYSIS, "data", "data", "processed", "pr01_02_data_v2")
PRODUCE = os.path.join(WORKSPACE, "乙_产出")
TABLES = os.path.join(PRODUCE, "tables")
FIGURES = os.path.join(PRODUCE, "figures")
OUT = os.path.join(ANALYSIS, "out")

for _d in (PRODUCE, TABLES, FIGURES, OUT):
    os.makedirs(_d, exist_ok=True)

# ---------------------------------------------------------------- 口径常量
TSS_INDEX_1BASED = 61          # 序列内 1-based TSS 位置
WINDOW_UPSTREAM = 60           # 相对 TSS 上游 60 bp
WINDOW_DOWNSTREAM = 20         # 相对 TSS 下游 20 bp
SEQ_LEN = 81

# 教科书已知元件（M2 用的是"已知答案"，不是从数据里发现的 motif）
KNOWN_ELEMENTS = {
    "-35 box (TTGACA)": "TTGACA",
    "-10 box (TATAAT)": "TATAAT",
}

SPLITS = ["discovery", "development", "holdout"]
SPLIT_CN = {"discovery": "发现集(训练)", "development": "开发集(调参)",
            "holdout": "留出集(验证)", "external_validation": "外部验证集"}

SOURCE_CN = {
    "regulondb_ecoli": "RegulonDB 大肠杆菌",
    "dbtbs_bsub": "DBTBS 枯草芽孢杆菌",
    "tjupan_bacillus_subtilis": "TJU 枯草芽孢杆菌",
    "tjupan_baumannii": "TJU 鲍曼不动杆菌",
    "tjupan_bradyrhizobium": "TJU 慢生根瘤菌",
    "tjupan_diphtheria": "TJU 白喉杆菌",
    "tjupan_escherichia_coli": "TJU 大肠杆菌",
    "tjupan_staphylococcus": "TJU 金黄色葡萄球菌",
}

# 甲_σ组可用性清单_20261006.md 的判定结果（乙不重新评级，直接沿用）
ELIGIBILITY_ORDER = {"confirmation_candidate": 0, "exploratory_only": 1, "descriptive_only": 2}
ELIGIBILITY_CN = {"confirmation_candidate": "可做正式结论",
                  "exploratory_only": "只能探索",
                  "descriptive_only": "只能描述"}
# 只有这两组够做正式结论；组间比较的"主比较"在它们之间，其余标注探索性
FORMAL_GROUPS = ["Sigma70", "SigA"]


# ---------------------------------------------------------------- 读取
def read_tsv(path):
    with io.open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield row


def load_sequence_master():
    return list(read_tsv(os.path.join(DATA, "sequence_master.tsv")))


def load_sigma_associations():
    return list(read_tsv(os.path.join(DATA, "sigma_associations.tsv")))


def load_sigma_group_summary():
    return list(read_tsv(os.path.join(DATA, "sigma_group_summary.tsv")))


def load_manifest():
    import json
    with io.open(os.path.join(DATA, "manifest.json"), "r", encoding="utf-8") as fh:
        return json.load(fh)


def read_fasta(path):
    """产出 (header, sequence)；序列统一大写。"""
    name, buf = None, []
    with io.open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name, buf = line[1:], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf).upper()


# ---------------------------------------------------------------- 派生
def build_sigma_index(assoc_rows):
    """
    返回 {sequence_id: {"sigma_raw": set, "sigma_id": set, "n_rows": int}}
    多标签一律全量保留，不硬选。
    """
    idx = collections.defaultdict(lambda: {"sigma_raw": set(), "sigma_id": set(), "n_rows": 0})
    for r in assoc_rows:
        sid = r["sequence_id"]
        idx[sid]["n_rows"] += 1
        idx[sid]["sigma_raw"].add((r.get("sigma_raw") or "").strip())
        idx[sid]["sigma_id"].add((r.get("sigma_id") or "").strip())
    return idx


def single_sigma_label(sigma_sets):
    """
    单标签判定：非 unknown 的 σ 名恰好 1 个 → 返回该名；否则返回 None。
    多标签或 unknown 都不强行归类。
    """
    named = {s for s in sigma_sets if s and s.lower() != "unknown"}
    if len(named) == 1:
        return sorted(named)[0]
    return None


def gc_percent(seq):
    if not seq:
        return float("nan")
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    return 100.0 * gc / len(s)


def offset_of_index(i_1based):
    """序列内 1-based 下标 → 相对 TSS 的偏移（TSS 本身为 0）。"""
    return i_1based - TSS_INDEX_1BASED


def write_tsv(path, header, rows):
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v) for v in r) + "\n")
    return path


def fmt(x, nd=3):
    if x is None:
        return ""
    try:
        f = float(x)
    except (TypeError, ValueError):
        return str(x)
    if f != f:  # NaN
        return "NA"
    return ("%%.%df" % nd) % f


def fmt_p(p, nd=4):
    """p 值格式化：小于 1e-4 时用科学计数法，避免出现无信息的 0.0000。"""
    if p is None:
        return ""
    try:
        f = float(p)
    except (TypeError, ValueError):
        return str(p)
    if f != f:
        return "NA"
    if f < 1e-4:
        return "%.2e" % f
    return ("%%.%df" % nd) % f
