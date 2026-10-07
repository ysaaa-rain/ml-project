# -*- coding: utf-8 -*-
"""
乙 - 交付包自检

用途：拿到这个包的人，不需要原始工作目录，就能验证：
  1. 乙的报告之间的内部链接是否都能打开；
  2. 报告引用的产出文件（表、图、JSON）是否都在包里；
  3. 各文件的 SHA256 是否与随包的清单一致。

用法（在本包的 乙_analysis/ 目录下运行）：
    python scripts\\07_check_package.py                  # 自检
    python scripts\\07_check_package.py --write-manifest # 生成/刷新哈希清单

判定规则：
  * **乙交付文件**（`乙*.md`、`README.md`）的链接**必须**全部有效，失效即判 FAIL；
  * **随包参考文档**（`M1M2*.md`、`甲*.md`）的链接只作提示，不影响结论——
    它们是别人写的文件，内含的失效引用不由乙负责，也不应被乙擅自修改。

注意：本脚本假定自己位于 `<包>/乙_analysis/scripts/`。在原始工作目录里运行没有意义
（那里没有随包的 README.md），请只在打包后的目录内运行。
"""
import hashlib
import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)          # <包>/乙_analysis
PACKAGE = os.path.dirname(ANALYSIS)       # <包>
OUT = os.path.join(ANALYSIS, "out")
MANIFEST_NAME = "package_manifest.json"
MANIFEST_REL = os.path.join("乙_analysis", "out", MANIFEST_NAME).replace(os.sep, "/")


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def own_report(name):
    """乙自己写的报告 → 它的链接必须有效。"""
    return name == "README.md" or name.startswith("乙")


def write_manifest():
    files = {}
    for dirpath, dirnames, filenames in os.walk(PACKAGE):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, PACKAGE).replace(os.sep, "/")
            if rel == MANIFEST_REL:
                continue
            files[rel] = sha256(full)
    payload = {
        "package": os.path.basename(PACKAGE),
        "created": "2026-10-07",
        "file_count": len(files),
        "files": dict(sorted(files.items())),
    }
    os.makedirs(OUT, exist_ok=True)
    with io.open(os.path.join(OUT, MANIFEST_NAME), "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    print("已写出哈希清单：%d 个文件 → %s" % (len(files), MANIFEST_REL))


def main():
    # ---- 0. 是不是在包里运行 ----
    if not os.path.exists(os.path.join(PACKAGE, "README.md")):
        print("警告：%s 下没有 README.md —— 这里不像打包后的交付目录。" % PACKAGE)
        print("      本脚本请只在交付包内运行。")

    # ---- 1. 报告内部链接 ----
    reports = [n for n in sorted(os.listdir(PACKAGE)) if n.endswith(".md")]
    own_fails, ref_fails, n_links = [], [], 0
    referenced = set()
    for name in reports:
        with io.open(os.path.join(PACKAGE, name), encoding="utf-8") as fh:
            content = fh.read().replace("\u2212", "-")
        for m in re.finditer(r"\]\(([^)#]+?)\)", content):
            target = m.group(1).strip()
            if target.startswith(("http://", "https://", "#")):
                continue
            n_links += 1
            referenced.add(target)
            if not os.path.exists(os.path.join(PACKAGE, target.replace("/", os.sep))):
                (own_fails if own_report(name) else ref_fails).append((name, target))

    # ---- 2. 引用但缺失的产出文件 ----
    missing_outputs = [t for t in sorted(referenced)
                       if (t.startswith("乙_产出") or t.startswith("乙_analysis"))
                       and not os.path.exists(os.path.join(PACKAGE, t.replace("/", os.sep)))]

    print("包目录            : %s" % PACKAGE)
    print("报告文件数        : %d（其中乙交付 %d，随包参考 %d）"
          % (len(reports), len([r for r in reports if own_report(r)]),
             len([r for r in reports if not own_report(r)])))
    print()
    print("内部链接          : %d 条" % n_links)
    print("  乙交付文件失效  : %d 条%s"
          % (len(own_fails), "" if not own_fails else "  <- 判为 FAIL"))
    for r, t in own_fails:
        print("     BROKEN  %-44s -> %s" % (r, t))
    print("  参考文档失效    : %d 条（提示，不计入结论）" % len(ref_fails))
    for r, t in ref_fails:
        print("     note    %-44s -> %s" % (r, t))
    print("引用但缺失的产出  : %d 个" % len(missing_outputs))
    for t in missing_outputs:
        print("   MISSING %s" % t)

    tdir = os.path.join(PACKAGE, "乙_产出", "tables")
    fdir = os.path.join(PACKAGE, "乙_产出", "figures")
    sdir = os.path.join(PACKAGE, "乙_analysis", "scripts")
    n_tables = len([f for f in os.listdir(tdir) if f.endswith(".tsv")]) if os.path.isdir(tdir) else 0
    n_figs = len([f for f in os.listdir(fdir) if f.endswith(".png")]) if os.path.isdir(fdir) else 0
    n_scripts = len([f for f in os.listdir(sdir) if f.endswith(".py")]) if os.path.isdir(sdir) else 0
    print()
    print("结果表            : %d 张" % n_tables)
    print("图                : %d 张" % n_figs)
    print("脚本              : %d 个" % n_scripts)

    # ---- 3. 哈希校验 ----
    man_path = os.path.join(OUT, MANIFEST_NAME)
    hash_fail = []
    if os.path.exists(man_path):
        with io.open(man_path, encoding="utf-8") as fh:
            manifest = json.load(fh)
        files = manifest.get("files", {})
        ok = 0
        for rel, expect in files.items():
            p = os.path.join(PACKAGE, rel.replace("/", os.sep))
            if not os.path.exists(p):
                hash_fail.append((rel, "缺失"))
            elif sha256(p) != expect:
                hash_fail.append((rel, "哈希不一致"))
            else:
                ok += 1
        print("哈希校验          : %d/%d 一致，%d 个问题" % (ok, len(files), len(hash_fail)))
        for rel, why in hash_fail:
            print("   %-12s %s" % (why, rel))
        if manifest.get("file_count") != len(files):
            print("   note 清单自称 %s 个文件，实际列出 %d 个"
                  % (manifest.get("file_count"), len(files)))
    else:
        print("哈希校验          : 未找到 %s（可加 --write-manifest 生成）" % MANIFEST_REL)

    bad = len(own_fails) + len(missing_outputs) + len(hash_fail)
    print()
    print("自检结论：%s" % ("通过" if bad == 0 else "发现 %d 个问题" % bad))
    if ref_fails:
        print("（另有 %d 条参考文档内的失效引用，属甲文件既有问题，见包 README 已知事项）"
              % len(ref_fails))
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    if "--write-manifest" in sys.argv:
        write_manifest()
    else:
        sys.exit(main())
