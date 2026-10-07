# -*- coding: utf-8 -*-
"""
乙 - 打包脚本（可复现地生成交付包）

干什么：把乙的全部交付物 + 报告引用的参考文档，装进一个文件夹，并写出哈希清单。

为什么要有这个脚本：
  手工复制踩过两次坑——
    (1) 只存在于包内的文件（README.md、07_check_package.py）在重建时被抹掉；
    (2) 包说明文件放在工作目录根层时，会被 `乙*.md` 规则连带复制进包，造成重复。
  所以把打包规则写成脚本，谁重建都不会再犯。

设计要点：
  * 报告一律**平铺**在包根层——乙的报告里用相对路径引用甲的报告与乙_产出/，
    平铺后这些链接在包内依然有效。
  * `乙_analysis/repo/` 只保留 `reports/M1M2报告框架.md`（乙的报告引用它），
    保持与原工作目录**相同的相对路径**，链接才不会断。
  * **不打包**：`data/`（派生数据，按甲的许可表不得再分发）、
    `repo/` 其余内容（可从公开仓库取得）、`__pycache__`。
  * 包首页 `README.md` 的来源是 `乙_analysis/package_readme.md`
    （故意不带 `乙` 前缀，避免被 `乙*.md` 规则复制成第二份）。

用法（工作目录 = 项目根）：
    python 乙_analysis\\scripts\\09_build_package.py
打包后自检：
    cd 乙交付_20261007\\乙_analysis && python scripts\\07_check_package.py
"""
import hashlib
import io
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)          # <根>/乙_analysis
ROOT = os.path.dirname(ANALYSIS)          # 项目根
PACKAGE_NAME = "乙交付_20261007"
PACKAGE = os.path.join(ROOT, PACKAGE_NAME)

README_SOURCE = os.path.join(ANALYSIS, "package_readme.md")
FRAMEWORK_REL = os.path.join("repo", "reports", "M1M2报告框架.md")
MANIFEST_REL = os.path.join("乙_analysis", "out", "package_manifest.json")


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def copytree_clean(src, dst):
    shutil.copytree(src, dst)


def main():
    if not os.path.exists(README_SOURCE):
        print("缺少包首页来源文件：%s" % README_SOURCE, file=sys.stderr)
        return 2

    if os.path.exists(PACKAGE):
        shutil.rmtree(PACKAGE)
    os.makedirs(os.path.join(PACKAGE, "乙_产出"))
    os.makedirs(os.path.join(PACKAGE, "乙_analysis"))

    n_report = 0
    # 1) 报告（平铺）
    for name in sorted(os.listdir(ROOT)):
        if not name.endswith(".md"):
            continue
        if name.startswith("乙") or name.startswith("甲") or name.startswith("M1M2"):
            shutil.copy2(os.path.join(ROOT, name), os.path.join(PACKAGE, name))
            n_report += 1

    # 2) 表格与图
    copytree_clean(os.path.join(ROOT, "乙_产出", "tables"),
                   os.path.join(PACKAGE, "乙_产出", "tables"))
    copytree_clean(os.path.join(ROOT, "乙_产出", "figures"),
                   os.path.join(PACKAGE, "乙_产出", "figures"))

    # 3) 脚本、机器可读结果、下载溯源
    copytree_clean(os.path.join(ANALYSIS, "scripts"),
                   os.path.join(PACKAGE, "乙_analysis", "scripts"))
    copytree_clean(os.path.join(ANALYSIS, "out"),
                   os.path.join(PACKAGE, "乙_analysis", "out"))
    shutil.copy2(os.path.join(ANALYSIS, "download_provenance.json"),
                 os.path.join(PACKAGE, "乙_analysis", "download_provenance.json"))

    # 4) 报告框架（保持相对路径不变）
    fw_src = os.path.join(ANALYSIS, FRAMEWORK_REL)
    fw_dst = os.path.join(PACKAGE, "乙_analysis", FRAMEWORK_REL)
    os.makedirs(os.path.dirname(fw_dst), exist_ok=True)
    shutil.copy2(fw_src, fw_dst)

    # 5) 包首页
    shutil.copy2(README_SOURCE, os.path.join(PACKAGE, "README.md"))

    # 6) 清理构建产物
    for dirpath, dirnames, _f in os.walk(PACKAGE):
        for d in list(dirnames):
            if d == "__pycache__":
                shutil.rmtree(os.path.join(dirpath, d))
                dirnames.remove(d)

    # 7) 哈希清单（不含清单自身）
    files = {}
    for dirpath, dirnames, filenames in os.walk(PACKAGE):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, PACKAGE).replace(os.sep, "/")
            if rel == MANIFEST_REL.replace(os.sep, "/"):
                continue
            files[rel] = sha256(full)
    payload = {"package": PACKAGE_NAME, "created": "2026-10-07",
               "file_count": len(files), "files": dict(sorted(files.items()))}
    with io.open(os.path.join(PACKAGE, MANIFEST_REL), "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)

    total = sum(len(f) for _d, _s, f in os.walk(PACKAGE))
    size = sum(os.path.getsize(os.path.join(d, f))
               for d, _s, fs in os.walk(PACKAGE) for f in fs)
    print("包目录        : %s" % PACKAGE)
    print("根层报告      : %d 份" % n_report)
    print("结果表        : %d 张"
          % len([f for f in os.listdir(os.path.join(PACKAGE, "乙_产出", "tables"))
                 if f.endswith(".tsv")]))
    print("图            : %d 张"
          % len([f for f in os.listdir(os.path.join(PACKAGE, "乙_产出", "figures"))
                 if f.endswith(".png")]))
    print("脚本          : %d 个"
          % len([f for f in os.listdir(os.path.join(PACKAGE, "乙_analysis", "scripts"))
                 if f.endswith(".py")]))
    print("文件总数      : %d（含哈希清单）" % total)
    print("体积          : %.2f MB" % (size / 1024.0 / 1024.0))
    print("哈希清单      : %d 个文件已登记" % len(files))
    print()
    print("下一步：")
    print("  cd %s\\乙_analysis" % PACKAGE_NAME)
    print("  python scripts\\07_check_package.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
