# -*- coding: utf-8 -*-
"""
乙 - 步骤 0：获取真实数据

为什么用 zip 而不是逐文件下载：
  - raw.githubusercontent.com 在本机不可达（连接超时）；
  - GitHub 未认证 API 限额 60 次/小时，而数据有 180+ 个文件，逐文件会超限；
  - codeload 的仓库 zip 一次请求即可拿到全部分支内容。

数据来源：https://github.com/ysaaa-rain/ml-project  main 分支
许可提示：data/processed/ 属派生数据，按 甲_许可状态表 不得提交回 Git，
          仅在本机用于 M2 分析。
"""
import hashlib
import io
import json
import os
import sys
import time
import urllib.request
import zipfile

REPO = "ysaaa-rain/ml-project"
BRANCH = "main"
ZIP_URL = "https://codeload.github.com/%s/zip/refs/heads/%s" % (REPO, BRANCH)

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
DATA_DIR = os.path.join(BASE, "data")
ZIP_PATH = os.path.join(BASE, "_repo_main.zip")


def fetch(url, tries=4, timeout=600):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 dsh"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = "%s: %s" % (type(e).__name__, e)
            print("  下载失败（第 %d 次）：%s" % (i + 1, last), file=sys.stderr)
            time.sleep(3)
    raise RuntimeError("下载失败：%s -> %s" % (url, last))


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    if not os.path.exists(ZIP_PATH):
        print("下载仓库 zip：%s" % ZIP_URL)
        blob = fetch(ZIP_URL)
        with open(ZIP_PATH, "wb") as fh:
            fh.write(blob)
        print("  已保存 %s（%d 字节）" % (ZIP_PATH, len(blob)))
    else:
        print("已存在 zip，跳过下载：%s" % ZIP_PATH)

    with open(ZIP_PATH, "rb") as fh:
        zip_sha = hashlib.sha256(fh.read()).hexdigest()
    print("  zip SHA256 = %s" % zip_sha)

    root = "ml-project-%s/" % BRANCH
    wanted_prefix = root + "data/processed/pr01_02_data_v2/"
    n_data = 0
    with zipfile.ZipFile(ZIP_PATH) as zf:
        names = [n for n in zf.namelist() if n.startswith(wanted_prefix) and not n.endswith("/")]
        print("zip 内 pr01_02_data_v2 文件数 = %d" % len(names))
        for name in names:
            rel = name[len(root):]
            dest = os.path.join(DATA_DIR, rel.replace("/", os.sep))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with zf.open(name) as src, open(dest, "wb") as out:
                out.write(src.read())
            n_data += 1
        # 顺带取出 manifest 里 cited 的文档，便于核对
        for extra in [root + "README.md", root + "AGENTS.md"]:
            try:
                with zf.open(extra) as src:
                    dest = os.path.join(DATA_DIR, "_repo", os.path.basename(extra))
                    os.makedirs(os.path.dirname(dest), exist_ok=True)
                    with open(dest, "wb") as out:
                        out.write(src.read())
            except KeyError:
                pass

    print("已解压 %d 个文件到 %s" % (n_data, DATA_DIR))

    prov = {
        "source_repo": "https://github.com/%s" % REPO,
        "branch": BRANCH,
        "download_url": ZIP_URL,
        "zip_sha256": zip_sha,
        "dest": DATA_DIR,
        "files_extracted": n_data,
    }
    with open(os.path.join(BASE, "download_provenance.json"), "w", encoding="utf-8") as fh:
        json.dump(prov, fh, ensure_ascii=False, indent=2)
    print("溯源写入 download_provenance.json")


if __name__ == "__main__":
    main()
