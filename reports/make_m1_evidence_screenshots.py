"""Create presentation-ready evidence boards from the real M1 demo outputs."""

from __future__ import annotations

import json
from pathlib import Path
import textwrap

from PIL import Image, ImageDraw, ImageFont
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "reports"
DEMO = ROOT / "tmp" / "m1_demo_configured"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = (
        Path("C:/Windows/Fonts/msyhbd.ttc") if bold else Path("C:/Windows/Fonts/msyh.ttc"),
        Path("C:/Windows/Fonts/simhei.ttf"),
        Path("C:/Windows/Fonts/arial.ttf"),
    )
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default()


F_TITLE = font(34, True)
F_SUBTITLE = font(20)
F_HEAD = font(22, True)
F_BODY = font(18)
F_SMALL = font(15)
F_MONO = font(15)

NAVY = "#17324D"
BLUE = "#2F75B5"
LIGHT_BLUE = "#EAF3FA"
GREEN = "#2E7D5B"
LIGHT_GREEN = "#EAF6EF"
ORANGE = "#C86B2C"
LIGHT_ORANGE = "#FFF1E7"
RED = "#B84444"
LIGHT_RED = "#FCECEC"
GRAY = "#607080"
LIGHT_GRAY = "#F3F6F8"
GRID = "#D7E0E8"
WHITE = "#FFFFFF"


def board(title: str, subtitle: str) -> tuple[Image.Image, ImageDraw.ImageDraw]:
    image = Image.new("RGB", (1600, 900), WHITE)
    draw = ImageDraw.Draw(image)
    draw.rectangle((0, 0, 1600, 112), fill=NAVY)
    draw.text((54, 26), title, fill=WHITE, font=F_TITLE)
    draw.text((56, 72), subtitle, fill="#C8D8E8", font=F_SUBTITLE)
    draw.text((1420, 34), "M1 / DEMO", fill="#C8D8E8", font=F_SMALL)
    return image, draw


def card(draw: ImageDraw.ImageDraw, xy: tuple[int, int, int, int], title: str, value: str, note: str, fill: str = LIGHT_BLUE) -> None:
    x1, y1, x2, y2 = xy
    draw.rounded_rectangle(xy, radius=16, fill=fill, outline=GRID, width=2)
    draw.text((x1 + 20, y1 + 16), title, fill=GRAY, font=F_SMALL)
    draw.text((x1 + 20, y1 + 42), value, fill=NAVY, font=font(27, True))
    draw.text((x1 + 20, y2 - 23), note, fill=GRAY, font=F_SMALL)


def table(draw: ImageDraw.ImageDraw, x: int, y: int, columns: list[tuple[str, int]], rows: list[list[str]], row_height: int = 34) -> None:
    current_x = x
    for label, width in columns:
        draw.rectangle((current_x, y, current_x + width, y + row_height), fill=NAVY)
        draw.text((current_x + 10, y + 7), label, fill=WHITE, font=F_SMALL)
        current_x += width
    for row_index, row in enumerate(rows, start=1):
        current_x = x
        fill = WHITE if row_index % 2 else LIGHT_GRAY
        for value, (_, width) in zip(row, columns):
            draw.rectangle((current_x, y + row_index * row_height, current_x + width, y + (row_index + 1) * row_height), fill=fill, outline=GRID, width=1)
            draw.text((current_x + 10, y + row_index * row_height + 7), str(value), fill=NAVY, font=F_SMALL)
            current_x += width


def make_clean_board() -> None:
    report = json.loads((DEMO / "quality_report.json").read_text(encoding="utf-8"))
    clean = pd.read_csv(DEMO / "promoters_clean.tsv", sep="\t")
    image, draw = board(
        "实际完成：数据清洗与分组",
        "真实运行产物：tmp/m1_demo_configured/quality_report.json + promoters_clean.tsv",
    )
    card(draw, (55, 145, 335, 245), "输入记录", str(report["input_rows"]), "demo smoke run")
    card(draw, (365, 145, 645, 245), "有效记录", str(report["output_rows"]), "进入后续分析", LIGHT_GREEN)
    card(draw, (675, 145, 955, 245), "过滤/去重", str(report["removed_rows"]), "不是静默删除", LIGHT_ORANGE)
    card(draw, (985, 145, 1265, 245), "发现 / 验证", "4 / 2", "按数据源标记", LIGHT_BLUE)
    card(draw, (1295, 145, 1545, 245), "长度范围", "30–31 bp", "GC 中位数 0.383", LIGHT_GRAY)

    draw.text((60, 285), "减少的 2 条记录具体原因", fill=NAVY, font=F_HEAD)
    removed_rows = [["高 N 含量", "1", ">10% N，进入过滤日志"], ["完全重复序列", "1", "保留第一条，避免重复计数"], ["其他质量问题", "0", "空值/非法字符/缺少身份字段"]]
    table(draw, 60, 322, [("原因", 210), ("数量", 100), ("解释", 600)], removed_rows)

    draw.text((60, 470), "清洗后保留的字段与分组（前 6 条全部展示）", fill=NAVY, font=F_HEAD)
    rows = []
    for row in clean.itertuples():
        rows.append([row.sequence_id, row.species, row.source_dataset, row.sigma_factor_type, row.split, f"{row.sequence_length} bp", f"{row.gc_fraction:.3f}"])
    table(draw, 60, 507, [("ID", 170), ("物种", 220), ("数据源", 130), ("σ 类型", 110), ("分组", 110), ("长度", 100), ("GC", 90)], rows, 32)
    draw.text((60, 830), "说明：以上是最小 demo，用来证明清洗、去重和分组流程能运行，不是正式生物学样本统计。", fill=RED, font=F_SMALL)
    image.save(OUTPUT / "M1证据_数据清洗与分组.png")


def make_b0_board() -> None:
    summary = pd.read_csv(DEMO / "known_elements" / "known_element_summary.tsv", sep="\t")
    hits = pd.read_csv(DEMO / "known_elements" / "known_element_hits.tsv", sep="\t").head(8)
    image, draw = board(
        "实际完成：B0 已知元件扫描",
        "共识 PWM：TTGACA（-35 box） / TATAAT（-10 box），支持正负链",
    )
    draw.text((60, 150), "demo 扫描汇总", fill=NAVY, font=F_HEAD)
    summary_rows = [[row.motif_name, str(row.sequence_count), str(row.hit_count)] for row in summary.itertuples()]
    table(draw, 60, 188, [("已知元件", 250), ("命中序列数", 180), ("命中次数", 160)], summary_rows)
    draw.rounded_rectangle((720, 145, 1535, 278), radius=14, fill=LIGHT_ORANGE, outline="#F1C9A8", width=2)
    draw.text((750, 166), "这一步实际验证了什么？", fill=ORANGE, font=F_HEAD)
    draw.text((750, 208), "程序可以读取清洗后的 FASTA/TSV，扫描已知短片段，输出位置、方向和分数。", fill=NAVY, font=F_BODY)
    draw.text((750, 238), "注意：当前 p_value_proxy 只是排序代理，不是正式 FIMO 显著性。", fill=RED, font=F_SMALL)

    draw.text((60, 330), "实际命中记录示例", fill=NAVY, font=F_HEAD)
    rows = []
    for row in hits.itertuples():
        rows.append([row.motif_name, row.sequence_id, f"{row.start}–{row.end}", row.strand, f"{row.score:.2f}", row.matched_sequence])
    table(draw, 60, 368, [("元件", 190), ("序列 ID", 170), ("位置", 120), ("方向", 80), ("分数", 100), ("实际匹配片段", 220)], rows, 38)
    draw.text((60, 820), "说明：命中数量来自 demo，不用于宣称真实元件富集或生物学显著。", fill=RED, font=F_SMALL)
    image.save(OUTPUT / "M1证据_B0已知元件扫描.png")


def make_repro_board() -> None:
    manifest = json.loads((DEMO / "m1_manifest.json").read_text(encoding="utf-8"))
    environment = manifest["environment"]
    report = manifest["quality_report"]
    image, draw = board(
        "实际完成：可复现性与环境检查",
        "每次运行保存输入哈希、配置哈希、随机种子和工具状态",
    )
    card(draw, (55, 145, 375, 245), "Python", environment["python"]["version"], "CPython", LIGHT_BLUE)
    card(draw, (405, 145, 725, 245), "随机种子", str(report["seed"]), "固定运行条件", LIGHT_GREEN)
    card(draw, (755, 145, 1075, 245), "输入 SHA256", report["input_sha256"][:12] + "…", "写入 quality_report", LIGHT_GRAY)
    card(draw, (1105, 145, 1425, 245), "配置 SHA256", report["config_sha256"][:12] + "…", "写入 quality_report", LIGHT_GRAY)

    draw.text((60, 300), "当前环境实测", fill=NAVY, font=F_HEAD)
    tool_rows = []
    for tool, info in environment["meme_suite"]["tools"].items():
        tool_rows.append([tool, "未检测到" if not info["available"] else "已检测", info["version"] or "—"])
    table(draw, 60, 338, [("工具", 180), ("状态", 180), ("版本", 430)], tool_rows)
    draw.rounded_rectangle((900, 320, 1535, 565), radius=16, fill=LIGHT_RED, outline="#E4B5B5", width=2)
    draw.text((935, 348), "当前真实阻塞项", fill=RED, font=F_HEAD)
    warning = "MEME Suite 主线工具尚未安装：\nmeme / streme / dreme / fimo / tomtom\n\n所以目前只能确认 Python 数据管线和 B0\n基线可运行，不能宣称正式 motif 实验已完成。"
    y = 400
    for line in warning.splitlines():
        draw.text((935, y), line, fill=NAVY, font=F_BODY)
        y += 30
    draw.text((60, 600), "关键 Python 包版本", fill=NAVY, font=F_HEAD)
    packages = environment["packages"]
    package_names = ["numpy", "pandas", "pyyaml", "scikit-learn", "scipy", "pytest"]
    package_rows = [[name, packages.get(name) or "未安装"] for name in package_names]
    table(draw, 60, 638, [("包", 220), ("版本", 260)], package_rows, 27)
    draw.text((60, 850), "这张图展示的是实验准备状态，不是模型效果图。", fill=RED, font=F_SMALL)
    image.save(OUTPUT / "M1证据_可复现环境检查.png")


if __name__ == "__main__":
    make_clean_board()
    make_b0_board()
    make_repro_board()
    print("generated=3")
