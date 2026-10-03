# MEME Suite 环境与 smoke 测试记录

更新日期：2026-09-30
目的：记录正式 M3 motif 发现前的工具安装、版本、回归测试、最小运行检查和当前正式运行状态。合成 smoke 序列只验证程序可执行，不构成项目实验数据或生物学证据。

## 安装来源与版本

- 软件：MEME Suite 5.5.9 官方源码包。
- 官方下载页：[MEME Suite releases](https://meme-suite.org/meme/doc/download.html)。源码包 URL：`https://meme-suite.org/meme/meme-software/5.5.9/meme-5.5.9.tar.gz`。
- 下载文件 SHA256：`0406fb7b1dc27f6aab3d6d3a29ecdf617bbdd946690cfa97ca2129bc210bfd12`。
- 本机源码、编译目录和安装前缀均位于 `tmp/`，不会进入普通 Git 提交；安装目录为 `tmp/meme-suite-5.5.9/`。
- 按官方 macOS/Linux 源码安装流程，配置时启用随包构建的 libxml2/libxslt；`make` 成功，官方 `make test` 汇总为 **192 PASS、0 FAIL、0 ERROR**，`make install` 成功。
- `meme`、`dreme`、`streme`、`fimo` 均报告版本 **5.5.9**。
- 官方安装说明见[安装指南](https://meme-suite.org/meme/doc/install.html)。

本机复现流程（示意，安装前缀应使用项目本地 `tmp/`，不覆盖个人目录）：

```bash
tar -xzf tmp/meme-suite-build/meme-5.5.9.tar.gz -C tmp/meme-suite-build
cd tmp/meme-suite-build/meme-5.5.9
./configure --prefix=../../meme-suite-5.5.9 \
  --enable-build-libxml2 --enable-build-libxslt
make -j4
make test
make install
```

## 合成数据 smoke test

使用 6 条人为构造的短 DNA 序列，输入与输出均留在本机被忽略的 `tmp/meme-suite-build/`。检查结果如下：

| 工具 | 检查 | 结果 | 解释边界 |
| --- | --- | --- | --- |
| MEME | classic/zoops 小样本发现 | 退出码 0；生成文本、XML、HTML 和 motif logo 文件 | 合成集上出现的任意 motif 不是项目发现 |
| DREME | 默认 `-e 0.05` | 退出码 0；无达到阈值的 motif | 小样本无显著 motif 是允许结果，不代表程序失败 |
| DREME | 仅用于验证结果输出的宽松 `-e 1` | 退出码 0；生成 1 个 motif | 此宽阈值只供 smoke 测试，不能用于正式结论 |
| FIMO | 对合成 MEME motif 文件扫描 | 退出码 0；生成 TSV、XML、GFF、HTML 等文件 | 只确认定位工具能运行，不是正式定位分析 |

MEME 本机没有 Ghostscript 或 ImageMagick 的 `convert`，所以 smoke 输出中的 EPS logo 没有自动转 PNG。MEME Suite 官方说明可在依赖不可用时通过 HTML 浏览器界面查看图形；正式展示图需要另行解决可复现的图像导出方式，不能把缺少 PNG 误报为 motif 搜索失败。

DREME 官方文档当前将其标记为 deprecated，并建议考虑 STREME；本项目仍保留 DREME 环境 smoke 记录。若 M3 将 DREME 作为正式互证方法，应在结果报告中说明工具维护状态并固定版本与参数；不能把 DREME 与 MEME 的不同搜索目标下的 E-value 直接当作同一种统计量比较。[DREME 官方说明](https://meme-suite.org/meme/doc/dreme.html)

## RegulonDB 辅助 M3 运行和结果（2026-09-30）

全体输入重建目录：`tmp/motif_inputs_20260930/`。输入序列 SHA256 与历史登记相符：

| 数据 | 条数 | SHA256 |
| --- | ---: | --- |
| 全体 RegulonDB TSS 对齐 discovery 序列 | 3,807 | `8514a6521df4a708d621aa06ca92573f0b561d0d091987360daaf82ec06df874` |
| 一一匹配的 N1 二核苷酸打乱对照 | 3,807 | `757f25fd64e6a29cac08a9a1f5c220f3532c6546d8830cb6fb27ad65f6961c78` |

FASTA 窗口为 `[-60,+20]`，共 81 bp；原始负链已反向互补到统一转录方向。每条 N1 对照都与对应真实序列不同，且逐条保持二核苷酸计数。未知 sigma 的 640 条序列留在全体发现中，不组成额外 sigma 分组。

### 运行策略调整记录

曾用 `-allw` 对 6–20 bp 内每个宽度穷举起始点。该探索运行约 4.6 小时后仍在搜索第 9 个 motif；虽然已写出 8 个暂存候选，但进程以退出码 130 中断，未正常结束。因此这些部分输出仅用于确认计算和文件格式正常，**不作正式结果、不引用 motif 数、E-value 或生物学解释**，目录保留在 `tmp/meme_runs/regulondb_tss_20260930/all_allw_interrupted_20260930/`。

为让全体及六组都能在同一资源约束下完成，正式运行保留同一 6–20 bp 宽度范围，改用 MEME 默认启发式 seed-width 策略（不传 `-allw`）。改变只涉及宽度起始点的搜索策略；版本、主/对照 FASTA、差异富集目标、`zoops`、双链、motif 上限和随机种子保持一致。启发式策略会评估代表性 seed widths，不保证遍历 6–20 的每一个整数宽度；这个取舍须随结果报告。

### RegulonDB 全体 MEME 搜索命令

```bash
tmp/meme-suite-5.5.9/bin/meme \
  tmp/motif_inputs_20260930/all.fasta \
  -dna -mod zoops -objfun de \
  -neg tmp/motif_inputs_20260930/all_n1.fasta \
  -nmotifs 10 -minw 6 -maxw 20 -revcomp \
  -seed 20260930 -brief 50 \
  -oc tmp/meme_runs/regulondb_tss_20260930/all
```

选择 `-objfun de` 是因为本轮有一一匹配、保持二核苷酸频率的 N1 对照；MEME 的差异富集目标用于比较 primary 与 control 序列中的 motif 位点富集。`-revcomp` 允许搜索反向互补方向；`zoops` 表示每条序列至多一个该 motif 位点，适合先做启动子核心模式发现，但结果解释仍需结合序列覆盖和位置分布。MEME 官方说明指出，差异富集模式使用 primary hold-out 子集评估 motif，因此其 E-value 的解释依赖相应搜索设计和输入对照。[MEME 官方命令说明](https://meme-suite.org/meme/doc/meme.html)

#### 全体 3,807 条运行结果（初步质检）

- 正式命令退出码：0；MEME XML 使用 `xmllint --noout` 检查通过。XML 中 primary/control 均为 3,807 条、各 308,367 bp；共 10 个 motif 矩阵；停止原因是达到 `-nmotifs 10` 上限。MEME 记录的运行时间约 2,557.54 秒（42 分 38 秒）。
- 结果文件 SHA256：`meme.txt`=`6ed1cbb08090e3f29f2f325059c3180243865735d535d20024daa577199e1957`；`meme.xml`=`707b40e1004383b41a97f2cc456701c6b33f8fd79e78798cdb80a32aa1da02b0`；`meme.html`=`e61b665392ecbd701c38cd319db405f747e51cd28f3c735afb1a6a759c013d0b`。
- MEME XML 显示 `searchsize=100000`（本次未显式设置，是 5.5.9 默认值），作为候选起始点搜索的字符数上限；primary/control 全长数据各 308,367 bp、全 3,807 条都仍在整体训练集/差异富集设计中。

| 序号 | MEME consensus | 宽度 | sites | E-value | 当前判读 |
| ---: | --- | ---: | ---: | ---: | --- |
| 1 | `CYGGCA` | 6 | 1,487 | `4.6e-18` | 差异富集候选 |
| 2 | `TTTTTYRNY` | 9 | 2,308 | `3.0e-8` | 差异富集候选；连续 T，须重点查低复杂度/组成偏差 |
| 3 | `CTGGTG` | 6 | 152 | `3.6e-1` | 未达到 0.05 |
| 4 | `YGCCGC` | 6 | 623 | `7.9e-2` | 未达到 0.05 |
| 5 | `BCCTGH` | 6 | 1,078 | `1.9e-5` | 差异富集候选 |
| 6 | `ATTSAYTAHATSAATRTATW` | 20 | 20 | `7.1e-1` | 未达到 0.05 |
| 7 | `TAATGCGC` | 8 | 37 | `1.1` | 未达到 0.05 |
| 8 | `ATTWATWWWTAATWTAT` | 17 | 31 | `1.4` | 未达到 0.05 |
| 9 | `GTTKACCCADYTTGSTYTCR` | 20 | 17 | `1.0e-1` | 未达到 0.05 |
| 10 | `ATYATTTGATSCWWAYSTT` | 19 | 31 | `1.0e-1` | 未达到 0.05 |

这里的“差异富集候选”仅指 MEME differential-enrichment 统计中 E-value≤0.05；不等于已验证的调控元件。特别是 motif 2 的连续 T 需要与局部碱基组成/低复杂度核对，并在分组及后续定位实验中复核。Motif 的位置、与 `-10/-35/UP` 重叠、唯一序列覆盖率和跨数据集重复性都还未由本轮 MEME 发现解决。由于 `-brief 50`，MEME XML 的 `contributing_sites` 为空；sites 数不是独立重算的唯一序列覆盖数，覆盖率留待 FIMO 阶段计算。

### 六个 sigma 分组搜索结果

六组均已按同一有效搜索设置完成。分组命令显式传入 `-searchsize 100000`（等于本版本默认值），并传 `-nostatus` 仅关闭进度输出；其余生物学分析参数与全体运行相同。每组 primary/control 条数相同，XML 已通过 `xmllint --noout`，每组输出 10 个 motif 并达到 `-nmotifs 10` 上限。

```bash
for group in Sigma70 Sigma24 Sigma32 Sigma38 Sigma28 Sigma54; do
  tmp/meme-suite-5.5.9/bin/meme "tmp/motif_inputs_20260930/$group.fasta" \
    -dna -mod zoops -objfun de \
    -neg "tmp/motif_inputs_20260930/${group}_n1.fasta" \
    -nmotifs 10 -minw 6 -maxw 20 -searchsize 100000 -revcomp \
    -seed 20260930 -brief 50 -nostatus \
    -oc "tmp/meme_runs/regulondb_tss_20260930/$group"
done
```

| 分组 | primary/control 条数 | 输出 motif 数 | E-value≤0.05 数 | 候选 consensus（sites；E-value） |
| --- | ---: | ---: | ---: | --- |
| Sigma70 | 1,942 / 1,942 | 10 | 5 | `WTDHTGCBSRTTTTTTBHYY`（374；3.7e-7）；`GCRGCR`（651；2.9e-4）；`GTGAYYTTNWTCAC`（48；6.7e-3）；`SCCTGA`（161；0.011）；`AATGATAA`（49；0.041） |
| Sigma24 | 513 / 513 | 10 | 4 | `KCWGGC`（352；8.7e-7）；`GKTCAG`（72；4.1e-5）；`TATCAG`（158；0.0021）；`TTYAYB`（342；0.031） |
| Sigma32 | 300 / 300 | 10 | 2 | `CTGGCG`（145；1.3e-5）；`CTGGWWWN`（133；0.0006） |
| Sigma38 | 176 / 176 | 10 | 1 | `AGTGTAGM`（55；0.032）；另有 `GMAAADBDSVMVVA` 的 E-value=0.053，未达到 0.05。 |
| Sigma28 | 140 / 140 | 10 | 1 | `TCGKCA`（66；3.1e-5） |
| Sigma54 | 96 / 96 | 10 | 1 | `TGCCAK`（58；0.0015） |

`E-value≤0.05` 是对各自单独搜索结果使用的候选筛选线，不等于 7 次搜索上的联合多重校正。`sites` 不是不同启动子计数；`-brief 50` 没有保存逐序列完整贡献位置。Sigma70 的 20 bp T-rich 模式和全体的 T-rich 模式都需重点防范低复杂度/碱基组成解释。

### `-10/-35/UP` 初步序列相似性比较

**2026-10-03 参考矩阵审计：** 发现参考文件中 σ70 `-10` 和 `-35` 两个 motif 名称对应的矩阵原先分别编码成 `ATATTC`、`TTTAGA`，而不是文件说明中的 `TATAAT`、`TTGACA`。现已修正、添加自动共识测试，并重跑全部七组/280 条 Tomtom 配对；下表及比较 manifest 均已替换为修正后结果。旧配对和旧 p/q 值不得引用。MEME de novo 搜索本身未受影响。

基线参考文件为 `baselines/known_promoter_elements.meme`：包含 σ70 `-10 TATAAT`、σ70 `-35 TTGACA`、UP 近端 `AAAAAARNR` 与远端 `AWWWWWTTTTT` 四个理想化矩阵。σ70 consensus 参考[细菌 σ 因子综述](https://academic.oup.com/femsre/article/22/3/127/643680)，UP 子位点参考[原始 UP 元件研究](https://genesdev.cshlp.org/content/13/16/2134)。Tomtom 设置为 Euclidean distance、最小重叠 6 nt、扫描正反向参考、`-thresh 1` 输出所有配对；每组 10×4，七组合计 280 条 query-reference 行。

Tomtom 对只有 4 个 target motif 的库发出警告：该规模不足以准确估计匹配 p-value；每个 query 只有 8 个正反向 target 比较，q 值估计同样很不稳定。故 p/q 只留作审计和排序线索，不作显著性或已知元件重合结论。按 MEME E-value≤0.05 的候选筛选后，各组排序靠前的线索为：全体 `CYGGCA`→σ70 `-35`（p=0.0108，q=0.0861，重叠 6 nt；T-rich motif→UP 近端 p=0.0603，q=0.4828）；Sigma70 长 T-rich motif→UP 近端（p=0.000943，q=0.00754，重叠 9 nt，低复杂度风险显著）；Sigma24 `TATCAG`→σ70 `-35`（p=0.0542，q=0.4338）；Sigma32 `CTGGCG`→σ70 `-35`（p=0.0169，q=0.1348）；Sigma38 `AGTGTAGM`→σ70 `-10`（p=0.1113，q=0.8903）；Sigma28 `TCGKCA`→σ70 `-35`（p=0.00349，q=0.0279）；Sigma54 `TGCCAK`→σ70 `-35`（p=0.00226，q=0.0181）。除 Sigma70 T-rich motif 以外，其余配对最多只有 6 nt 重叠。这些只是四个模板内部的初步排序；尤其 σ70 核心共识不应被当成其他 sigma 的专属参考。相似性不代表 motif 已位于 `-10/-35/UP` 的基因组坐标位置；位置重合要等 FIMO。

汇总脚本会逐项核对 70 个 MEME motif 与 280 条 Tomtom 配对，并把输入 FASTA、MEME TXT/XML、Tomtom TSV、参考库、矩阵与摘要的 SHA256 写入 `results/motif/m3_20260930/comparison_manifest.json`。可在输入和输出路径仍存在时复现摘要：

```bash
./.venv/bin/python -m experiments.summarize_motif_comparison \
  --meme-root tmp/meme_runs/regulondb_tss_20260930 \
  --tomtom-root tmp/meme_runs/known_element_comparison_20261003/regulondb \
  --input-dir tmp/motif_inputs_20260930 \
  --reference baselines/known_promoter_elements.meme \
  --output-dir results/motif/m3_20260930
```

经 `meme2meme` 整理的 7 个矩阵文件及全部配对/最佳匹配摘要位于 `results/motif/m3_20260930/`。原始 MEME XML/HTML、FASTA、Tomtom HTML/XML 和日志仍放在 Git 忽略的 `tmp/`。RegulonDB 衍生结果的再分发许可尚未确认，因此虽然小型结果已本地整理，**目前没有提交或推送这些结果**。RegulonDB 辅助分析链和报告已形成；TJU Pan 主数据的 M3 仍在运行，logo 图与中期汇报材料尚待整理，FIMO 和其余 M4 分析未开始。

## 当前后续验收

RegulonDB 辅助分析的全体与六个 sigma 组搜索及参考矩阵初比已完成，详细结果见辅助分析报告。当前正在推进的是 TJU Pan 六物种主数据：先依次完成 MEME 搜索，再逐物种核验退出状态、样本量、XML、候选模式、背景和低复杂度风险。主数据各组未验收前，不将 M3 标记为完成。

原始 FASTA、序列 ID、MEME 原生目录和日志留在 Git 忽略路径。经许可审查后，只将可公开的矩阵、摘要、图与 manifest 整理进 `results/motif/`。
