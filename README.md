# PR01-02：多物种启动子调控模式发现

本项目围绕课程题目 PR01-02，研究细菌启动子 DNA 序列中反复出现的短调控模式（motif），分析 motif 与已知元件的关系、不同 σ 因子组的差异，以及规律能否跨物种、跨数据集重现。

## 研究问题

- **RQ1：** 从启动子中发现的 motif 与已知 −10/−35 box、UP 元件等有多大重合？
- **RQ2：** 不同 σ 因子类型的启动子在 motif 组成、出现频率、间距和排列上有什么差异？
- **RQ3：** 发现的规律能否在其他物种或独立数据集中重现？哪些是共同规律，哪些只出现在特定数据中？

## 数据角色

| 数据 | 在项目中的用途 | 当前边界 |
| --- | --- | --- |
| 天津大学课程云盘 reg_and_gen 六物种 | 主数据；正类启动子分别做 MEME/STREME、FIMO、位置及物种比较；负类作差异富集对照 | CSV 无 σ 标签或链方向列；对 TSS 方向的位置解释须先审计 |
| RegulonDB E. coli K-12 | σ 因子分组分析；与云盘 E. coli 去除序列重叠后的跨数据集比较 | 原始 3,807 条 motif 与六个 σ 组的 MEME 结果已存在，当前作为辅助数据使用 |
| DBTBS release 4.1 B. subtilis | 有余力时进行增强验证 | 当前坐标与方向解析仍需复核 |
| 云盘 strenth | 实验 5 的候选强度数据 | 来源、测量单位、样本映射和许可通过审计后再用 |

云盘 reg_and_gen 共 18,370 条记录、9,724 条正类启动子，覆盖 B. subtilis、A. baumannii、Bradyrhizobium、C. diphtheriae、E. coli 和 Staphylococcus。每个物种只使用一份 Dataset.csv 建立输入；train/dev/test 是它的拆分，不能当成额外独立数据。上游项目说明窗口为 81 bp、相对 TSS −60 到 +20 bp；原始 CSV 没有证明负链是否统一到转录方向。

原始云盘文件保存在本地 data/raw/tju_pan_promoter/，由 SHA256 清单校验并被 Git 忽略。正式结果引用到的表、图和运行 manifest 应进入仓库；原始数据、完整序列、运行缓存和大型中间文件不进入普通 Git。

## 当前进度

- **M1：** PR01-02 的目标、问题、数据分工与审计规则已固化。
- **M2：** RegulonDB 辅助数据和 TJU Pan 六物种的描述性 EDA 均已生成；云盘方向未核实，所以其位置统计只按序列内坐标描述。主数据报告见[reports/M2_TJUPan阶段报告.md](reports/M2_TJUPan阶段报告.md)。
- **M3：** TJU Pan 六物种全量 MEME、60 个 motif 矩阵/summary、统一 logo 和 240 条 −10/−35/UP Tomtom 初步比较已完成；这是矩阵相似性探索，不是 motif 基因组位置重叠的证明。完整证据见[六物种 M3 阶段报告](reports/M3_TJUPan阶段报告.md)。
- **M4：** 六物种第一轮 de novo FIMO 已完成：每 motif q≤0.05 共 353 个位点；五个“物种×motif”正负覆盖差异在 60 项 BH 后 q≤0.05，均属同批发现数据上的探索性线索。另已扫描四个理想化大肠杆菌 −10/−35/UP 参考矩阵；site q≤0.05 下没有参考位点通过，因此当前阈值下不能据此评估坐标重叠，不能把“零重叠”解释成这些元件不存在。motif pair 首轮共现/间距描述已完成，但 540 项检验无一经 BH 显著，且没有位置随机化；方向坐标、敏感性分析、跨物种和跨数据集复核仍待完成。详见[M4 FIMO 阶段报告](reports/M4_TJUPan_FIMO阶段报告.md)。
- **扩展实验：** 小组选定实验 6（motif 组合）和实验 5（motif 与启动子强度关联），均纳入计划；之后推进调控语法和自定义 RQ4。实验 7 为有余力时的扩展。

当前每个实验的准确状态、限制和下一步见 docs/当前进度与协作交接.md。

## 环境与运行

需要 Python 3.11+、requirements.txt 中的依赖以及 MEME Suite 5.5.9。MEME Suite 官方源码构建记录见 docs/MEME环境与Smoke测试_20260930.md。

~~~bash
python -m pip install -r requirements.txt
python -m experiments.prepare_tjupan_motif_inputs
python -m experiments.run_tjupan_meme
python -m experiments.render_meme_logos --meme-root tmp/meme_runs/tjupan_promloop_20261003 --output-dir results/motif/tjupan_m3_20261003/logos --species bacillus_subtilis baumannii bradyrhizobium diphtheria escherichia_coli staphylococcus
python -m experiments.summarize_tjupan_meme --meme-root tmp/meme_runs/tjupan_promloop_20261003 --run-manifest results/motif/tjupan_m3_20261003/run_manifest.json --output-dir results/motif/tjupan_m3_20261003/summary
python -m experiments.run_tjupan_fimo
python -m experiments.summarize_tjupan_fimo
python -m experiments.run_tjupan_reference_fimo
python -m experiments.summarize_tjupan_reference_overlap
python -m experiments.analyze_tjupan_motif_pairs
~~~

输入准备命令会逐条核验六个物种的正类文件、序列长度、碱基字符、拆分关系和 SHA256，并生成本地 FASTA。MEME 命令按物种逐个搜索，使用 positive_samples.csv 作 primary、Dataset.csv 的 label=0 作真实负类对照。参数和软件版本写入运行 manifest。当前正类共 9,724 条，四条不等长负类已排除。复现说明见 docs/运行说明.md。

## 项目目录

~~~text
data/            本地原始与处理中数据（不提交）
preprocessing/   输入验证、清洗、来源适配和序列窗口处理
models/          可解释 PWM 基线
experiments/     数据审计、motif 发现和统计分析入口
results/         审核过的结果表、图和 manifest
baselines/       已知启动子元件参考矩阵
configs/         固定的实验参数
docs/            任务定义、数据、方法、复现与进度记录
reports/         阶段报告和最终报告
~~~

## 协作与记录

每项正式分析都要记录输入来源、样本数、软件版本、完整参数、随机种子、命令、输出位置、SHA256、解释边界和下一步。Git 提交应表达实际改动，例如：feat: prepare six-species promoter motif inputs。数据与结果许可未确认前，不推送原始记录或可还原原始序列的文件。
