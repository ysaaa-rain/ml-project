# PR01-02 启动子调控模式发现

本仓库用于《机器学习综合实践》PR01-02 项目：以 DNA 预训练大模型嵌入为主要方法，从多物种启动子序列中发现调控结构，同时保留一套紧凑的课程 motif/间距线进行校验。

## 项目边界

本项目是纯模式发现任务，当前范围只包括：

- 多物种启动子序列与元数据构建；
- 冻结 DNA 预训练模型的序列级和局部表示提取；
- 嵌入空间中的聚类、近邻检索和分组差异分析；
- 从大模型表示中定位重要序列区域，并与已知 `-10/-35 box` 对照；
- 跨物种、跨数据集的一致性验证；
- 一套紧凑的 PWM + MEME/DREME/STREME + FIMO 课程线，用于全序列/分组 motif、间距分析和大模型结果的交叉验证；
- 资源允许时进行轻量微调或 motif/embedding 融合对照。

启动子分类、启动子强度预测、启动子生成/优化不属于本仓库范围。

## 当前状态

当前已完成真实来源下载、字段适配、第一轮清洗与 EDA、L0-L2 大模型运行骨架、自动化测试，以及多轮真实数据向量分析。已按老师要求落地 RegulonDB（E. coli）和 DBTBS release 4.1（B. subtilis）：合并 4,740 条，清洗后保留 4,489 条。2026-09-22 在来源原生窗口上完成 64 条 pilot 和全部 4,489 条 embedding；2026-09-23 完成 28 条 RegulonDB TSS 对齐代表序列的遮挡探索；2026-09-24 对 3,807 条 TSS 对齐序列完成分组向量诊断。当前 sigma 与聚类、近邻的一致性接近随机基线，尚无 sigma 特异结构证据。由于 DBTBS 仍缺少可确认的统一 TSS/链方向，现有结果不能作为最终跨来源生物学结论。详细状态见 [`docs/GateA复核_20260922.md`](docs/GateA复核_20260922.md)、[`docs/分组向量诊断_20260924.md`](docs/分组向量诊断_20260924.md)、[`docs/数据下载记录_20260916.md`](docs/数据下载记录_20260916.md) 和 [`docs/进度日志.md`](docs/进度日志.md)。

## 快速开始

```bash
git clone --recurse-submodules https://github.com/ysaaa-rain/ml-project.git
cd ml-project
python -m pip install -r requirements.txt
```

完整实验命令和各类产物说明见 [`docs/运行说明.md`](docs/运行说明.md)，并与当前代码同步维护。

## M1/M2 最小运行示例

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python -m experiments.run_m1 --config configs/m1_demo.yaml
python -m experiments.run_embedding --config configs/l0_embedding_hf_demo.yaml
```

前两条入口会生成清洗后的 TSV/FASTA、随机背景、质量报告和 EDA；第三条会用真实预训练权重生成 embedding、聚类、近邻和遮挡窗口结果。demo 仅用于工程验证。

问题建模的正式定义、符号、判定规则和实验矩阵见 [`docs/03_问题建模.md`](docs/03_问题建模.md)。

真实实验按 [`docs/04_两周实验推进计划.md`](docs/04_两周实验推进计划.md) 逐节点推进。`preflight` 已确认本机可运行 PyTorch、Transformers 和 Accelerate；MEME Suite 可执行工具当前未发现，课程 motif/FIMO 线仍不能写成已完成。模型许可证和数据再利用边界也仍需按实验记录继续确认。

DNA 大模型嵌入仍是主线：已锁定 50M 多物种 Nucleotide Transformer，并完成来源原生窗口的聚类、近邻和偏差诊断；下一步是固定规则的局部遮挡、TSS 可比子集/独立来源复核，以及课程要求的全序列/分组 motif、FIMO 和间距分析。传统线不扩展成第二条算法主线。具体设计见 [`docs/05_DNA大模型嵌入主线方案.md`](docs/05_DNA大模型嵌入主线方案.md)。

## 目录约定

```text
data/            原始、处理中和标准化后的数据（默认不提交大文件）
preprocessing/   数据清洗、标准化、分组和背景构造
models/          DNA embedding 后端、无监督分析和小型 PWM 基线
experiments/     可复现实验入口和参数配置
results/         表格、图形、日志和报告中间产物
baselines/       锁定的公开基线依赖
docs/            需求、研究设计、数据说明、实验记录和进度文档
reports/         阶段报告和最终报告
```

## 协作与版本记录

每次代码、配置或实验结论变更都必须有清晰的 Git 提交信息，并同步更新相关文档。提交信息应说明“做了什么”和“为什么做”，例如：

```text
docs: 固化 PR01-02 研究问题与验收矩阵
data: 增加多物种启动子数据字段规范
feat: 增加固定版本 DNA 模型的 embedding 入口
analysis: 记录跨来源 embedding 稳定性与偏差诊断
```

不提交原始大数据、临时文件、个人密钥和无关对话内容。

目录整理记录：2026-09-21 已将项目源文件集中到仓库目录，详细实验状态以 [docs/进度日志.md](docs/进度日志.md) 为准。
