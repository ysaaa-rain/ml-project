# PR01-02 启动子调控模式发现

本仓库用于《机器学习综合实践》PR01-02 项目：从多物种启动子序列中发现调控 motif 和间距模式，与已知元件对比，并验证结果在不同物种和数据集上的一致性。

## 项目边界

本项目是纯模式发现任务，当前范围只包括：

- 多物种启动子序列与元数据构建；
- \`-10/-35 box\`、UP element 等已知元件的分布分析；
- 全体与分组 de novo motif 发现；
- motif 定位、元件间距和排列顺序分析；
- 跨物种、跨数据集的一致性验证；
- Motif 与启动子强度的关联分析；
- Motif 组合共现模式分析。

启动子分类、启动子强度预测、启动子生成/优化不属于本仓库范围。

## 当前状态

当前已完成项目初始化、课程要求固化、候选基线锁定、Git 子模块接入，以及可复现的数据清洗和 M2 EDA 基础管线。仓库仍未提交真实外部数据，当前 demo 只用于验证代码链路，不作为生物学实验结论。详细状态见 [\`docs/进度日志.md\`](docs/进度日志.md)。

## 快速开始

\`\`\`bash
git clone --recurse-submodules https://github.com/ysaaa-rain/ml-project.git
cd ml-project
python -m pip install -r requirements.txt
\`\`\`

完整实验命令将在 M1 数据构建和 M2 预处理完成后补充到 [\`docs/运行说明.md\`](docs/运行说明.md)，并与当前代码同步维护。

## M1/M2 最小运行示例

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python -m experiments.run_m1 --config configs/m1_demo.yaml
```

该命令会生成清洗后的 TSV/FASTA、序列打乱背景、质量报告、分组统计、位置碱基组成和 PNG 图形。将配置中的 `input` 替换为经来源审计的真实 metadata 表后，可复用同一流程。

问题建模的正式定义、符号、判定规则和实验矩阵见 [`docs/03_问题建模.md`](docs/03_问题建模.md)。

第一阶段课堂汇报统一使用 [`reports/M1建模阶段汇报材料_20260911.md`](reports/M1建模阶段汇报材料_20260911.md)，每一页同时包含页面内容、讲解提示、展示建议和备查问题，不另列演讲稿。

后续真实实验按 [`docs/04_两周实验推进计划.md`](docs/04_两周实验推进计划.md) 逐节点推进。运行 `python -m experiments.preflight --output tmp/preflight.json` 可检查当前环境是否具备 MEME Suite 主线工具。

## 目录约定

\`\`\`text
data/            原始、处理中和标准化后的数据（默认不提交大文件）
preprocessing/   数据清洗、标准化、分组和背景构造
models/          motif 表示、相似性和组合分析所需的模型/算法封装
experiments/     可复现实验入口和参数配置
results/         表格、图形、日志和报告中间产物
baselines/       锁定的公开基线依赖
docs/            需求、研究设计、数据说明、实验记录和进度文档
reports/         阶段报告和最终报告
\`\`\`

## 协作与版本记录

每次代码、配置或实验结论变更都必须有清晰的 Git 提交信息，并同步更新相关文档。提交信息应说明“做了什么”和“为什么做”，例如：

\`\`\`text
docs: 固化 PR01-02 研究问题与验收矩阵
data: 增加多物种启动子数据字段规范
feat: 增加按 sigma 因子分组的 motif 运行入口
exp: 记录 MEME 与 STREME 参数对比结果
\`\`\`

不提交原始大数据、临时文件、个人密钥和无关对话内容。
 
"������Ϊ���Դ��뽡���Ⱥ����֤��¼��" 
