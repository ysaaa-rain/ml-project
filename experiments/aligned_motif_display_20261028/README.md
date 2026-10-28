# 机器学习课程实验：横向对齐与结果展示

2026-10-28 已完成：MEME补跑两背景共20个PWM；五方法主比较、真实Logo、统一Development位置图及对齐验收表已生成并目视检查。24组集合、384次独立标量扫描、旧14组结果及24个图中候选核验通过，339个历史文件哈希不变。天然背景首次HTML导出故障已单独恢复，原退出码和日志保留。8项相关测试通过。

范围固定为枯草芽孢杆菌五方法完整案例；其他物种列覆盖清单。已有BERT、One-hot和CNN不重训，仅补齐固定Discovery上的MEME天然及Shuffle两背景。当前状态以 [核验](summary/verification.json) 与 [报告](summary/REPORT.md) 为准，脚本存在不代表已运行。

公共规则先冻结于 [protocol.json](protocol.json)：同一分区、Discovery排序与90%分位阈值、双链扫描、TOP5/ALL、天然与配对Shuffle分开；CNN的标签监督和不同方法的候选生成条件明确注明。完整五方法比较只适用于枯草芽孢杆菌。10月3日旧MEME矩阵来自全量数据，仅保留历史材料。

在隔离项目副本中依次执行（已有目录时拒绝覆盖）：

```powershell
python -B -m experiments.aligned_motif_display_20261028.pipeline audit
python -B -m experiments.aligned_motif_display_20261028.pipeline discover
python -B -m experiments.aligned_motif_display_20261028.pipeline score
python -B -m experiments.aligned_motif_display_20261028.report
python -B -m experiments.aligned_motif_display_20261028.verify
python -B -m experiments.aligned_motif_display_20261028.course_stats
```

`discover` 使用项目本地WSL环境 `mlproject-meme-20261028` 中的 `/opt/meme-5.5.9/bin/meme`。环境文件位于忽略的 `tmp/aligned_motif_runtime_20261028/`，原生输出、日志位于忽略的 `artifacts/`。MEME官方源码包SHA256为 `0406fb7b1dc27f6aab3d6d3a29ecdf617bbdd946690cfa97ca2129bc210bfd12`；运行参数、输入输出哈希和工具哈希见 `summary/meme_runs.json`。公开文件仅包括PWM矩阵、聚合表、图表和核验记录，原始序列及逐序列位点留在本地。

环境安装参照 [MEME官方安装说明](https://meme-suite.org/meme/doc/install.html) 与 [Microsoft WSL导入说明](https://learn.microsoft.com/en-us/windows/wsl/use-custom-distro)。环境准备失败不允许以旧全量PWM冒充新Discovery结果。

M1/M2统计是只读课程补充，输出在 `results/eda/m1m2_course_20261028/`；固定模板扫描仅Discovery/Development，全split只做数量、GC与长度质检。保留历史阶段报告，统一统计另附。
