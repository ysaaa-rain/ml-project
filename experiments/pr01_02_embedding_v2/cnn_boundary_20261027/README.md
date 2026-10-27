# CNN参考无关窗口边界与方向消融（2026-10-27）

依据[上一轮报告](../cnn_features_20261026/REPORT.md)的评分后诊断，检验候选窗口存在1–2bp偏移时，是否可只依据Discovery序列使同一filter的实例整体平移；另比较双向共享和单向训练/提取。冻结协议在 `protocol.json`，新数据种子、独立版本，不改上一轮模型、结果或真实数据。

每filter先按固定激活阈值选择Discovery正序列，筛去任一±2bp候选窗口会越界的实例；RAW和SHIFT共享这批序列。SHIFT只比较5个统一偏移的PWM相对Discovery空白背景的信息量，按固定平局规则取一个；没有输入参考motif或植入位置。此处不是逐实例自由对齐，信息量选择可能偏爱随机波动，必须比较新Test。

3条件×5数据种子×2算法种子，预定60次CNN训练、30次OH48聚类、150组PWM池。每个模型固定末轮，同一组模型RAW和SHIFT仅重建/选择候选，不计为新训练；方向消融共享初始化与批次流，但同时改变训练及提取方式。FULL、CAP5、匹配边界前后及匹配OH48共420组评价集合，全部在Dev/Test评分前锁定。原始非内缩CNN并未复跑，边界敏感序列剔除规模在结果中报告。

主要终点是30%退化、匹配候选预算下SHIFT−RAW的参考辅助全12bp定位F1；先在各数据种子平均两次算法运行。全部候选的完整结果、误报、成本和零候选一起报告；次要方向对照和OH48对照不是相同监督量的纯表征对比。真实Holdout关闭，未进行真实训练、SAE或功能验证。

```powershell
python -m pytest experiments/pr01_02_embedding_v2/cnn_boundary_20261027/test_pipeline.py experiments/pr01_02_embedding_v2/cnn_features_20261026/test_pipeline.py -q
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.pipeline fit
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.pipeline score
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.pipeline verify
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.report
```

`fit`拒绝覆盖既有artifacts。人工序列、实例和拟合对象留本地；公开代码、协议、训练曲线、摘要、诊断与哈希。状态以 `summary/manifest.json` 为准。
