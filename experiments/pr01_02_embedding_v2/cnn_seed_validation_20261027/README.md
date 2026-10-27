# CNN边界处理：多初始化独立复核（2026-10-27）

已完成：[结果报告](REPORT.md)。30%退化下RAW/SHIFT F1约0.506/0.539，5数据2胜1负2平；4个新算法种子中1个平均差为正，未通过预设进一步合成研究门槛。30次算法运行中26次F1不变、3次提高、1次下降；同批OH48约0.510但监督条件不同。算法种子同时控制权重初始化与训练批次，真实Holdout未使用。

依据[上一轮结果](../cnn_boundary_20261027/REPORT.md)：30%退化时边界处理的收益集中于两个初始化中的一个。新方案冻结在 `protocol.json`，优先验证收益是否跨初始化重现，不更换模型或按旧Test改算法。

强模式、30%退化、混合反向链15%退化各5个**全新数据种子**，每组使用6个模型/聚类随机种子（上一轮2个初始化加4个全新初始化）。RAW/SHIFT共享同一个CNN、同一批边界有效实例；完整One-hot48每数据/算法种子单独拟合。计划90次CNN训练、90次OH48聚类、270组PWM池和900组预算评价；RAW/SHIFT重用模型，不计额外训练。单向CNN上轮消融结果不一致，本次只复核双向共享方法。

主要终点：30%退化、匹配候选数的SHIFT−RAW定位F1；先对每个数据种子平均6次初始化，再比较5个数据种子。同时报告4个**新初始化**各自跨5数据的配对平均差，不能把共享数据的初始化当独立样本。比较范围、预设推进门槛和误报约束见协议。全部Discovery拟合、PWM、平移、排序和预算在任何新Development/Test评分前锁定；真实Holdout关闭，未进行真实数据拟合或SAE。

沿用上一轮实现：同一个CNN算法种子同时决定模型初始权重和训练批次顺序；因此本轮观察的是这两项共同变化的稳定性，不能将差异只归因于初始权重。OH48使用相同编号各自作为聚类随机种子。

```powershell
python -m pytest experiments/pr01_02_embedding_v2/cnn_seed_validation_20261027/test_pipeline.py experiments/pr01_02_embedding_v2/cnn_boundary_20261027/test_pipeline.py -q
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.pipeline fit
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.pipeline score
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.pipeline qa
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.qa_cross_version
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.report --figure-only
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.pipeline verify
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.report
```

拟合拒绝覆盖已有artifacts。序列、实例、模型与聚类对象留本地；公开冻结协议、代码、分数据/初始化统计、图表及哈希。状态以 `summary/manifest.json` 为准。
