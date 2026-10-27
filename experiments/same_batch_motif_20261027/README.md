# 同批 Motif 候选比较

2026-10-27 完成同批比较：[合成数据报告](../../reports/same_batch_motif_equal_20261027/REPORT.md)｜[六物种真实数据报告](../../reports/same_batch_real_equal_20261027/REPORT.md)。每个比较单元使用同一输入划分、同一评价及同一候选数；两份报告覆盖的方法不同，不能把两份结果合并成一个总排名。

合成：复用既有BERT、One-hot候选，在同一批15组数据上增加45次CNN训练；每组和算法种子按五种表示中共同可用候选数取最多5个，重新运行参考辅助片段定位。测试集先前已观察，本轮是事后探索；CNN使用序列级监督。

真实：复用STREME既有主候选、One-hot和ProkBERT各3个种子，按同一Discovery天然对照校准阈值、同一Discovery规则排序，再对同一Development集进行双向扫描。每个物种使用三方法共同可用的候选数。STREME候选早已筛选，Embedding候选来自各运行的PWM池，产生条件不同；CNN和MEME不在该比较中。未使用真实Holdout。

代码、[协议](protocol.json)、拟合哈希锁及输出明细保留；训练产生的序列相关文件、模型和中间实例只保留在本地忽略的 `artifacts/`。本次拟合已完成，`fit` 对已有目录拒绝覆盖。评估结果和报告不能在已有输出目录上直接重建，避免覆盖既有结果。重算时应在隔离副本下使用：

```powershell
python -m experiments.same_batch_motif_20261027.pipeline fit
python -m experiments.same_batch_motif_20261027.equal_budget
python -m experiments.same_batch_motif_20261027.report
python -m experiments.same_batch_motif_20261027.real_compare
python -m experiments.same_batch_motif_20261027.real_report
```
