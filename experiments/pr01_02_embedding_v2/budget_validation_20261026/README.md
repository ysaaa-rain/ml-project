# 背景拒绝的独立复核与候选预算匹配（2026-10-26）

冻结协议：`protocol.json`。本轮聚焦30%退化和混合反向链15%退化，保留强模式控制；每条件5个新数据种子、3个聚类种子。比较Frozen ProkBERT局部128维、One-hot PCA24和不降维One-hot48，RAW/BG_ONLY均执行。实际135次聚类、270组PWM池；预算视图复用这些池，共900组评价集合，不增加拟合。

候选仅按Discovery正类减空白命中率排序，平局按信息量、支持序列数、ID排序，禁止参考PWM/植入真值/Test参与排名。先锁定全部拟合、排名和候选集合，再评分Development/Test。

- FULL：各自全部候选，反映完整流程。
- CAP5：各自最多5个，候选不足保留并报告，不伪造占位PWM。
- MATCH_OH24 / MATCH_OH48：分别让BERT和对应基线各取相同m个，m=min(5,双方可用数量)。m=0保留为0分，不删除；必须结合FULL/CAP5，避免共同裁零掩盖失败。

主要终点预设为30%退化、BG_ONLY、MATCH_OH48的BERT减One-hot48定位F1。先在数据种子内平均3个聚类运行，再汇总5组数据。报告逐种子差值、胜负与候选集误报，配对符号翻转检验只有32种符号组合，5组数据的双侧p最小0.0625，不能宣称p<.05。门槛只决定是否值得进一步研究，不是功能/总体优势证明。

```powershell
python -m pytest experiments/pr01_02_embedding_v2/budget_validation_20261026/test_pipeline.py -q
python -m experiments.pr01_02_embedding_v2.budget_validation_20261026.pipeline fit
python -m experiments.pr01_02_embedding_v2.budget_validation_20261026.pipeline score
python -m experiments.pr01_02_embedding_v2.budget_validation_20261026.pipeline verify
```

旧实验、冻结v2与真实Holdout只读；不微调、不换编码器、不叠加对齐。完整人工序列、特征、实例和拟合对象留本地；公开代码、配置、排名/预算统计、图表与哈希。`fit`拒绝覆盖artifacts。最终结果见完成后生成的`REPORT.md`。
