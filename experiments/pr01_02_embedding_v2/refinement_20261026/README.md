# 候选片段方向/偏移对齐与背景拒绝（2026-10-26）

本轮是上一轮之后的新版本小试，冻结协议见 `protocol.json`。4类合成条件（强模式、15%/30%退化、混合反向链15%退化）×3个新数据种子×2个聚类种子，One-hot24与Frozen ProkBERT局部有序128维均运行相同5种候选后处理。实际48次聚类，240组发现输出；后处理复用聚类，不重复计作新拟合。

- RAW：每簇每序列最近窗口直接构建PWM。
- RC_ONLY：只允许候选片段反向互补。
- RC_SHIFT：允许反向互补和相对原窗口±2bp偏移。
- BG_ONLY：用RAW的Discovery空白最大分95百分位拒绝低分实例，至少20条重建PWM。
- RC_SHIFT_BG：对RC_SHIFT应用同样的拒绝规则。

对齐模板只由实际Discovery候选片段学习，初始种子是距簇中心最近的窗口，固定3轮，无参考PWM/植入真值参与。偏移始终相对原窗口，不累计。这里只调整候选片段，不翻转整条TSS标准化序列。BG为**聚类后的实例拒绝**，不宣称更换了聚类算法。

所有输出先锁定，再统一评分Development/Test。主指标沿用全宽参考cosine≥.8、方向正确、±2bp定位F1；另报告所有候选的单PWM和整个候选集的空白误报率，不受参考门槛筛除。先在每个数据种子内平均2个聚类运行，再汇总3组数据；仍是小试，不作总体优劣/稳定性或真实功能结论。

```powershell
python -m pytest experiments/pr01_02_embedding_v2/refinement_20261026/test_pipeline.py -q
python -m experiments.pr01_02_embedding_v2.refinement_20261026.pipeline fit
python -m experiments.pr01_02_embedding_v2.refinement_20261026.pipeline score
python -m experiments.pr01_02_embedding_v2.refinement_20261026.pipeline verify
```

拟合拒绝覆盖已有artifacts。模型/特征/人工序列/逐实例留本地；公开源码、协议、聚合统计、图表和哈希。旧结果与冻结v2只读，真实Holdout不启用。最终结果见完成后生成的`REPORT.md`。
