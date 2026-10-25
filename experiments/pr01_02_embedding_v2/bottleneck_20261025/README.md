# Embedding 定位瓶颈小规模实验（2026-10-25）

本目录是首轮和 `next_round/` 之后的新版本，不覆盖历史运行。冻结协议见 `protocol.json`。本轮只做合成机制诊断，不使用真实 Holdout，不新增生物学 motif 结论。

比较七条流程：One-hot；历史 A0 的全上下文碱基/窗口均值；ProkBERT-mini 全上下文与 24bp 局部重编码各自的均值和保留位置拼接；ProkBERT-mini-c 单碱基表示的全上下文均值。四组 MEAN/ORDER 使用同一批完全落在核心窗口内的 token；A0 另列，不能混同支持范围。

强模式、15%/30% 退化模式均随机位置植入，两个新数据种子，共 6 个条件、42 个拟合。固定算法种子，12bp 模式与 12bp 窗口一致，避免把部分锚点当完整模式恢复。此抽象连续模式不代表真实 −35/−10 间距结构。二维因素比较限定于 mini 的四组 MEAN/ORDER。局部输入长度变化仍有分布偏移；mini-c 与 mini 的权重/结构也不完全相同，不能据比较作纯分词因果结论。

发现只用 Discovery 序列，PCA、24 簇聚类和每簇每序列最近窗口规则统一。全部拟合后写入 `fit_lock.json`，再评价 Development/Test；不根据 Test 更改方法或选阈值。阈值来自独立合成校准空白序列最大分第95百分位，严格大于。参考辅助选择一个候选并要求全宽 cosine≥0.8，定位要求方向正确、起点误差≤2bp；该评价不是无监督发现 precision。所有候选另报不受参考门槛过滤的空白误报率。

```powershell
python -m experiments.pr01_02_embedding_v2.bottleneck_20261025.pilot prepare
python -m pytest experiments/pr01_02_embedding_v2/bottleneck_20261025/test_pilot.py -q
python -m experiments.pr01_02_embedding_v2.bottleneck_20261025.pilot run
python -m experiments.pr01_02_embedding_v2.bottleneck_20261025.pilot verify
```

`run` 拒绝覆盖已有 `artifacts/`。本地保留人工合成序列、实例、模型拟合；权重/缓存不入 Git。公开 `summary/` 聚合统计、逐候选空白误报、来源/输出哈希和加载报告，以及最终 `REPORT.md`。两个数据种子仅提供描述性诊断；多算法种子、实例过滤/对齐、反向链植入、高GC/Markov、真实数据和合成 STREME 均留待下一轮。
