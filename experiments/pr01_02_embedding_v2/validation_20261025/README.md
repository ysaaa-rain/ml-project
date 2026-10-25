# 有序局部表示复核与瓶颈拆分（2026-10-25）

本轮使用全新合成种子，不覆盖 `bottleneck_20261025/`。冻结参数见 `protocol.json`。

本轮已完成，见[阶段报告](REPORT.md)。原流程多种子15%退化F1约0.695，未稳定复现上一轮约0.900；扩大PCA到128维约0.790，仍低于One-hot约0.900。监督探针显示30%退化的局部表示仍可读出约73.0%的位置，提示不能把无监督失败直接归为无信息。固定距离筛选无一致提升，反向链混合仍是短板；未进入真实数据。

主实验：3难度×5数据种子×3聚类种子，比较One-hot与局部token拼接的PCA24/64/128维及固定距离筛选。压力测试：高GC、Markov相关背景、混合反向链，各2数据种子×3聚类种子×3固定方法。合计234次实际聚类拟合、324组PWM输出；筛选对照复用同一聚类，不重复计作新拟合。

先完成并锁定所有无监督发现，再评价新Development/Test；最后用30%退化的5个数据条件执行25个监督定位探针。探针只回答已知标签训练后位点信息是否可读出，不算无监督motif发现或真实功能证据。

```powershell
python -m pytest experiments/pr01_02_embedding_v2/validation_20261025/test_pipeline.py -q
python -m experiments.pr01_02_embedding_v2.validation_20261025.pipeline fit
python -m experiments.pr01_02_embedding_v2.validation_20261025.pipeline score
python -m experiments.pr01_02_embedding_v2.validation_20261025.pipeline probes
python -m experiments.pr01_02_embedding_v2.validation_20261025.pipeline verify
```

模型版本和哈希复用上一轮固定元数据，优先只读旧缓存；缺缓存时下载到本轮目录，不改旧审计。每条81bp序列的70个核心窗口只需编码58个不同24bp裁剪，边缘共享裁剪，token和核心映射不变。

PCA固定一个Discovery基底，24/64/128为其前缀，固定PCA随机种子；只改变KMeans种子。CORE50在每簇内按距离中位数保留窗口，再每序列取最近窗口，至少20条序列形成候选。所有门槛仅用Discovery，无真值参与发现。

主统计先在同一数据种子内平均3个聚类种子，再汇总5个数据种子。保留每次运行及配对差值；不把15次运行说成15个独立生物学重复。压力测试仅两个数据种子，结论更有限。真实Holdout继续禁用，STREME当前运行环境不可用时如实记录，不能宣称总体胜过传统方法。

原流程复核另按 `replication_protocol.json` 执行90次拟合：直接复用未修改的上一轮fit函数，每个算法种子重拟合PCA24/KMeans。其输出与共用PCA基底的消融分开，避免计算路径变化混入“原流程复现”。主矩阵拟合完成后、任何Test评分前运行 `python -m experiments.pr01_02_embedding_v2.validation_20261025.replicate fit`，之后分别score和verify。

原始人工序列、局部特征、拟合对象、实例和模型缓存留本地。公开协议、代码、聚合统计、误报和诊断表、来源及输出哈希。`fit`拒绝覆盖现有artifacts。最终结果见运行后生成的`REPORT.md`。
