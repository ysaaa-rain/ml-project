# 同批枯草芽孢杆菌候选比较

已完成：[报告](summary/REPORT.md)｜[候选比较图](summary/comparison.png)｜[CNN位置图](summary/cnn_positions.png)｜[核验](summary/verification.json)。

统一复用传统方法的冻结真实数据，以 Discovery 排序并固定90%分位 PWM 扫描阈值，在同一 Development 正类、天然负类及配对打乱背景上评价。包括原 STREME 自然背景运行全部10个候选、已有两个种子的CNN、One-hot、ProkBERT候选。CNN使用序列标签，STREME只运行一次；已有Development结果已观察，为探索。MEME完整PWM尚未找到，暂不列入比较。

在隔离副本中依次运行 `python -m experiments.pr01_02_embedding_v2.real_bsub_comparison_20261027.run`、`python -m experiments.pr01_02_embedding_v2.real_bsub_comparison_20261027.report`、`python -m experiments.pr01_02_embedding_v2.real_bsub_comparison_20261027.verify`。已有 summary 目录时拒绝重算覆盖；原始模型、实例和数据不上传。
