# 扩展 M3 聚合证据（已完成）

54 套结果、540 个按运行编号的 PWM，344 次 FIMO 全部成功。详见 run_evidence.json 与 ../../.. 下 reports/M3_扩展双方法与sigma实验_20261009.md。

8 次辅助 SEA 因 development 仅 1 条正样本失败，FIMO 不受影响；DBTBS 六组 discovery<10 登记跳过双方法。

cluster_tests.tsv 包含簇级出现率差、95% 保守区间、原始 p 与按 split/family 的 BH；新增 holdout 全为探索性。原 17 主候选与旧结果保持不变。原始序列、逐位点 matched_sequence、受限来源 PWM 不发布。

validation.json 记录冻结数据与原结果哈希核验、工具运行计数；publication_sha256.json 记录本目录公开文件哈希。完整原生输出及补充 Tomtom 执行收据仅本地保存。
