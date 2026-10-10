# 后续研究材料入口

先读 [完整研究报告](REPORT.md)，再按其相对路径查看证据。本目录为首轮实验的新增研究轮次，不覆盖首轮运行或传统主线。

| 目录/文件 | 内容 |
| --- | --- |
| protocol.json / configs/ | 本轮协议、严格评价规则、历史实现快照 |
| stability_diagnostics/ | 家族匹配、位置容差、偏移、中心描述 |
| gc_diagnostics/ | GC 分布、CMH 探索、Shuffle 配对效应 |
| pooling_ablation/ | A0/A2 权重等价性 |
| synthetic_benchmark/ | 27 数据条件、132 方法拟合、严格评价、总体间距分母 |
| real_data_comparison/ | 18 个 A2 与 6 个 A0 重实现真实运行及汇总 |
| run_manifests/ | 输入/输出/代码哈希、耗时、环境；最终状态看 run_inventory_final.tsv |
| figures/ | 12 张真实数据图或明确标记的原理图 |
| verification.json | 5626 项哈希及新增片段/坐标/真值核查 |
| STATUS.json | 已完成项、未执行项、三页汇报约束 |

结论证据表与图索引保留在 `../presentation/`。旧课堂 PPT、逐页讲稿、摘要、问答及制作预览已按要求移至回收站；实验代码、运行结果、研究图表和完整报告保留。

`verify_results` 是证据核查；`figures` 与 `write_materials` 是派生材料更新；原始拟合脚本拒绝覆盖已有结果。若新增方法或使用本轮 Test 结果改方案，创建独立版本输出并用新的合成测试种子。真实 Holdout 不应在未冻结最终协议前启用。
