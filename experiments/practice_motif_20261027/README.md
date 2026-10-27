# PR01-02 课程实践候选总览

本版本将真实序列传统方法、真实序列Embedding历史结果和已有合成实验整理到统一入口，面向候选查看与方法学习。保留所有候选，数字来自已完成运行；不改变历史结果或原始数据。配置见 `protocol.json`，输出位于 `reports/practice_motif_20261027/`。

合成试验中，从每个已锁定候选集合出发，按已知合成参考寻找重叠至少8/12bp的最佳局部匹配，允许正反向与偏移，将参考起点投影回输入序列并按±3bp计算参考辅助定位F1；不按相似度设额外入选截止值。该评价使用合成真值，是对已看过的Development/Test的事后描述，不是新的盲验证。匹配长度、相似度、候选集合空白误报及空候选均一并保留，局部相似不等于功能motif。

```powershell
python -m pytest experiments/practice_motif_20261027/test_overview.py -q
python -m experiments.practice_motif_20261027.build
```

输出目录已存在时拒绝覆盖。历史主线和先前Embedding报告仍为数据来源与历史口径；搜索结果是候选记录，运行间重复不算新的独立motif。
