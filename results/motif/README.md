# motif 结果目录约定

本目录只保留当前仍有明确证据边界的 motif 派生产物。原始 FASTA、下载数据、完整 MEME HTML/XML、临时日志、大型缓存等不作为本目录正式结果；每个正式 run 必须记录输入哈希、工具版本、参数、样本量、输出哈希和证据限制。

## 2026-10-05 清理状态

已确认旧 RegulonDB 处理对 1,985 条已按转录方向提供的负链序列再次 reverse-complement，旧记录级 holdout 还存在近同源跨 split 泄漏。因此以下旧派生结果已经从当前工作树删除，不再作为课程报告或后续分析输入：

- `m3_20260930/`：基于旧 3,807 条 RegulonDB 输入的 MEME/Tomtom 汇总与矩阵；
- `regulondb_sigma_fimo_20261004/`：基于修复前六 sigma 输入的 FIMO/coverage/specificity/spacing；
- `regulondb_sigma_holdout_20261004/`：基于修复前输入和旧记录级 80/20 切分的 holdout 复核。

删除的是**确定由错误输入产生的派生结果**；原始 RegulonDB 下载快照、原始 GFF3/基因组及新的 G0 修复数据不删除。Git 历史仍可用于审计旧文件，但当前分支不得继续引用这些结果作为科学证据。

当前数据角色为：TJU Pan / PromLoop 六物种是 RQ1/RQ3 主序列数据；RegulonDB/DBTBS 主要用于 sigma annotation、机制解释和同物种跨数据集验证。在新数据方案验收前不运行新的正式 MEME/FIMO。

现存 `tjupan_m3_20261003/` 属于 TJU 主数据的**历史探索结果**。其输入不是上述 RegulonDB 方向错误数据，因此本轮不因 RegulonDB 清理而删除；但正式新实验仍需先完成标签一致性、exact/RC/近同源隔离和背景敏感性验收，不能自动把历史结果升级为最终结果。
