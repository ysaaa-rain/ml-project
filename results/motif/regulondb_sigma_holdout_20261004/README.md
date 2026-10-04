# RegulonDB 六 σ 组留出复核

本目录汇总同一 RegulonDB 数据集内的固定训练/留出复核。每组用约 80% 配对启动子及 N1 对照发现 10 个 MEME motif；剩余 20% 只用于 FIMO 检验。此结果用于内部可重复性检查，不代表跨数据集验证。

## 主要结果

- 6 组共留出 636 对序列；其中 15 个训练 motif 的 MEME E-value≤0.05。
- 来源组正类与 N1 的命中差异：raw site p≤0.01 时无 motif 通过 BH；p≤0.001 时 3/15 通过；p≤0.05 时命中率接近饱和。
- 来源 σ 组与其余组比较：Sigma70 MEME-10 和 Sigma24 MEME-3 在 p≤0.01、0.001 两档通过 BH。
- 组内 270 个 motif pair：p≤0.01 时 13 个正类共现高于 N1 且通过 BH；p≤0.05 和 p≤0.001 时均为 0 个。
- 严格 FIMO site q≤0.05 下候选命中稀少；结果不支持普遍、稳定的 σ 特异调控语法。

## 文件

- `sigma_motif_coverage.tsv`：360 行，motif 在来源组和目标组的正类/N1 命中率及配对检验。
- `sigma_motif_group_specificity.tsv`：1,080 行，来源 σ 组与其余五组的命中率比较。
- `sigma_motif_pair_spacing.tsv`：810 行，共现检验、最近有符号间距、重叠和排列描述。
- `matrices/`：六个 σ 组各自从训练集发现的 MEME PWM 矩阵。
- `sigma_fimo_manifest.json`：样本拆分、随机种子、MEME/FIMO 参数、版本及输出哈希。

原始输入、拆分 FASTA、MEME 原生文件和逐位点 FIMO 文件保存在项目内 `tmp/regulondb_sigma_holdout_20261004/`，不纳入 Git。复现命令：

```bash
.venv/bin/python -m experiments.validate_regulondb_sigma_holdout \
  --work-dir tmp/regulondb_sigma_holdout_reproduction \
  --result-dir results/motif/regulondb_sigma_holdout_reproduction
```
