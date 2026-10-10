# PR01-02 DNA Embedding 独立实验路线

最新成果入口：[2026-10-10 实验与结果归档](SYNC_20261010.md)。

本目录是传统 Motif 主线之外的新实验分支。EMB 编号只描述本分支工作包，原 M1/M2/M3 状态和冻结数据不改写。实际成果以 `audit/` 和 `runs/*/manifest.json` 为准。

## 命名和初始协议

- 模型名统一为 **Frozen ProkBERT-mini**，模型仓库 `neuralbioinfo/prokbert-mini`；禁止使用 promoter 分类微调模型代替预训练编码器。
- 模型和实现代码各自固定完整 commit，见 `config.json`。本地下载代码、词表、权重全部记 SHA256，不使用随机初始化回退。
- 方法名 `EMB-PKM`：contextual token → base → window → PCA → MiniBatch K-Means → PWM；对照 `OH-PKM` 使用相同 PCA/聚类设置的窗口 One-hot。
- 诊断方法 `EMB-TW-PKM` / `tokenwindow` 只平均完整落在窗口内的token，保留完整启动子上下文；不是独立局部模型编码。`multiscale` 将五尺度按相同中心对齐后拼接，在Discovery拟合PCA；恢复实例使用16bp外包络。
- 运行名 `{source}__{representation}__L{layer}W{width}__s{seed}`。同名运行拒绝覆盖；失败运行保留。
- 初始模型层 3/6、窗口 6/8/10/12/16 均可缓存；初始发现固定 L6/W10，PCA=24、K=24，三种子 20261005/06/07。**缓存多个尺度不等于已经完成多尺度比较。**

## 数据与坐标

唯一主入口 `data/processed/pr01_02_data_v2/`。使用冻结 sequence master 和 FASTA，按 sequence_id 排序，不重拆分。81 bp 中 TSS 第61位/Python索引60，相对位置0。Token i覆盖 `[i,i+6)`；CLS/SEP不映射到碱基。碱基向量为覆盖该碱基的token向量平均；窗口再对碱基平均，因而带有窗口外侧和完整启动子的上下文，不等同于独立局部编码。

Discovery正类用于拟合PCA/聚类/PWM，无Motif标签。正样本入口是标签筛选后的启动子集合，不能称为对全部无标签序列的盲发现。每簇每启动子只保留距中心最近窗口，避免重叠窗口重复加权；固定同长、正向窗口对齐。相似去重要求完整窗口两方向中心化cosine≥0.95，保留合并映射。该起步方案尚未实现可变边界和gapped alignment，可能使Motif漂移或混合。

## 独立验证口径

PWM固定后扫描Development正类、天然label0和配对Shuffle。两链零阶log-odds扫描，每条序列只取最高分位点。背景为Discovery天然对照的互补对称组成，阈值为该Discovery对照最大分数95%分位点，严格 `score > threshold`。阈值不看Development或Holdout；此经验阈值不是FIMO site p/q。

天然对照按冻结leakage_group聚合为是否命中；正负两侧共有组从Fisher两侧剔除，单侧富集检验并报OR/区间。Shuffle继承父启动子group并做配对单侧二项检验，避免将配对背景当独立样本。两套检验各自BH，初步单运行家族含所有Embedding与STREME候选；跨物种汇总时重新统一BH，保留原p。多种子用于稳定性而不是新增独立样本。

现有STREME固定PWM使用同一自定义扫描和阈值规则作方法对照，原FIMO结果不修改，也不把两种扫描的q值等同。PWM相似只做描述，不能把余弦相似当Tomtom显著性。Grammar目前仅输出位点顺序/edge gap，不宣称置换显著或功能协同。

Holdout默认关闭。传统项目已查看过该留出结果，未来即使锁定Embedding配置也只能如实注明探索性再评估；不能声称新盲确认。RegulonDB外部候选需另验与全部TJU split的关联；DBTBS一条独立留出不做统计。TJU上游TSS对齐未逐记录重建。

## 可执行验收标准

| 工作包 | 完成判定 |
|---|---|
| EMB-01 | raw/output全哈希一致；81bp/TSS/冻结split正确；Shuffle父ID可追溯且二核苷酸组成核验；物种/多σ/辅助来源审计；环境和模型访问已记录 |
| EMB-02 | 精确版本权重全部加载；官方token与显式6mer坐标逐输入一致；坐标/边界测试通过；缓存含层/尺度/ID/文件哈希；覆盖六物种才记全量完成 |
| EMB-03 | Discovery拟合；实例可逐条还原DNA和坐标；PFM/PWM/logo；固定扫描及Development天然/Shuffle分别评价；阴性结果保留 |
| EMB-04 | 同一数据与评价下实际运行局部/上下文、尺度、层、One-hot或kmer对照；多尺度必须实际融合及消融，单独缓存不计完成 |
| EMB-05 | 顺序/位置/gap可复算；独立组样本充分后才做冻结的置换检验；σ证据不足时描述；必须区分描述与统计 |
| EMB-06 | 基础可解释候选经稳定性和独立富集验证后，才考虑SAE；三类SAE每类有重构/稀疏/死亡特征/定位/消融证据；无前提则记录未启动原因 |
| EMB-07 | 全物种公平对照、合成发现恢复（不仅扫描）、多种子匹配、混杂审计、外部独立性、成本、锁定配置后留出及失败报告齐备 |

工程通过不等于生物学成立。候选提升需Development中两类对照分别有BH支持，三种子有PWM/位点稳定性，且相对One-hot/传统方法有明确增益；数据库无匹配不是新功能证据。基础结果不足时保留负结果，暂不训练SAE。

## 运行

从仓库根运行（现有Python环境；固定实际版本见 `audit/environment.json`）：

```powershell
python -m experiments.pr01_02_embedding_v2.audit
python -m pytest experiments/pr01_02_embedding_v2/test_pipeline.py -q
python -m experiments.pr01_02_embedding_v2.run --sources tjupan_bacillus_subtilis --representations prokbert onehot --seeds 20261005 20261006 20261007
python -m experiments.pr01_02_embedding_v2.ablate
python -m experiments.pr01_02_embedding_v2.synthetic --tag synthetic_tokenwindow
python -m experiments.pr01_02_embedding_v2.summarize
python -m experiments.pr01_02_embedding_v2.verify
```

模型文件、缓存均放本分支内并忽略Git；参考资料文本和审计保留。本分支不自动commit/push。CPU推理固定线程，运行耗时计入Manifest；没有CUDA时不报告GPU显存。失败原因写Manifest，模型载入阶段失败另见环境/最终报告，不生成伪Embedding结果。

已存在同名输出时不能直接重跑；复制本分支到新的项目内运行目录或使用独立checkout再复现，不删除原始结果。完整六物种命令及本次结果见 `REPORT.md`。
