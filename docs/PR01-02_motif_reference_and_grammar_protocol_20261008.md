# PR01-02：已知 motif 参考与空间语法检验预注册（2026-10-08）

状态：**分析规则冻结；真实参考 PWM 文件仍待按本规范构建与 SHA256 验收**。此文件不宣称已有可审计的 reference_library_v1.meme，也不代表已运行 MEME/STREME/FIMO。

## 1. 项目数据角色与主程序

- TJU 六物种 v2 discovery 为 RQ1/跨物种主输入；holdout 不参与 motif 发现或参数选择。
- **六物种主发现工具为 STREME**。每物种正类使用 `discovery.positive.fasta`，唯一显式对照为同物种 `discovery.natural_control.fasta`（发布的 label=0）。固定 `--dna --p <positive> --n <natural_control> --minw 5 --maxw 15 --nmotifs 10 --seed 20261005 --thresh 0.05`，按实际版本与软件帮助核验参数。`--nmotifs 10` 为搜索上限，不代表十个 motif 均显著；只按软件给出的显著性及事前筛选规则解释。
- 不新造 GC-matched 负样本，不按 GC 分层或聚类，不把已生成 `dinucleotide_null` 纳入本轮正式主比较；这些旧文件可留作可追溯历史，不得当成实际使用的对照。**限制**：发现的区分力只能声称相对于原发布 label=0；正负 GC 差异可能解释部分 A/T 偏好，不能宣称已排除此混杂。
- MEME 用于选定组（如 TJU E. coli 和 B. subtilis）的验证性方法交叉检查，不重新从验证集发现；DREME 非必做。RegulonDB 为 sigma 专项，DBTBS 主要作 B. subtilis sigma/TSS/位点注释，不作为有效规模的独立验证集（432 核心，外部 holdout 仅 1 条）。
- 所有命令只读取 v2 冻结路径，不能继续读取 20260930/20261003 旧默认输入。正式启动前固定工具版本、参数、输入 hash 与输出目录。

## 2. 已知 motif 参考库——已冻结的来源选择和构建规范

**主参考库：按物种+sigma+功能分别建 PWM，不建立“六物种通用 E. coli −10/−35”的假参考。**

A. E. coli：以 RegulonDB 实验支持的 promoter 核心元件注释与可追溯论文位点为第一来源。按可靠 sigma/−10/−35/适用的 UP element 分开；不得把 GFF3 整个 81bp `Sequence` 冒充已知 −10/−35 的实测短位点。发布或文献已给出的经审计 PWM 也可直接纳入。

B. B. subtilis：以 DBTBS 4.1 页面的 sigma/相对 TSS 信息和红色标注的 cis-elements 对应文献为来源；要确认具体短位点归属和实验支持，不能将整个 cis-element 序列自动当作单一 PWM。独立来源 DNA 与 TJU 高度重叠，不按独立复现实验计数。

C. 其他四种：只纳入物种/菌株、sigma 家族、生物学功能及文献可核对的已知元件；缺少合格参考的物种只做 de novo 发现/位置描述，明确记为“无可用已知参考”，不强行比对 E. coli σ70 共识。

D. 额外的 TFBS 可使用 CollecTF 的实验验证位点；这些属于 TFBS，不能替代 −10/−35 的 sigma-promoter 参考。RegPrecise 的比较基因组预测仅可另列探索参照，不能标记为“实验验证”。

参考记录至少有 `motif_id, species, strain, sigma_family, element_type, source_database, release_or_date, source_record_id, PMID_or_DOI, experimental_evidence, site_sequence_or_PWM, n_sites, genome_accession_if_coordinate_based, orientation, evidence_tier`。从多个位点自己构建 PWM 前，冻结对齐、方向、去重复、伪计数规则；少于 10 个独立实测位点的自建 PWM 降为描述性参照。理想共识（E. coli σ70 −10 TATAAT、−35 TTGACA）单列 `consensus_only`，不能冒充实验 PWM。

输出固定为 `references/reference_library_v1.meme` 和 `references/reference_library_v1.tsv`（或明确指定的相对路径），保存来源下载快照/引用、文件 SHA256 和 motif 数量。**这些文件尚未产生前，参考来源和构建规范已定，但参考库本身尚未冻结**。

Tomtom 对照按适用物种/家族做；输出全部匹配、offset、overlap、orientation、p/E/q。非常小的目标库 p/q 容易不稳，若工具警告样本库过小，只作相似度排序/描述，不作“已验证已知 motif”推断。参考位点与 discovery 重合不得充当独立验证。

## 3. FIMO 与统计分析冻结

- 使用固定 discovery PWM 扫描 discovery、development、holdout 和对应 natural_control；阈值只在 discovery/development 上定，holdout 不调。正式 site 过滤为 **FIMO q ≤ 0.05**。另将 raw p < 1e-4 定义为探索候选，不能混同。FIMO 检验家族、`--max-stored-scores`、q 值算法和输出完整性写入 manifest。不要用截断 raw p 表自己伪算全局 q。
- FIMO start/stop 是本 81bp 输入序列内 1-based inclusive；TSS 在第 61 位。motif 中心相对 TSS = `(start+stop)/2−61`。motif hit strand 是相对输入串的方向，不代表原始 genomic strand。
- 选模仅用 discovery：每物种最多 5 个通过 STREME 事前显著性条件的非完全重复 motif（相似去重规则和选择顺序要随参数文件保存）；固定列表进入下游统计，绝不使用 holdout 调整其身份。
- 每个 promoter 对每个 motif 只取 q 最小、其次 score 最大、再次 start 最小的一个主 site；另保留所有 site 供完整审计。没有 site 记为 absent，不把位置缺失填成 0。对同一 motif-pair，使用两主 site；排序为转录方向位置靠前/靠后，`gap=downstream.start−upstream.stop−1`，负值为重叠，另存 center distance、二者 strand、原始两个 site。
- 模型固定后，holdout 中按相同方法扫描正样本及天然负类；统计单位为序列/同源簇，不将一个 promoter 的多个 hit 当独立样本。

## 4. BH 与 grammar 的主检验家族

全部基于事先选出的每物种最多 5 个 motif，以同物种 **holdout positive vs holdout natural_control** 为主对比。各家族内跨六个物种一起做一次 BH 校正，统一报告 raw p、BH q、效应量、95% CI；阈值 **BH q≤0.05**。若因 motif 数少而没有相应检验，清楚报告分母实际数量。

1. **单 motif 出现率**：每个物种每个 motif 一次 Fisher exact（至多 6×5=30 项），比较两类启动子至少一个严格 site 的比例。记录 OR 和 CI。
2. **双 motif 共现率**：每个物种每个无序 motif 对一次 Fisher exact（至多 6×C(5,2)=60 项），比较同一启动子同时命中两个 motif 的比例。这只能说明共同出现的组间差异，不能自动声称功能协同作用或 sigma 机制。
3. **相对 TSS 的位置分布**：每物种每 motif 比较含位点样本的中心坐标分布，采用固定的 cluster-aware permutation two-sample KS 或相同冻结检验（至多 30 项）；每组严格命中不足 20 个独立簇时只描述不做确认统计。
4. **pair 间距分布**：每物种每 motif 对，在两类样本均有至少 20 个独立 promoter/簇双命中时，用 5000 次按完整簇重排类别标签的 KS 检验比较 gap 分布（至多 60 项）；不足则记 NA、报描述图及原因，不临时合并不同 motif 或改阈值。

四个家族分别 BH，不能把每家族单独挑的最小 p 当作全局已校正。执行若需修改预注册设计，必须在看到 holdout 数据结果之前另存版本，并解释变更。各类阴性/零命中结果必须照常报告。

**机制专项预设**：仅对具有合格 E. coli σ70（或相应 B. subtilis sigma）−35/−10 元件归属证据的 motif 对，报告核心盒位置、转录方向顺序与 edge-to-edge gap 分布；E. coli σ70 可另以 15–19bp 作为事先列明的生物学描述范围。不向所有物种/σ 强推 σ70 架构；σ54 的 −24/−12 单列且不混入 −35/−10。

SpaMo 仅作方法对照：81bp 短窗要特别说明 margin/range、扫描集合和丢弃数，不能把它与自编 edge-gap 结果机械视为同一种检验或独立生物验证。

## 5. 验收后顺序

1) 重新运行 `preprocessing.verify_tjupan_positive_subset` 并归档输出（不动 raw）。
2) 确认 v2 manifest/hash 完整，保存本预注册规则和主 STREME 配置。
3) TJU 六物种 discovery STREME；其后固定 motif，再进行必要的 MEME 方法对照、已知参考库构建与 Tomtom、FIMO、position/spacing、holdout。
4) RQ2 以 RegulonDB E. coli sigma 为主；DBTBS 仅注释/探索。严格跨来源测试只报有效样本规模，不用重合的数据拼出独立性。
