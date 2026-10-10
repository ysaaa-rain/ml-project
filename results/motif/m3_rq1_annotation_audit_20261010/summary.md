# M3 RQ1 本地只读补证（2026-10-10）

参考库位点：462；按来源 ID、assembly、TSS、strand 和精确片段核验后可投影 462 条。投影失败：0。

RegulonDB/DBTBS 的重合表按来源记录 ID 关联，并要求唯一 TSS/strand、参考版本一致、6bp 元件片段与冻结 81bp 窗口一致。一个预测位点与多个真实区间重叠时按唯一预测区间和唯一注释区间分别去重；≥1bp 为主定义。

TJU 坐标覆盖逐条为 0/26980；原因是 PromLoop CSV 未提供逐条基因组坐标，本报告只统计 FIMO 输入中的相对 TSS 位点，不借用同序列外部坐标。

## RQ1 初步结论（已观察 / 未证实 / 无法评估）

**已观察：** 两种发现方法在六物种发现部分均有原生显著 motif；同源 PWM 比对产生 61 对 q≤0.05，其中两侧方法各自原生显著的 9 对。已知元件库 Tomtom 的 320 对中 15 对 q≤0.05，且同时满足 discovery 原生显著条件的 2 对均为 TJU B. subtilis SigA −10 相似候选。RegulonDB 与 DBTBS 注释元件坐标按来源 ID 成功投影，在发现序列内，RegulonDB −10/−35 的 18/13 个注释位点及 DBTBS SigA −10/−35 的 141/124 个注释位点中，native-significant strict FIMO 均未出现 ≥1bp 重叠（各运行/候选分母见逐 motif 表；注释位点重复展示不跨 motif 相加）。

**未证实：** PWM 相似或 FIMO 区间相交不单独证明功能；三条预先报告的探索性 holdout 阳性不构成新的确认结论，原有确认性统计未被改写。B. subtilis 两条 natural-control holdout motif 的阳性集合 GC=32.6%、AT=67.4%，负样本 GC=46.4%、AT=53.6%；候选命中位点中心中位数分别为 −35.5bp 与 −24bp，说明组成差异可能参与自然负对照关联，需按探索性解释。位置窗口内命中率是描述统计。

**无法评估：** TJU 缺逐条 genomic coordinate / strand，不能计算真实已知元件的坐标重合率；当前合格参考库没有合格 UP PWM；RegulonDB 本地 operon 注释含 1 条 guaBp UP 区间叙述（−59 至 −38），但底层 citation record `RDBECOLIPRC22477` 在快照中未解析出 PMID/证据对象，因此暂不纳入正式 recovery。结构化坐标化 TFBS 仍无合格覆盖，不计算恢复率。

全输出 motif 分析范围：540 条报告 PWM；其中 MEME 原生显著 50，STREME 29。位置与重合结果均在 `tss_window_distribution.tsv` / `element_overlap.tsv` 中按 `all_reported` 与 `native_significant_only`、split 和 FIMO 门槛区分。

位置表的 −10/−35 中心窗口仅用于来源/σ分组确有统一参考架构的子集；TJU 缺逐条 σ、RegulonDB/DBTBS overall 是混合组，因此只给相对 TSS 中心分布并将架构窗口标 NA。DBTBS SigW development discovery 输入为空，相关 120 行标记为 FIMO 未运行，不记作零命中。

## 证据边界与目录

旧首次 FIMO 诊断的 SEA 34 项（33 项 native E≤0.05）、raw p<1e−4 的 1,494 个位点、原 strict 汇总 45 与完整 FIMO 表 40 分开保留。Staphylococcus development 的 19 与 14 是不同 q-value 估计路径/零分布规模下的两个证据口径，不合并。

逐条映射及注释记录仅保存在本机 `tmp/m3_rq1_annotation_audit_20261010/`，本次公开文件只有不可逆汇总、短 motif 共识和统计图。输入/输出校验和见 `input_output_sha256.json`。
