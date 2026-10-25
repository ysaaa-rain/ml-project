"""Build a reviewable Chinese report from completed, verified pilot outputs."""
import pandas as pd

from .pilot import HERE, PROTOCOL, SUMMARY, read

LABELS = {
    'OH': 'One-hot', 'MINI_FULL_A0': 'mini 完整上下文 A0',
    'MINI_FULL_MEAN': 'mini 完整上下文 token 均值',
    'MINI_FULL_ORDER': 'mini 完整上下文 token 拼接',
    'MINI_LOCAL_MEAN': 'mini 局部重编码 token 均值',
    'MINI_LOCAL_ORDER': 'mini 局部重编码 token 拼接',
    'MINIC_FULL_MEAN': 'mini-c 完整上下文碱基均值',
}


def markdown_table(frame):
    header = '| ' + ' | '.join(map(str, frame.columns)) + ' |'
    separator = '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |'
    rows = ['| ' + ' | '.join(map(str, row)) + ' |' for row in frame.itertuples(index=False, name=None)]
    return '\n'.join([header, separator, *rows])


def main():
    cfg = read(PROTOCOL); manifest = read(SUMMARY / 'run_manifest.json')
    verified = read(SUMMARY / 'verification.json')
    assert manifest['status'] == 'COMPLETED' and verified['status'] == 'PASSED'
    summary = pd.read_csv(SUMMARY / 'test_summary.tsv', sep='\t')
    metrics = pd.read_csv(SUMMARY / 'metrics.tsv', sep='\t')
    primary = metrics[(metrics.split == 'test') & (metrics.tolerance_bp == 2)]
    table = pd.DataFrame({'方法': [LABELS[m] for m in cfg['methods']]})
    for scenario in [s['id'] for s in cfg['scenarios']]:
        indexed = summary[summary.scenario == scenario].set_index('method')
        table[scenario] = [f'{indexed.loc[m, "f1_mean"]:.3f} [{indexed.loc[m, "f1_min"]:.3f}, {indexed.loc[m, "f1_max"]:.3f}]'
                           for m in cfg['methods']]
    seed_table = primary.pivot(index=['scenario', 'data_seed'], columns='method', values='f1').reset_index()
    for method in cfg['methods']:
        seed_table[method] = seed_table[method].map(lambda v: f'{v:.3f}')
    nulls = pd.read_csv(SUMMARY / 'candidate_null_fpr.tsv', sep='\t')
    null_summary = nulls[nulls.split == 'test'].groupby('method').raw_blank_fpr.agg(['mean', 'min', 'max']).reindex(cfg['methods'])
    null_summary = null_summary.map(lambda value: f'{value:.3f}').reset_index()
    null_summary['method'] = null_summary.method.map(LABELS)
    fits = pd.DataFrame(manifest['runs'])
    pca = fits.groupby('method').agg(dimensions=('input_dimensions', 'first'),
        pca_variance_min=('pca_explained_variance', 'min'), pca_variance_max=('pca_explained_variance', 'max')).reindex(cfg['methods'])
    for key in ['pca_variance_min', 'pca_variance_max']:
        pca[key] = pca[key].map(lambda value: f'{value:.3f}')
    pca = pca.reset_index(); pca['method'] = pca.method.map(LABELS)
    text = f'''# Embedding 定位瓶颈与 mini-c 对照：阶段实验报告

日期：2026-10-25。版本：{cfg['version']}。状态：本轮小规模合成实验完成，整体生物学发现仍未完成。

## 本轮完成范围

完成三个合成条件（强模式、15%/30%逐位退化），每个条件两个新数据种子，共6个条件、{len(fits)}个实际拟合。加入真实执行的局部重新编码、保留 token 位置顺序与 Frozen ProkBERT-mini-c 对照。全部拟合先完成并锁定 PWM，之后统一评价 Development/Test；没有依据 Test 调参。人工合成的连续12bp模式不是生物学 −35/−10 间距模型。

真实 v2 不重建、不重拆分；本轮没有真实数据拟合，也没有使用真实 Holdout。原首轮和 next_round 结果不覆盖。

## 协议与可比性

- 新数据种子为 {cfg['data_seeds']}，算法种子固定 {cfg['algorithm_seed']}；两个数据重复只提供描述性证据，不能把不同方法或跨难度条件当成额外独立重复。
- 每条件 Discovery 256、Development 128、Test 128 条植入序列，每个 split 有同数空白序列；另有256条阈值校准空白。各分区使用独立随机流，生成时核验全部分区内/间无 exact/RC 重复。没有重新执行全量近同源审计。
- 序列81bp，核心窗口12bp，与完整植入模式等长。Discovery仅用正序列拟合 PCA24维、MiniBatchKMeans24簇、n_init=5；每簇每序列取最近窗口，至少20条序列建 PWM。所有方法使用同一下游流程，未加入实例过滤、对齐或微调。
- mini 的四组 token MEAN/ORDER 均使用完整落在核心窗口内的7个6-mer hidden states；均值和拼接的 token 支持范围相同。A0使用历史 token→base→window 均值另作对照，不能混同这两种聚合定义。
- 局部编码是真正重新运行编码器；输入24bp，通常核心两侧各6bp，边缘裁剪向内移动，70个窗口的偏移见 `summary/crop_map.tsv`。每条件17920个局部片段，六条件共107520个；输入变短仍带来分布变化，不能把差异归于单一上下文机制。
- mini-c使用单碱基 token，但其权重和模型结构与mini也有差异。本轮是完整编码器对照，不是纯分词因果实验。

## 评价规则

使用 Discovery 空白序列计算对称零阶背景。每候选双链最大 log-odds 分数的阈值固定为独立校准空白的95百分位，严格 score>threshold，不是FIMO位点p/q。

参考 PWM 只在评价阶段匹配候选，匹配要求全12列 cosine≥0.8；只对最佳匹配家族计算参考辅助定位。TP要求方向正确且起点误差≤2bp；位置错误计FP，所有植入位点进入召回分母，空白预测也进入precision分母。该规则含参考辅助家族选择，不能称为无监督发现的绝对precision。门槛未通过会得到0分，不证明hidden state完全没有信息。

所有候选另行报告**不受参考匹配门槛过滤**的空白误报率，避免0分掩盖原始误报。`metrics.tsv`保留0/1/2/3bp容差、cosine、门槛、TP/预测/真值分母。

## 新合成 Test 结果

表中为两个数据种子的F1均值 [最小值, 最大值]，不是置信区间。

{markdown_table(table)}

逐数据种子结果：

{markdown_table(seed_table)}

每个种子的相对One-hot与A0配对差值分别在 `summary/paired_delta_vs_OH.tsv` 与 `summary/paired_delta_vs_MINI_FULL_A0.tsv`。不根据候选数或只展示最佳种子评价方法。

## 本轮判断与下一步选择

1. **当前优先继续改表示与抽取流程，不直接更换主模型。** 在强模式条件，完整上下文拼接、局部拼接与One-hot的F1均约0.977；A0与完整上下文均值未通过参考门槛。这说明当前模型配合不同表示构造可以恢复模式，不能把原流程失败直接归为模型完全没有局部信息。
2. **局部重编码与保留顺序的组合有继续验证价值。** 15%退化条件中，局部拼接两个数据种子分别约0.898/0.901，均值约0.900；One-hot约0.919。完整上下文均值/拼接各只有一个数据种子过门槛；局部均值两个种子都失败。单独裁短或单独拼接都不足以解释本轮全部结果。这是两个样本重复下的描述性比较，尚未建立普遍机制或稳定性结论。
3. **重度退化仍是未解决的失败模式。** 30%退化中所有Embedding流程未达到参考匹配门槛，而One-hot约0.553。保留0分与原始cosine，不降低门槛重新包装成功；下一轮应检查PCA降维、窗口聚类和实例混合，而不是先增加模型规模。
4. **mini-c均值方案暂不升级。** 本轮三个条件、两个数据种子均未过参考门槛；结论仅适用于该编码器配合当前均值/PCA/聚类流程，不否定mini-c的其他用法，也不能作纯分词因果判断。
5. **下一步优先做有序局部表示的复核与瓶颈拆分。** 新版本用新的合成Test种子，将数据重复扩到至少5个，单独加入多个算法种子；先区分PCA24维是否丢失定位信息，再检验实例距离过滤/偏移与方向对齐。补齐高GC、Markov、反向链和STREME后，再决定是否进入真实数据小试。SAE与微调继续后置。

本轮窗口12bp、算法种子与植入位置方案和历史10bp实验不同；应比较本轮同条件方法，不能把跨轮结果拼成单一提升曲线。

## 全部候选的空白误报

下表汇总每方法六条件、每条件24个候选在Test空白上的FPR，未按参考匹配过滤。每个候选独立校准，不能把这里的均值解释为整个候选集的家族误报率。

{markdown_table(null_summary)}

## 表示维度与 PCA

各方法统一降到24维，但原始维度和保留方差不同。这是当前发现流程的受控比较，仍需后续分离PCA瓶颈；较高解释方差本身不是motif信息充分的证明。

{markdown_table(pca)}

## 资源、版本与核验

本地Python {manifest['python']}、Torch {manifest['dependencies']['torch']} CPU、Transformers {manifest['dependencies']['transformers']}，4线程，无GPU。`run`端到端约{manifest['elapsed_seconds']:.1f}秒，包含载入、共享编码、拟合与统一评价，不包含模型下载和独立核验。逐方法 `fit_seconds` 仅是降维/聚类/PWM阶段，不能当作该方法独立端到端成本。

mini固定revision `feb2520a43cd9cdb5b3d8477e47209dbcb55d1dc`；mini-c固定revision `417ebe73d355bc68b01d4fbe9dc2b0cb274b86e0`。两者均使用固定共享实现revision `366d9336b6b89d7c1feddf56d9b01be8065ccc43`，编码器加载均无missing/unexpected/mismatched/error keys。发布者LFS权重SHA256、配置与全部模型文件哈希见 `model_info.json`。

核验：{verified['hashes_checked']}项产物哈希、{verified['instances_checked']}条建模片段、{verified['truth_sites_checked']}条植入真值、{verified['crop_core_checks']}次裁剪/核心还原通过。冻结v2输出、53份历史运行manifest、旧报告等{verified['protected_files_unchanged']}个保护文件哈希不变。哈希一致不等于独立重跑全部历史实验。

相关测试13项通过（新边界/评价测试6项及原分支相关测试7项）。

## 未完成与下一步边界

本轮是协议固定的小规模机制诊断，不构成总体方法排名或新增功能motif证据。未执行：合成STREME、新真实数据、多个算法随机种子、反向链植入、不同模式长度/出现比例、高GC/Markov背景、局部侧翼长度梯度、实例过滤与对齐、监督定位探针、SAE、微调。合成STREME仍需可用的MEME Suite运行环境补齐。

下一步按本轮失败模式设计独立版本，优先检查PCA与实例抽取，不默认扩大模型数量。若依据本Test结果更改方法，必须换新的合成Test种子；本轮Test不能继续称为未查看的最终集。真实Development此前已反复使用，项目Holdout已被传统主线查看，未来不能包装成全新盲验证。

## 复现入口

```powershell
python -m experiments.pr01_02_embedding_v2.bottleneck_20261025.pilot prepare
python -m pytest experiments/pr01_02_embedding_v2/bottleneck_20261025/test_pilot.py experiments/pr01_02_embedding_v2/test_pipeline.py experiments/pr01_02_embedding_v2/next_round/test_methods.py -q
python -m experiments.pr01_02_embedding_v2.bottleneck_20261025.pilot run
python -m experiments.pr01_02_embedding_v2.bottleneck_20261025.pilot verify
python -m experiments.pr01_02_embedding_v2.bottleneck_20261025.report
```

`run`拒绝覆盖已有artifacts；独立重跑请使用新checkout或新版本目录。Git只发布代码、协议、聚合表、逐候选空白统计、哈希、模型加载与核验证据；权重、局部缓存、人工序列、逐实例与拟合对象留本地并忽略。不会提交无关删除项或历史PPT。
'''
    (HERE / 'REPORT.md').write_text(text, encoding='utf-8')
    print('REPORT.md generated from completed verified outputs')


if __name__ == '__main__':
    main()
