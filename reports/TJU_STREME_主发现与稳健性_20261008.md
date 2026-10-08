# TJU 六物种 STREME 主发现与背景稳健性（2026-10-08）

已执行 STREME 5.5.9 十二组（六物种×两背景），全部退出码0。冻结数据v2不变，仅discovery输入；未扫描或分析项目development/holdout。原始CSV一致性核验六组通过，新增7项测试通过，全套142项通过（16项既有依赖警告）。原始输入和210项v2产出哈希全部匹配。完整原生PWM/XML/HTML/site表/日志在`results/motif/tju_streme_v2_20261008/`，实际命令和输入输出SHA256见manifest.json。

## 固定参数与软件纠错

```text
--dna --p <discovery.positive.fasta> --n <discovery.natural_control.fasta 或 discovery.dinucleotide_null.fasta> --order 2 --minw 5 --maxw 15 --thresh 0.05 --nmotifs 10 --seed 20261005 --o <新的独立目录>
```

参数顺序必须保留：本地5.5.9源码`tmp/meme-suite-build/meme-5.5.9/src/streme.c:350–356`证明后解析`--thresh`会把nmotifs清零。最初照示例顺序执行的12组保留于`tmp/streme_v2_20261008_parameter_order_diagnostic/`，不作为正式结果。本次先thresh后nmotifs，XML确认每组10个输出。命令与预注册统计目标相同；纠错时未读取项目holdout。完整初次/正式输出均保留，不比较择优。

显式negative下order2用于估计Markov背景，不重新shuffle输入。STREME默认内部10%测试样本来自discovery，不能称项目外部holdout；最终PWM会使用全部discovery正类。nmotifs上限输出包括不显著候选。参见[官方STREME说明](https://meme-suite.org/meme/doc/streme.html)。

## 实际结果

各组输出10个motif；表中只计达到相应软件显著性标准的数量。p与E分别保留；这里的p不是FIMO site q，也不是后续BH q。

| species subset | natural p≤0.05 | dinucleotide p≤0.05 | natural E≤0.05 | dinucleotide E≤0.05 |
| --- | ---: | ---: | ---: | ---: |
| bacillus_subtilis | 3 | 0 | 1 | 0 |
| baumannii | 2 | 1 | 1 | 0 |
| bradyrhizobium | 4 | 3 | 2 | 0 |
| diphtheria | 3 | 2 | 2 | 1 |
| escherichia_coli | 3 | 1 | 1 | 0 |
| staphylococcus | 2 | 2 | 2 | 2 |

natural共17个通过事前p≤0.05规则，按p升序、ID次序，每物种最多5个，exact/RC PWM去重复后仍17个。`selected_main_motifs.tsv`保存主名单，每物种`selected_main.meme`供后续固定扫描。这是发现候选，不是已确认生物功能；主列表未根据稳健性结果改变。

稳健性组通过p规则9个，其中B.subtilis为0；E规则总计3个。不能用不同背景的p值大小直接衡量同一个PWM的富集变化，两组会重新发现不同PWM。已完成六组Tomtom（Pearson、minimum overlap5、threshold1保留匹配），完整offset/overlap/orientation/p/E/q见各物种`robustness_tomtom/tomtom.tsv`。目标库每组仅10个，匹配只作PWM相似度描述，不证明功能，也不把匹配不到解释为没有富集。要回答固定主PWM在两背景的富集差异，还需后续固定PWM扫描。

## 下一步与边界

DBTBS只承担TJU B.subtilis可靠映射的sigma/坐标/文献注释，不将1条外部样本构成正式统计分支。Dataset.csv是正式序列/标签唯一原始入口，positive_samples保留raw归档；核验程序读取它仅用于独立原始一致性审计，不参与实验输入构建。natural为主背景、dinucleotide为稳健性，不建立GC-matched。

真实reference_library_v1.meme仍未生成。下一阶段按已冻结规范构建真实短位点/PWM参考，再进行已知元件Tomtom；FIMO、位置/间距及四组BH分析尚未执行，不能将本报告当作已完成RQ1–RQ3。

复现入口：`python -m experiments.run_tju_streme_v2`（拒绝覆盖既有输出）；`python -m experiments.summarize_tju_streme_v2`。固定运行配置见实际manifest。此前数据修复文档中的“未运行STREME”仅描述2026-10-05数据阶段；本报告为当前实验状态。
