# CNN 在标准真实数据上的小试

已完成：[结果报告](summary/REPORT.md)｜[候选图](summary/candidate_profiles.png)｜[来源与核验](summary/verification.json)。两种子分别生成15和13个候选；前3候选 Development 正类组命中率0.848/0.818，天然对照组0.325/0.247，配对 Shuffle 0.667/0.652。CNN 使用序列级标签监督，结果为已观察数据上的探索。

使用 MEME/STREME 同一套冻结 TJU B. subtilis 数据及划分，两个模型种子、天然负类监督，生成原窗口和 Discovery 信息量平移后的 PWM。Discovery 天然负类90%分位阈值、最少10条支持序列用于探索性候选；统一双链 PWM 扫描及 Discovery 90%分位阈值，分别报告 Development 正样本、天然负类、配对 Shuffle。具体见 `protocol.json`；不使用 Holdout。结果将以实际产物为准。

拟合只用 Discovery，随后锁定模型、实例、PWM 与排序再扫描；Development 之前已被项目观察，本轮属于探索。CNN 使用序列级标签，候选与既有同源数据 STREME PWM 的相似是描述性比较，不代表功能确认。模型、序列实例与中间文件保留在本地忽略的 artifacts 中。

```powershell
python -m experiments.pr01_02_embedding_v2.cnn_real_pilot_20261027.pipeline fit
python -m experiments.pr01_02_embedding_v2.cnn_real_pilot_20261027.pipeline score
python -m experiments.pr01_02_embedding_v2.cnn_real_pilot_20261027.pipeline verify
python -m experiments.pr01_02_embedding_v2.cnn_real_pilot_20261027.report
```
