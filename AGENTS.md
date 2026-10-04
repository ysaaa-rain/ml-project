# Repository instructions

## 后续实验路线（用户指定，2026-10-04）

开始任何实验设计、数据处理、运行、结果解释或报告前，阅读 `docs/07_后续实验执行规范_20261004.md` 和 `docs/参考资料/PR01-02_后续实验指导_20261004.txt`。后续实验以用户新指导为准；与旧路线冲突时按执行规范处理，保留旧结果的原始证据边界。

采用经审计 RegulonDB/DBTBS 发现 → 传统新实验1–5 → 新6 sigma预测消融 → 新7 CNN/归因/TF-MoDISco独立发现的路线。embedding/SAE为有明确motif提取机制的选做扩展。不能把pooled embedding、attention或分类性能直接称为motif/因果发现。

不编造TSS/链/缺失窗口；按同源簇隔离训练与测试，motif发现和特征选择仅用训练数据。验证集不重发现、不调阈值。严格site q与raw p候选分开；阴性结果须保留。课程原编号与新规范编号分开标注，课程必做消融不得由ML消融替代。

保留原始资料与用户已有修改。同步远端前检查工作区；只提交明确授权路径，不自动提交/推送。每次新运行单独保存配置、结果、失败/限制和manifest。

## 已确认的审计纠错要求

先读 `docs/08_全方位核查审计_20261004.md`。旧RegulonDB坐标Gate已撤回：Sequence已经转录定向，不得再按genomic strand翻转；重新适配原始GFF3才可窗口处理。保留所有sigma标签，多标签不能只取首行。旧80/20留出存在跨split近重复，仅限显式历史复现；正式划分须按全局同源簇。软件测试及hash通过不等于旧生物学实验重跑完成。

## PowerPoint default workflow

For any task that creates, edits, redoes, polishes, lays out, reproduces, or visually reviews PowerPoint, PPT, PPTX, slides, or decks, default to loading and following `$ppt-quality-loop`.

Unless the user explicitly requests file-only output or explicitly declines visual review, when PowerPoint MCP is available do not declare PowerPoint work complete without real PowerPoint screenshot visual verification.
