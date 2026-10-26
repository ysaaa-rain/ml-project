# 老师材料核对与下一条路线

2026-10-26。阅读本地原PDF及其已保存文本，并渲染核对关键页面。原材料含部分未解析的参考占位符，方法选择以文档内容和可核对的原始论文共同约束，不把材料中的所有书目信息视为已经证实。

| 文件 | SHA256 | 本轮核对位置 |
| --- | --- | --- |
| `docs/参考资料/motif-embedding-methods.pdf` | `0a9ed1b7ce927933f384cf05f791e2bc6f7e2e3b84f67db5a887392c0364a8a5` | pp.1–4 CNN，pp.7–12 SAE，pp.22–23比较与分层 |
| `docs/参考资料/motif-embedding-overview.pdf` | `68a79e8bfb775d5aed77562e6c1cd74668bd1a15a03878f5e1ba67c8b55750a2` | p.1 CNN，pp.2–5路线与比较 |

老师材料把 representation、extraction、reconstruction、validation 分开：CNN/GNN是架构，PLM是表示模型，SAE是分解模型，attention/归因是提取信号。材料列出六类，未规定必须依次运行；methods p.7特别强调SAE值得展开，并不意味着已有可直接用于本项目的DNA模型。

| 文档方法 | 当前安排 | 理由与边界 |
| --- | --- | --- |
| CNN/DeepBind式局部卷积表示→高激活实例→PWM | **本轮先执行** | 直接针对DNA局部模式，可在CPU训练；有序局部感受野、提取链清楚；需要序列级标签，与旧无监督聚类的监督预算不同 |
| SAE/MotifAE/InterPLM | 后续候选，未训练 | 优先级高，但代表研究主要为蛋白质；需明确DNA输入、重构/稀疏/死特征/定位验收，不能把蛋白权重直接当DNA模型 |
| PLM/DNABERT/NT/GPN/HyenaDNA | 当前不切换大模型 | 更换backbone不会自动解决提取机制；已收束的是ProkBERT聚类实例链，不否定整个家族 |
| 注意力解释 | 暂缓 | 延续Transformer依赖，attention不是因果重要性 |
| ISM/归因/TF-MoDISco | CNN有可解释结果后再评估 | 要有已冻结、任务有效的预测模型；分类准确不能替代位点恢复 |
| GNN/MMGAT | 暂缓 | 图构建与节点/边先验尚未定义，当前先检验单个局部模式 |
| 位置保留k-mer | 可作简单对照，非本轮优先路线 | 前轮建议是候选假设；按用户最新要求，先落地老师文档直接列举的CNN路线 |

论文核对：[DeepBind原文](https://www.nature.com/articles/nbt.3300)介绍从实验标签学习DNA/RNA结合特异性；[CNN架构研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4908339/)讨论卷积层作为motif扫描器与高响应片段重建。当前实现只借鉴该范式，不声称使用其发布模型或完全复现其超参数/数据。

[InterPLM原文](https://doi.org/10.1038/s41592-025-02836-7)使用ESM-2蛋白表示进行SAE分解；[MotifAE原始稿](https://pmc.ncbi.nlm.nih.gov/articles/PMC12637460/)同样研究蛋白语言模型。两者给出机制启发，不能作为本项目DNA功能验证。

新小试使用新合成数据、固定12bp窗口和随机位置，不接触真实Holdout。CNN只见Discovery序列级正负标签，不见植入坐标、方向、pattern或Test；OH48保持无监督发现加背景筛选。结果只能比较各自完整流程，不能单独归因于“embedding优劣”。
