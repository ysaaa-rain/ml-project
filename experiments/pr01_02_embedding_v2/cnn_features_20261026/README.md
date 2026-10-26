# 老师文档推荐路线：CNN局部特征与卷积核激活提取

本轮已完成：[结果报告](REPORT.md)。36次CNN训练、18随机模型、18次OH48聚类；弱模式匹配5候选F1为0.277/0.547，未通过小试门槛；混合链约0.894/0.871，有局部正信号。评分后边界诊断提示1–2bp错位，未据此改模型或分数。真实Holdout未用。

新方法家族，非BERT调参。依据 `motif-embedding-methods.pdf` pp.1–4及 `motif-embedding-overview.pdf` p.1，采用One-hot→局部卷积24维特征→高激活实例→PWM→固定扫描评价；按文档分层区分表示、提取和验证。方法选择与SAE适用性核对见 [老师材料记录](../closeout_20261026/teacher_methods.md)。这是借鉴DeepBind思想的自建小模型，未使用DeepBind发布权重，不声称完整复现。

冻结方案见 `protocol.json`：3个新数据种子×2算法种子×3条件。CNN_LABEL、打乱标签CNN、未训练随机CNN、完整One-hot48背景拒绝四对照。实际计划36个CNN训练、18个随机初始化模型、18次KMeans，72组PWM池，FULL/CAP5/匹配预算共180组评价集合。

CNN_LABEL使用Discovery序列级正负标签，未使用植入位置/方向/参考；与OH48的无监督发现不同，任何提升不能单独归为embedding优势。两个方向共享卷积权重，额外引入方向归纳偏好，也不等同于只替换特征。分类AUROC单列，不能冒充motif定位。全部拟合与候选锁定后再统一评估，真实Holdout关闭。固定最后80轮模型，不按Development/Test挑checkpoint。零候选、预算不足与阴性结果保留。

```powershell
python -m pytest experiments/pr01_02_embedding_v2/cnn_features_20261026/test_pipeline.py -q
python -m experiments.pr01_02_embedding_v2.cnn_features_20261026.pipeline fit
python -m experiments.pr01_02_embedding_v2.cnn_features_20261026.pipeline score
python -m experiments.pr01_02_embedding_v2.cnn_features_20261026.qa
python -m experiments.pr01_02_embedding_v2.cnn_features_20261026.diagnose
python -m experiments.pr01_02_embedding_v2.cnn_features_20261026.verify_artifacts
python -m experiments.pr01_02_embedding_v2.cnn_features_20261026.report
```

拟合拒绝覆盖artifacts。人工序列、神经权重、聚类对象、激活和实例留本地；公开源码、协议、训练曲线/汇总、候选排名/预算、图表、核验和哈希。实际状态以summary/manifest.json和REPORT.md为准。预注册的pipeline.verify曾因排名来自原精度PWM而文本保存8位小数导致浮点字段检查失败；verify_artifacts从实例重建原精度排名，保留差异审计，不改变冻结方法和评分，详见报告。
