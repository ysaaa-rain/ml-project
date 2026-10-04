# 2026-10-04 审计证据

对应 docs/08_全方位核查审计_20261004.md，基线远端8aa70af。

- public_audit_summary.json：方向、历史split重建、泄漏计数、代码SHA和核验边界。
- data_profiles.json：原始来源缺失/标签/证据统计，不含序列。
- artifact_hash_checks.tsv：80项manifest引用已匹配；跨run矩阵引用按来源目录定位。
- legacy_holdout_near_duplicates.tsv：4对跨split近重复，2条留出记录；全长81bp任一方向≥73个匹配，无shift/indel认证。

完整文件清单、机器环境和原始对照保存在本地忽略目录，不提交原始序列。可复现只读审计：python -B -m experiments.audit_repository --output-dir results/data_audit/<new_run>。本轮未运行新的motif/ML实验。
