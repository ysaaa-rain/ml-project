# 甲的工作（数据验收与说明）

> 本文件是任务总纲：上半部分说明数据长什么样，下半部分是八步工作计划，每步标明状态并链到对应产出文件。

## 交付物索引

| # | 文件 | 内容 | 对应步骤 |
| --- | --- | --- | --- |
| 1 | [甲_验收记录_20261006.md](甲_验收记录_20261006.md) | 文件完整性与数量核对证据（含 48 项明细） | 第 1、2 步 |
| 2 | [甲_σ组可用性清单_20261006.md](甲_σ组可用性清单_20261006.md) | 19 个 σ 组哪些能做正式分析 | 第 3 步 |
| 3 | [甲_排除记录汇总_20261006.md](甲_排除记录汇总_20261006.md) | 395 条排除记录，7 类原因 | 第 4 步 |
| 4 | [甲_数据字典_20261006.md](甲_数据字典_20261006.md) | 字段说明、口径说明、5 条已知限制 | 第 5 步 |
| 5 | [甲_许可状态表_20261006.md](甲_许可状态表_20261006.md) | 各数据来源的许可状态 | 第 6 步 |

目录位置（绝对路径）：

```text
D:\AI\ml-project\docs\数据三人小小组结果
```

五个产出文件的文件名：

```text
甲_验收记录_20261006.md
甲_σ组可用性清单_20261006.md
甲_排除记录汇总_20261006.md
甲_数据字典_20261006.md
甲_许可状态表_20261006.md
```

## 进度

| 步骤 | 状态 |
| --- | --- |
| 第 0 步 准备 | ✅ 已完成 |
| 第 1 步 文件完整性核对 | ✅ 已完成 |
| 第 2 步 数量核对 | ✅ 已完成 |
| 第 3 步 σ 分组可用性清单 | ✅ 已完成 |
| 第 4 步 排除记录汇总 | ✅ 已完成 |
| 第 5 步 数据字典 | ✅ 已完成 |
| 第 6 步 许可状态表 | ✅ 已完成 |

每步的详细说明见下方「三、我的工作（甲）」。

---

## 一、data 里面有什么

data （根目录下的data文件夹）分三层，就像**原料 → 半成品 → 成品**：

| 目录 | 是什么 | 有多少 |
| --- | --- | --- |
| `data/raw/` | **原料**：从官网和云盘直接下载的文件，一个字节都没改过 | 1,105 个文件 |
| `data/interim/` | **半成品**：从原料里解析出来的中间表 | 3 个文件 |
| `data/processed/` | **成品**：能直接拿去做分析的数据 | 315 个文件 |

**原料（data/raw）里有三样东西：**

| 内容 | 说明 |
| --- | --- |
| `tju_pan_promoter/` | 云盘那套，六个物种 × 5 个 CSV，还有强度数据，36 个文件 |
| `regulondb_*.gff3` + `regulondb_*.fna` | 大肠杆菌的启动子清单 + 完整基因组 |
| `dbtbs_v4.1_20261005/` | 枯草芽孢杆菌的 1,062 个网页 + 参考基因组 |

**成品（data/processed）里有两个目录，但只有一个算数：**

| 目录 | 状态 |
| --- | --- |
| **`pr01_02_data_v2/`** | ✅ **现在的正式数据，我们要用的就是它** |
| `g0_rebuild_20261004/` | ⚠️ 旧的 RegulonDB 版本，现在降级成"辅助资料" |

## 二、`pr01_02_data_v2` 里面到底是什么
### 一、整体地图
一句话：8 个来源文件夹 + 8 个总账文件。

#### 8 个来源文件夹
```text
pr01_02_data_v2\
│
├─ 8 个全局文件（直接放在最外层）
│
├─ tjupan_bacillus_subtilis\main\      ← 枯草芽孢杆菌（TJU 云盘）
├─ tjupan_baumannii\main\              ← 鲍曼不动杆菌
├─ tjupan_bradyrhizobium\main\         ← 慢生根瘤菌
├─ tjupan_diphtheria\main\             ← 白喉杆菌
├─ tjupan_escherichia_coli\main\       ← 大肠杆菌
├─ tjupan_staphylococcus\main\         ← 金黄色葡萄球菌
├─ regulondb_ecoli\                    ← RegulonDB 数据库
│   ├─ core\                            （核心子集）
│   ├─ extended\                        （扩展子集）
│   └─ core\sigma\Sigma70\ 等            （按 sigma 再分）
└─ dbtbs_bsub\core\                    ← DBTBS 数据库
    └─ core\sigma\SigA\ 等              （按 sigma 再分）
```
#### 最外层那8 个总账文件
这 8 个是**全局汇总表**，把下面所有来源的信息合并在一起：

| 文件 | 是什么 | 大小 |
| --- | --- | ---: |
| `manifest.json` | **总台账**：所有文件的哈希、数量、参数、已知限制 | 166 KB |
| `sequence_master.tsv` | 所有序列的身世表 | 7.2 MB |
| `source_record_mapping.tsv` | 所有记录追溯到原始文件 | 7.9 MB |
| `sigma_associations.tsv` | 所有 sigma 关联关系 | 3.4 MB |
| `tju_sigma_annotation_candidates.tsv` | TJU 的 sigma 注释候选 | 130 KB |
| `leakage_links.tsv` | 同源关系表（哪些序列太像） | 468 KB |
| `exclusions.tsv` | 被排除的记录及原因 | 122 KB |
| `sigma_group_summary.tsv` | 每个 σ 组多少条、能不能用 | 977 B |

**平时你只需要看最后两个**，其他是备查用的。

### 二、第 2 层：`main` / `core` 里面有什么

打开任意一个层文件夹，你会看到 12 个文件（TJU）或 11 个（RegulonDB/DBTBS）。命名规律是：

**`<划分>.<类型>.fasta`**

- **划分**：`discovery`（练习题）/ `development`（模拟考）/ `holdout`（期末考）/ `external_validation`（外部验证）
- **类型**：`positive`（启动子）/ `natural_control`（天然对照）/ `dinucleotide_null`（打乱对照）

拿大肠杆菌举例，实际就是这个：

```text
tjupan_escherichia_coli\main\
├─ discovery.positive.fasta           1,143 条  ← 找规律用这个
├─ discovery.natural_control.fasta             ← 对照
├─ discovery.dinucleotide_null.fasta           ← 对照
├─ development.positive.fasta           163 条  ← 调参数用这个
├─ development.natural_control.fasta
├─ development.dinucleotide_null.fasta
├─ holdout.positive.fasta               329 条  ← 最后验证用这个
├─ holdout.natural_control.fasta
├─ holdout.dinucleotide_null.fasta
├─ sequence_master.tsv                          ← 这个来源的身世表
├─ sigma_associations.tsv                       ← 这个来源的 sigma 关联
└─ source_record_mapping.tsv                    ← 这个来源的溯源表
```

**为什么 RegulonDB/DBTBS 只有 11 个文件？** 因为它们少了 `natural_control` 那一类——RegulonDB 和 DBTBS 数据库里**只有启动子，没有"非启动子"的对照序列**。TJU 云盘里带了 label=0 的负样本，所以有天然对照。

### 三、第 3 层：`sigma` 里面

只有 RegulonDB 和 DBTBS 有这一层，因为**TJU 数据没有 sigma 标签**。

```text
regulondb_ecoli\core\sigma\
├─ Sigma70\   （6 个文件）
├─ Sigma24\
├─ Sigma32\
├─ Sigma38\
├─ Sigma54\
└─ Sigma28\

dbtbs_bsub\core\sigma\
├─ SigA\  SigB\  SigD\  SigE\  SigF\  SigG\  SigH\
└─ SigI\  SigK\  SigL\  SigM\  SigW\  SigX\
```

每个 sigma 文件夹里是 6 个文件（少了 external_validation 和 natural_control）：

```text
Sigma70\
├─ discovery.positive.fasta           583 条
├─ discovery.dinucleotide_null.fasta
├─ development.positive.fasta          67 条
├─ development.dinucleotide_null.fasta
├─ holdout.positive.fasta             167 条
└─ holdout.dinucleotide_null.fasta
```

**有些 sigma 文件夹只有 2 个或 4 个文件**（比如 SigI 只有 discovery），那是因为**那一组样本太少，根本不够分三份**——文件少本身就是个信号：这组别指望做正式分析。

## 三、我的工作（甲）

**职责：数据验收与说明；数据处理由队长完成。**

下面每一步都按 **干什么 → 为什么 → 怎么做 → 做到什么程度** 写。

### 第 0 步：准备（30 分钟）✅ 已完成

**干什么**：确认自己在正确的数据上工作。

**为什么**：如果搞错版本，后面全部白做。仓库这两天换过三次数据口径，一个目录名不对就全盘错。

**怎么做**：

```powershell
git log --oneline -3          # 确认仓库是最新的
git status --short            # 确认没有意外改动
Get-ChildItem data\processed  # 应该看到 pr01_02_data_v2 和 g0_rebuild_20261004 两个目录
```

**做到什么程度**：能明确回答三个问题——当前数据版本叫什么、主数据是哪个来源、这轮要分析哪批。

答案：版本是 `pr01_02_data_v2`；主数据是 **TJU 云盘六物种**；RegulonDB/DBTBS 用于 sigma 标注和独立来源验证。

### 第 1 步：文件完整性核对 ✅ 已完成

**干什么**：确认 210 个文件一个不少、一个没坏。

**为什么**：数据是通过网盘或 U 盘传的，传输过程可能丢文件、可能被压缩工具改坏。`manifest.json` 里存了每个文件的 SHA256，这是唯一能证明"文件和我拿到的原版一模一样"的办法。

**怎么做**：跑下面这段脚本，它会逐个文件算哈希比对。

```powershell
$root = "data\processed\pr01_02_data_v2"
$m = Get-Content -Raw "$root\manifest.json" | ConvertFrom-Json
$ok = 0; $bad = @(); $missing = @()
foreach ($p in $m.output_sha256.PSObject.Properties) {
    $f = Join-Path $root $p.Name
    if (-not (Test-Path $f)) { $missing += $p.Name; continue }
    $actual = (Get-FileHash $f -Algorithm SHA256).Hash.ToLower()
    if ($actual -eq $p.Value.ToLower()) { $ok++ } else { $bad += $p.Name }
}
Write-Host "一致: $ok  缺失: $($missing.Count)  不匹配: $($bad.Count)"
$missing | ForEach-Object { "  缺失 $_" }
$bad     | ForEach-Object { "  不匹配 $_" }
```

**已经查过的结果**：

- manifest 里登记了 **210 个文件**
- 磁盘上实际有 **211 个文件**
- 差的这 1 个是 `manifest.json` 自己（它没法给自己算哈希），**属于正常情况**
- manifest 里登记的 210 个文件，逐个算了 SHA256 比对，一个不差、一个没坏。磁盘上多出来的那 1 个是 manifest.json 自己（它没法给自己算哈希），属于正常。

### 第 2 步：数量核对 ✅ 已完成

**干什么**：数每个 FASTA 里有多少条序列，把 manifest 里登记的每个数字，跟磁盘上文件里实际的序列条数对一遍

**为什么**：哈希只能证明"文件没被动过"，**不能证明"这版数据是对的"**。数量是判断版本对不对最直接的办法；数量对不上，说明拿到的可能是另一版。

**怎么做**：

```powershell
Get-ChildItem data\processed\pr01_02_data_v2 -Recurse -Filter '*positive.fasta' |
  Sort-Object FullName |
  ForEach-Object {
    $n = (Select-String -Path $_.FullName -Pattern '^>' -AllMatches).Count
    "{0}{1}{2}" -f $_.FullName.Replace('D:\AI\ml-project\data\processed\pr01_02_data_v2\',''), "`t", $n
  }
```


### 第 3 步：σ 分组可用性清单（1 小时）✅ 已完成

**干什么**：把 `sigma_group_summary.tsv` 翻译成"哪些组能做正式分析"。

**为什么**：样本量不足的组做统计，结果不可靠。

**怎么做**：打开 `data\processed\pr01_02_data_v2\sigma_group_summary.tsv`，只有 19 行。判定规则：

- **总数 ≥50 且 discovery ≥30 且 holdout ≥20** → 可做正式结论
- **20–49 条** → 只能探索
- **<20 条** → 只能描述

整理结果：

| 等级 | σ 组 | 含义 |
| --- | --- | --- |
| ✅ 可正式分析 | **Sigma70**（817 条）、**SigA**（249 条） | 能下结论、能写进报告 |
| ⚠️ 只能探索 | Sigma24、Sigma32、Sigma38、Sigma54（RegulonDB）；SigB、SigE、SigG、SigK（DBTBS） | 可以做，但结论要加"探索性" |
| ❌ 只能描述 | Sigma28、SigD、SigF、SigH、SigI、SigL、SigM、SigW、SigX | 只报样本量，不做统计 |

**做到什么程度**：19 组逐组列出总数、三个划分数量和可用性判定。

**产出**：一份 σ 组可用性清单（含每组的三个划分数量）。

**已完成，清单见 → [甲_σ组可用性清单_20261006.md](甲_σ组可用性清单_20261006.md)**

结论：19 个 σ 组中，**只有 Sigma70（817 条）和 SigA（249 条）够做正式结论**；8 组只能探索；9 组只能描述。

### 第 4 步：排除记录汇总（1 小时）✅ 已完成

**干什么**：统计有哪些记录被排除了、为什么。

**为什么**：验收明确要求"排除原因可追溯"。报告里必须写清楚"总共排除了多少条、分别是什么原因"，**这一节不能省**。而且这直接关系到结论可信度——排除很多说明数据质量有限。

**怎么做**：打开 `data\processed\pr01_02_data_v2\exclusions.tsv`（122 KB），按 `reason` 列统计。

已数好的结果：

| 原因 | 条数 | 这意味着什么 |
| --- | ---: | --- |
| 多个/冲突的 TSS 标记 | 169 | 网页上有好几个可能的起点，无法确定用哪个 |
| Location=ND 或坐标无法解析 | 112 | 网页没给可用的位置信息 |
| 缺序列或无效 TSS | 92 | 原始记录本身就残缺 |
| Location 端点不一致 | 15 | 相对位置和绝对位置对不上 |
| 非 81nt / 非 ACGT | 4 | 长度不对或含非法字母 |
| 划分记录冲突 | 2 | 同一批数据的标签自相矛盾 |
| 参考窗口不匹配 | 1 | 从参考基因组切出来的和记录对不上 |

**做到什么程度**：能说出"总共排除多少条、每类多少条、每类为什么"，并能解释**为什么这些排除是合理的、不是随便删的**。

**产出**：一张排除原因表 + 一段说明文字。

**已完成，汇总见 → [甲_排除记录汇总_20261006.md](甲_排除记录汇总_20261006.md)**

结论：共排除 **395 条**、7 类原因；其中 DBTBS 占 296 条（全部是坐标类问题）。

### 第 5 步：写数据字典（半天到一天，**最核心的活**）✅ 已完成

**干什么**：给 `pr01_02_data_v2` 写一份 README。

**为什么**：这个目录有 211 个文件、多张表和十几个字段，原本没有任何字段说明；验收清单明确要求提供字段说明。

**怎么做**：新建 `data\processed\pr01_02_data_v2\README.md`，必须包含这七节。

**① 数据来源与版本**

- 8 个来源分别是什么数据库/数据集
- 版本号、访问日期、SHA256 在哪查

**② 目录结构地图**

- 用树状图画出 `来源/层/文件` 的规律
- 说明 `main` / `core` / `extended` 的区别

**③ 三个划分的含义**

| 文件 | 用途 | 禁令 |
| --- | --- | --- |
| discovery | 找规律 | 只能用它发现 motif |
| development | 调参数 | 只能用它定阈值 |
| holdout | 最终验证 | **只能用一次，不许调参** |

**④ 三类序列的含义**

| 类型 | 是什么 |
| --- | --- |
| positive | 启动子正样本 |
| natural_control | 天然负样本（只有 TJU 有） |
| dinucleotide_null | 电脑打乱生成的对照 |

**⑤ 每张表的每一列**（重头戏）

以 `sequence_master.tsv` 为例，17 列都要写：

| 列名 | 含义 |
| --- | --- |
| sequence_id | 唯一编号 |
| species | 物种 |
| source | 来自哪个数据库 |
| sequence | DNA 序列（81bp） |
| length | 长度 |
| local_tss_index | TSS 在序列里的位置（TJU 固定是 61） |
| split | discovery / development / holdout |
| tier | core / extended |
| label | 正样本还是对照 |
| orientation | 基因组方向 |
| 其余列 | 同样逐列写清楚 |

**⑥ 五个容易踩的坑**（必须写，否则一定有人搞错）

1. `split` 三种划分不能混用，尤其 holdout 不许调参
2. `tier` 的 core 和 extended 证据等级不同，不能混着算
3. `orientation`（基因组方向）和 `input_orientation`（输入序列方向）**是两回事**
4. 三个窗口不能当三倍样本
5. 多 sigma 的记录不能只取第一个标签

**⑦ 已知限制**（`manifest.json` 里现成的 5 条，翻译成人话）

- TJU 的 TSS 对齐是**采信上游声明**，不是逐条独立核验的
- DBTBS 用的参考基因组版本和 2005 年原版未必一致
- promoter 证据强 ≠ sigma 关联证据强
- 小 σ 组必须先查样本量再分析
- 相似度阈值是操作规则，不等于生物学上的独立

**做到什么程度**：字段含义、取值、口径说明和已知限制全部写明。

**已完成，字典见 → [甲_数据字典_20261006.md](甲_数据字典_20261006.md)**

内容：数据来源、目录结构、三个划分、三类序列、8 张表逐列说明、字段口径说明、5 条已知限制、关键数字速查。

### 第 6 步：许可状态表（1 小时）✅ 已完成

**干什么**：列出每个数据来源的许可状态。

**为什么**：这是 M1 三项目前**仍未通过**的要求之一。原始数据和含序列的文件现在不能推 GitHub，但报告和图表能不能公开、公开到什么程度，取决于这个。

**怎么做**：查三个来源各自的许可条款，做成表：

| 来源 | 许可状态 | 能不能公开 |
| --- | --- | --- |
| TJU / PromLoop | Apache-2.0（已核验 30/30 文件） | 汇总结果可以 |
| RegulonDB | **未确认**（软件仓库是 Apache，不等于数据） | 待定 |
| DBTBS | **未确认**（页面有版权声明） | 待定 |

**做到什么程度**：**不要写"已确认"**——查不到就写"未确认，需要进一步核实"。诚实记录空缺比编一个答案安全得多。

**已完成，状态表见 → [甲_许可状态表_20261006.md](甲_许可状态表_20261006.md)**

结论：TJU/PromLoop 许可明确（Apache-2.0）；RegulonDB 与 DBTBS 的数据再分发条款**未确认**。

