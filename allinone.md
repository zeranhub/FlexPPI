# FlexPPI apo–holo 生物学验证：文献检索与方案审计

检索日期：2026-09-13。范围：只读代码核查和公开文献检索；未启动批量预测、MD 或训练。下述规模、阈值和图表均为建议，不是已经完成的结果。本文保留 paper-search 六源补充检索的全部返回项；它不是系统性文献综述。主要研究建议以另行核对的原始论文为依据。

## 全部六源检索结果

检索词：`apo holo protein flexibility`；年份：2000–2026；每源最多 1 条。限额仅用于检验相关候选，不代表召回完整。

### Semantic Scholar（1 条）

| 题目 | 年份 | 期刊 | 接口引用数 |
|---|---|---|---|
| [Inherent versus induced protein flexibility: Comparisons within and between apo and holo structures](https://www.semanticscholar.org/paper/aab3483059f4888834d26a9de5ede9f1d4a3cfce) | 2019 | PLoS Computational Biology | 83 |

作者：Jordan J. Clark、Mark L. Benson、Richard D. Smith、Heather A. Carlson。原文已通过 [PLOS](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006705) 核对。研究的是不同静态结构间的变化，不是逐残基 MD RMSF/Hχ 真值。其 apo–apo/holo–holo 基线设计值得借鉴。

### OpenAlex（0 条；接口不可用）

运行环境兼容错误；不能解释为没有相关论文。原始报错见文末。

### arXiv（1 条）

| 题目 | 年份 | 载体 | 接口引用数 |
|---|---|---|---|
| [AlphaFold predicts the most complex protein knot and composite protein knots](http://arxiv.org/abs/2207.07410v1) | 2022 | arXiv | 0 |

作者输出：Maarten A. Brems、Robert Runkel、Todd O. Yeates 等。与当前结合响应问题不直接相关，未用作支持性证据。接口返回的 0 不保证现实引用数为零。

### OpenReview（0 条）

本次检索没有返回匹配结果，不代表该领域没有论文。

### Crossref（1 条）

| 题目 | 年份 | 载体 | 接口引用数 |
|---|---|---|---|
| [Faculty Opinions recommendation of Inherent versus induced protein flexibility: Comparisons within and between apo and holo structures](https://doi.org/10.3410/f.734951582.793560475) | 2019 | Faculty Opinions – Post-Publication Peer Review of the Biomedical Literature | 0 |

作者：Reinhard Sterner、Rainer Merkl。这是对 Clark 论文的推荐/评论，不是一项独立研究，不能重复计作证据。

### DBLP（0 条；连接错误）

TLS 连接错误，不能解释为没有相关论文。

### Model knowledge（0 条未核实新增）

本文不添加仅凭记忆、未核实的论文。下文的额外文献均通过原作者、期刊或 PubMed/PMC 页面核对，属于定向网络检索。

## 检索总结

### 1. Overview

六源补充检索共返回 3 条记录，其中 1 篇相关原始研究、1 条针对同一研究的评论和 1 篇不相关预印本。OpenAlex 和 DBLP 失败，部分接口限流，因此不能把此结果作为完整研究版图。另通过 agent-reach/Exa 与原始文献页面查证了下列关键研究。

### 2. Trends

样本太少且检索不完整，不适合推断发表量、热门方向或作者趋势。可明确区分的是：静态构象差异、MD 涨落统计和 NMR 动态代理属于不同证据层次，不能互相替代。

### 3. Key themes

PPI 结合后的局部柔性重分配：Grünberg 2006 是最直接的先例。侧链动力学与主链运动不必同步：Frederick 2007 提供蛋白–肽 NMR 证据。远端响应：Fuentes 2004 的 PDZ2 体系提供可检验位置。界面预组织：Cole 2002 采用计算侧链熵，但其抽出单体不能当成实验 apo。配对分析质控：Clark 2019 的同态结构变化基线与 Wankowicz 2022 的严格结构配对可借鉴。

### 4. Keywords frequency

仅统计 3 条接口记录题目中的不区分大小写关键词，不代表领域热度；评论与原文重复计数已在此披露。

| 关键词 | 题目数 |
|---|---:|
| protein | 3 |
| flexibility | 2 |
| apo | 2 |
| holo | 2 |
| AlphaFold | 1 |

### 5. Most cited by accepted paper

| 排名 | 论文 | 年份 | 接口引用数 |
|---|---|---|---:|
| 1 | Clark et al., Inherent versus induced protein flexibility | 2019 | 83 |

只有这一项确认属于相关、已发表原始研究；不编造其余名次。

### 6. Most cited by first author

| 作者 | 返回记录数 | 接口引用数合计 |
|---|---:|---:|
| Jordan J. Clark | 1 | 83 |
| Reinhard Sterner | 1 | 0 |
| Maarten A. Brems | 1 | 0 |

该表只是技能要求的检索元数据摘要，包含评论和不相关预印本，不宜用于研究者排名。

### 7. Recommendations for reading

| 优先阅读 | 原始研究 | 对 FlexPPI 的用途与边界 |
|---|---|---|
| 1 | [Grünberg et al., Structure 2006, Flexibility and conformational entropy in protein-protein binding](https://pubmed.ncbi.nlm.nih.gov/16615910/) | 17 个 PPI 及未结合组分，提示结合引起柔性重分配而非统一刚化。原指标和采样方案不同，复现趋势不等于复现原数值。 |
| 2 | [Cole & Warwicker, Protein Science 2002](https://pmc.ncbi.nlm.nih.gov/articles/PMC2373749/) | 计算结果提示界面在抽出单体中较预组织；21/25 同二聚体及 8/9 大型异二聚体组分呈较低界面侧链熵/表面积。抽出单体不是实验 apo，且不是 MD Hχ。 |
| 3 | [Frederick et al., Nature 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC4156320/) | CaM 与六种靶肽；跨复合物主链快速运动较相近，侧链运动依赖伙伴。对照必须都是钙饱和状态；不要混入去钙变化。 |
| 4 | [Fuentes et al., JMB 2004](https://pubmed.ncbi.nlm.nih.gov/14698303/) | PDZ2 与肽结合后的远端侧链动态响应；适合作为位置级外部实验检验。NMR 序参量不是 Hχ，须明确时间尺度和探针差异。 |
| 5 | [Wankowicz et al., eLife 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9084896/) | 743 对严格配对晶体数据，界面局部刚化与远端增柔趋势；主要为小分子配体，当前模型不能直接编码这些配体，因此是设计先例，不是可直接整套跑的 PPI 真值集。 |

## 当前模型核查

新五折检查点路径为 `artifacts/performance_first_v2_atlas_hchi_v1`。各折为结构输入绝对状态模型：`enable_response_head=false`、`use_condition_features=false`，输出 RMSF 和 chi_entropy。当前不能沿用旧响应头流程/成绩声称新版差值已验证。应同一检查点分别预测两个状态、还原到原尺度后按相同残基作差。代码依据：`src/flexppi/model/network.py`、`model/geometry.py`、`model/ipa.py`、`data/structures.py`、`inference.py`。

模型使用氨基酸类别和 N/CA/C/O 骨架，包括跨链几何；不直接编码侧链坐标、任意小分子/核酸或实验温度。仅改变静态侧链而骨架与序列相同，模型没有区分信号。RMSF 是链独立对齐后的 Cα 涨落；Hχ 是 rotamer 联合状态熵，不是完整蛋白构象熵。与 NMR、HDX、B 因子不能无条件一一等同。

旧 ITrace 244 个固定坐标伙伴消融体系是上下文敏感性测试，不是新检查点实验 apo–holo 有效性的证据。五折模型一致性只是模型稳健性，不是五次独立生物实验。

## 推荐的研究设计

先固定论文提出的可证伪假设，纳入增柔、刚化、不变和失败案例，不按模型结果筛选“符合文献”的样本。先建立少量独立配对真值集，验证差值本身；再扩展结构预测图谱。可从 [Docking Benchmark 5.5 官方表](https://zlab.wenglab.org/benchmark/) 中已有 bound/unbound 配对起步，随后用 [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) 精确映射同一残基。

每对保留三个输入：真实实验未结合状态 A；从 holo 抽出且坐标不变的目标组分 B；完整 holo C。比较 C−A 为总预测响应；C−B 是固定坐标伙伴上下文消融；B−A 是目标骨架/自身环境差异。三者构成代数拆分，不能解释成严格物理因果或自由能分解。若未结合状态本身为生理多聚体，应保留该状态所需的其他链。

先用约 20–30 个高质量配对作可行性筛查，可在预算允许时对代表性体系补充双状态各 3 条独立 MD；轨迹长度由收敛与目标时间尺度决定，不把固定 100 ns 当成充足保证。定量主证据使用一致处理协议的 ΔRMSF、ΔHχ，NMR 用作独立探针层次的支持。先比较零变化、同一位点 AA-only（差值必为零）、局部接触/埋藏度简单模型，以及适用的结构柔性方法。

严格控制实验构建体、物种/异构体、突变、缺失段、配体/辅因子、寡聚状态、结构解析质量。完整状态独立建模，仅在共有可靠残基上比较；不可为求匹配而预先删除各状态不同的真实结构环境。真正无序到有序、apo 无坐标的片段不应擅自补值作为验证真值。训练集排除不只检查 PDB ID，还应检查双方蛋白家族与同源关系。若使用五折集成，测试集须对全部五折训练集干净；若做 OOF 配对，每一对两端要由同一个未见该家族的模型预测。

统计以蛋白/家族为单位先汇总再比较；需要置信区间时按体系或家族重采样，不能把相关残基或五个折当成独立重复。报告差值 MAE、体系等权相关、在测量/MD 噪声之外的方向准确率、动态变化位置识别和失败比例。增加 apo–apo、同伙伴 holo–holo 的背景变异对照。

## 建议主图

一张四面板可分别给出：界面/近界面/远端的 ΔRMSF 与 ΔHχ 双指标总体效应；每体系界面 ΔHχ 对远端 ΔHχ 的象限图；PDZ 等独立动态实验响应位置与预测的比对；大规模图谱按家族组织的主链/侧链联合变化。MD 定量验证及基线必须有独立主结果图或在主图中占据明确位置，不能只有模型预测之间的相关性。

论文主张上限：无独立动态参照时是“发现与文献相容的预测规律”；差值通过独立 MD 后可称“支持跨状态响应预测”；再对外部实验定位/趋势验证才有更强生物学证据。不能由大规模预测本身宣称完整变构通路、动力学机制或结合自由能已被重建。

## 六源原始输出与错误

为避免隐藏检索失败，保留本次补充查询的完整控制台输出。第一次更宽查询同样只作为探索，没有据其未完整展示的输出新增证据；上文逐项表格对应下面这一次完整可追溯结果。

```text
[open_alex] unavailable (import failed: unsupported operand type(s) for |: 'type' and 'NoneType'); skipping this source.
Rate limited. Waiting 3 seconds...
[dblp] Error: HTTPSConnectionPool(host='dblp.org', port=443): Max retries exceeded with url: /search/publ/api?q=apo+holo+protein+flexibility&format=json&h=1&f=0 (Caused by SSLError(SSLEOFError(8, 'EOF occurred in violation of protocol (_ssl.c:1129)')))
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...

============================================================
  open_alex: 0 papers found
============================================================

============================================================
  dblp: 0 papers found
============================================================

============================================================
  crossref: 1 papers found
============================================================
  [1] Faculty Opinions recommendation of Inherent versus induced protein flexibility: Comparisons within and between apo and holo structures.
      Authors: Reinhard Sterner, Rainer Merkl
      Year: 2019  Citations: 0  Venue: Faculty Opinions – Post-Publication Peer Review of the Biomedical Literature
      URL: https://doi.org/10.3410/f.734951582.793560475


============================================================
  arxiv: 1 papers found
============================================================
  [1] AlphaFold predicts the most complex protein knot and composite protein knots
      Authors: Maarten A. Brems, Robert Runkel, Todd O. Yeates...
      Year: 2022  Citations: 0  Venue: arXiv
      URL: http://arxiv.org/abs/2207.07410v1


============================================================
  openreview: 0 papers found
============================================================

============================================================
  semantic_scholar: 1 papers found
============================================================
  [1] Inherent versus induced protein flexibility: Comparisons within and between apo and holo structures
      Authors: J. Clark, Mark L. Benson, Richard D. Smith...
      Year: 2019  Citations: 83  Venue: PLoS Comput. Biol.
      URL: https://www.semanticscholar.org/paper/aab3483059f4888834d26a9de5ede9f1d4a3cfce


Total: 3 papers from 6 sources.
```

## 2026-09-13 补充检索：Hidden dynamic allostery

检索窗口 2000–2026，每源最多 2 条。以下完整保留六库本次实际返回；不是全部相关文献，也不是系统综述。定向追踪的 13 篇原始论文/预印本及机制分析另见 artifacts/apo_holo_mechanistic_literature_v2_20260913/README_机制验证候选.md。

### Semantic Scholar（2 条）

| 编号 | 题目 | 日期 | 刊物 | 引用数 |
|---|---|---|---|---|
| [S1](https://doi.org/10.1073/pnas.0904492106) | Hidden dynamic allostery in a PDZ domain | 2009-10-27 | PNAS | 342 |
| [S2](https://www.semanticscholar.org/paper/3e4f7d4b31323bd38e757ea29a72a40b283281c5) | Markov State Models and Perturbation-Based Approaches Reveal Distinct Dynamic Signatures and Hidden Allosteric Pockets in the Emerging SARS-Cov-2 Spike Omicron Variant Complexes with the Host Receptor: The Interplay of Dynamics and Convergent Evolution Modulates Allostery and Functional Mechanisms | 2023-08-07 | Journal of Chemical Information and Modeling | 18 |

后端提示原文 `Rate limited. Waiting 3 seconds...` 共 13 次，最终返回以上两条。引用数为此次 API 读数，非质量评分。

### OpenAlex（0 条，后端失败）

`[open_alex] unavailable (import failed: unsupported operand type(s) for |: 'type' and 'NoneType'); skipping this source.`

### arXiv（0 条，后端失败）

`[arxiv] Error: The read operation timed out`

### OpenReview（0 条）

本次未返回记录，不代表所有相关论文都不存在。

### Crossref（2 条）

| 编号 | 题目 | 日期 | 刊物 | 引用数 |
|---|---|---|---|---|
| [C1](https://doi.org/10.1016/j.bpj.2025.11.1496) | BPS2026 – Mutual allostery: A hidden mechanism shaping RNA biology | 2026-02 | Biophysical Journal | 0 |
| [C2](https://doi.org/10.7554/elife.79915.sa0) | Editor's evaluation: Dynamic allostery in substrate binding by human thymidylate synthase | 2022-06-08 | API 未提供 | 0 |

C1 是会议相关条目，C2 是编辑评价；均不当作已筛中的 MD 加实验原创研究。

### DBLP（0 条，后端失败）

`[dblp] Error: Expecting value: line 1 column 1 (char 0)`

### Model Knowledge（0 条）

不补入仅凭记忆、未核实的引用。定向搜索并经原始来源核验的论文不冒充模型记忆来源。

### Overview

查询 Hidden dynamic allostery，2000–2026，六源共返回 4 条记录，其中 2 条是此小样本中的研究论文，另 2 条是会议/编辑评价条目。定向全文追踪另外构建了 PDZ3、TCR、凝血酶、KIX、GIV–Gαi、BAP1 等证据链。

### Trends

样本仅 4 条，且多源失败，不能据此推断年份热度、主流团队或方法增长趋势。

### Key themes

本次返回涵盖 PDZ 动态变构（S1）、病毒变体构象/变构分析（S2）、RNA 变构会议主题（C1）、底物结合动态变构的编辑评价（C2）。每个主题仅一条，无法提供 2–3 篇独立代表文献。

### Keywords frequency

按“标题加已返回摘要中含该词干的记录数”计算，不是词的总出现次数。

| 词干/概念 | 记录数 |
|---|---|
| alloster / 变构 | 4 |
| dynamic / 动态 | 3 |
| binding / 结合 | 2 |
| conformation / 构象 | 1 |
| mutation / 突变 | 1 |

### Most cited by accepted paper

只有两条明确研究论文，故不凑满五条。

| 排名 | 题目 | 年份 | 引用数 |
|---|---|---|---|
| 1 | S1 — Hidden dynamic allostery in a PDZ domain | 2009 | 342 |
| 2 | S2 — Markov State Models and Perturbation-Based Approaches… | 2023 | 18 |

### Most cited by first author

这里只统计上述两篇研究论文，不将评价作者或会议条目与原创研究混算。

| 排名 | 作者 | 集合内论文数 | 总引用数 |
|---|---|---|---|
| 1 | C. Petit | 1 | 342 |
| 2 | Sian Xiao | 1 | 18 |

### Recommendations for reading

先读 [Petit 2009 的实验](https://doi.org/10.1073/pnas.0904492106)，理解四状态侧链动态与量热结果；再读定向补充检索的 [Mostarda 2012 MD](https://doi.org/10.1371/journal.pcbi.1002429)，理解主链时间尺度及长肽作用；接着读 [Murciano-Calles 2014 MD/NMR/ITC](https://doi.org/10.1371/journal.pone.0090030)，引入相反方向突变与功能对照；最后读 [Rangarajan 2018 TCR](https://doi.org/10.1074/jbc.RA118.003832)，建立跨结构域响应案例。这里后三篇来自额外的原始论文定向检索，不是上述六库本次返回。

原始 JSON、来源、结构审计和适用边界均随候选报告保存。此前另一个聚合搜索进程的会话已不可恢复，没有把其未取回的结果记入统计。

## 2026-09-13 机制案例作图补充检索（本次运行，完整记录）

查询：`Ensemble-Based Analysis Dynamic Allostery PSD-95 PDZ3`，2000–2026年，每库最多1篇。本次运行使用paper-search的六库CLI；API未取回有效结果，不把它解释成没有相关论文。另经Exa、PubMed和期刊原站核对的论文已列在机制案例图README中，不混入本次API结果计数。

### Semantic Scholar（0篇）

`[semantic_scholar] unavailable (import failed: No module named 'requests'); skipping this source.`

### OpenAlex（0篇）

`[open_alex] unavailable (import failed: No module named 'requests'); skipping this source.`

### arXiv（0篇）

`[arxiv] Error: The read operation timed out`

### OpenReview（0篇）

`[openreview] Error: openreview not installed. pip install openreview-py`

### Crossref（0篇）

`[crossref] unavailable (import failed: No module named 'requests'); skipping this source.`

### DBLP（0篇）

`[dblp] unavailable (import failed: No module named 'requests'); skipping this source.`

### Model Knowledge（0篇）

没有额外加入仅依赖记忆且未核实的新题名。前述定向原站核实来源另计。

### Overview

本次六库查询为上述PSD-95动态变构主题，2000–2026年，取回0篇；运行失败与空命中不同。任务主体的图形设计仍有已核对的原始论文依据。

### Trends

这次API集合为空，不能推断年份、期刊或作者团队趋势。

### Key themes

空集合不能产生主题聚类。实际作图借鉴来自另外核对的三条原始证据链：局部相互作用与残基动态对应、NMR位点与结构对应、不同伙伴下的局部柔性与实验热稳定性对照。

### Keywords frequency

没有可计算的标题/摘要语料，关键词频次及前5名均不可用，未造数。

### Most cited by accepted paper

本次API集合没有论文，不能列引用排名。

### Most cited by first author

本次API集合没有论文，不能列第一作者引用排名。

### Recommendations for reading

额外定向核对的阅读顺序为[Mostarda2012图2/4](https://doi.org/10.1371/journal.pcbi.1002429)、[Rangarajan2018图4/7/8](https://pubmed.ncbi.nlm.nih.gov/30135211/)、[Peiffer2022图6–8](https://doi.org/10.1371/journal.pcbi.1009977)。分别用于组织局部结构、实验扰动映射、伙伴选择性对照；不是本次失败API返回的命中。

## 2026-09-14：PDZ3 α3 聚焦检索的完整日志与报告

本节为追加记录，不替换前面的检索结果。当前任务已完成 429 个无配体模型输入、五 checkpoint 的实际测试，输出见 output/flexppi_pdz3_alpha3_focus_20260914。下面的失败仅指通用六源 API，不等于没有相关论文；原站定向核对和公开实验表格读取另行成功。

### Semantic Scholar（0 篇）

`[semantic_scholar] unavailable (import failed: No module named 'requests'); skipping this source.`

### OpenAlex（0 篇）

`[open_alex] unavailable (import failed: No module named 'requests'); skipping this source.`

### arXiv（0 篇）

`Rate limited, retrying in 3s...`

`[arxiv] Error: The read operation timed out`

### OpenReview（0 篇）

`[openreview] Error: openreview not installed. pip install openreview-py`

### Crossref（0 篇）

`[crossref] unavailable (import failed: No module named 'requests'); skipping this source.`

### DBLP（0 篇）

`[dblp] unavailable (import failed: No module named 'requests'); skipping this source.`

### Model Knowledge（0 篇）

没有补入仅凭记忆的未核对题名。这里列的是通道返回值，不把另外成功核对的原始来源伪装成 API 命中。

### Overview

查询 `PSD95 PDZ3 alpha3 allostery`，年份 2009-2026，每 API 源最多 2 篇。返回集合为空，原因是依赖/连接失败，不是确认没有匹配文献。

### Trends

空集合不能推断发表趋势、主要期刊或研究团队趋势。

### Key themes

API 空集合不做主题聚类。另行原站核验的工作覆盖 NMR 动态变构、MD 的局部环区及突变响应、实验推断能量景观三个互补主题；这些是定向证据链，不是本次 API 检索统计结果。

### Keywords frequency

没有可计算的标题摘要集合，关键词及前 5 名频次不可用。

### Most cited by accepted paper

没有 API 返回论文，不进行引用量排名，不填补缺失引用数。

### Most cited by first author

没有 API 返回论文，不进行第一作者总引用排名。

### Recommendations for reading

另行核对的阅读路径为 [Petit 等 2009：Hidden dynamic allostery in a PDZ domain](https://doi.org/10.1073/pnas.0904492106)，用于理解最初的动态变构命题；[Mostarda 等 2012：Beyond the Binding Site](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002429)，用于明确局部环区与 V328 方向检验；[Hidalgo-Carcedo 等 2026：Allosteric and energetic remodeling of a PDZ domain by protein domain extensions](https://www.nature.com/articles/s41467-026-69673-w)，用于核对 F337 的 α3 依赖性及公开实验数值。2009 DOI 原站本轮直连报 Internal Error，相关书目信息通过成功打开的原始论文参考文献交叉核对；2012、2026 原站及 2026 Supplementary Data 2/3 已读取。没有将这些原站核对结果称为六库返回。

Agent Reach / Exa 通道也遇到限流；未盲目重复重试，改用原始论文网页和原始公开数据。本轮没有安装新的检索依赖或更新工具，所有失败均不影响已完成的本地模型推断。


# 2026-09-19：FlexPPI 4.2 信息来源与原图检索记录

最终精选图形与实验设计见 [调研记录](docs/flexppi_information_sources_figure_review_2026-09-19.md)。六库通用查询为 protein flexibility sequence structure，2024–2026，每源上限 3。下列命中相关性有限，均未直接用作精选图形证据；精选原图来自另行定向检索与原文核对。未把零返回接口当作证明不存在相关论文。

## semantic_scholar（3 条）

| 标题 | 年份 | 来源提供的会议/期刊 | 来源提供的引用量 |
|---|---|---|---|
| [RNA-EFM: energy-based flow matching for protein-conditioned RNA sequence-structure co-design](https://www.semanticscholar.org/paper/f9f3b794159dc5cf44b7e7760e7f323e8b4f8e3f) | 2025 | bioRxiv | 6 |
| [Extending the classical sequence-structure-function paradigm through protein dynamics and context-dependent behavior.](https://www.semanticscholar.org/paper/dcfd266c94ddedab72c277e36bdb30aeabf4173d) | 2026 | FEBS Letters | 1 |
| [Predicting high-risk clinical missense variants of SMARCB1 in rare neurogenetic disorder schwannomatosis (nerve tumor) through sequence, structure, and molecular dynamics analyses](https://www.semanticscholar.org/paper/15e6c78c67933c095ac0bced0dd8bd0d2c3bacc9) | 2025 | Neurogenetics | 2 |

## open_alex（0 条）

无返回；具体接口失败见下。

## arxiv（3 条）

| 标题 | 年份 | 来源提供的会议/期刊 | 来源提供的引用量 |
|---|---|---|---|
| [Mass Balance Approximation of Unfolding Improves Potential-Like Methods for Protein Stability Predictions](http://arxiv.org/abs/2504.06806v1) | 2025 | arXiv | 0 |
| [Quantum Simulation of Protein Fragment Electronic Structure Using Moment-based Adaptive Variational Quantum Algorithms](http://arxiv.org/abs/2601.00656v1) | 2026 | arXiv | 0 |
| [idMotif: An Interactive Motif Identification in Protein Sequences](http://arxiv.org/abs/2402.05953v1) | 2024 | arXiv | 0 |

## openreview（0 条）

无返回；具体接口失败见下。

## crossref（3 条）

| 标题 | 年份 | 来源提供的会议/期刊 | 来源提供的引用量 |
|---|---|---|---|
| [Reviewer #2 (Public review): Interpretable protein-DNA interactions captured by structure-sequence optimization](https://doi.org/10.7554/elife.105565.3.sa1) | 2025 | not provided | 0 |
| [Reviewer #3 (Public review): Interpretable protein-DNA interactions captured by structure-sequence optimization](https://doi.org/10.7554/elife.105565.3.sa2) | 2025 | not provided | 0 |
| [Reviewer #3 (Public review): Interpretable Protein-DNA Interactions Captured by Structure-Sequence Optimization](https://doi.org/10.7554/elife.105565.1.sa0) | 2025 | not provided | 0 |

## dblp（0 条）

无返回；具体接口失败见下。

## Model knowledge

未补入未在线核验的记忆条目。

## 汇总

总计 9 条原始命中；没有形成主题一致的论文集合，不从这批结果计算趋势、关键词频次或被引排名，以免将通用检索噪声表述为领域结论。阅读路径为 Flexpert 图 5/6 → BackFlip 图 A.2/表 A.1 → Zhang 与 Liu 图 2/3 → FastProtFlex 图 7/8。精确链接、实看状态和局限见上述调研记录。

## 接口失败原文

```text
[dblp] unavailable (import failed: No module named 'requests'); skipping this source.
[openreview] Error: openreview not installed. pip install openreview-py
[semantic_scholar] unavailable (import failed: No module named 'requests'); skipping this source.
[crossref] unavailable (import failed: No module named 'requests'); skipping this source.
[open_alex] unavailable (import failed: No module named 'requests'); skipping this source.
[open_alex] unavailable (import failed: unsupported operand type(s) for |: 'type' and 'NoneType'); skipping this source.
Rate limited. Waiting 3 seconds...
[dblp] Error: Expecting value: line 1 column 1 (char 0)
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Rate limited. Waiting 3 seconds...
Jina Reader: AuthenticationRequiredError: You have been blocked from performing anonymous queries due to bad network reputation (AS30058). Please authenticate.
Agent Reach check-update: 无法检查更新（GitHub API 速率限制，已重试 3 次）
```
