# FlexPPI 候选分析标的调研：除 TCR-pMHC 外有公开 MD 数据的蛋白-蛋白复合物体系

> 调研日期：2026-08（本轮）
> 目的：为 FlexPPI 寻找 TCR-pMHC 之外、具备公开可下载 MD 轨迹、且与 FlexPPI 的 apo/holo 结合响应框架兼容的蛋白-蛋白相互作用（PPI）复合物标的。
> 调研方式：agent-reach Exa 搜索 + 多路 web 检索 + Zenodo/Harvard Dataverse 官方记录核实 + 论文原文交叉验证；两个后台子代理分别深挖"大型数据集"与"Zenodo 单项轨迹"，全部下载链接均来自官方记录页实测。

## 下载状态（2026-08-25 更新）

**第一梯队 Zenodo 6638504（8 复合物，20.2 GB 压缩）已全部下载到本机并解压**：

- 位置：`E:\data\protein_md_datasets\zenodo_6638504_ppi_md\`
  - 8 个 `*.tar.gz`（MD5 全部校验通过，与 Zenodo 官方 checksum 一致）
  - `extracted\<PDB>\` 解压目录：每个复合物含 `A/`（游离链 A）、`B/`（游离链 B）、`complex/without_water/`（去水复合物）、`complex/water/`（含水复合物），均为 `start.gro` + `md.xtc`
- 体积：压缩包 20.2 GB，解压后约 40.4 GB（含全部文件）
- 注意：轨迹为 gro+xtc 格式（无 tpr 拓扑文件），对接 FlexPPI 时需用 PDB 结构做输入、xtc 算标签（与 DynaRepo 管线相似但需适配 gro 拓扑）

**第二梯队 ThermoPCD 2.0 已部分下载（策略：全量 550 GB 超出 E 盘剩余空间，先取"纯标签"部分）**：

- 位置：`E:\data\protein_md_datasets\ThermoPCD_labels\`
- **235 个纯标签数据集条目全部下载完成**（1.5 GB，3787 个文件，229 个目录无空目录）：
  - 每条目含 0/50/100(200/500)ns 的 PDB 快照 + `rmsf.xvg`、`rmsd.xvg`、`sasa.xvg`、`rg.xvg`、`hbond.xvg` 及 `.tab` 版标签（共 238 个 rmsf.xvg）
  - 覆盖 4ZQK(PD-1/PD-L1)、6UMT(PD-1/PD-L2)、5DK3(pembrolizumab)、7DV4(CTLA-4)、4G3Y(TNFα/infliximab)、5FB8(IL-16)、1ZTX(西尼罗)、2R29(登革)、6OC7、6OY4、5HDQ 等体系
  - 数据源：Harvard Dataverse ThermoPCD 分馆（[链接](https://dataverse.harvard.edu/dataverse/ThermoPCD)），全部免注册下载，每文件独立 DOI
- **未下载（超出 E 盘空间，需另行决策）**：
  - 551 个"轨迹条目"（含 GROMACS tpr/trr/xtc 原始文件，约 523 GB）——其中 nvt.trr/npt.trr 平衡中间文件占大头（约 1.1 GB/条目），若只取生产轨迹 `dynamic-nopbc.xtc` + PDB 快照可降至约 400-500 MB/条目
  - 原始 .xtc 全轨迹还需联系作者（官网说明："To obtain the .xtc files etc, use the contact details"）
- 索引文件：`E:\data\protein_md_datasets\thermopcd_index.csv`（1061 条目 × doi/name/files/bytes）

---

## 0. 结论速览（TL;DR）

| 优先级 | 标的 | 体系 | MD 数据 | apo/holo | 可用性 |
|---|---|---|---|---|---|
| ★★★ | **本地 DynaRepo 复合物** | PDBbind 406 + DynaBench 210 + test53 39 复合物（1B2S、1AK4、1CGI、1WEJ、2PCC…） | 每体系 3×500 ns 全原子，本地已有 structure.pdb + trajectory.xtc | ✅ test53 已建 37 对 apo/holo 配对 | 本机 D:\data\dynarepo_non_nucleic_full，**零下载成本，立即可用** |
| ★★★ | **Zenodo 6638504 八复合物集** | barnase-barstar(1BRS)、RNaseSa-barstar(1AY7)、colicin-免疫蛋白(1EMV)、LIF-gp130(1PVH)、Vav/GRB2 SH3(1GCQ)、cyclophilin-HIV衣壳(1AK4)、proteinaseB-抑制剂(3SGB)、ubiquitin-连接酶(2OOB) | 每复合物 unbound(A/B)+complex，含 water/无 water 两套 | ✅ **唯一严格 apo+holo 配对公开集** | [Zenodo 6638504](https://zenodo.org/records/6638504)，CC-BY-4.0 免费 |
| ★★★ | **ThermoPCD 2.0（抗体-抗原）** | ~85–90 个唯一 PDB 体系：PD-1/PD-L1(4ZQK)、PD-1/PD-L2(6UMT)、CTLA-4、HER2、pembrolizumab(5DK3) 等免疫复合物 | 310–313 K × 100/200/500 ns × 1–3 种子，GROMACS 原始文件+每 50 ns PDB 快照，**自带 RMSF/SASA/Rg/RMSD/H-bond 标签** | ⚠️ 部分（~6–7 个游离抗体/抗原条目，无系统性配对） | [Harvard Dataverse](https://dataverse.harvard.edu/dataverse/ThermoPCD)，免注册，1061 个数据集条目 |
| ★★ | **AIFM1–AK2** | AIFM1 二聚体（PDB 4bur/9FL7）+ AK2 | 全原子 CHARMM36，2×2 μs × 2 体系 = 8 μs，0.5 ns 快照 | ⚠️ 状态对比：FAD vs NAD+FAD（NADH 结合响应） | [Zenodo 15077557](https://zenodo.org/records/15077557)，2×2.3 GB |
| ★★ | **RIPK1–VHL（PROTAC）** | RIPK1 激酶域(4NEU) + VHL(5NVW) | 全原子 Amber14SB，~420 μs（Folding@Home 31 种子）+ 4 个 PROTAC 三元状态 | ⚠️ 无 apo；含 encounter 复合物与三元状态 | [Zenodo 18155129](https://zenodo.org/records/18155129)，11.8+4.5 GB |
| ★★ | **HADDOCK MD-scoring（25 复合物）** | Docking Benchmark 5 的 25 个复合物（1BJ1…3F1P 等） | 全原子 CHARMM36m，每复合物 100 ns；8 HADDOCK 模型×2 重复+晶体×4 重复 | ⚠️ native vs non-native 模型对比（无 apo） | [Zenodo 4629895](https://zenodo.org/records/4629895)，11.3 GB |
| ★ | **Protein A–IgG** | SpA 五域(链E) + IgG1 Fc(链B/D) | 全原子 CHARMM36m，2×400 ns | ⚠️ pH 7.3 vs 5.6 条件对比 | [Zenodo 17039261](https://zenodo.org/records/17039261)，21 MB |
| ★ | **SARS-CoV-2 RBD–ACE2（SIRAH）** | RBD+ACE2-B0AT1（PDB 6M17） | 粗粒化 SIRAH，4 μs | ❌ 无 | [Zenodo 3942566](https://zenodo.org/records/3942566)（CG，需 backmap） |
| ☆ | **DynaPPI / DynPro** | 大规模 PPI / 单链 MD 提案 | 未发布 | — | NeurIPS 2025 数据提案，**跟踪发布**（[arXiv 2608.10435](https://arxiv.org/html/2608.10435)） |

---

## 1. 先纠正一个前提：本地数据并非只有 TCR-pMHC

用户前提是"目前只用 TCR-pMHC 做具体分析"。核对了 `artifacts/data/dynarepo_states.jsonl` 与 `artifacts/data/dynarepo_pairs.jsonl`：

- **本地 DynaRepo 共 922 条带 MD 轨迹的状态**，分组为：
  - `pdbbind_complex`：406 个 PPI 复合物（如 1AK4 = cyclophilin/HIV 衣壳，1B2S = 蛋白-蛋白界面突变响应体系，1PVH = LIF/gp130）
  - `dynabench_complex`：210 个（如 1CGI = 胰凝乳蛋白酶-抑制剂，2PCC = 细胞色素 c 过氧化物酶复合物，1WEJ = IgG1 Fab-细胞色素 c 抗体-抗原）
  - `antigen_single`：215 个抗原单体（1WEJ、2DD8、3HI6、5ZUF、6HER…，来自 SAbDab）
  - `test53_single`：52 个、`test53_complex`：39 个（transient 结合基准）
- 每条含 `structure.pdb` + `trajectory.xtc`（约 86 MB/条）+ `topology.tpr`，轨迹按 DynaRepo 协议为 **3×500 ns、310 K、CHARMM36m**，标签已用 100 ns 窗口计算。
- `dynarepo_pairs.jsonl` 已有 **37 对通过 0.90 一致度/0.80 覆盖度过滤的 apo/holo 配对**（1F6M、1I2M、1NCA、1WDW、1ZHH、2AYO、2C0L…），只是这些配对目前只用于训练/微调，尚未像 TCR-pMHC 那样产出深度案例分析报告。
- 本地 ATLAS 1,938 条为单链（无伙伴环境），不在本次"PPI 复合物标的"范围内。

**含义**：若要快速扩大深度分析覆盖面，第一个动作不是下载新数据，而是从本地 DynaRepo 的 test53/DynaBench/PDBbind 复合物里挑 3–5 个代表体系做与 7L1D/7RRG 同规格的 PyMOL 着色 + 响应分析。它们自带真值轨迹（TCR-pMHC 深度分析反而没有 MD 真值），可以直接与 FlexPPI 预测 RMSF 对照。

---

## 2. 严格 apo/holo 配对：Zenodo 6638504（首选新数据源）

Martin & Frezza, *A dynamical view of protein-protein complexes: studies by molecular dynamics simulations*, **Front. Mol. Biosci. 2022**, DOI [10.3389/fmolb.2022.970109](https://doi.org/10.3389/fmolb.2022.970109)；数据 [Zenodo 6638504](https://zenodo.org/records/6638504)（CC-BY-4.0）。

- **8 个二元复合物**，每个包含 **A/B（unbound 游离链）+ complex（结合态）** 轨迹，且 complex 再分 `without_water` 与 `water` 两套 —— 这是本轮调研中**唯一严格满足 FlexPPI "apo 受体 + holo 复合物" 成对输入要求**的公开数据集。
- 体系列表（均为经典 docking/affinity benchmark 成员，结构/生物功能注释齐全）：

| PDB | 复合物 | 功能类别 |
|---|---|---|
| 1BRS | barnase/barstar | 酶-抑制剂（RNase/抑制剂，高亲和经典体系） |
| 1AY7 | RNase Sa/barstar | 酶-抑制剂 |
| 1EMV | colicin E9 DNase/免疫蛋白 Im9 | 细菌毒素-免疫 |
| 1PVH | LIF/gp130 | 细胞因子-受体 |
| 1GCQ | Vav/GRB2 SH3 域 | 信号传导（SH3-配体） |
| 1AK4 | cyclophilin A/HIV-1 衣壳 | 病毒-宿主 |
| 3SGB | proteinase B/抑制剂 | 丝氨酸蛋白酶-抑制剂 |
| 2OOB | ubiquitin/泛素连接酶（UBA 域） | 泛素化 |

- **契合度**：完全对齐 FlexPPI 的 `--apo-pdb` + `--holo-pdb` 双输入接口与"结合响应 ΔRMSF/ΔSASA/Δ侧链熵"任务；8 个体系覆盖酶-抑制剂、细胞因子、SH3 信号、病毒-宿主、泛素化多种 PPI 类别，可检验模型泛化而非单一类别。
- **注意事项**：2022 年轨迹（力场/时长需查记录页；该数据集不含温度声明但为全原子）；下载体积需核对记录页文件列表。

---

## 3. 大规模免疫复合物库：ThermoPCD 2.0（抗体-抗原，与 TCR-pMHC 同属免疫识别）

- **旧版**：Singh & Stan, *ThermoPCD: a database of molecular dynamics trajectories of antibody–antigen complexes at physiologic and fever-range temperatures*, **Database (Oxford) 2024**, DOI [10.1093/database/baae015](https://doi.org/10.1093/database/baae015)，[PMC10950042](https://pmc.ncbi.nlm.nih.gov/articles/PMC10950042/)。
- **2.0 版**：*ThermoPCD 2.0: an updated database…*, **Journal of Big Data 2025**, DOI [10.1186/s40537-025-01335-y](https://link.springer.com/article/10.1186/s40537-025-01335-y)。
- **规模**：约 50 个免疫复合物（v1）+ 41 个新结构（v2：病毒 10/细菌 12/原虫 6/过敏 7/肿瘤自免 7 + 分离抗体）+ 11 项旧条目扩展；唯一 PDB 体系约 85–90；HDR 共 **1,061 个数据集条目**（每复合物 × 温度 × 种子独立 DOI）。
- **轨迹**：310–313 K（3–4 个温度点）× 100/200/500 ns × 1–3 随机种子，GROMACS 原始文件 + 每 50 ns PDB 快照；**数据库自带逐残基 RMSF、SASA、Rg、RMSD、H-bond 处理标签**（论文 + 官网双重确认）。
- **体系构成**：PD-1/PD-L1（PDB 4ZQK）、PD-1/PD-L2（6UMT）、CTLA-4、HER2、pembrolizumab（5DK3）、trastuzumab/pertuzumab 等治疗性抗体与肿瘤/自免抗原复合物；包含 PD-1、PD-L1、CTLA-4 等免疫检查点。
- **下载**：[https://dataverse.harvard.edu/dataverse/ThermoPCD](https://dataverse.harvard.edu/dataverse/ThermoPCD)，免注册，示例条目 DOI 10.7910/DVN/WS83QK。
- **apo/holo 情况**：⚠️ 仅部分覆盖 —— 约 6–7 个游离抗体/抗原条目（如 5DK3 全长 pembrolizumab、1OAY），无系统性 apo/holo 配对；但同一复合物多温度 × 多种子天然构成"温度响应"标签，且可用游离条目近似 apo。
- **契合度**：与 TCR-pMHC 同属免疫识别（抗体-抗原），论文叙事连贯；自带处理标签可免去自行算轨迹；PD-1/PD-L1 体系本身就是热门免疫检查点 PPI。

---

## 4. Zenodo 单项轨迹（状态对比型，非严格 apo/holo）

子代理逐一核实（全部经 Zenodo 记录页 + 论文原文交叉验证）：

### 4.1 AIFM1–AK2（首选单项）★
- 体系：AIFM1 二聚体（PDB 4bur / 新 cryo-EM 9FL7）+ AK2（腺苷酸激酶 2）
- 轨迹：全原子 CHARMM36，GROMACS 2020.4，两个体系各 2×2 μs（共 8 μs），0.5 ns 快照 → 每复合物 2×4001 构象
- 状态对比：AIFM1(FAD):AK2 vs AIFM1(NAD+FAD):AK2 —— **NADH 结合引起结合响应变化**，是"结合响应"型分析的自然对照；二聚体中未结合 AK2 的 protomer 可作内部参照
- 下载：[Zenodo 15077557](https://zenodo.org/records/15077557)，2×2.3 GB，CC-BY-4.0
- 论文：Schildhauer et al., *An NADH-controlled gatekeeper of ATP synthase*, **Molecular Cell 2025**, DOI [10.1016/j.molcel.2025.06.007](https://doi.org/10.1016/j.molcel.2025.06.007)（bioRxiv 2024.08.22.609182）

### 4.2 RIPK1–VHL（PROTAC 体系）★
- 体系：RIPK1 激酶域（PDB 4NEU）+ VHL E3 连接酶（PDB 5NVW）+ PROTAC 三元复合物（DiffDock 对接）
- 轨迹：全原子 Amber14SB + TIP3P + GAFF2，**约 420 μs**（Folding@Home 31 个种子结构，仅溶质轨迹，1 ns 存帧）+ 4 个三元状态（Site2/Site3 × Linker State3/4/5，各 5 条轨迹）
- 下载：[Zenodo 18155129](https://zenodo.org/records/18155129)，Raw_Trajs 11.8 GB + Docking_Trajs 4.5 GB，CC-BY-4.0
- 论文：Wu et al., *Metastable Protein–Protein Interactions as a Design Principle for PROTACs: Insights from the RIPK1–VHL System*, **JACS Au 2026**, DOI [10.1021/jacsau.6c00260](https://doi.org/10.1021/jacsau.6c00260)
- 注意：无 apo；轨迹不含溶剂/离子（SASA 的水介导指标受限）

### 4.3 HADDOCK MD-scoring（25 复合物）★
- 体系：Docking Benchmark 5 的 25 个 PPI 复合物（训练 15 + 额外 5 + 独立测试 5：1BJ1, 1BUH, 1E4K, 1E6J, 1EAW, 1JPS, 1KAC, 1NSN, 1OC0, 1PXV, 1RV6, 1VFB, 2A5T, 2NZ8, 3EO1, 1GPW, 1KTZ, 2OOB, 3S9D, 4G6M, 1YVB, 2HRK, 2O8V, 2Z0E, 3F1P）
- 轨迹：全原子 CHARMM36m，GROMACS 2019，每复合物 100 ns；8 个 HADDOCK 模型 × 2 重复 + 参考晶体 × 4 重复（测试集无晶体重复）
- 对比类型：**native（晶体）vs non-native（HADDOCK 错误模型）** —— 天然"结合稳定性/结合响应"对照
- 下载：[Zenodo 4629895](https://zenodo.org/records/4629895)，11.3 GB；脚本 [github.com/haddocking/MD-scoring](https://github.com/haddocking/MD-scoring)
- 论文：Jandova, Vargiu, Bonvin, **JCTC 2021**, DOI [10.1021/acs.jctc.1c00336](https://doi.org/10.1021/acs.jctc.1c00336)

### 4.4 Protein A–IgG（pH 响应）☆
- 体系：金黄色葡萄球菌 Protein A（链 E，含 E/D/A/B/C 五域）+ 人 IgG1 Fc（链 B/D 同源二聚体），子代理实测下载解析 PDB 确认
- 轨迹：全原子 CHARMM36m，pH 7.3 与 pH 5.6 各 400 ns（累计 ~1 μs）
- 下载：[Zenodo 17039261](https://zenodo.org/records/17039261)，2×10.6 MB，极易下载
- ⚠️ 无关联论文；初始结构为 docking 产物、无明确 PDB ID

### 4.5 SARS-CoV-2 RBD–ACE2（SIRAH 粗粒化）☆
- 体系：Spike RBD + ACE2-B0AT1 六聚（PDB 6M17）
- 轨迹：粗粒化 SIRAH 2.2 + Amber18，4 μs（含二硫键更新版）
- 下载：[Zenodo 3942566](https://zenodo.org/records/3942566)；论文 [Front. Med. Technol. 2021](https://doi.org/10.3389/fmedt.2021.644039)
- ⚠️ 粗粒化无法直接算全原子逐残基 RMSF/侧链熵/SASA，需 backmap，仅作粗粒化扩展方案

### 4.6 RBD–ACE2 全原子（Zenodo 17780782）✗
- 35 GB 9 分卷；元数据极简、时长/apo/论文均未确认 —— 不推荐优先。

---

## 5. 尚未发布但值得跟踪的提案

- **DynaPPI**（NeurIPS 2025 AI4Science 数据提案竞赛，[arXiv 2608.10435](https://arxiv.org/html/2608.10435)）：设计为 1 万复合物的"解离 → 结合"重组 MD 轨迹（40% 蛋白-蛋白，含 barnase-barstar、HIV-1 蛋白酶二聚体等），计划含 ΔSASA/接触图标签，论文声明将托管于 Zenodo/Dryad。**数据尚未发布**，可周期性检查 Zenodo/Dryad。
- **DynPro**（NeurIPS 2025 提案）：数万单链体系 × ≥100 μs 自适应采样，未提及复合物/apo-holo，暂不适用。

---

## 6. 推荐行动清单（按投入产出排序）

1. **零成本第一步**：从本地 DynaRepo test53 配对（已有 37 对 apo/holo + 标签）和 DynaBench/PDBbind 复合物中选 3–5 个代表体系（建议：1BRS 或 1CGI 酶-抑制剂、1WEJ 抗体-抗原、2PCC 电子传递复合物、1GCQ SH3 信号），按 7L1D/7RRG 同规格产出 PyMOL 着色 + apo/holo 响应报告 —— 与 TCR-pMHC 案例形成跨类别泛化证据，且有本地 MD 真值可直接验证预测 RMSF。
2. **首个新数据源**：下载 Zenodo 6638504（8 复合物，严格 apo/holo 配对），扩展现有 `index-dynarepo`/`compute-dynarepo-labels` 管线处理这批轨迹，直接进入配对微调与响应评估。
3. **免疫方向扩面**：接入 ThermoPCD 2.0（Harvard Dataverse 免注册），利用自带 RMSF/SASA 标签 + 多温度 × 多种子，做 PD-1/PD-L1（4ZQK）、PD-1/PD-L2（6UMT）、CTLA-4 等免疫检查点的"温度响应 + 结合响应"分析，叙事上与 TCR-pMHC 免疫识别主线一致。
4. **状态对比型单项**：AIFM1–AK2（NADH 结合响应，Molecular Cell 2025）或 RIPK1–VHL（PROTAC 多状态，JACS Au 2026）作为"非免疫 PPI 泛化"案例；HADDOCK MD-scoring 的 native/non-native 对比可做"结合稳定性"标签。
5. **跟踪**：DynaPPI 发布后评估其解离→结合轨迹对 FlexPPI 响应头训练的价值。

---

## 7. 检索说明与证据边界

- 主要检索入口：agent-reach（Exa）、web_search 多路查询、Zenodo Records API、Harvard Dataverse API、RCSB/Semantic Scholar/PMC 交叉验证、两个后台子代理（分别深挖大型数据集与 Zenodo 单项）的独立核实。
- 下载链接均来自官方记录页原文，未臆造；子代理明确标注了无法确认的字段（如 Protein A–IgG 的来源 PDB ID、Zenodo 17780782 的时长/论文）。
- ThermoPCD 旧版期刊为 *Database (Oxford)* 2024（非 NAR），已按 PMC 记录更正。
- DynaPPI/DynPro 为"数据提案"，本报告不将其计为可用数据；引用其 arXiv/OpenReview 页面仅为跟踪用途。
- 本地 DynaRepo 数据来自公开库 [dynarepo.inria.fr](https://dynarepo.inria.fr/)（MDDB 法国节点，bioRxiv 2025.08.14.670260，NAR Database 2025），本地目录为 2026 年已下载快照。
