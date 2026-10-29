# TCR–pMHC 真实状态输入核查（2026-09-15）

本笔记只记录文献、结构文件和几何核查。未运行新的 FlexPPI 推理、训练或 MD。结构源文件保存在 `artifacts/tcr_pmhc_reassessment_20260915/sources`；19 个原始 mmCIF 的下载 URL、SHA256、作者链号、逐位主链原子、缺失原子和实体序列保存在 `structure_sequence_audit.json`。可重复脚本为本目录 `audit_notti_structures.py`。

后续核查更新：本页1G4表按每对结构各自共同位点计算，不同配对的位点数不同；主报告新增四结构统一366位点的结果，三者RMSD为8.117、7.641和0.719 Å，跨条件展示优先使用该统一版本。另外，Knapp2019 A6实际首帧β链仅209位点，缺失多个远端片段，不能直接作为完整A6的H3/FG机制参考。综合优先级与新证据以[主报告](TCR_pMHC_机制重评与验证方案.md)为准。

## 1. 首选：1G4(LY) 完整 TCR–CD3 的真实纳米盘 apo / holo

Notti et al., Nature Communications (2025), DOI [10.1038/s41467-025-66939-7](https://www.nature.com/articles/s41467-025-66939-7)。真实无配体闭合/压缩构象为 [8TW4](https://www.rcsb.org/structure/8TW4) (ND-I) 与 [8TW6](https://www.rcsb.org/structure/8TW6) (ND-II)；结合态为 [9C3E](https://www.rcsb.org/structure/9C3E)；[9BBC](https://www.rcsb.org/structure/9BBC) 为去垢剂中无配体的开放对照。这里的 apo 是独立实验重建的无配体样品，既没有从 holo 删除 pMHC，也没有将结构形变人工设为零。

4 个条目都用作者链 A=α，B=β，D=CD3δ，E/F=两份CD3ε，G=CD3γ，X/Y=CD3ζ。9C3E 另有 H=单链 NY-ESO-1(C9V)–β2M–HLA 融合配体；H 是作者链号，mmCIF label 链为 I。结构建模中的 pMHC 并非普通三条独立链，不能把编号1–593连续当成一个无缺口天然蛋白。初次预测可保留所有有完整主链的已观察蛋白链，具体解析需保留缺口和未知残基掩码。

| PDB / 真实实验状态 | 全局分辨率 | α已建模Cα / 实体序列长度 | β已建模Cα / 实体序列长度 | 与9C3E共同α/β位点 | 共同αβ整体拟合Cα RMSD |
|---|---:|---:|---:|---:|---:|
| 8TW4 / 纳米盘无配体 ND-I | 3.3 Å | 186 / 274 | 241 / 556 | 186 / 239 | 8.464 Å (425位点) |
| 8TW6 / 纳米盘无配体 ND-II | 3.1 Å | 196 / 274 | 219 / 556 | 196 / 219 | 7.636 Å (415位点) |
| 9C3E / 纳米盘HLA结合态 | 3.5 Å | 252 / 272 | 281 / 557 | 252 / 281 | — |
| 9BBC / 去垢剂无配体 | 3.3 Å | 253 / 273 | 284 / 305 | 252 / 281 | 0.829 Å (533位点) |

上述共同位点按作者链号和残基编号匹配，在所有 TCR/CD3 共同位点上均为完全相同的氨基酸；这些已观察位点均有 N、CA、C、O。匹配全集、左右独有位点及差异已输出到 `notti_matched_residue_comparisons.json`。然而完整表达构建体不能说完全相同：apo和holo的β实体分别注释mCherry/EGFP融合，其他实体端部长度也不同。实体长度不是可见胞外结构长度。共同可见核心序列一致支持比较，但标签、膜环境及重建条件需记录。

`notti_geometry_audit.json` 给出逐残基对应和拟合结果。预先选取的粗略结构域窗口 A21–130、A140–230、B22–140、B150–270 各自独立拟合后，ND-I→holo的RMSD约1.23–1.76 Å，ND-II约1.48–2.91 Å；这些是探索窗口而非经IMGT核实的精确结构域边界。与整体7.6–8.5 Å相比，说明很大一部分差异来自域间姿态重新排列。这些数值是新计算的**静态结构差异RMSD**，不可改名为动力学波动RMSF。

建议主要比较 ND-I/ND-II 各自相对9C3E，在每个真实状态完整输入后只评分共同观察位点；9BBC作为膜环境/构象对照单列。必须另查缺失位点附近的敏感性，因为apo中未观察残基更多，边界残基可能受到缺口影响。不能将apo缺失位点置零或解释成刚性，也不能把“密度更清楚”直接转成逐残基实验RMSF。此组可检验模型是否为真正不同构象给出与论文区域性稳定化一致的结果；没有同条件轨迹时仍不能量化模拟替代性能。锁定突变只有细胞功能数据，未在本轮得到对应独立突变实验结构。

## 2. A6：此前输入有真实apo，但文献MD起点与静态删伙伴是不同问题

Rangarajan et al. (2018) [原文Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6187629/)明确：unbound MD以1QRN D/E抽出的TCR作初始结构；bound以4FTV为模板回突变，再加水/离子、最小化、逐步放松，使用300K、NAMD/CHARMM27+CMAP。原文时长表述有内部冲突：Results写1200 ns，并在Results和Methods都写舍弃前200 ns、用最后1000 ns计算RMSF；Methods另有“300 ns production”一句。因此不能把300 ns作为已经核实的总时长，应核查轨迹后再确定。四个结构域分别拟合，去除了域间运动。其结论来自之后可运动的轨迹，不能用未放松的删除结构替代。原文列出的变体对应顺序/重复133位点存在歧义，应以坐标实际序列和补充材料核对，而非逐字直接生成回突变。

| A6输入 | 作者目标链 | 主链Cα覆盖 | 与3QH3观察序列差异 |
|---|---|---|---|
| 3QH3 独立无配体晶体 | A/B；C/D为另一个晶体副本 | 194/245；副本D少4位 | apo本身有稳定化工程构建背景 |
| 1AO7 原始Tax复合物 | D/E | 115/209 | β共同序列2处不同；α恒定区大段未建模 |
| 1QRN Tax-P6A复合物 | D/E | 200/243 | α2处、β4处共同位点不同，另有端部/编号差异 |
| 4FTV A6c134高亲和变体 | D/E | 200/245 | α1处、β5处共同位点不同，另有端部/编号差异 |

`classical_apo_holo_sequence_mapping.json`保存了以观察序列全局比对得到的每个对应位置与不匹配。实际3QH3→4FTV差异是α158C/164T，β98A/99M、99G/100S、100G/101A、101R/102Q、172C/173S（左为3QH3，右为4FTV编号）。这些差异说明“同名A6”不自动等于同构建体。1AO7虽有天然Tax序列，却因缺少大部分Cα结构而不适合恒定区远端主张。

A6可以继续研究，但应将“结构状态预测”与“远端动态传播能力”分开报告。优先寻找原论文平衡后的apo/bound起始状态和逐残基RMSF数据；若只能取得不同晶体，应使用多apo/多holo副本并报告构建体差异，不把总差异先验视为几何伪影，也不能据此证明结合导致的因果效应。

## 3. 其他真实apo/holo路线：优先B4.2.3与DMF5，保留经典组作扩展

| 体系 | 真实apo→holo（作者链） | 本次观察序列核对 | 适合验证内容及限制 |
|---|---|---|---|
| B4.2.3 | 5IW1 A/B→5IVX E/F | α189、β232相同对应位点，零替换 | 同研究的NMR+功能证据；apo另有C/D、E/F副本，第三副本占有率0.6；apo3.0Å且孪晶，holo2.1Å，须测副本敏感性 |
| DMF5 | 3QEU D/E→3QDG D/E | α199、β242相同对应位点，零替换 | 2016MD研究的真实无配体/结合结构，可比较其CDR动态与A6；主链动力学幅度可测，协同运动不能从单残基标量推出 |
| LC13 | 1KGC D/E→1MI5 D/E | α201相同，β240相同+β192 C/A一处差异 | 2019多重复MD体系之一；须保留突变和晶体差异，不称严格同构建体 |
| 2C | 1TCR A/B→2CKB A/B | α202、β237相同，零替换 | 2CKB另有C/D复合物副本；可做经典构象可塑性外部结构检验；当前未核实同条件原始动态数据 |
| 1G4可溶胞外结构 | 2BNU A/B→2BNR或2BNQ D/E | α202相同+128 R/D一处差异，β241相同 | 2BNR含天然9C肽，2BNQ为9V肽；不能和2025的1G4(LY)全长构建体混为一组 |

B4.2.3原文与数据：[Natarajan2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5440810/)，PDB5IW1/5IVX，BMRB26977为β主链/甲基化学位移归属，其余数据声明向作者合理请求；BMRB归属本身不是完整结合前后RMSF真值。DMF5/A6原文：[Ayres2016](https://www.nature.com/articles/srep25070)，其中自由DMF5选3QEU第二个分子；其从bound删伙伴后再跑500ns的反向实验，反映构象弛豫与能垒，不能替换为静态删除一次推理。2C结构原文：[RCSB1TCR](https://www.rcsb.org/structure/1TCR)、[RCSB2CKB](https://www.rcsb.org/structure/2CKB)。

## 4. 可公开取得的MD/结构数据入口

2019 [Knapp等四TCR多重复研究](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007338)共26μs，LC13每状态100×100ns、A6/JM22/1G4每状态各10×100ns。更新核查：Data Availability明确给出四个Figshare归档，API已确认[A6](https://figshare.com/articles/A6_tar_gz/8067827)约399MB、[1G4](https://figshare.com/articles/1G4_tar_gz/8067743)约420MB、[LC13](https://figshare.com/articles/LC13_tar_gz/8067746)约4.24GB、[JM22](https://figshare.com/articles/JM22_tar_gz/8067803)约424MB，准确字节数、校验值和下载地址见literature/knapp2019_*_figshare.json。其主要分歧是“不存在跨四种TCR保守的动态图样”，并不是否认所有远端变化；A6的AB环反而出现约50%的柔性增加，而另三种没有相同模式。适合检验体系特异性，不能强求全部TCR同向。

2024 [McMaster结构数据集](https://github.com/benjiemc/tcr-pmhc-interface-analysis)开放391结构、22TCR、19MHC等位基因、79肽的整理代码与结构，用于筛选真正apo/holo配对和构象差异；它不是MD标签集。2024 [框架工程论文](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1345368/full)给出A6/LC13/868的MD和功能验证，并指向[Zenodo10685013](https://zenodo.org/records/10685013)，值得下一步精确查看文件清单。2025 [TCRDab归档17552276](https://zenodo.org/records/17552276)有3.9GB的TCR Fv数据，但采用CALVADOS3粗粒化再cg2all回映，不能直接当全原子侧链熵真值；其[作者仓库](https://github.com/oxpig/CDR_MD_simulations)另指向[17525665](https://zenodo.org/records/17525665)中一例TCR的全原子MD。以上大型归档未下载。

## 优先顺序

首先以1G4(LY)纳米盘4结构建立真实状态与环境对照，区域性主链柔性结果对应晶体/冷冻电镜支持及锁定突变的功能定位。其次以B4.2.3真实apo/holo和NMR位点做远端区域富集验证，同时清楚区分化学环境改变与柔性改变。第三以A6–DMF5对照及可取得的MD数据测试同为TCR却使用不同柔性模式的机制；A6原2018结论需要原轨迹或条件匹配的MD读出，不能由未放松的结构删伙伴试验作成败判定。
