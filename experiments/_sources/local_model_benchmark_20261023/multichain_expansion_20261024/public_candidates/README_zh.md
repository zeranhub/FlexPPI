# 公开多链同源多聚体候选：MDposit 筛查

2026-10-24 接续「protein 组不含单链」的扩充要求，本轮把公开数据这一步做完。结论是：**公开库里确实还有一批可用的同源多聚体 MD 体系，全部在 MDposit，共 23 条记录（21 个 PDB 条目、约 19 个不同蛋白），且没有一个落在现有 3,063 条队列里。** 这批记录的结构与逐原子涨落数据已经下载并核对；但还没有算成可用标签、也没有跑模型，不能计入任何成绩。

## 一个会漏掉 632 条记录的接口陷阱

MDposit 的 `/projects` 元数据里有两个序列字段，含义不同：`metadata.PROTSEQ` 按实体去重（977 个多链项目里有 632 个是 `chains=2` 而 `PROTSEQ=1`，正是同源二聚体），`metadata.SEQUENCES` 才是**逐链**序列。只按 `PROTSEQ` 统计多链，同源二聚体会被整体判成「蛋白链不足两条」；实测按 `PROTSEQ` 只能得到 11 个同源候选，改按 `chains` 顺序取 `SEQUENCES` 逐链比对后是 41 个。这个差异不是口径争论，而是字段语义问题，必须记录在案。

## 筛查口径与结果

数据源为 `https://mdposit.mddbr.eu/api/rest/v1/projects`，查询条件 `chains.1 $exists`，每页 100 条共 10 页，`filteredCount=977`。逐项目按 `chains` 顺序取 `SEQUENCES`，两两做序列比对（`difflib.SequenceMatcher` 比值），再按下列顺序过滤：

| 过滤步骤 | 剩余 |
|---|---:|
| 多链项目总数 | 977 |
| 逐链序列两两一致（identical） | 37 |
| 逐链一致性 ≥95% 且 <100% | 4 |
| 以上 41 条中不含核酸链 | 41 |
| 且总残基数 ≤1,024（Full 输入上限） | 23 |
| 且 `trajectory.xtc` 存在 | 23 |

23 条的规模为 74–912 残基、2 或 4 条链，全部来自 IRB/MMB 节点，总轨迹约 55 GB（不需要全下，见下）。PDB 条目去重后 21 个；其中 SARS-CoV-2 主蛋白酶（6LU7×3、6Y84、6Y2E）是同一蛋白的五条不同模拟，按蛋白去重约 19 个。与现有 3,063 条队列的 `pdb_id` 交集为 **0**。

明细见 `local_inventory/mdposit_homooligomer_candidates.csv`（筛查结果，含 node/local/规模/轨迹大小）与 `public_candidates/structure_verification.csv`（下载后的实测链数与残基数）。

## 已完成的核对

结构用 `/projects/{local}/structure` 逐条下载（23/23 成功），涨落用 `/projects/{local}/analyses/fluctuation` 下载（23/23 成功），请求清单与 SHA 在 `public_candidates/fetch_manifest.json`。解析 PDB 的 ATOM 记录后：

- 23 条全部是真实多链拓扑，链数 2 或 4，与元数据一致；
- 20 条全部链序列完全相同；1SV4（75/76）、1C0E（151/150）、1G1U（221/219；四个链中一个短 2 个残基）是端部缺口的近同源多聚体，不是异源复合物；
- 逐原子涨落数组与结构原子数在 22/23 条上 1:1 对齐，可直接按原子序号映射到残基；MCV1900316 的拓扑比分析多 220 个原子，暂不可用。

## 单位问题已用轨迹实测解决

涨落数据的单位在文档里没有写。用最小的 A015M（1WFA，74 残基、902 原子、10,000 帧、108 MB 原始 float32 坐标）实测：自行做 Kabsch 叠合后重算逐原子 RMSF，与 MDposit 的 `fluctuation.rmsf.data` 逐原子相比，**相关系数 0.91，比值中位数 0.0987**。因此该数组的单位是**纳米**，乘 10 即 Å；换算后 1AAP 的残基平均 RMSF 为 1.17 Å，与 ATLAS 的 0.99 Å 量级一致。这也说明不能直接拿原始数值当 Å 用——错了就是整体差一个数量级。

`/trajectory` 端点返回的不是 XTC，而是无头的 float32 XYZ 数组（`trajectory.bin`），形状 10,000×902×3；要用轨迹必须自己 reshape，MDAnalysis 直接按 XTC 读会报 `XDR read error = magic`。

## 标签定稿

按 GROMACS/ATLAS 的残基内 RMS 聚合（`sqrt(mean(rmsf_atom²))`）换算成 Å，输出 `labels/per_residue_rmsf_angstrom.csv`。22 条可用（MCV1900316 因原子数不匹配跳过）。逐记录平均 RMSF 为 0.74–5.95 Å，22 条均值 2.02 Å。

## 一个会静默丢残基的问题

MDposit 的结构文件用 Amber 残基名（CYX、HID/HIE、GLH 等）。直接喂给 `flexppi predict` 时，1AAP 的 112 个残基只返回 100 个——**12 个 CYX（二硫键半胱氨酸）被解析器整段丢掉**，不报错。把残基名改回 CYS/HIS/ASP/GLU/LYS 之后（只改第 18–20 列的名字，不动原子顺序、编号和坐标），覆盖率回到 100%。这一步必须写进后续任何脚本，否则丢的正好是双硫键位点。

## 五折推理与首批成绩

用 `artifacts/performance_first_v2_atlas_hchi_v1/fold_{0..4}/epoch_040.pt` 做五折等权平均，22 条 × 5 折共 110 次推理，单次约 3 秒，全部在本机 CPU 上跑完。结果在 `scores_flexppi_new.csv`：

**mean ρ 0.735 · AP 0.566 · MAE 0.942 Å · bias −0.704 Å（预测整体偏低）。**

ρ 0.735 至少说明这批标签没有错位：若残基对不上，相关性会塌到 0 附近。但 bias −0.70 Å 比 ATLAS 上的 −0.33 Å 大一倍，在把它解释成「模型在复合物上更保守」之前，得先分清是模型偏低还是 MoDEL 短轨迹的标签本身偏大。

## 仍未完成

1. BackFlip 1.0/2.1、DynaProt-M、FastProtFlex、ANM/GNM、Static/AA 尚未在这 22 条上运行；
2. 主图仍然只用原来的 10 条多链蛋白，没有混入新记录；
3. 训练暴露（序列同源）未查，不能声称零样本；
4. MoDEL 条目多为单副本 10 ns，ATLAS 为 3×100 ns，参考标签质量不同，将来进图必须写明。

## 复现

`public_sources/screen_by_chain_sequences.csv` 记录全部 977 条的逐链一致性与规模。`public_candidates/` 下：`structures/`、`structures_std/`（改过残基名的推理副本）、`fluctuation/`、`predictions/`（110 个五折输出）、`labels/per_residue_rmsf_angstrom.csv`、`scores_flexppi_new.csv`、`fetch_manifest.json`、`structure_verification.csv`、`inference_manifest.csv`。本目录不修改 3,063 条冻结成员、原横评图或 v63 论文。
