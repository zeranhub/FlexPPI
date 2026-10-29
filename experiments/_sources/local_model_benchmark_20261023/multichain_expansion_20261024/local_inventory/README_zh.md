# 本机真实多链蛋白扩充清点

2026-10-24按用户新澄清的“protein组实际排除单链”要求，本清点只寻找真实MD输入至少包含两条蛋白链、且可支持同一种蛋白身份的候选。没有把单链轨迹拼成多链，也没有把异源伙伴改名为protein。本报告不更改此前3,063条冻结分类，专门记录可供新实验方案继续核查的数据。

本机DynaRepo下载清单位于`D:/data/dynarepo_non_nucleic_manifest.json`，有922个实际下载项目；原3,063条队列只使用890个DynaRepo项目，差额32个。对这32个项目仅读取小型`structure.pdb`中第一个模型的蛋白Cα残基，保留实际链、完整序列、序列哈希及链长，并只对大型轨迹做文件存在性与大小统计，没有打开、解码、遍历或计算其坐标。差额32条中，29条具有至少两条蛋白链、3条为单链；29条多链中，28条含明显不同的序列伙伴，不能在没有生物身份依据的情况下充当同源聚体。只有 **A00AJ / 6B6U** 的两条链是强同蛋白候选：A链525 aa、B链519 aa，较短链全部519个残基精确匹配较长链，长度覆盖98.86%。其真实结构1,309,643 bytes、XTC 307,813,948 bytes、TPR 4,674,204 bytes，均位于`D:/data/dynarepo_non_nucleic_full/A00AJ/`。

6B6U目前仍标作“待官方蛋白身份确认的同源聚体候选”，没有冒充已经认证的新测试集。它共1,044个残基，超过原Full输入上限1,024；此前训练冻结也明确因长度排除了该记录。是否能按新实验方案完整输入、是否与训练集存在同源关联，以及RMSF标签如何与原任务严格一致，都需要分别核查。项目目录中未发现独立的预计算RMSF标签；真实轨迹存在不等于本轮已经生成可直接评分标签。本清点不执行标签计算或模型推理，也不把“不在原3,063条精确记录中”写成“已证明外部无同源泄漏”。

完整差额表在`dynarepo_32_not_in_original_cohort.csv`和同名JSON，逐链实际序列在`dynarepo_32_actual_chain_sequences.csv`。表中每个项目都保存原来源描述、所有真实链长、相对第一链的序列匹配比例、总残基数、完整目录，以及结构/轨迹/拓扑文件的绝对路径与大小；重复模拟或重复副本没有被算成不同蛋白体系。

原队列已经确认的10条多链protein的本机资产另外列于`existing_10_multichain_protein_local_assets.csv`：A000J/1E3U、A0046/2X89、A004V/3F4Y、A007Q/4LLO、A0088/4Q96、A009H/5JJD、A00C7/3H6G、A00CS/4ZRJ、A00K8/7CQC、A00PW/1RKE。十条均有原始真实结构与XTC；A00PW目录未找到TPR，但有冻结的校正输入结构。其标签可定位到原远程NPZ路径、已存在的本机神经网络导出容器`target[:,0]`与`valid[:,0]`，以及共同参考残基表的sample_id筛选；本清点只定位这些容器，不读取预测或计算成绩。它们已经处于原五折队列，每条在另外四个模型的训练中出现，不能当作新增独立外测。对应留出模型可提供原记录OOF，仍不能因此宣称新家族隔离外测。

这10条也不能全部称为天然同源多聚体：其中1E3U、2X89、4Q96、3H6G为同蛋白/近相同变体的多链输入；3F4Y是gp41派生设计片段体系，4LLO、5JJD、4ZRJ、1RKE是同一成熟蛋白拆开的结构域，7CQC是同一工程融合蛋白在输入中拆成多链。新实验若要求天然同源聚体，应预先把这些含义分开，而不能以既有protein名称替代真实组成审查。

其他指定位置已做目录与摘要级初筛：`E:/Data/protein_md_datasets`中的ATLAS、DPF、chameleon归档按PDB_链组织；存在多链原始PDB编号不证明其模拟包含多链，不能把独立链轨迹合并扩充。BioEmu MegaSim/CATH已下载卡片明确将multi-chain conformational dynamics列为不覆盖用途，不能按ZIP数量当作多链测试体系。`E:/Project/Data/db_DynaRepo`主要是相同来源的静态PDB；`db_AlphaFlow`当前可见大量逐链PDB也不构成真实MD多链标签证据。`E:/Project/FlexPPI/data/structural_discovery`有静态结构但没有本轮确认的匹配真实MD标签；`E:/Project/FlexPPI-archive`查到的清单以旧ITrace异源复合物为主，不能计作新同源聚体。ThermoPCD另由并行审查专门核对`D:/data/ThermoPCD_trajectories/manifest.jsonl`、`data/`及`ThermoPCD_labels/`，本报告不重复扫描或提前计数其候选。

当前严格结论是：已找到 **1个有真实多链结构与MD轨迹的新增同源强候选，尚未确认可独立测试的新增同源体系**；本机DynaRepo差额本身不能支持“新增几十个同源蛋白”。审计与可复现脚本为`local_inventory_audit.json`、`inventory_dynarepo.py`。仅32个小拓扑PDB及清单计算SHA256，大型轨迹不计算哈希、不读取内容；全部派生输出仅写本目录，原冻结成员、输入、标签和成绩均未修改。
