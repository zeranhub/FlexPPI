> 当前执行策略更新（2026-10-10）：用户明确要求“MD跑一次也够了，先聚焦于论文”。原计划的新12条三重复生产与5项后续分析已在实际运行父队列中取消，不会启动；现有模型训练及CPU对照保留。本地历史四状态每状态一条100ns已重新验证全部400段轨迹和质量表并整合进v18。下文三重复制备与运行方案仅为此前计划的保留记录，不能作为当前仍将运行的任务或已完成结果。现行策略见 `pipeline_control.json`，运行时取消证明见 `pipeline_md_cancellation_audit_20261010.json`，单轨迹完整核查见 `kix_existing_reuse_20261010/reuse_audit.json`。

# KIX 四状态 MD 完成审计与新三重复批次方案

审计日期：2026-10-10。权威资产为 `E:\Project\FlexPPI`，并与旧工作树4666内的可复现脚本核对。本任务没有运行GPU、改变原脚本或原数据、安装本机依赖，也没有修改稿件。下文的生产命令由根代理在 hit-gpu 的新环境调度；准备参数化和通过烟测不等于完成MD实验。

## 结论与已完成证据

历史四状态 apo/myb/mll/both 各有一条真实完成的 `rep01`：每条10 ns带重原子位置约束平衡与100 ns无约束生产，55,000,000个2 fs积分步，生产阶段100个连续1 ns分段、每段1,000帧、间隔1 ps。历史队列于2026-09-16北京时间完成400 ns生产后正常退出；`PILOT_POLICY.json`与`progress.json`明确 deferred_repeats=[2,3,4,5]，没有已完成的多重复结果。2026-10-10的新授权支持补充重复，不能将旧的“仅rep01”策略误读为当前仍需再向用户请求许可。

本次从本地逐个重新计算全部400个DCD的SHA256，全部与历史分析摘要登记一致；文件字节数、连续文件名0000–0099以及DCD头部每段1,000帧也全部一致。合计8,780,110,400字节，约8.78 GB。四份逐残基CSV及protein.pdb的SHA与各自摘要及prepared元数据一致。完整400段证据、每状态完成标记、坐标/指标哈希及缺失系统文件记录在 `md_local_input_audit.json`；这次执行是只读文件验证，没有重新计算MD指标，也没有把相邻帧当独立重复。

| 状态 | 实际链及保留作者编号 | 历史生产帧数 | 本次DCD/指标哈希 | 新批次独立重复 |
|---|---|---:|---|---|
| apo | KIX B586–672，87残基 | 100,000 | 全部一致 | rep01/02/03 |
| myb | cMyb A291–315 + KIX B | 100,000 | 全部一致 | rep01/02/03 |
| mll | KIX B + MLL C839–857 | 100,000 | 全部一致 | rep01/02/03 |
| both | A+B+C，131残基 | 100,000 | 全部一致 | rep01/02/03 |

## 原样迁移的缺口与最终决定

本地保存了四份完整蛋白拓扑坐标、四份元数据、`prepared_READY.json`、修复审计和400份protein-only DCD，但没有四套完整的 `system.xml`、`solvated.pdb`、完整溶剂二进制检查点。旧工作树4666也未找到这些文件。蛋白DCD没有水和离子坐标，不能据此重建历史溶剂体系。旧服务器目录为 `/ranze/runs/flexppi_kix_md_20260914`，环境为 `/ranze/envs/flexppi-kix-md-20260914`；SSH配置中历史候选 `alias=ranze@172.16.196.191:22` 的只读BatchMode探测连接超时；hit-gpu 当前没有该 `/ranze` 根目录。

根代理已决定新建隔离的重建批次：四状态各重新运行rep01、rep02、rep03，共12条轨迹，1,200 ns生产、1,320 ns含平衡；旧rep01只保留为单独历史敏感性参照，绝不计入新批次三重复均值、标准差或状态差值。若后续找回原XML和全溶剂PDB，可以保留并比较，但本批次已经明确标为重新参数化，不以事后系统相似性重新混合两批轨迹。OpenMM的二进制检查点还要求相同硬件和OpenMM版本，跨机器迁移不宜直接当作可靠续跑；本批次每个新重复都从新制备坐标独立最小化和重新生成速度。[OpenMM官方说明](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)

实际payload为 `kix_md_triplicate_rebuild_20261010/`，已经准备冻结的四状态重原子输入、protocol、三重复MD程序、按状态分片的安全队列、重建前后输入审计、独立逐重复分析、完整delta分析脚本和环境安装脚本。较早的 `kix_md_triplicate_20261010/` 是“恢复原系统后只补rep02/03”的候选包，已被新重建方案替代，不能上传调度该旧包。

## 物理方案、输入和随机种子

新 `run_md.py` 与原脚本仅改变重复次数5→3、队列范围1–5→1–3，以及任务和总时长计数20/2,000/2,200→12/1,200/1,320。原物理参数及seed规则保持不变：原版CHARMM36（不是CHARMM36m）、配套CHARMM TIP3P水/离子、100 mM NaCl并中和、298 K、NVT、不加压力耦合、Langevin摩擦5 ps⁻¹、2 fs、PME、12 Å非键截止、10–12 Å势能切换、HBonds约束、刚性水，10 ns平衡的重原子约束为 U=k r²、k=4,184 kJ mol⁻¹ nm⁻²；生产约束置零。两阶段最小化均为OpenMM L-BFGS，每阶段最多2,000步。1 ns提交一次检查点和JSON；生产每1 ps保存全部蛋白原子、每10 ps保存热力学日志。上述方案与Peiffer原文有最小化器、LJ切换实现等差异，不写成其方法的逐项复刻。

从2AGH模型1保留原作者编号；MLL严格用839–857，不能套用旧840–858构建体的预测。先对完整三元复合物补齐D841的CG/OD1/OD2及新末端T857的OXT，再从同一修复重原子骨架删除伙伴派生四状态。D841按实际PDB存档保留Asp。KIX五个His为HID、链端带电未封端，属于明确的制备假设，不将NMR样本pH5.5视作模拟pH的证据。模型预测依赖的N/CA/C/O骨架在原冻结输入中完整。

速度及Langevin随机种子采用 `20260914 + repeat*100 + state_index`；state_index按apo/myb/mll/both为0/1/2/3。新三重复各自从起始坐标重新最小化并生成速度，彼此不延长已有轨迹，不复制第一条轨迹片段；相同构建体内的三个不同种子是独立启动。沿用历史数值seed规则不意味着可以将历史rep01并入新批次。

| 状态 | rep01 | rep02 | rep03 |
|---|---:|---:|---:|
| apo | 20261014 | 20261114 | 20261214 |
| myb | 20261015 | 20261115 | 20261215 |
| mll | 20261016 | 20261116 | 20261216 |
| both | 20261017 | 20261117 | 20261217 |

## 环境与可执行命令

独立MD环境固定为 `/root/flexppi_completion_20261010/md_env`，使用Python3.11，尽量精确重现历史 numpy2.2.6、scipy1.15.3、OpenMM8.4.0.post2、OpenMM-CUDA-12 8.4.0.post2、MDAnalysis2.9.0、pdbfixer1.12.0。它与主Torch环境隔离。`install_md_env.sh` 仅向该新目录安装，pip缓存和官方源码存档也在新实验目录。pdbfixer官方v1.12 tag的setup.py确认版本为1.12.0，固定commit为 `94cfa4c0ca551cdc5f13320f9a658efd59f2b881`，安装官方tar.gz并保存实际下载SHA；不用master，也不自动升级主环境。[官方release](https://github.com/openmm/pdbfixer/releases/tag/v1.12)、[固定版本setup.py](https://raw.githubusercontent.com/openmm/pdbfixer/v1.12/setup.py)、[OpenMM固定发行](https://pypi.org/project/OpenMM/8.4.0.post2/)

上传payload至下例目录后按顺序执行；PY311替换为根代理已建的绝对Python3.11路径。这里提供命令，没有在本审计中执行安装或GPU动作。

```bash
TASK=/root/flexppi_completion_20261010
MDROOT="$TASK/kix_md_triplicate_rebuild_20261010"
MDPY="$TASK/md_env/bin/python"
bash "$MDROOT/install_md_env.sh" "$PY311"
"$MDPY" "$MDROOT/code/run_md.py" prepare --root "$MDROOT"
"$MDPY" "$MDROOT/code/audit_rebuilt_prepared.py" --root "$MDROOT"
CUDA_VISIBLE_DEVICES=0 "$MDPY" "$MDROOT/code/run_md.py" smoke --root "$MDROOT"
# 在生产前生成下述Full+Static冻结预测及FROZEN.json。
CUDA_VISIBLE_DEVICES=0 "$MDPY" "$MDROOT/code/queue_triplicate.py" --root "$MDROOT" --verify-only
CUDA_VISIBLE_DEVICES=0 "$MDPY" "$MDROOT/code/queue_triplicate.py" --root "$MDROOT"
# 所有12条已完整生产、各自分析及QA后：
"$MDPY" "$MDROOT/code/analyze_triplicate_delta.py" --root "$MDROOT"
```

若根代理按GPU资源选择两个并行状态分片，可以分别指定 `--states apo myb` 与 `--states mll both`，每个进程绑定不同CUDA_VISIBLE_DEVICES；不在同一GPU同时开多个上下文。队列对每个 `runs/<state>/repXX/writer.lock` 取得独占锁，重复调度同一作业会失败。MD内部的根级progress.json被不同分片交替更新，只能当最近消息；完整进度以各作业分段/DONE及 `extension_worker_<states>.json` 为准。后台启动、GPU分配及资源占用由根代理处理。新根中的STOP_AFTER_SEGMENT用于1 ns段结束后安全停止；原任务的停止文件不复制、不删除。

## 制备、生产、分析的真实完成判据

进入生产前必须有新prepared/READY.json；`audit_rebuilt_prepared.py` 核查力场XML哈希与历史固定版本一致、四状态全部输入哈希相同、重原子只发生整体平移及不超过0.002 Å的PDB舍入差、骨架和作者编号完整、五His为HID、电荷中性、周期重原子间隔下界至少2 nm、非键开关及截止与方案一致，并记录新solvated.pdb/system.xml的SHA。必须在当前GPU四状态各通过40 ps技术烟测，包括1 ps DCD往返、40帧、温度、共价几何和当前平台二进制检查点恢复。烟测不计入100 ns、不证明平衡或收敛。

每个重复必须有10个平衡提交JSON和100个连续生产提交JSON、对应完整检查点，生产DCD恰好100,000帧，时间间隔1 ps、蛋白原子数与拓扑一致，无不有限坐标或断链/PBC巨跳。DONE中必须production_ns=100、equilibration_ns=10、final_step=55,000,000。`analyze_job.py`逐条计算RMSF/Hχ，保留每段SHA、全部帧几何检查和前后半程敏感性。最后 `analyze_triplicate_delta.py` 仅接受12条新批次重复全部完成，再次核查全部DCD哈希、CSV哈希、作者编号和统一有效位点，产生完整绝对/delta/重复敏感性源表；没有完整12条时会报错，不能用缺失状态或前几条先宣布实验完成。

质量解释还须检查热力学日志、生产过程KIX RMSD、伙伴界面接触保持情况、每条前后50 ns和累积10 ns窗口变化，并展示三条独立轨迹差异。旧quality_pilot.py写死rep01，已在新目录制作 `quality_one_repeat.py`，显式接受 `--repeat rep01/rep02/rep03`，验证源重复编号并给全部输出写入实际重复名称；原全帧数学分析保持不变，改动追溯见QUALITY_CODE_PROVENANCE.json。每个重复可执行 `"$MDPY" "$MDROOT/code/quality_one_repeat.py" --root "$MDROOT" --repeat rep01`，默认分别保存quality/rep01等目录；三个重复都要保留。需要区分“生产已完成”和“慢运动收敛充分”：三条100 ns仍可能没有覆盖所有构象，质量不佳或重复间差异大时保留结果并报告限制，不因符号与预期不一致剔除轨迹或追加至出现预期方向。

## 已固定的绝对及delta分析

新方案在看到历史pilot之后、但新生产开始前固定。因此可称“追加实验预先固定分析”，不可写成从未看过任何KIX结果的盲态预注册。详细机器可读规则在payload的ANALYSIS_CONTRACT.json。主模型使用当前chemical_identity_corrected_v1的固定fold1 Full；Static也使用匹配训练范围的现存fold1静态模型，不重新用MD标签拟合。旧frozen_predictions/对应performance_first模型，不能与本次修正模型混称同一冻结预测。新模型及Static需要在新生产前生成四状态全部B586–672逐残基输出，权重、代码、四份输入和CSV全部记录SHA；Intra及AA可同时冻结为辅助对照。冻结预测CSV字段为 `state,method,chain,author_residue_number,amino_acid,predicted_rmsf_angstrom,predicted_hchi_nats`；method采用full/static/intra/aa。FROZEN.json需有prediction_csv_sha256及四状态input_sha256。安全队列在缺少Full或Static、完整87位点或哈希不匹配时拒绝生产。

RMSF按每条轨迹KIX自身87个Cα对齐第一生产帧，再计算围绕对齐后轨迹平均位置的波动；每条重复独立计算，再等权平均三个逐残基读数。不把三条轨迹串接来计算所谓三重复结果，不转成Peiffer论文优势簇参考的RMSF。Hχ沿用本项目三个120°χ角箱的联合状态频率与自然对数Shannon熵，每条100 ns独立计算后再重复等权平均，不先合并状态计数。主侧链分析纳入所有实际有χ角的位点，包括Pro，KIX为77位点、L12为6位点；Ala/Gly为不适用，不能填零。旧AGP排除的74/4位点结果仅保留为明确标注的历史敏感性。L12仍固定作者编号614–621，RMSF八个位点；区域算术均值为主，sqrt(mean(RMSF²))另名保存，不混用。

绝对分析保存每状态、每重复、每方法的全KIX和L12原单位MAE、RMSE、偏差、Pearson/Spearman及均值，再报告相对三条MD均值的分数。RMSF Å和Hχ nats分别分析。状态内百分位是辅助排序读数：在每个状态/重复内，用全部有效残基的平均秩计算100×(rank−1)/(n−1)，再对三个MD重复等权平均；每个冻结模型在其自身该状态的全部有效位点上独立转秩。百分位差的单位是百分点，不能继续标Å或nats，也不能把幅度归一化后的一致性替代原单位delta。

状态变化保留所有六个两两对比：myb−apo、mll−apo、both−apo、mll−myb、both−myb、both−mll。逐残基delta是后状态的三重复平均减前状态的三重复平均；Full/Static的delta用相应冻结预测相减。主要端点是全KIX原单位RMSF delta相对零响应及Static的MAE收益，以及L12平均RMSF的mll−myb和both−apo；其余状态、Hχ、百分位及both−myb−mll+apo为预先规定的辅助结果。非加和读数只是模拟观察量的四状态差，不解释成结合自由能协同或实验亲和力改变。

零响应对照对所有位点预测delta=0。Static信息对照使用同一四状态实际骨架和接触摘要预测的delta。对每项比较保存 MAE(零响应)−MAE(Full) 与 MAE(Static)−MAE(Full)，正值表示Full相应误差较低；不拟合缩放、平移或选择检查点。不按MD变化大小、符号、局部相关或后验最优位点筛选，全87/77位点及固定L12全部保留。还保存各状态三条逐残基标准差/最小最大值和全部3×3跨状态参考组合的delta及评分，这九个组合相互相关，不是九个独立MD重复；跨状态编号相同并不表示物理配对，不能做配对检验。三重复主要用于审计方向与变化幅度的可重复性，当前不报告bootstrap置信区间或p值，不以残基、帧数、模型折数代替重复数，也不将所有九种组合方向相同写成显著性证明。

## 功能注释及可支持的结论

2AGH一手结构研究支持KIX在两个不同表面结合cMyb和MLL，并显示MLL结合时局部环和末端螺旋的结构变化，可用来交代识别背景和伙伴界面。[De Guzman等，J Mol Biol 355:1005–1013，2006；DOI10.1016/j.jmb.2005.09.059](https://www.rcsb.org/structure/2AGH) Peiffer等针对四种组成的研究明确定位L12–G2作者编号614–621，并结合MD构象微状态、甲基序参量和DSF讨论其功能联系，支持将这八个位点事先定义为生物学关注区域。[Peiffer等，PLOS Comput Biol 18(4):e1009977，2022；DOI10.1371/journal.pcbi.1009977](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009977) 该文的构象熵/甲基序参量不等于本文Hχ，本批次只有MD比较，不直接测亲和力、催化、热力学熵或转录功能。功能区注释成立与模型伙伴诱导变化预测正确是两项不同证据。

历史pilot显示全KIX绝对柔性排序较好，但预定义L12伙伴诱导变化未获得支持：旧五检查点模型对mll−myb的L12 RMSF平均变化为−0.413 Å，旧MD为+0.014 Å；both−apo分别−0.898和+1.017 Å。当前稿件fold1修正模型已有另一组绝对RMSF排序结果，不能继续直接援引旧五模型数值。本次新三重复应完整报告绝对定位和伙伴变化，包括方向冲突、零响应不弱于Full及Static更优的情形。即使新结果支持部分变化，也仍是已知近同源KIX上的模拟案例，不是完全去家族同源的外部泛化或功能性变构机制验证。

## 文件与验证记录

`md_local_input_audit.json` 已实际生成并确认400段一致；新payload的Python脚本经过本地语法检查，真实新MD尚未准备、烟测或生产。delta分析另通过12作业合成数据的数值及文件校验测试，包括正确恢复已知delta、Full相对零响应/Static的收益、平均秩与百分位、77/6位点掩码及故意破坏DCD文件后的拒绝；这些仅为明确标注的合成代数/元数据测试，不是MD结果，见ANALYSIS_QA.json。`PAYLOAD_MANIFEST.json`明确prepared=false、production_started=false，并保存输入及代码SHA。`ANALYSIS_CONTRACT.json`保存固定端点/掩码/重复单位；`install_md_env.sh`保存隔离环境及官方PDBFixercommit；新目录的所有文件均为实验自己的准备产物。后续根代理完成运行或修改自己的新代码时应更新manifest，保留本审计与历史数据，不改原worktree或既有E盘资产。
