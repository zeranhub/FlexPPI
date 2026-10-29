# 全位点结构特征与强非线性基线

本目录是新增实验；原始结构、标签、模型、划分和 Word 文档不做修改。输入是冻结 corrected pipeline 的全部 **3,063 条记录**，不是只用 83 条共同窗口记录训练。原静态特征共 961,321 行，RMSF/Hχ 有效位置分别为 961,321/824,988。83 条优先表仅供另一分支的配对与构象占据分析；它的 OOF 值和全体表完全相同。

## 运行

```powershell
& 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe' -u 'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\feature_baseline\run_feature_baselines.py' all --workers 3 --threads 3
```

可单独使用 prepare/extract/assemble/train/summarize 阶段；每条记录的 SASA 缓存具有输入 SHA256 和球面点数检查，可中断后恢复。3 个提取进程各用一个线程；拟合顺序执行，每次限制 3 线程。原 Static、Full、Intra 预测直接复用其真实 OOF，不重新拟合或用五模型平均替换。

## 结构特征口径

主表键为 `sample_id + residue_index`。后者是冻结解析器/OOF 中 0 起始的记录内位置，不是作者残基编号；作者位置保留为 `chain + author_residue_id`，包括插入码，空链按原解析器写成 `_`。完整复合物和逐链独立输入使用**同一组坐标**；删除其他链后计算的 SASA 不代表真实游离态结构或游离态模拟。

SASA 采用 MDTraj Shrake–Rupley，1.4 Å 水探针，每原子 960 球面点，所有被冻结解析器认可的蛋白残基重原子参加遮蔽；排除氢、其他配体、水和离子。侧链排除 N/CA/C/O/OXT/OT1/OT2。`sidechain_sasa_complex_A2` / `residue_sasa_complex_A2` 为完整结构的侧链/整残基面积；`*_isolated_A2` 为原链在同坐标下单独出现时的面积；`delta_sidechain_sasa_A2` / `delta_residue_sasa_A2` 为 isolated−complex，故正值代表伙伴遮蔽。单链输入直接复制面积，差值严格为零。缓存保留负差值截断前的最小值和计数；小于 −0.001 Å² 会失败，轻微负值才容许截断。

原 47 个特征按原名保留，接触数的原始特征是 log(1+count)。另有解码后的 `intra_ca_contacts_6A/8A/10A/12A` 与 `partner_ca_contacts_6A/8A/10A/12A`。链内接触排除序列相隔两个以内的位置；伙伴接触只计其他链。邻域距离使用原 Static 的严格 `<8 Å` 等定义；`interface_ca8` 为伙伴 Cα 8 Å 接触至少一个。它不区分同源链和异源链的生物身份。链长、离端位置、主链 φ/ψ 正余弦和 N–Cα–C 夹角等也随表保留。

SASA 使用完整侧链原子，而 FlexPPI 与原 47 特征 Static 主要只使用主链 N/Cα/C。新增基线是**获得更多结构输入的解释性强基线**，不能称作同信息预算的神经网络消融。来源、模拟温度、record ID、生物系统 ID、fold、MD 或其他预测值不作为基线输入。

## 完整性审计

`sidechain_complete` 按每种标准氨基酸的全部预期侧链重原子集合检查；`chi_input_complete` 单独检查定义所有 χ 所需原子。另保存预期/实际侧链重原子数、非标准残基、丢弃 altloc 原子计数。SASA 为零时可结合完整性判断它是否确为完整侧链的埋藏状态；不能将缺原子得到的小面积解释为埋藏。逐记录 QA 在 `quality_priority83.csv` / `quality_all.csv`。敏感性评价额外只保留完整侧链，不依靠筛选结果改变模型配置。

83 条优先记录已核验：38,599 行，33,493 个有效 Hχ；完整侧链缺失、有效 Hχ 输入缺角、非标准残基、altloc 丢弃和负 ΔSASA 计数均为零。1,453 个有效 Hχ 位点初始侧链面积为零且侧链完整。

## 基线拟合

预先固定三个输入集：AA+SASA（中心氨基酸种类和六个 SASA/遮蔽面积特征）、AA+SASA+contacts（再加原 Static 的链内和伙伴 Cα 接触特征）、Static47+SASA（原 47 特征加六个面积特征）。各自对 RMSF 与 Hχ、五个原始 OOF 留出折独立拟合 Histogram Gradient Boosting。配置与原 Static 一致：120 轮、learning_rate=0.07、max_leaf_nodes=31、min_samples_leaf=50、l2=1、max_bins=63、seed=7，禁用 early stopping；未用测试表现选择特征或配置。每次只在其余四折训练，按与原 Static 相同的来源/系统/记录/有效残基平衡权重拟合。类别分箱、缺失值处理等由训练数据内完成。

统计保存逐记录 MAE、RMSE、Spearman、偏差；总体使用记录等权汇总，并按来源、界面和完整侧链子集拆开。需要 PDB/家族整簇置信区间和同环境配对的分支可直接读取 `oof_predictions_all.parquet`，不要把所有残基当作独立生物重复。这里没有额外拟合家族划分版本，因为现有家族 Full/Intra 逐位置缓存尚不可用，无法公平对照。

## 文件

- `cohort_records.csv/json`：3063 记录、原始路径、SHA256、fold、表内位置区间和 83 优先标记。
- `feature_contract.json`：列和测量口径；`base_residues.parquet`：身份、真值、原 OOF 和原 47 特征。
- `features_priority83.parquet`：83 条优先全位点表；`features_all.parquet`：3063 条全位点表。
- `cache/*.npz`：单记录 SASA、完整性和遮蔽差原始审计。
- `models/*`：逐特征集/折/目标模型与留出预测；`*_oof.npz`：新基线全位点 OOF。
- `oof_predictions_all.parquet`：真值及全部比较方法的统一键预测表。
- `per_record_metrics.csv`、`record_macro_*.csv`：逐记录与记录等权汇总。
- `baseline_contract.json`、`training_progress.json`、`sasa_progress.json`、`run.log`：配置、进度和可恢复执行记录。
