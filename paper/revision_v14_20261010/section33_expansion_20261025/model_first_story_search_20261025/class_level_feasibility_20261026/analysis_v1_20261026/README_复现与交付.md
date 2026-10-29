# 类别分析复现与交付

本文件夹保存2026-10-26功能肽识别类别分析。先读[正式结果与论文取舍](类别分析_正式结果与论文取舍.md)，再查源表。这里没有改写论文Word，没有重新运行神经模型或MD生产模拟；本轮使用已冻结的家族整组留出神经预测，重新训练了结构对照，并从已有轨迹重新计算统一窗口参考。

## 直接检查结果

`biological_annotation/all54_PDB_biological_annotations.csv`保存54PDB的分类、核心成员标志和原始依据；`formal_results/cohorts.json`保存全部敏感性队列。`formal_results/per_record_all_methods_references.csv`保存全部模型、参考、组分和区域的记录级指标；`summary_all_scopes.csv`为汇总，`paired_MAE_gain_cluster_CI.csv`为成对增益及区间。`sameAA_match_type_equal_cluster_CI.csv`包括严格匹配阈值、覆盖数和阴性结果。`figures/figure_all480_values.csv`及`figure_macro_means_CI.csv`与图逐点对应。

`controls/training_progress.json`保存20次新训练的完成状态；`baseline_contract_family.json`保存参数、特征和输入哈希；`models/`同时包含训练模型和预测缓存。`uniform_md/compute_summary.json`保存53条成功记录、3条时间不足记录和掩码检查；`record_audit/`及`per_record/`保存逐体系对齐和计算记录。`independent_statistics_review/`为独立核算，`audit/`为最终交付检查。`PROTOCOL_类别分析.md`说明固定口径及探索性质。

## 从交付数据重算统计和图

使用Python 3.12，安装numpy、pandas、scipy、pyarrow、matplotlib、Pillow。用于原始轨迹提取的环境另含MDAnalysis、mdtraj；重训结构对照另需scikit-learn、joblib、psutil、threadpoolctl。确切scikit-learn版本见训练合同。图字体为Arial。

PowerShell示例；把`$analysisDir`改为实际解压目录即可。变量只是给这一分析文件夹起一个简短名字，`& $pythonExe`表示调用该Python解释器。

```powershell
$analysisDir = 'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\class_level_feasibility_20261026\analysis_v1_20261026'
$pythonExe = 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe'
& $pythonExe "$analysisDir\formal_root_evaluation.py"
& $pythonExe "$analysisDir\figures\plot_class_model_comparison.py"
& $pythonExe "$analysisDir\audit_delivery.py"
```

`formal_root_evaluation.py`使用已交付的`controls/all56_features_and_family_predictions.parquet`、生物注释和统一标签表，重算主表与严格匹配结果，因而不必重新处理轨迹或训练。绘图脚本只读正式记录表。绘图后QA标志会重置为pending，因为新渲染文件应重新查看；原交付版本的人工审核与哈希见`figures/QA_人工审核.md`及`audit/final_delivery_audit.json`。

## 从外部原始输入重做前处理与训练

`annotate_cohort_from_primary_metadata.py`使用文件夹内官方元数据和原始摘要缓存重建分类。`uniform_md/preflight.py`、`extract_uniform_labels.py`使用本地原始轨迹；`controls/train_family_baselines.py`使用项目原有全量结构特征和3063条family划分；`evaluate_controls.py`组装候选残基表。这些脚本保留原机器绝对路径及外部数据依赖，移动机器时需按合同和轨迹清单改路径。解压包本身不是完整外部训练/MD数据集的镜像。

原始神经预测来源为`C:\Users\DELL\.codex\worktrees\2a52\FlexPPI\output\non_tcr_case_screen_20261025\archive_extract\family_folds\evaluation\oof_per_residue.csv.gz`，SHA256为`0852469b8fb6dd2977796d040ce8bc7fe04f3de29f217e7fc563213c064ad7f3`；候选部分另保存于`shared/oof_candidates56.parquet`。原始结构特征目录为`E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\feature_baseline`，family划分目录为`E:\Project\FlexPPI\experiments\_sources\family_folds`。外部原始轨迹及时间轴来源、大小、修改时间和哈希在`uniform_md/trajectory_inventory.csv`及各记录审计中。

压缩包包含统计复现所需的本轮源表、脚本、模型、缓存、图和审计；不包含数GB的原始轨迹或全量神经训练数据，未包含临时_staging和__pycache__目录。不能把压缩包称为从原始序列起的独立端到端复现包。MD标签是计算参考而非新实验测量，分块时窗不是独立重复，Full/Intra不是游离态/结合态的状态对照。
