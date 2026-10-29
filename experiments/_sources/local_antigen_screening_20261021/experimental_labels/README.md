# 本地抗原排序实验：真实SPR标签

本目录于任何新面板FlexPPI推断结果被读取前，依据原始实验文献和HERMES公开数据整理。只有两个独立TCR系列：1G4–NY-ESO-1为14个唯一肽，A6–Tax为16个唯一肽，共30个候选。1G4来自Pettmann2021的9个肽和Aleksic2010的11个肽，重合6个，因此有20条研究级记录；加A6的16条，共36条已发表的“肽–研究均值”。这些不是36个独立候选，也不是原始技术重复。

## 使用文件

- **SPR_candidates.csv**：30个候选的无实验标签结构输入清单。`candidate_id`稳定唯一，包含肽、TCR、HLA和标签无关的最近序列模板。
- **SPR_labels_candidates_primary.csv**：主终点30行。`delta_pKD_vs_parent`是先在每个研究内扣除其自身测得的母肽pKD，再取跨研究中位数；两个研究时等于等权均值。正值表示比母肽亲和力更高。`study_delta_pKD_json`保留各研究结果，应逐研究作敏感性分析。
- **SPR_labels_observations.csv**：主终点的36条研究级记录，含单位、37°C、误差类型/大小、n、母肽、源表位置和哈希。`parent_same_measurement_batch=not_reported`：只能确认同研究、同协议和温度，不能声称逐次SPR实验严格配对。
- **SPR_labels_observations_reported_hybrid.csv / SPR_labels_candidates_reported_hybrid.csv**：主文件的明确版本名。
- **SPR_labels_observations_all_constrained.csv / SPR_labels_candidates_all_constrained.csv**：HERMES兼容的完整敏感性。
- **SPR_labels_observations_all_fitted.csv / SPR_labels_candidates_all_fitted.csv**：传统Bmax拟合的完整敏感性。三种拟合版本不能扩大样本量。
- **SPR_labels_HERMES_aggregate_secondary.csv**与旧兼容别名**SPR_labels_candidates.csv**：仅复现HERMES把不同研究的绝对-logKD均值合并的次分析；不要用于主终点。`SPR_labels_study_level.csv`是原HERMES全约束列的36条独立来源记录。
- **SPR_labels.json**：以上所有明确命名版本；`observations`为主版本，`HERMES_aggregate_secondary`为次版本。**templates.json**记录实际PDB所有链的序列及作者残基编号。

## 预先固定的实验定义

主版本遵循Pettmann2021原文Methods中的报告规则：约束Bmax得到的KD大于20 μM时使用该约束值，否则使用传统Bmax拟合值。该选择依据文献，不依据模型预测，并经评估协议代理在看到新预测之前固定。HERMES仓库则统一使用约束列，因此另存全约束敏感性。Aleksic2010的常规SPR均值在三版本中不变。规则及原文位置见`label_rule_freeze.json`。

两个1G4研究均有SLLMWITQV（9V）母肽，A6有LLFGYPVYV母肽。`pKD_M = -log10(KD/M) = 6 - log10(KD/μM)`。HERMES图中的`-log10(KD/μM)`与pKD相差6，不能混用。更重要的是，HERMES的closest/averaged CSV在均值处理后仍保留第一研究的KD列，该列可能不等于合并后的-logKD；本目录复现次表已经重算几何平均KD，主分析按研究内母肽差值处理。

原表全部给出数值点估计，没有`>`或`<`删失项，因此`censoring=none_reported`，检测上限留空。不能将最大的1309或2321 μM擅自改成检测界限。Pettmann表中部分A6肽仅n=1，原表SD为N/A；本目录保留缺失，未伪造方差。Aleksic表为至少3次独立实验的均值±SEM，不能与Pettmann的SD直接混合。

## 源文献与模板

HERMES正式版：doi:10.1073/pnas.2504783122，PMC12557503；作者仓库`StatPhysBio/tcr_antigen_design`固定提交`cfd85f3fefa652f0f00d3fd1e58cb2b4256f5ce2`。已保存正式全文XML、39页补充PDF、Dataset S1/S2和原CSV。首次PMC普通HTML请求遭遇reCAPTCHA，其文件标记为blocked且未用于任何证据；正式XML及补充通过Europe PMC公开API正常取得。

Pettmann et al., eLife2021，doi:10.7554/eLife.67092。直接来源为`Figure 1—source data 2`（1G4）和`source data 3`（A6），公开v3 CSV第7–10列为约束Bmax，第3–6列为传统拟合；分别保留均值、SD和n。源表逐行与HERMES原始CSV核对36条（含下述11条），全部数值一致。

Aleksic et al., Immunity2010，doi:10.1016/j.immuni.2009.11.013。来源为补充`Table S1`、PDF第16页上表，KD单位μM，均值±SEM，至少3次实验；11条固定WT HLA的肽序列及数值已逐项核对。原表另有6种HLA替换，已完整保存在`Aleksic2010_TableS1_complete_pMHC_panel.csv`，不混入固定HLA抗原面板。该PDF页面已实际打开目视核验。

真实模板序列经坐标核对：**2BNQ=SLLMWITQV，2BNR=SLLMWITQC**；PNAS主图2A文字将两编号反写，但其补充表S1、仓库constants和PDB实际序列一致，本目录采用后者。A6主模板为**1AO7=LLFGYPVYV、1QSF=LLFGYPVAV**。最近模板按肽Hamming距离选择，同距按预先列定的顺序，完全不看SPR值；主固定骨架测试仍可统一使用2BNQ和1AO7。完整PDB可能有多个晶体复合物，建模方须明确选用一个生物组装，不把所有晶体复制体当作伙伴。

## 核验与复现

`source_manifest.json`记录保存的原文件SHA-256和来源，`label_validation_report.json`记录36项原始数值匹配与30项模板/合并复核，`label_rule_freeze.json`记录主规则和三版本哈希。已打开原文两个关键PDF页面检查原单位、表脚注和PDB序列。所有新写入均限于本目录，未运行模型推断、训练或MD，也未读取模型输出挑选候选。

完整重新取得与重建顺序：`retrieve_sources.py`→`fetch_panel_files.py`→`fetch_primary_data.py`→`build_labels.py`→`build_label_variants.py`。已有raw文件时仅需后两步；`build_labels.py`直接核验原表。联网获取使用agent-reach的GitHub/API及原始网页流程；无网页登录或联系作者。

只有两个独立TCR系列，且均为有数值亲和力的已研究配体；这里未设置任意二分类“命中”阈值。候选间排序和主终点可用于探索性验证，不能把技术重复、两种拟合或36条研究均值当作更大的独立筛选验证集。
