from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'label_sensitivity/runtime'))
import pandas as pd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
model=json.loads((ROOT/'model_audit/verified_summary.json').read_text(encoding='utf8'))
infer=json.loads((ROOT/'fresh_inference_audit.json').read_text())
numeric=json.loads((ROOT/'numeric_rescreen/independent_numeric_QA.json').read_text())
md=json.loads((ROOT/'replica_recheck/all_rechecked_complexes_audit.json').read_text(encoding='utf8'))
rep=pd.read_csv(ROOT/'promising_cases_replica_validation/recomputed_replica_complete_region_summary.csv')
assert model['ITrace_train_overlap']==[] and model['strict_ITrace_overlap']==[]
assert model['all_ITrace_in_fold0_test'] and infer['forward_passes']==84
assert md['status']=='passed' and md['existing_MD_trajectories']==57
main=[]
for name,a,b,scope,metric in [
 ('KRAS G12D：未替换肽位点主链','7ow5','7ow6','whole_peptide_unchanged_AA','rmsf'),
 ('DMF4：完整pMHC界面侧链','3qdm','3qeq','natural_pMHC_interface','hchi'),
 ('DMF5：同组肽变化的受体对照','3qdg','3qdj','natural_pMHC_interface','hchi'),
 ('MEL5：IMP2→BST2完整CDR3主链','7q99','7q9a','CDR3','rmsf'),
 ('MEL5：Melan-A→IMP2对照','4qok','7q99','CDR3','rmsf'),
 ('MEL5：Melan-A→BST2对照','4qok','7q9a','CDR3','rmsf'),
 ('PRAME：HLA-A3→A11界面侧链','8rym','8ryn','natural_pMHC_interface','hchi')]:
    f=rep[rep.a.eq(a)&rep.b.eq(b)&rep.scope.eq(scope)&rep.metric.eq(metric)].copy()
    assert len(f)==3,(name,len(f))
    f.insert(0,'candidate',name);main.append(f)
selected=pd.concat(main,ignore_index=True)
selected.to_csv(ROOT/'优先候选_三重复复核与双模型结果.csv',index=False,encoding='utf-8-sig')

paper=ROOT.parent.parent/'paper_v40_20261017_中文修订_TCR识别环境变化实例_三线表.docx'
paper_sha=sha(paper)
assert paper_sha=='bd9851a0f6dabdeb6f4a7c2fbbe9fe0dd20ec0820a156578db7f624cd990f5ca'
fwd_errors=[r.get('max_error_vs_existing_cache',0) for r in infer['predictions'] if r['model']=='family_full_no_ITrace']

text='''# ITrace内部实例搜索：未用ITrace训练的模型复测（2026-10-17）

有合适的ITrace内部候选。本轮最值得优先推进的是KRAS G12D（7OW5/7OW6）的完整肽段主链响应；侧链方向优先考虑DMF4（3QDM/3QEQ）的pMHC接触面状态变化。它们的优势是生物学对照更直接，而且两套从未用ITrace训练的已有模型都能在完整预定区域上降低变化预测误差。MEL5识别IMP2/BST2（7Q99/7Q9A）也有可重复的CDR3主链证据。以上都是候选发现结果，不改论文，也没有删改原有SB27及相关位点内容。

## 使用了什么模型，是否真的未训练ITrace

已逐条审计实际训练记录和远程检查点。主模型为联合RMSF/χ联合熵Full模型，训练数据为ATLAS 1866条与DynaRepo 564条，ITrace为0条，验证集为空；244条ITrace完整留作测试。更早的strict Full模型训练2431条，同样没有ITrace。两者均为已有模型，没有重新训练、微调或按本次结果更换权重。原先普通交叉验证模型中有195/196条ITrace参与训练的版本未用于本轮。模型和预测的SHA、实际清单比较见model_audit目录；主模型为ff297440…15ec5，旧strict为b3eed060…2b39f。

本轮先按原文中的明确生物差异列出12组系列、17个真实对照，再固定完整肽段、全部六条识别环及界面等区域统一测试。额外加入系列内对照和保留实例，共重新推断42个结构、两模型84次，整理67对不重复结构比较。已从原始文件重读19个结构的57条现有MD轨迹，以原标签方法复算均值、方差并核对；所有源拓扑、轨迹哈希保持不变，没有新启动MD。更早结果中已有ITrace候选，这次的新增工作是补齐生物学关系、实际训练来源审计，并用旧strict模型平行复测。

## 优先候选及能论证的能力

**第一优先：KRAS G12D，7OW5→7OW6。** 同一工程TCR、同一HLA-A11，分别结合野生型VVVGAGGVGK与G12D肽VVVGADGVGK，只改变肽第6位；原文有明确的识别/亲和差异。因此这里可问一个具体的模型问题：一个肽位点的改变，能否使模型读出其它未替换位点的主链运动变化。[原始研究](https://www.nature.com/articles/s41467-022-32811-1)。

完整10肽的主链变化相关在主模型/旧strict分别为0.806/0.879，相对于全部预测不变的误差改善38.9%/60.3%。更关键的是，去掉G6/D6本位后，剩余全部9个未改变氨基酸的位置仍有相关0.733/0.833、误差改善40.2%/62.8%；主模型8/9、旧strict9/9个变化方向正确。9位的MD平均RMSF下降0.303 Å，模型分别下降0.135与0.220 Å。两种状态各3次MD中，有7/9位的取值范围完全分开。这里的正面证据覆盖其余整条肽，不是由突变本位或一个较好的残基支撑。Intrachain控制在同9位的相关为−0.083，误差比预测不变增加66.1%；它虽将方向都判为下降，却明显高估幅度，说明“方向都对”也不足以算准确。

这个实例适合论证“模型能从复合物结构读出单点抗原变化伴随的间接主链响应”。它尚不能把侧链高相关一起作为强证据：未变肽位点只有5个有效熵位置，变化均未超过两倍三重复标准误，而且旧strict的侧链相关不为正。原文的亲和力差异提供生物学背景，模型预测的是局部运动幅度，不是结合亲和力数值。7PB2亲本受体在突变肽上的结构作为另一条工程对照保留，未与固定TCR的WT/突变肽比较混合。

**第二优先：DMF4识别MART-1十肽/九肽，3QDM→3QEQ。** 同一DMF4和HLA-A2识别ELAGIGILTV与AAGIGILTV，原文观察到不同识别几何和结合强弱，适合检验未替换的pMHC界面侧链状态重排。[原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC3166883/)。完整pMHC一侧界面19个可自然对应、氨基酸未替换且可计算χ熵的位置，主模型/旧strict的变化相关为0.779/0.553，变化误差改善25.2%/27.1%；7个位点的变化超过两倍三重复标准误，方向分别为7/7与5/7正确。肽长不同，源表明确保留序列对应和未对齐位置，不能将二者逐编号生硬相减。并行DMF5的同组十肽/九肽对照3QDG→3QDJ未显示同样的预测准确性，完整保留，故此例论证限定于DMF4的侧链环境响应。

**第三候选：MEL5识别IMP2/BST2，7Q99→7Q9A。** TCR和HLA相同，抗原不同；两套模型对全部22个CDR3位点的主链变化相关为0.713/0.553，误差改善42.9%/41.8%，主模型20/22、旧strict17/22个方向正确。MD平均变化为+0.112 Å，模型为+0.048/+0.104 Å；它并不是“识别更强必然变得更硬”的简单故事。Intrachain的相关仅0.226、误差改善18.2%。但Melan-A第三抗原4QOK的两项对照不能省略：4QOK→7Q99两模型都较差，4QOK→7Q9A主模型虽然相关0.824，误差却比预测不变增加11.0%。因此这一对值得作为具体主链案例，尚不足以写成MEL5跨三种肿瘤抗原都准确。[原始研究](https://pubmed.ncbi.nlm.nih.gov/37490916/)。

**PRAME跨HLA呈递环境可作后备侧链例。** 固定S2受体和ELFSYLIEK肽、由HLA-A3换到A11（8RYM→8RYN），完整同AA pMHC界面15个熵位置相关为0.486/0.686、误差改善14.3%/24.1%；它有环境响应信息，但Intrachain也表现相近，未显示完整复合物模型明显更强。S2-198工程对照8RYO一并保留。[原始研究](https://www.nature.com/articles/s41598-024-75140-7)。

MAGE-A10亲本/单点变体、AHIII12.2三种单MHC变化、1G4 C9V、NRAS K/R和KRAS G12V，虽然生物学问题清楚，当前模型没有显示更有力的完整区域结果。1G4 CDR2工程在主模型上表现较好，但旧strict未复现，且环外构建差异可能影响按整链拟合的RMSF，暂不列为优先成功实例。完整负结果与边界问题均在源表保留。7QPJ不是7PDW的第二代工程体、7Q9B是MEL8而非MEL5，这些对应关系已按原文纠正。

## 统计口径与文件

相关衡量模型能否排出各位点变化大小；变化平均绝对误差同时约束幅度，基线为对所有位点均预测不变。主链单位Å，侧链联合熵单位nats。替换本位会改变侧链化学类型，因此同AA侧链另算，KRAS主链也额外报告全部9个未变位点。每状态3条独立MD是重复单位；残基或帧不作为独立实验重复。“大于2SEM”及“三重复范围分开”仅用于描述重复波动，未据此作显著性检验。以上属于有标签的探索性实例搜索，不将它当成另一个事先完全锁定的泛化性能测试。

完整候选文献与生物关系见biology_first/biology_first_ITrace_candidate_series.csv；12组17个真实对照的五个完整区域、两模型结果见numeric_rescreen/biology_first_12_series_17_contrasts_all5regions_2models.csv；全部8160项区域/指标/模型统计与逐位点数据在numeric_rescreen/。优先候选及三重复源数据见“优先候选_三重复复核与双模型结果.csv”和promising_cases_replica_validation/，原始MD重读审计在replica_recheck/。论文文件SHA与本轮状态保存在search_only_completion_status.json。
'''
report=ROOT/'ITrace内部候选_双模型复测结论_不改论文.md'
report.write_text(text,encoding='utf8')
status={'status':'candidate_search_and_existing_model_measurements_complete','manuscript_modified':False,
 'SB27_existing_content_preserved':True,'manuscript_path':str(paper),'manuscript_SHA256':paper_sha,
 'all_ITrace_in_test':True,'primary_ITrace_training_records':0,'old_strict_ITrace_training_records':0,
 'biological_series':12,'direct_biological_contrasts':17,'structures_fresh_inferred':42,'existing_models_fresh_inferred':2,
 'new_forward_passes':84,'all_unique_PDB_pair_comparisons':67,'region_metric_model_rows':8160,
 'existing_MD_trajectories_rechecked':md['existing_MD_trajectories'],'new_training':0,'new_MD_production':0,
 'max_primary_fresh_prediction_error_vs_previous_cache':max(fwd_errors),
 'prioritized_backbone_case':['7OW5','7OW6'],'prioritized_sidechain_case':['3QDM','3QEQ'],
 'secondary_backbone_case':['7Q99','7Q9A'],'report':report.name,
 'conclusion':'ITrace contains informative biological contrasts. KRAS G12D backbone and DMF4 pMHC-interface side-chain responses reproduce positive change-MAE improvements in both models; neither is a uniformly strong two-output result across every region.'}
(ROOT/'search_only_completion_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8')
files=[report,ROOT/'优先候选_三重复复核与双模型结果.csv',ROOT/'search_only_completion_status.json',ROOT/'analysis_plan.json',ROOT/'fresh_inference_audit.json',ROOT/'model_audit/verified_summary.json',ROOT/'promising_cases_replica_validation/replica_comparison_audit.json',ROOT/'numeric_rescreen/independent_numeric_QA.json']
(ROOT/'report_source_manifest.json').write_text(json.dumps({str(p.relative_to(ROOT)):{'SHA256':sha(p),'bytes':p.stat().st_size} for p in files},ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(status,ensure_ascii=False,indent=2))
