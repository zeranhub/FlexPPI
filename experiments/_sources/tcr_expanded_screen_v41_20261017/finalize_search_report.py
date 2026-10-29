from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012'/'plot_runtime'))
import numpy as np
import pandas as pd
read=lambda p:pd.read_csv(p,float_precision='round_trip')
OLD=ROOT.parent/'tcr_case_v40_20261017'
old=read(OLD/'change_screen'/'all_candidate_change_metrics.csv')
new=read(ROOT/'construct_tolerant_screen'/'expanded_candidate_pairs.csv')
pairs=lambda d:set(tuple(sorted([a,b])) for a,b in zip(d.a,d.b))
counts={'old_distinct_pairs':len(pairs(old)),'expanded_distinct_pairs':len(pairs(new)),
        'new_pairs_relative_to_all_old_candidates':len(pairs(new)-pairs(old)),
        'union_distinct_pairs':len(pairs(new)|pairs(old)),
        'held_out_ITrace_complexes':244,'newly_reprocessed_existing_MD_trajectories':38,
        'new_model_forward_passes_this_search':8,'new_training_runs':0,'new_MD_production_runs':0}
assert counts['old_distinct_pairs']==251 and counts['union_distinct_pairs']==343
rows=[]
def append(case,a,b,region,metric,n,rho,mae,zero,gain,source,note=''):
    rows.append(dict(case=case,PDB_before=a,PDB_after=b,region=region,metric=metric,n_positions=int(n),
                     change_spearman=float(rho),change_MAE=float(mae),zero_change_MAE=float(zero),
                     relative_MAE_improvement=float(gain),source_csv=source,note=note))

mel=read(ROOT/'mel5_validation'/'corrected_predefined_region_change_metrics.csv')
mel_scopes=['TCR_CDR_all_aligned','TCR_CDR_common_AA','CDR3_all_aligned','CDR3_common_AA',
            'CDR3α_all_aligned','CDR3β_all_aligned','MHC_plus_peptide_interface_common_AA',
            'MHC_interface_common_AA','all_interface_common_AA','TCR_interface_common_AA']
for _,r in mel[mel.scope.isin(mel_scopes)].iterrows():
    append('MEL5 → α24β17','3HG1','4JFF',r.scope,r.metric,r.n,r.rho_delta,r.MAE_delta,
           r.zero_change_MAE,r.relative_MAE_gain,'mel5_validation/corrected_predefined_region_change_metrics.csv',
           'All-aligned scopes retain substitutions; same-AA entropy scopes isolate unchanged residue chemistry.')
pik=read(ROOT/'pik3ca_validation'/'predefined_region_change_metrics.csv')
for _,r in pik[pik.scope.isin(['all_interface_common_AA','MHC_plus_peptide_interface_common_AA','peptide_common_AA','CDR3_common_AA'])].iterrows():
    append('PIK3CA two TCRs','7L1D','7RRG',r.scope,r.metric,r.n,r.rho_delta,r.MAE_delta,
           r.zero_change_MAE,r.relative_MAE_gain,'pik3ca_validation/predefined_region_change_metrics.csv',
           'MHC/peptide natural positions are primary; cross-TCR matching contributes strongly to backbone improvement.')
jm=read(ROOT/'jm22_f5l_validation'/'model_MD_change_summary.csv')
for _,r in jm[jm.scope.isin(['all9','common8'])].iterrows():
    append('JM22 F5L','1OGA','5HHM',r.scope,r.metric,r.n,r.rho_delta,r.MAE_delta,
           r.zero_change_MAE,r.relative_MAE_gain,'jm22_f5l_validation/model_MD_change_summary.csv',
           'all9 includes F5L replacement; common8 excludes it. Most unchanged-site changes are underpredicted.')
catch=read(ROOT/'catchbond_validation'/'regional_change_metrics.csv')
for _,r in catch[catch.input_kind.eq('experimental_matched')&catch.region.isin(['all_unchanged_aa','peptide_all9'])].iterrows():
    append('TCR156 WT → S32M','9NMU','9NMX',r.region,r.metric,r.n,r.rho,r.change_MAE,
           r.zero_change_MAE,r.relative_MAE_reduction,'catchbond_validation/regional_change_metrics.csv',
           '10 independent existing MD replicas per state. Not a stronger illustrative case.')
for _,r in old[old.a.eq('2ak4')&old.b.eq('4jry')&old.scope.eq('peptide_unchanged')&old.method.eq('full')].iterrows():
    append('SB27 → SB47 (retained)','2AK4','4JRY',r.scope,r.metric,r.n,r.rho_delta,r.mae_delta,
           r.mean_abs_md_delta,1-r.mae_delta/r.mean_abs_md_delta,'../tcr_case_v40_20261017/change_screen/all_candidate_change_metrics.csv',
           'Existing manuscript case retained; Q7 local evidence is stronger than the whole-peptide entropy result.')
metrics=pd.DataFrame(rows)
metrics.to_csv(ROOT/'verified_candidate_comparison.csv',index=False,encoding='utf-8-sig')
(ROOT/'expanded_search_counts.json').write_text(json.dumps(counts,ensure_ascii=False,indent=2),encoding='utf8')

# Preserve the superseded narrative as an audit trail; current interpretation uses corrected alignment.
judgment=ROOT/'mel5_validation'/'复算判断.md'
historical=ROOT/'mel5_validation'/'旧对齐结果_已被修正取代.md'
if not historical.exists():
    historical.write_text('此文件保存早期对齐下的历史判断，数值已经被corrected_*和当前复算判断.md取代，不得用于论文。\n\n'+judgment.read_text(encoding='utf8'),encoding='utf8')
judgment.write_text('''# MEL5 / α24β17：修正对齐后的复算判断

本判断使用corrected_pair_all_aligned_positions.csv、corrected_predefined_region_change_metrics.csv和plot_sources，取代旧SequenceMatcher结果。4JFF的两条TCR链相对MEL5各多一个N端残基；六环统一采用WT边界并通过全链全局对齐映射，避免β链CDR编号错一位及工程替换片段被跳过。全长对应位置824个；TCR共有19个替换，其中六环内17个，与原始研究报告的19个序列差异一致。全部六环45个主链位置保留，完整CDR3含α链9位、β链13位，合计22位，其中8位为工程替换。先前14位ρ=0.938的结果已废止。

全部六环45位的主链变化ρ=0.440，平均绝对误差0.211 Å，相对于预测不变的0.287 Å改善26.6%；完整22位CDR3的ρ=0.779，误差由0.292降至0.196 Å，改善32.8%。CDR3平均MD波动降低0.292 Å，模型降低0.097 Å，幅度仍偏小。两个CDR3并不均衡：α链9位ρ=0.467，误差改善6.3%；β链13位ρ=0.786，误差改善42.4%，13个位置的下降方向均被预测正确，其中7个位置两组各三次MD的取值范围完全分开。CDR2α的主链变化预测较差，图中必须保留该完整区域，不能把合并CDR3的较高相关写成六环均准确。

去掉CDR3的工程替换位点，同AA主链14位ρ=0.829，误差改善40.3%；但同AA侧链仅8位，ρ=−0.048，误差比预测不变增加8.7%。包含替换的全部CDR3侧链ρ=0.802，明显受到侧链化学类型改变的影响，不能用于证明环境敏感的侧链预测。TCR界面同AA侧链ρ=−0.467、误差增加33.1%，也不支持这一表述。

另一个较强结果位于未改变序列的pMHC一侧。按任一结构中Cα距离≤8 Å定义界面并取并集，保留全部30个自然对应位置，其中23个有有效χ角。侧链熵变化ρ=0.755，平均绝对误差0.202 nats，相对于预测不变的0.269 nats改善24.8%；19/23个变化方向正确。MHC单独18个熵位置ρ=0.705、误差改善23.4%。K68的MD熵由3.282降至3.047 nats（−0.235），模型由3.017降至2.769（−0.247）；三条工程TCR轨迹的K68熵均低于三条WT轨迹。R65下降幅度被低估，K66和Y159的增加被明显低估，均保留在图中。整个肽没有重复支持充分的熵变化，不把它包装为主要变化区域。

新增模型输入检查：保持每个结构的pMHC坐标不变，仅移去TCR后，23个界面侧链的变化相关由0.755降至0.057，误差变为0.271 nats，与预测不变的0.269相当。K68预测变化由−0.247变为+0.036 nats。完整复合物预测与原缓存最大差1.99×10⁻⁵，重复推断可复现。该检查说明这组侧链差异预测使用了TCR上下文；它没有模拟真实解离，也不证明熵变化导致亲和提高。

作为论文实例，MEL5比SB27/SB47有更明确的工程前后关系，并将同AA侧链预测证据从单个位点扩展到完整pMHC接触面，因此值得作为新增实例。它尚未显示所有区域、两项输出都比SB27更准确；应保留SB27现有实例，让两例承担不同论证。MEL5可以论证“模型识别到工程后TCR主链约束和未改变序列的pMHC侧链状态重排，并在后一项上使用结合伙伴信息”，不把局部构象熵等同于结合熵、自由能或亲和力。

每组均使用已有3条200 ns MD的所有保存帧，无新增MD或训练。全蛋白原始复算均值、方差与缓存误差<5×10⁻⁷，源拓扑和轨迹SHA未改变。独立复核脚本从原缓存、逐重复NPZ和PDB坐标重建824个对应位置、六环及8 Å界面，通过全部核对。重复数指独立MD轨迹，不将残基或帧当作独立重复。
''',encoding='utf8')

report='''# TCR–pMHC新增实例搜索与复核（2026-10-17）

本次搜索找到最值得新增的候选是MEL5→α24β17（3HG1/4JFF）。它更有价值的地方是：具有清楚的工程前后关系，而且能在未改变序列的完整pMHC接触面上检验侧链变化，新增输入检查也直接显示模型使用了TCR信息。它并非在所有区域和所有指标上都超过SB27/SB47。因此保留论文原有SB27及其位点内容、图6和正文，新增结果单独归档。本次没有修改v40稿件，也未调整正式交付文件指针。

## 搜索范围与方法

本次扩展筛查使用当前联合RMSF/侧链熵模型，现有244个ITrace复合物均在该模型的家族划分测试集中。原有候选包含251个不同PDB对；放宽恒定区和末端构建差异、核对可变区及MHC沟槽后，扩展筛查得到313对，其中92对是原有251对之外的候选，两轮合并为343个不同结构对。这是候选发现范围，不是343项独立生物实验，也不是预先锁定的测试集。另盘点ATLAS、DynaRepo、ThermoPCD、DynaTPH及本地CatchBond数据。深度复算JM22、PIK3CA、MEL5和CatchBond，共重新分析38条已有MD轨迹；新增8次模型前向推断，没有新增训练或MD生产模拟。

比较采用变化量，即后状态减前状态。ρ衡量模型能否分出哪些位点变化较大、哪些较小；平均绝对误差还检查预测幅度是否接近MD。“误差改善”以所有位点均预测不变为基线，例如32.8%表示模型的变化预测误差比一律预测不变减少32.8%。侧链指标是多个χ角的联合旋转状态熵，单位nats。替换成另一种氨基酸本身会改变可用旋转状态，因此同AA位置另行比较，避免将化学类型差异当成环境预测能力。各区域及全部位置数据均保留，候选之间区域不同，不能把下表当成公平的模型排名。

## 候选比较

| 实例 | 完整比较区域 | 主链变化ρ / 误差改善 | 侧链变化ρ / 误差改善 | 判断 |
| --- | --- | --- | --- | --- |
| MEL5→α24β17 | 全六环45个主链位置 | 0.440 / 26.6% | 同AA六环侧链ρ=−0.038，未改善 | 工程关系清楚，主链效果有区域差异 |
| MEL5→α24β17 | 完整CDR3：22个主链位置，含8替换 | 0.779 / 32.8% | 同AA8位ρ=−0.048，未改善 | 主链证据主要来自β环 |
| MEL5→α24β17 | 同序列pMHC界面：30个主链、23个熵位置 | 0.421 / −28.7% | 0.755 / 24.8% | 最有价值的环境敏感侧链证据 |
| PIK3CA两种TCR | 同序列pMHC界面：33/26位 | 0.545 / 5.4% | 0.724 / 15.9% | 较宽界面，但主链误差改善小、重复波动较大 |
| JM22 F5L | 排除替换位点后的完整肽：8/6位 | 0.714 / 22.0% | 0.657 / 12.3% | 单点替换关系清楚，但多数变化幅度被低估 |
| CatchBond WT/S32M | 全部未替换位点：407/343位 | 0.022 / −57.6% | −0.020 / −2.8% | 不能作为更强实例 |

MEL5的原始实验报道了MART1–HLA-A2结合亲和力的显著增强，以及新增接触主要落在MHC上的结构证据；其pMHC序列相同，适合区分工程改造后的主链约束和侧链状态变化。[原始研究，JBC 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3696650/)。本次复算中，完整CDR3主链平均MD变化为−0.292 Å、模型为−0.097 Å，仍低估幅度。CDR3β的13个位点全部预测为波动下降，与MD方向一致；CDR3α仅轻度改善，CDR2α则存在反方向预测。图中保留全部六环45位，使这一差别可以直接检查。全CDR3侧链的较高相关受到工程替换影响，去掉替换后不成立，因此MEL5的侧链论证放在未改变序列的pMHC一侧。

这一侧的23个有效侧链位置中，19个变化方向正确；平均绝对误差0.202 nats，预测不变的误差为0.269 nats。最接近的示例之一是MHC K68：MD变化−0.235、模型−0.247 nats，三条工程状态MD的取值均低于三条亲本状态。模型移去TCR、保留每个结构的相同pMHC坐标后，完整23位的变化相关由0.755降为0.057，误差由0.202增为0.271 nats；K68预测变化转为+0.036 nats。由此能够推进的模型论证是：模型利用结合伙伴信息识别相同pMHC中的侧链状态重排。图中同时保留R65幅度低估、K66/Y159熵增加未能充分预测等结果。局部构象状态变化并不等于结合熵，也不能单独解释整个亲和力变化；原研究还显示溶剂贡献，因此论文将功能背景与模型实际预测的量分别说明。

PIK3CA比较21LT2-2与0606T1-2对同一新抗原–HLA的识别。其完整界面双输出ρ约0.68/0.67，但去掉跨TCR配对位置，保留自然逐位置对应的MHC和肽后，主链误差仅改善5.4%；全界面主链误差改善的90.6%来自跨TCR匹配点。因此它提供较宽的侧链变化对照，尚不足以称为更强的整体预测实例。[原始研究，Nature Medicine 2022](https://www.nature.com/articles/s41591-022-01786-3)。JM22的抗原单点变体F5L具有明确的结合实验背景；排除替换位点后两项变化相关仍为正，但F7联合熵MD增加0.495 nats、模型仅增加0.038 nats，V6/L9亦明显低估，不能只根据相关较高认定效果更好。[原始研究，PNAS 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4843436/)。CatchBond使用每种状态10次独立已有MD、共20条轨迹，在全体未替换位点和完整肽上均未优于预测不变；准备结构和实验结构两套输入结果一致，故不推进为论文实例。[数据来源](https://zenodo.org/records/17065125)。

## 可用于论文的新增论证与交付文件

SB27现有内容保留。MEL5作为另一种任务关系补充：SB27比较同一pMHC在不同TCR识别环境下的差异，MEL5比较有亲本关系的工程前后状态。后者的证据链为“工程前后结构→完整识别环主链变化→同序列pMHC接触面侧链重排→移去TCR后侧链变化预测减弱”。可采用的表述是：“在MEL5亲和成熟实例中，模型辨认出CDR3β整体主链波动降低，并在未改变序列的pMHC界面上捕捉侧链状态重排；后者依赖输入中的TCR上下文。”完整区域的吻合程度和幅度误差决定了这句话的范围，不扩展成亲和力预测。

图采用英文，保留六环45个主链位置及完整pMHC界面30个位置，其中23个可计算χ熵；无χ角位置单独标记。灰色为3条独立MD的均值，误差棒为两状态均值差的标准误，蓝色为完整复合物预测，橙色为移去TCR后的模型预测。所有位置可回溯至逐重复源表；残基或保存帧不作为独立MD重复，不报告将它们当作独立样本得到的显著性。

主要文件：`mel5_review_figure/mel5_complete_region_review.pdf`、同名PNG/SVG；`verified_candidate_comparison.csv`；`mel5_validation/plot_sources/`；`mel5_context/regional_context_comparison.csv`；`independent_MEL5_review_audit.json`；`expanded_search_counts.json`。全量筛查位于`construct_tolerant_screen/`和`existing_pair_rank/`，深度复算位于各候选的`*_validation/`目录。MEL5旧对齐曾跳过工程替换并导致CDR3统计虚高，已修正并独立核查；旧14位ρ=0.938结果作废，当前完整CDR3采用22位ρ=0.779。
'''
(ROOT/'TCR更多实例_搜索复核结论.md').write_text(report,encoding='utf8')

paper=ROOT.parent.parent/'paper_v40_20261017_中文修订_TCR识别环境变化实例_三线表.docx'
sha=hashlib.sha256(paper.read_bytes()).hexdigest()
assert sha=='bd9851a0f6dabdeb6f4a7c2fbbe9fe0dd20ec0820a156578db7f624cd990f5ca'
status={'status':'search_and_validation_complete','manuscript_modified':False,'existing_SB27_case_preserved':True,
        'manuscript_path':str(paper),'manuscript_SHA256':sha,'recommended_additional_case':'MEL5 / α24β17',
        'conclusion':'Stronger complementary same-pMHC side-chain and engineering-context evidence; not uniformly better at both outputs across all regions.',
        'report':'TCR更多实例_搜索复核结论.md','comparison':'verified_candidate_comparison.csv','counts':counts,
        'MEL5_context_test':'23 same-AA interface Hchi positions: rho 0.75494 full vs 0.05731 TCR removed; MAE 0.20241 vs 0.27067 nats.'}
(ROOT/'search_completion_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(status,ensure_ascii=False,indent=2))
print('Comparison rows:',len(metrics))
