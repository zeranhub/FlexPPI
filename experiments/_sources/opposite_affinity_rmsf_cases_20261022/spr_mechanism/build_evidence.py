from pathlib import Path
import csv, json, hashlib, math

BASE=Path(r'E:\Project\FlexPPI\experiments\_sources')
OUT=BASE/'opposite_affinity_rmsf_cases_20261022'/'spr_mechanism'
OLD=BASE/'local_antigen_screening_20261021'
DIRECT=BASE/'local_peptide_max_rmsf_drop_20261022'/'direct_affinity'

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
labels=read(OLD/'experimental_labels'/'SPR_labels_observations_reported_hybrid.csv')
ranks=read(DIRECT/'primary_26_single_mutants_ranked.csv')
focus=[r for r in ranks if r['mutation'] in ['G4A','L3Y','L3A']]
cases=[]
for r in focus:
    observations=[]
    for o in labels:
        if o['candidate_id']!=r['candidate_id']:continue
        p=next(p for p in labels if p['observation_id']==o['parent_sample'])
        observations.append({
          'study_id':o['study_id'],'assay':'SPR','temperature_C':37,
          'KD_uM':float(o['KD_uM']),'KD_error':float(o['error_value']) if o['error_value'] else None,
          'error_type':o['error_type'],'n':int(o['n']) if o['n'] else None,
          'n_min':int(o['n_min']) if o['n_min'] else None,
          'parent_KD_uM':float(p['KD_uM']),'parent_KD_error':float(p['error_value']),
          'parent_error_type':p['error_type'],'parent_n':int(p['n']) if p['n'] else None,
          'KD_ratio_mutant_over_parent':float(o['KD_uM'])/float(p['KD_uM']),
          'affinity_ratio_mutant_over_parent':float(p['KD_uM'])/float(o['KD_uM']),
          'delta_pKD':float(o['delta_pKD_vs_parent']),'doi':o['primary_doi'],
          'source_url':o['primary_source_url'],'source_locator':o['primary_locator'],
          'same_individual_assay_batch':'not_reported','selected_fit':o['selected_fit'],
          'source_sha256':o['primary_source_sha256']})
    cases.append({'panel':r['panel_id'],'mutation':r['mutation'],'parent_sequence':r['parent_peptide'],
      'mutant_sequence':r['peptide_sequence'],'parent_crystal_template':r['template'],
      'max_drop_A':float(r['score_max_rmsf_drop_A']),'max_drop_positions':r['max_drop_positions'],
      'rmsf_all9_decrease':all(float(r[f'drop_RMSF_P{i}_A'])>0 for i in range(1,10)),
      'profile':[{'peptide_position':i,'parent_RMSF_A':float(r[f'parent_RMSF_P{i}_A']),
          'mutant_RMSF_A':float(r[f'mutant_RMSF_P{i}_A']),
          'drop_A':float(r[f'drop_RMSF_P{i}_A'])} for i in range(1,10)],
      'SPR_observations':observations,
      'input_scope':'all N/CA/C/O coordinates fixed to parent/template; amino-acid identity changed; no mutant conformational relaxation',
      'dynamic_prediction_experimental_validation':'No matched mutant MD/NMR/residue-resolved RMSF validation identified in the retrieved primary sources.'})

contacts=read(OUT/'parent_crystal_peptide_contacts.csv')
data={
 'analysis_id':'opposite_affinity_rmsf_SPR_mechanism_20261022',
 'knowledge_scope':'Exploratory mechanism audit of already-seen predictions and affinity labels; not blinded selection or new experimental validation.',
 'cases':cases,
 'all26_primary_single_mutants_positive_delta_pKD':[
     {'candidate':r['candidate_id'],'delta_pKD':float(r['delta_pKD_vs_parent'])}
     for r in ranks if float(r['delta_pKD_vs_parent'])>0],
 'A6_G4A_mechanistic_evidence':{
   'primary_paper':'Piepenbrink KH et al. Nature Communications 2013;4:1948',
   'doi':'10.1038/ncomms2948','url':'https://www.nature.com/articles/ncomms2948',
   'fulltext_url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC3708045/',
   'fulltext_xml_url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3708045/fullTextXML',
   'evidence_locators':['Results P13 / Fig3a','Discussion P28 / Fig7 caption','Methods Surface plasmon resonance'],
   'experiment':'SPR double-mutant cycles, 25 C; inter-residue energetic coupling, not direct G4A dynamics.',
   'measured_coupling_kcal_mol':{'T98alpha_R65HLA':-2.8,'D99alpha_R65HLA':[-3.4,-3.0],'D99Nalpha_R65HLA':-2.5},
   'source_structural_explanation':'Tax peptide Gly4 lacks a Cbeta, allowing the bound CDR3alpha conformation that forms strong T98alpha/D99alpha interactions with HLA Arg65. Adding a position-4 Cbeta is shown to create steric crowding against the Thr98alpha carbonyl in Fig7.',
   'important_distinction':'Fig7 added-Cbeta illustration is a structural rationale; no experimentally solved G4A-bound structure or G4A-specific MD/NMR validation is claimed.',
   'no_cross_temperature_energy_conversion':True,
   'checked_parent_1AO7_polar_heavy_atom_distances_A':[
      {'atom1':'D:THR98:OG1','atom2':'A:ARG65:NH2','distance':2.942307710647583},
      {'atom1':'D:ASP99:OD1','atom2':'A:ARG65:NE','distance':2.850537061691284},
      {'atom1':'D:ASP99:OD1','atom2':'A:ARG65:NH2','distance':3.1347994804382324}],
   'geometry_definition':'Atom-atom distances in original parent PDB 1AO7; N/O pairs <=3.5 A are polar-contact-compatible geometry only, without hydrogen or angle test. No hydrogen-bond occupancy, free energy or mutant loss is inferred from distances alone.',
   'verbatim_short_quote':'formation of the strong interactions with between CDR3α and Arg65 is dependent on the presence of glycine at peptide position 4',
   'quote_note':'Original PMC Fig7 caption contains duplicated words with between; preserve verbatim. 25 words by whitespace.'},
 '1G4_parent_structure_observations':{
   'actual_sequences':{'2BNQ':'SLLMWITQV','2BNR':'SLLMWITQC'},
   'peptide_P1_min_TCR_heavy_atom_distance_A':6.984382152557373,
   'peptide_P3_min_TCR_heavy_atom_distance_A':5.04121208190918,
   'P1_and_P3_contacts_within_4_5A':'Only HLA in parent 2BNQ; full contacts CSV retained.',
   'interpretation':'Largest P1 RMSF decreases for L3Y and L3A concern the N-terminal MHC-groove region; not direct evidence that TCR-facing contacts improve.',
   'L3Y_mutant_structure_or_dynamics':'Not identified in the retrieved primary sources; cannot claim a newly observed specific lost hydrogen bond, steric clash, or entropy mechanism.',
   'L3A_affinity_interpretation':'Approximately preserved affinity; about 1.11-fold point-estimate improvement, not a strong enhancement and no significance inferred.'},
 'dynamic_background_reference':{
   'citation':'Scott DR et al. J Mol Biol 2011;414:385-400','doi':'10.1016/j.jmb.2011.10.006',
   'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC3230710/',
   'methods':'Free A6 crystal structure; time-resolved fluorescence anisotropy; MD on free A6.',
   'applies_to':'A6 CDR3 loop conformational flexibility, not G4A mutant peptide RMSF.',
   'access':'Primary PMC text retrieved via agent-reach Exa search. EuropePMC fullTextXML endpoint returned JSON HTTP500 twice; archived errors, not accepted as XML.'},
 'reserve_9C9V_background_not_main_case':{
   'citation':'Chen JL et al. J Exp Med 2005;201:1243-1255','doi':'10.1084/jem.20042323',
   'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC2213140/',
   'SPR_25C':{'9C_KD_uM':13.3,'9V_KD_uM':5.0,'affinity_improvement_fold':2.66},
   'real_paired_crystals':{'9C':'2BNR','9V':'2BNQ'},
   'reason_not_main_opposite_group':'Only moderate 2.66-fold affinity improvement; native-template predicted P9 RMSF decline is only 0.0723 A, not the proposed large-drop positive arm. Different geometry inputs also differ from fixed-backbone mutant analysis.'},
 'manuscript_use_verdict':{
   'positive_model_validation':False,
   'suitable_role':'Exploratory counterexample / mechanistically informed negative control for interpreting predicted flexibility; can motivate context-aware residue-level interpretation.',
   'cannot_claim':['RMSF decrease experimentally verified for A6 G4A or 1G4 L3Y/L3A','Model predicted the G4A new bound conformation or loss of HLA Arg65 interactions','Reduced bound RMSF establishes a binding entropy benefit','A6/1G4 provide matched large-drop strong-affinity-up and strong-affinity-down groups'],
   'pairing_with_real_structure_MD_cases':'Do not pool fixed-parent virtual mutant predictions with real paired-structure/MD cases as equal evidence levels.'},
 'source_files':{}
}
for p in [DIRECT/'primary_26_single_mutants_ranked.csv',OLD/'experimental_labels'/'SPR_labels_observations_reported_hybrid.csv',OUT/'piepenbrink2013_fullText.xml',OUT/'chen2005_fullText.xml',OUT/'parent_crystal_peptide_contacts.csv']:
    data['source_files'][str(p)]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
(OUT/'evidence.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
md='''# 1G4/A6“局部变硬但亲和力不同”的原始证据审计

本次只解释既有结果，不修改论文，不重新训练、推理或运行MD。**现有26个单突变并不能组成“同样RMSF大降、一个强增亲和、一个强降亲和”的完整正反组**：1G4只有L3A的点估计比9V母肽略强，幅度约1.11倍；A6所有单突变都弱于母肽。1G4相关系数较高表示候选之间有排序信息，不表示其高分候选相对母肽大幅增强。

## 精确亲和力与模型读数

| 体系/变体 | 同研究母肽KD，μM | 变体KD，μM | 结合变化 | 最大预测RMSF下降 |
|---|---:|---:|---|---|
| A6 Tax G4A，Pettmann2021 | 2.01±0.437，SD，n=4 | 118±18.9，SD，n=3 | 弱58.71倍 | P4，0.942835→0.761021 Å，下降0.181814 Å |
| 1G4 L3Y，Pettmann2021 | 7.97±2.55，SD，n=10 | 47.3±8.58，SD，n=3 | 弱5.93倍 | P1，1.454331→1.003170 Å，下降0.451160 Å |
| 1G4 L3Y，Aleksic2010 | 7.2±0.5，SEM，n≥3 | 30±1，SEM，n≥3 | 弱4.17倍 | 同一份模型预测，不能当独立动力学验证 |
| 1G4 L3A，Pettmann2021 | 7.97±2.55，SD，n=10 | 7.04±1.49，SD，n=5 | 点估计强1.13倍，近于保持 | P1，1.454331→1.291882 Å，下降0.162449 Å |
| 1G4 L3A，Aleksic2010 | 7.2±0.5，SEM，n≥3 | 6.6±0.5，SEM，n≥3 | 点估计强1.09倍，近于保持 | 同一份模型预测 |

以上SPR均为37℃、同一研究与协议内比较的已发表均值，逐次实验是否同批配对未报告。误差类型不可混用；未用汇总误差推断显著性。Pettmann主拟合规则遵循原文：约束Bmax得到KD>20 μM时用约束值，否则用传统拟合；并非按预测效果选来源。其原始[1G4数据](https://cdn.elifesciences.org/articles/67092/elife-67092-fig1-data2-v3.csv)、[A6数据](https://cdn.elifesciences.org/articles/67092/elife-67092-fig1-data3-v3.csv)和[eLife正文](https://elifesciences.org/articles/67092)可复核。Aleksic来源为[原论文](https://doi.org/10.1016/j.immuni.2009.11.013)补充Table S1、PDF第16页上表。structured evidence.json保存逐行来源、原文件哈希和每个肽位点读数。

## A6 G4A：物理机制有依据，模型是否预测对动力学尚未验证

Tax是HTLV-1 Tax蛋白11–19位的九肽LLFGYPVYV；肽P4为原蛋白Gly14。**G4A的最大预测下降就在P4，不是远端；而且九个位点的预测RMSF全部下降。** 模型读数是将同一1AO7母本骨架上的P4氨基酸身份由Gly改为Ala之后的条件预测，未让骨架重新寻找突变体构象。

[Piepenbrink等2013年原始研究](https://www.nature.com/articles/ncomms2948)用25℃SPR双突变循环测量界面残基间的能量耦合：CDR3α Thr98/Asp99与HLA Arg65贡献很强。Fig7结合母本和游离A6结构提出，P4 Gly没有Cβ，使CDR3α能够采用形成这些接触的结合构象；P4换成Ala多出的Cβ造成拥挤，可能妨碍该构象。**这为“局部运动受限却结合变差”提供明确的结构理由，但不是G4A新构象已被直接测得。** 不能把2013年的25℃耦合能与2021年的37℃KD合成同一次热力学测量。

本次重新读取原始1AO7，链A为HLA、链C为肽、链D为TCRα：Thr98 OG1—Arg65 NH2为2.942 Å，Asp99 OD1—Arg65 NE为2.851 Å，Asp99 OD1—Arg65 NH2为3.135 Å，支持母本中确有相应极性接触的几何条件。这些是母本静态距离，不是突变体氢键丢失或占据率测量。未用简单距离直接计算相互作用能，也未声称G4A的9个位点真实动力学都下降。

[Scott等2011年原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC3230710/)确实做过游离A6晶体结构、时间分辨荧光各向异性和MD，说明CDR3α、CDR3β柔性和可达构象不同；该工作没有直接验证这里G4A肽逐残基RMSF。因而现有链条是：**G4A亲和减弱为直接SPR事实；关键结合几何与强耦合相互作用有原研究支持；模型预测G4A变硬是否真实，以及是否由该变化导致亲和减弱，尚无匹配动力学证据。**

## 1G4 L3Y/L3A：稳定的是HLA槽内区域，不能直接改写为TCR接触收益

1G4识别HLA-A*02:01呈递的NY-ESO-1 157–165九肽。此处母肽是优化过的9V序列SLLMWITQV，不能把它误写成天然9C。L3为原蛋白Leu159。L3Y与L3A最大的RMSF下降都在P1 Ser157；L3Y为0.4512 Å，L3A为0.1624 Å，突变P3本身的下降仅0.0196 Å和上升0.0301 Å。

读取真实母本2BNQ并采用任意重原子最短距离≤4.5 Å定义接触，P1、P3都只接触HLA，未直接接触TCR；其与TCR最近重原子距离分别为6.984 Å、5.041 Å。因此最大值主要突出了HLA槽内N端区域，而非TCR识别面上的直接增强。L3Y的大侧链改变可能调整HLA内堆积与相邻识别面的几何，属于值得核查的机制假设；**本次未找到L3Y/L3A匹配共晶或MD/NMR数据，不能声称发现了某条特定氢键断裂、侧链取向改变或结合熵改善。** L3A只能作为亲和力近于保持的对照，不宜写作高亲和力设计成功。

2BNQ=9V，2BNR=9C的真实结构序列与[Chen等2005原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2213140/)Accession numbers一致。该文有真实9C→9V优化案例，25℃SPR为13.3→5 μM，约2.66倍增强，可作为背景库存；但当前真实模板预测P9仅降0.0723 Å，且P7/P8升高，不符合本题“同样大幅下降”的正向组，暂不扩展为主案例。

## 对论文的用途判断

A6 G4A适合充当有原始结构机制支撑的解释对照：一个残基更少运动，并不回答它是否仍能形成正确的识别构象。1G4 L3Y/L3A说明最大下降值还会突出HLA槽内非直接TCR接触区域。两者能推进“柔性变化需要结合位置和界面几何解读”的论证，但**不能独立用来证明FlexPPI已准确解释了真实突变动力学，更不能组成强亲和上升/强下降的对称模型能力验证图**。如果主文已有真实成对结构及MD案例，应将这里固定母本的虚拟突变作为不同证据层级单列，不能混入同一组当同等验证。
'''
(OUT/'evidence.md').write_text(md,encoding='utf-8')
print('Wrote evidence.json and evidence.md')
print('Only primary positive',data['all26_primary_single_mutants_positive_delta_pKD'])
