import csv
import gzip
import json
from pathlib import Path
from collections import Counter

B = Path(r'E:\Project\FlexPPI\experiments\_sources')
O = B / 'candidate_case_atlas_20261021/catalog'
I = B / 'itrace_case_search_only_20261017'
T = B / 'tcr_case_v40_20261017'
X = B / 'tcr_expanded_screen_v41_20261017'
J = B / 'tcr_joint_case_search_20261019'

def rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

meta = {r['pdb']: r for r in rows(I/'all244_ITrace_metadata_and_sequences.csv')}
bio = rows(I/'biology_first/biology_first_ITrace_candidate_series.csv')
regional = rows(I/'numeric_rescreen/all_prespecified_region_model_metrics.csv')
reading = rows(I/'numeric_rescreen/biology_first_17_contrasts_reading_table.csv')
entries = []

def entry(cid, name, pdbs, contrasts, tier, background, supports, limits, sources, extras=None):
    pset = [p.lower() for p in pdbs.split()]
    pairs = [c.lower().split('>') for c in contrasts.split()]
    e = dict(case_id=cid, name=name, tier=tier, pdbs=pset, contrasts=pairs,
             background=background, model_can_support=supports, model_cannot_support_or_limitations=limits,
             local_primary_sources=[str(p) for p in sources], source_exists=[Path(p).exists() for p in sources],
             figures_to_make=['实验对照与真实结构定位', '预定完整区域的两状态模型绝对RMSF',
                              '同氨基酸位置的两状态模型Hχ', '完整区域变化与重复MD的辅助核对/模型信息控制'],
             raw_predictions={}, raw_structures={}, md_labels={}, complete_pair_tables=[],
             directly_recomputed_MD=False)
    for p in pset:
        cands = [I/'family_full_no_ITrace'/f'itrace_{p}.npz', T/'family_predictions/full'/f'itrace_{p}.npz']
        e['raw_predictions'][p] = [str(z) for z in cands if z.exists()]
        e['raw_structures'][p] = dict(model_input_pdb=[str(z) for z in B.rglob(f'{p}_model_input.pdb')],
            MD_topology=str(Path('E:/data/ITrace/trajectories')/f'{p}_run1_topology.pdb'),
            MD_topology_exists=(Path('E:/data/ITrace/trajectories')/f'{p}_run1_topology.pdb').exists(),
            relaxed_structure_archive='E:/data/ITrace/structures/relaxed_structures.tar.gz')
        labels = [z for z in I.glob(f'replica_recheck/{p}/*all_residue_per_replica_labels.csv')]
        labels += list(X.glob(f'*_validation/{p}_all_residue_per_replica_labels.csv'))
        e['md_labels'][p] = [str(z) for z in labels]
    for a,b in pairs:
        cands = [I/'promising_cases_replica_validation'/f'{a}_{b}_all_positions_three_replicas.csv',
                 I/'promising_cases_replica_validation'/f'{b}_{a}_all_positions_three_replicas.csv',
                 T/'change_screen'/f'pair_{a}_{b}_all_positions.csv', T/'change_screen'/f'pair_{b}_{a}_all_positions.csv',
                 X/'existing_pair_rank'/f'candidate_{a}_{b}_all_aligned_positions.csv',
                 X/'existing_pair_rank'/f'candidate_{b}_{a}_all_aligned_positions.csv']
        e['complete_pair_tables'] += [str(z) for z in cands if z.exists()]
    e['all_position_global_source'] = str(I/'numeric_rescreen/all_aligned_position_changes.csv.gz')
    e['expanded_screen_all_position_source'] = str(X/'construct_tolerant_screen/expanded_pair_all_positions.csv.gz')
    e['full_region_statistics'] = [r for r in regional if [r['a'],r['b']] in pairs or [r['b'],r['a']] in pairs]
    e['prespecified_reading_results'] = [r for r in reading if [r['a'],r['b']] in pairs]
    e['directly_recomputed_MD'] = bool(any(e['md_labels'].values()))
    if extras: e.update(extras)
    entries.append(e)
    return e

entry('C01','MEL5 → α24β17亲和成熟','3hg1 4jff 4jfd 4jfe','3hg1>4jff 4jff>4jfd 4jff>4jfe','A_深度复核主候选',
 '天然MEL5与亲和成熟α24β17识别同一MART-1 A27L十肽ELAGIGILTV–HLA-A2。原文SPR Kd为18 μM与0.6 nM，晶体显示新增接触主要位于MHC一侧。4JFD/4JFE是工程TCR的肽替换特异性控制，不能视作工程中间代。',
 '定位βCDR2/βCDR3主链约束增强，并在同序列pMHC界面读出侧链状态重排。Q72 Hχ 1.792→1.375 nats与未替换βGly54 RMSF 1.118→0.747 Å，结合Q72酰胺两端接触β53/54骨架的几何，可构成同一界面的双输出机制推断。',
 'αCDR2方向与MD不一致；全六环同AA Hχ相关弱。结合态Hχ不等于结合熵或亲和力；接触网络改变可共同约束两读数，不能把Hχ下降写成造成RMSF下降的力。',
 [X/'biology_sources/mel5_thermodynamics_2013.md',J/'literature/verified_primary_evidence.json'],
 {'complete_pair_tables':[str(X/'mel5_validation/corrected_pair_all_aligned_positions.csv')],
  'plot_sources':[str(X/'mel5_validation/plot_sources/mel5_all6_CDR_per_position_source.csv'),str(X/'mel5_validation/plot_sources/mel5_natural_pMHC_interface_source.csv')],
  'targeted_residue_evidence':str(J/'residue_story/trajectory_link/three_replica_summary.csv'),
  'context_control':str(X/'mel5_context/all_natural_pMHC_context_comparison.csv'),
  'numbering_warning':'WT TCRβ53 Ile→工程输入β54 Phe；论文共同编号β54 Gly→工程输入β55 Gly。必须沿修正对齐，旧CDR3 ρ=0.938作废。',
  'directly_recomputed_MD':True})
entry('C02','SB27 / SB47识别足迹差异','2ak4 4jry 4jrx','2ak4>4jry 2ak4>4jrx','A_深度复核主候选',
 'EBV BZLF1的13肽LPEPLPQGQLTAY由HLA-B*35:08呈递。SB27覆盖中央突起，SB47偏N端并绕过突起；CA5(4JRX)是同研究的第三受体控制。',
 '同肽Q7主链RMSF 0.850→1.449 Å、Hχ 1.726→2.825 nats，说明识别足迹变化伴随同一侧链和骨架的约束释放。去TCR输入后Q7 Hχ差异基本消失。',
 '全13肽RMSF变化ρ=0.945，但全11有效Hχρ=0.321，Q9不佳；不同天然受体并非单点干预，不能直接归因为一个残基。历史称谓“SB27&R7”未找到独立受体依据，不能增造R7案例。',
 [T/'biology_sources/sb27_sb47_primary_jbc_2013.md'],
 {'plot_sources':[str(T/'figure6/all13_model_and_MD_changes.csv')],
  'targeted_residue_evidence':str(T/'sb27_sb47/Q7_per_replica_physical_summary.csv'),
  'context_control':str(T/'sb27_sb47/context_inference.csv'), 'directly_recomputed_MD':True})
entry('C03','KRAS G12D单点新抗原','7ow5 7ow6 7pb2','7ow5>7ow6 7pb2>7ow6','A_深度复核主候选',
 '同一工程JDIa41b1与HLA-A11识别WT/G12D十肽，只改变P6。原文SPR、结构、MD及功能实验建立突变特异性，Kd约3 μM→743 pM。7PB2是亲本JDI结合G12D的另一工程对照。',
 '完整10肽及排除P6后的全部9个未替换位置均有主链响应。9位变化ρ=0.733、误差较零变化基线降低40.2%；可解释单点抗原替换伴随的邻位骨架约束与HLA介导识别。',
 'Hx仅5个未变肽位置，重复变化弱、旧模型未复现。不能以此声称双输出共同解释亲和力，也不能把P6自身化学类型改变当环境熵证据。',
 [T/'biology_sources/kras_g12d.md',I/'biology_first/sources/kras_g12d_2022_exa.txt'],
 {'plot_sources':[str(B/'kras_affinity_v48_20261019/figure4')], 'directly_recomputed_MD':True})
entry('C04','JM22识别流感M1 F5L','1oga 5hhm 5hho','1oga>5hhm','A_深度复核后备',
 'JM22识别HLA-A2–GILGFVFTL和P5 F→L变体；原文SPR WT 1.79±0.5 μM，F5L>200 μM，结合态晶体显示诱导契合。',
 '排除替换本位后，完整8个主链/6个Hχ位置变化排序ρ=0.714/0.657，可检验单点肽变体对未替换邻位的双读数响应。',
 '多数Hχ幅度明显偏小：F7 MD增0.495而模型增0.038 nats。游离到结合态诱导契合不能直接用两个结合态Hχ复述；相关高不足以证明功能机制。',
 [X/'biology_sources/flu_M1_crossreactivity_2016.md',X/'jm22_f5l_validation/复算判断.md'],
 {'complete_pair_tables':[str(X/'jm22_f5l_validation/all9_model_MD_physical_changes.csv')],
  'plot_sources':[str(X/'jm22_f5l_validation/common8_amino_acid_changes.csv')], 'directly_recomputed_MD':True})
entry('C05','PIK3CA H1047L的两受体识别','7l1d 7rrg','7l1d>7rrg','A_深度复核后备',
 '21LT2-2和0606T1-2识别同一PIK3CA H1047L新抗原ALHGGWTTK–HLA-A3，原研究有结构、免疫原性和治疗靶向背景。',
 '自然对应的同序列pMHC界面33个RMSF、26个Hχ位置可比较两类环境响应；侧链ρ约0.724，误差改善15.9%。',
 '自然pMHC主链误差仅改善5.4%，且重复波动大。全界面的好成绩很大部分来自跨不同TCR配对，不能把非同源CDR位置当严格同位，更不能声称整个界面双输出均准确。',
 [X/'biology_sources/pik3ca_public_neoantigen_2022.md',X/'pik3ca_validation/复算判断.md'],
 {'complete_pair_tables':[str(X/'pik3ca_validation/recomputed_pair_all_aligned_positions.csv')], 'directly_recomputed_MD':True})

bio_names = {'B01':'MAGE-A10 c728→c796','B02':'AHIII12.2识别三种MHC单点替换','B03':'1G4识别NY-ESO-1 C9V锚定变体','B04':'N17.1.2识别NRAS Q61K/Q61R','B05':'PRAME S2跨HLA-A3/A11','B06':'KRAS G12V 4TCR2→MH','B08':'PRAME S2→S2-198亲和工程','B09':'MEL5识别MART-1十肽/九肽','B10':'DMF4/DMF5识别MART-1十肽/九肽','B11':'MEL5识别三种肿瘤抗原','B12':'1G4 CDR2亲和成熟'}
bio_support = {
 'B01':'可检验受体单点变化对同序列抗原界面侧链环境的影响；主模型整肽7个Hχ位置ρ=0.643，误差改善15.7%。',
 'B02':'固定TCR/肽而改变MHC，是主链跨界面响应的直接问题；K66A全六环48位ρ=0.406、误差改善10.7%，需同时呈现W167A/T163A。',
 'B03':'可检验埋藏锚位替换是否经pMHC环境传到未变识别位点；目前完整肽RMSF未优于预测不变。',
 'B04':'同一受体识别两个新抗原，可检验侧链改变对相同骨架识别环境的影响；目前预定完整区域未提供有力模型支持。',
 'B05':'15个同AA界面Hχ位置两模型ρ=0.486/0.686，误差改善14.3%/24.1%，可支持呈递环境响应。',
 'B06':'真实亲本/工程对照可检验同一新抗原环境的识别环响应；当前全六环RMSF结果为反向，不能作为成功机制例。',
 'B08':'完整23个CDR3位置RMSF主模型ρ=0.559、误差改善37.2%，可提出工程后骨架变化的候选解释。',
 'B09':'NMR、晶体、热力学背景直接涉及肽锚定可塑性，适合检验模型是否识别环境变化；当前整环结果弱。',
 'B10':'DMF4自然pMHC界面19个Hχ位点ρ=0.779、误差改善25.2%，旧模型亦复现；可展示未改化学类型的界面侧链重排。',
 'B11':'IMP2→BST2全部22个CDR3主链位置两模型ρ=0.713/0.553、误差改善42.9%/41.8%；可支持特定抗原环境下的识别环响应。',
 'B12':'全六环49位主模型RMSFρ=0.680、误差改善59.7%；43个未替换位置ρ=0.710，可作为结构来源复核后的骨架响应候选。'}
bio_limits = {
 'B01':'旧模型ρ=0，7位仅1位变化超过描述性2SEM；不能称跨模型稳定机制。7QPJ c756并非c728衍生体。',
 'B02':'K66A旧模型未复现，W167A/T163A弱或失败。三种对照不能删去不利者；没有apo动态，Hχ不能替代ITC结合熵。',
 'B03':'完整肽RMSF相关−0.117、误差增加16.7%，无位点超过描述性2SEM；不能用P9化学类型差解释模型环境能力。',
 'B04':'全六环主链ρ=−0.365，误差增加38.4%。两种都被识别的新抗原不是WT免疫逃逸对照。',
 'B05':'Intrachain同样好，不能据此证明完整复合物信息提供增量；HLA等位基因有多处差异，不能仅归因163位。',
 'B06':'主模型全六环ρ=−0.507、误差增加60.7%；构建体/恒定区缺失也需关注。',
 'B08':'旧模型误差增加4.1%，同AA Hχ变化弱；不能以此解释亲和提升数量或更强A3偏好。S8的8RYP/8RYQ同时改HLA与TCR。',
 'B09':'肽长与编号改变，静态结构可能只含一种占据；主模型全六环ρ=−0.149、误差增加13.1%。6EQB是工程TCR结合态，不是apo。',
 'B10':'并行DMF5对照失败，必须一起保留；肽长度和锚身份同时改变，禁止按原编号生硬相减或泛化到所有受体。',
 'B11':'第三抗原Melan-A不能省：4QOK→7Q99失败，4QOK→7Q9A相关高但主模型误差更大；7Q9B为MEL8，抗原敲除杀伤证据主要不能移到MEL5。',
 'B12':'旧模型未复现；环外缺残基/构建差异可能影响按整链拟合的RMSF，需结构对应复核，不能包装成双模型一致。'}
cid = 6
for r in bio:
    k = r['case_id']
    if k=='B07': continue
    sources = []
    for v in r['source_local'].split(';'):
        p = (I/'biology_first'/v).resolve()
        if not p.exists():
            candidates = list(B.rglob(Path(v).name))
            if candidates: p=candidates[0]
        sources.append(p)
    pids=r['pdbs'].replace(';',' ')
    if k=='B01':pids+=' 7qpj'
    if k=='B08':pids+=' 8ryp 8ryq'
    if k=='B09':pids+=' 6eqb'
    if k=='B11':pids+=' 7q9b'
    if k=='B12':pids+=' 2p5e 2p5w 2pye'
    e=entry(f'C{cid:02}',bio_names[k],pids,r['contrasts'].replace(';',' '),
       'B_生物关系已核_含弱结果及负结果',r['biological_change']+'；'+r['controlled'],
       bio_support[k],bio_limits[k]+' '+r['caveat'],sources,
       {'original_biology_registry':r,'primary_doi':r['source_doi'],'primary_url':r['source_url']})
    cid+=1

entry('C17','TCR156 WT/S32M CatchBond','9nmu 9nmx','9nmu>9nmx','C_明确否决但应保留的负例',
 'TCR156识别PAP九肽TLMSAMTNL–HLA-A2。已有两状态各10条独立MD，Science2026原研究涉及力依赖识别与catch bond。',
 '可展示静态模型在该场景的适用边界；两种输入制备方式结果一致，不能将弱结果简单归因准备结构。',
 '407个同AA主链/343个Hχ位置变化ρ=0.022/−0.020，两项均不优于预测不变；静态平衡输出不能证明力依赖键寿命或信号，原数据水氧参数问题已披露。',
 [X/'other_data_inventory/结论与可执行路径.md'],
 {'complete_pair_tables':[str(X/'catchbond_validation/all408_MD_model_changes_both_inputs.csv')],
  'plot_sources':[str(X/'catchbond_validation/all408_ten_replica_mean_SD_SEM.csv')],
  'raw_predictions':{'9nmu':[str(X/'catchbond_validation/catchbond_WT_experimental_matched_prediction.npz'),str(X/'catchbond_validation/catchbond_WT_prepared_prediction.npz')],
                     '9nmx':[str(X/'catchbond_validation/catchbond_S32M_experimental_matched_prediction.npz'),str(X/'catchbond_validation/catchbond_S32M_prepared_prediction.npz')]},
  'raw_MD_directories':['E:/data/CatchBond/WT/WT','E:/data/CatchBond/S32M/S32M'],
  'directly_recomputed_MD':True})

extras = [
 ('C18','流感M1 JM22/LS01识别谱','2vlk 2vlr 5isz','2vlk>5isz 2vlr>5isz','少量同AA CDR3双读数初筛较好，9/5位RMSF/Hχρ=0.583/1.000。','只有5个熵点；不同TCR、不同环长的SequenceMatcher对应未经严格IMGT/结构复核，不能凭ρ=1写机制。','flu_M1_crossreactivity_2016.md'),
 ('C19','p53 R175H不同TCR识别','6vqo 7rm4 6vrm 6vrn','6vqo>7rm4 6vrm>6vrn','全界面25/19位双输出有弱正收益，可检验同一癌新抗原不同受体识别环境。','跨不同TCR少数同AA CDR对应待复核；不能把CDR3 8/5位高相关推广到整界面或单一因果接触。',None),
 ('C20','HCMV NLV的C7/C25','5d2l 5d2n','5d2l>5d2n','同一HCMV NLV肽不同TCR的界面主链变化ρ=0.927，侧链ρ=0.345，可作为识别方式差异候选。','初筛仅10个可配TCR界面位置；跨受体同位关系与原文功能需深核，不能称整体侧链精确。',None),
 ('C21','A6识别Tax、变体及交叉抗原','1ao7 1qrn 1qse 1qsf 3h9s 3pwp','1ao7>1qrn 1ao7>1qse 1ao7>1qsf 1ao7>3h9s 1ao7>3pwp','同一A6受体的点变体与Tel1p/HuD交叉识别可检验远端主链/侧链环境响应。','原先变体筛查未显示比保留实例更强；跨抗原多替换，不能用不同肽化学复杂度Hχ直接解释信号；Tel1p的跨来源无TCR辅助比较也较差。',None),
 ('C22','HLA-E呈递HBV变体','8rlt 8rlu 8rlv','8rlt>8rlu 8rlt>8rlv 8rlu>8rlv','同TCR识别HBV envelope371–379及S3N/L6I可检验非经典HLA环境的局部变化。','原先候选未达更强案例；需展示整肽和全部既定区域，不能选单一好位点；HLA-E特有结合模式不能类推所有HLA。',None),
 ('C23','gp100锚定变体的受体依赖效应','6vma 6vm9 6vmc 6vm7','6vma>6vm9 6vma>6vmc','T4H2对gp100 P2 T/M/L变体的识别，原研究展示结构近似但动态变构、受体依赖的识别后果。','原文生物问题非常贴近模型但当前未深复算；SILv44(6VM7)与T4H2不同，不得当成同受体工程前后。','gp100_dynamic_allostery_2021.md'),
 ('C24','HIV TCR55信号解耦与Pep20','6bj2 6bj3 6bj8','6bj3>6bj8 6bj2>6bj3','同TCR55的原肽IPLTEEAEL与Pep20 VPLTEDAEL，原研究显示结合可与信号解耦，力依赖行为不同。','同AA区域大多未获模型误差收益；不能把化学替换驱动的全肽Hχ收益当成功，也不能以平衡RMSF/Hχ证明catch bond或CD45排斥。','uncoupled_signaling_2018.md'),
 ('C25','MAGE-A3/Titin交叉反应','5brz 5bs0','5brz>5bs0','亲和增强TCR对肿瘤MAGE-A3与心肌Titin自体肽的分子模拟，是临床毒性相关的真实交叉识别背景。','多处肽替换、当前仅候选筛查未有完整机制复核；模型不能预测临床毒性或人体安全性，Hχ差异不能替代特异性实验。','mage_titin_2016.md'),
 ('C26','868识别HIV SL9逃逸系列','5nme 5nmf 5nmg','5nme>5nmf 5nme>5nmg','同868识别SL9、V6I和三突变；原文指向肽–HLA稳定性/呈递密度下降而非单纯TCR不结合。','两结合态局部动态不能解释细胞表面呈递量；三突变需剔除化学替换本位后检验，不能强行写成TCR亲和力逃逸。','hiv_sl9_escape_2017.md'),
 ('C27','1E6自身/病原多配体交叉识别','3uts 3utt 5c07 5c08 5c09 5c0a 5c0b 5c0c 5hyj','3uts>5c0a 3uts>5hyj','自身胰岛素前体反应1E6通过保守识别热点容忍大量抗原差异，可检验固定受体跨抗原的动态模式。','多处替换与不同晶型；当前只有初筛未有残基机制核实。5HYJ元数据DOI与原文主体不同，5C07 PDB序列与文献命名需逐一核查，不能直接照录。','autoimmune_crossreactivity_2016.md'),
 ('C28','公共/私人TCR及寡聚状态','9rup 9rxm','9rup>9rxm','同FTSDYYQLY抗原的TCRpub/pMHC dimer与TCRpriv/pMHC，初筛可变区154/134个同AA位置双输出有收益。','公开条目题名与寡聚体/构建差异尚未核实，原文当时To Be Published；只能作为待确认数据候选，不作功能机制结论。',None),
 ('C29','P14识别病毒/自体肽','5til 6g9q 5tje 5m00 5m01 5m02','5til>6g9q','P14对LCMV GP33变体与多巴胺单加氧酶自体肽的交叉识别；5TIL/6G9Q整肽双相关约0.767/0.750。','四处替换，主链误差仅改善6.9%；原条目缺清楚对应论文DOI，变体功能与序列映射待核，不宜按两相关包装机制。',None),
 ('C30','DMF5/α24β17识别同一MART-1','3qdg 4jff','3qdg>4jff','同一MART-1 A27L十肽由无亲本关系的两个TCR识别，初筛少数界面位点双输出改善。','3QDG是DMF5，4JFF是MEL5衍生α24β17，绝不能写成亲和成熟前后；11/6位对应与完整界面结果需复核。','dmf5_crossreactivity_2011.md'),
]
for cid,name,ps,cs,support,limit,src in extras:
    ss=[X/'existing_pair_rank/候选复核说明.md']
    if src:ss.append(X/'biology_sources'/src)
    if cid=='C21':ss.append(T/'biology_sources/a6_tax_tel1p.md')
    if cid=='C22':ss.append(T/'biology_sources/hbv_hlae.md')
    entry(cid,name,ps,cs,'D_已点名的发现性/文献候选_未完机制复核',support,support,limit,ss)

# Discovery universe: retain every pair without inventing biological stories.
old=rows(T/'local_pair_screen/biological_candidate_pairs.csv')
old_metrics=rows(T/'change_screen/all_candidate_change_metrics.csv')
new=rows(X/'construct_tolerant_screen/expanded_candidate_pairs.csv')
universe={}
for stage,rs in [('initial',old),('initial_metrics',old_metrics),('expanded',new)]:
    for r in rs:
        a=r.get('a') or r.get('pdb_a');b=r.get('b') or r.get('pdb_b')
        if not a or not b: continue
        key='|'.join(sorted([a.lower(),b.lower()]))
        z=universe.setdefault(key,dict(a=a.lower(),b=b.lower(),stages=[],metadata=r))
        if stage not in z['stages']:z['stages'].append(stage)
for e in entries:
    for a,b in e['contrasts']:
        key='|'.join(sorted([a,b]))
        if key in universe:universe[key].setdefault('catalog_case_ids',[]).append(e['case_id'])

for e in entries:
    e['PDB_metadata']=[{k:meta[p].get(k) for k in ['pdb','title','DOI','peptide_sequence','mhc_chain','peptide_chain','tcr_alpha_chain','tcr_beta_chain']} for p in e['pdbs'] if p in meta]
    e['prediction_coverage']=sum(bool(v) for v in e['raw_predictions'].values())
    if e['case_id']=='C17':e['prediction_coverage']=2
    e['metadata_coverage']=sum(p in meta for p in e['pdbs'])

result=dict(date='2026-10-21',scope='Prior explicitly nominated TCR-pMHC candidate series plus full discovery-pair index; no new training/inference/MD',
 total_case_series=len(entries),tiers=dict(Counter(e['tier'] for e in entries)),
 discovery_pair_count=len(universe),known_search_union_count=343,
 caution='30 series are not 30 successes; broad pair universe is discovery, not independent validation. R7 is unresolved historical naming and not a separate receptor case.',
 model=dict(name='existing joint Full family fold0',checkpoint_sha256='ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',training_ATLAS=1866,training_DynaRepo=564,training_ITrace=0),
 shared_sources=dict(all244_metadata=str(I/'all244_ITrace_metadata_and_sequences.csv'),biology_registry=str(I/'biology_first/biology_first_ITrace_candidate_series.csv'),all_position_source=str(I/'numeric_rescreen/all_aligned_position_changes.csv.gz'),expanded_all_position_source=str(X/'construct_tolerant_screen/expanded_pair_all_positions.csv.gz'),full_region_statistics=str(I/'numeric_rescreen/all_prespecified_region_model_metrics.csv')),
 aliases=[dict(historical='SB27 & R7',resolution='No separate R7 receptor or data found. Prior actual manuscript case is SB27/SB47, Q7 is the peptide glutamine. 4JRX=CA5, not R7.',case_id='C02')],
 auxiliary_cross_source_checks=[dict(name='DynaRepo SB27 3KXF(A005I) versus ITrace SB47 4JRY',rmsf_rho=0.8901,Hchi_rho=-0.0826,role='Different MD protocol corroboration only, not independent PDB case'),dict(name='DynaRepo 3H7B Tel1p without TCR versus ITrace3H9S A6',rmsf_rho=0.233,role='weak input-context auxiliary comparison'),dict(name='DynaRepo3HAE antibody-NYESO versus TCR2BNQ',role='no stable model advantage')],
 cases=entries)
(O/'candidate_inventory.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'all_discovery_pair_index.json').write_text(json.dumps(list(universe.values()),ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 既有TCR–pMHC潜在案例盘点','',f'本目录整理{len(entries)}个此前实际点名、推荐或否决的候选系列，含负结果；另保留{len(universe)}个结构配对的发现性索引。所有候选不等于预测成功。R7未发现独立受体证据，旧实际案例为SB27/SB47及肽Q7。','', '模型固定为未用ITrace训练的联合Full第0折。完整配对、预定区域和负结果均保留；RMSF与Hχ不能直接代替亲和力、结合自由能、力依赖键寿命或临床功能。','']
for e in entries:
    lines += [f"## {e['case_id']} {e['name']}",'',f"层级：{e['tier']}。结构：{' / '.join(p.upper() for p in e['pdbs'])}。",'',e['background'],'','模型能支持：'+e['model_can_support'],'','模型不能支持或当前缺陷：'+e['model_cannot_support_or_limitations'],'',f"预测NPZ覆盖：{e['prediction_coverage']}/{len(e['pdbs'])}；轨迹已深复算：{e['directly_recomputed_MD']}。",'']
    lines += ['原始文献/复核：'+ '；'.join(e['local_primary_sources']), '完整位置表：'+ '；'.join(e['complete_pair_tables'] or [e['all_position_global_source'],e['expanded_screen_all_position_source']]), '']
lines += ['## 全量索引与范围','', 'all_discovery_pair_index.json保留343个结构配对。未形成明确生物问题或未经原文核实的配对只应作为待确认数据索引，不可为每一对补造因果故事。30系列中扩展构建体、肽控制和第三状态列在同一案例，不重复计作独立成功。','']
(O/'candidate_inventory.md').write_text('\n'.join(lines),encoding='utf-8')
print(json.dumps({k:result[k] for k in ['total_case_series','tiers','discovery_pair_count']},ensure_ascii=False,indent=2))
