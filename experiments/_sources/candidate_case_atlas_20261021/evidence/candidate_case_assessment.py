from pathlib import Path
import json, hashlib, shutil
BASE=Path(r'E:\Project\FlexPPI\experiments\_sources')
OUT=BASE/'candidate_case_atlas_20261021/evidence'
p=OUT/'narratives.json'
d=json.loads(p.read_text(encoding='utf-8'))
rows=[
['C01','MEL5 → α24β17','残基机制复核','β环与对侧侧链的共同约束','多处替换；未分离单接触因果'],
['C02','SB27 / SB47','残基机制复核','同肽中央骨架与Q7取向差异','整肽侧链较弱；受体不同'],
['C03','KRAS WT / G12D','完整区域复核','突变邻近未替换肽骨架响应','侧链增量弱；不能预测亲和力'],
['C04','JM22原肽 / F5L','完整区域复核','未替换肽位点的部分双读数响应','多数侧链变化幅度不足'],
['C05','PIK3CA两种TCR','完整区域复核','相同pMHC侧链环境差异','自然界面主链改善较小'],
['C06','MAGE-A10 c728 → c796','完整区域筛查','整肽侧链存在单模型线索','较大变化位点少；模型间不一'],
['C07','AHIII12.2四种HLA','完整区域筛查','K66A有弱主链响应','其余对照弱；整体未复现'],
['C08','1G4原肽 / C9V','完整区域筛查','当前未形成有效变化证据','未优于预测不变；MD变化小'],
['C09','NRAS Q61K / Q61R','完整区域筛查','当前未形成有效环变化证据','完整识别环结果较弱'],
['C10','PRAME S2：A3 / A11','完整区域筛查','同类型界面侧链的部分响应','未显示链间信息额外优势'],
['C11','KRAS G12V：4TCR2 → MH','完整区域筛查','当前未形成有效环变化证据','构建体差异；主链变化小'],
['C12','PRAME S2 → S2-198','完整区域筛查','CDR3有单模型主链线索','模型间不一；侧链变化弱'],
['C13','MEL5：MART-1十肽 / 九肽','完整区域筛查','当前未恢复长短肽环响应','长度及锚定同时改变'],
['C14','DMF4 / DMF5长短肽系列','完整区域筛查','DMF4界面侧链变化','DMF5平行对照未成功'],
['C15','MEL5三种肿瘤抗原','完整区域筛查','IMP2 / BST2的CDR3主链差异','其余抗原对照较弱'],
['C16','1G4 → c49c50','待构建体复核','未替换环的单模型主链线索','构建和拟合影响；模型间不一'],
['C17','TCR156 WT / S32M','完整区域复核','当前未恢复整体动态差异','不计算力依赖结合寿命'],
['C18','流感M1：JM22 / LS01','发现阶段','少量配对位点有双读数线索','仅5个熵点；跨环对应待核'],
['C19','p53 R175H不同TCR','发现阶段','界面变化有弱正向线索','局部小样本强；生物关系待深核'],
['C20','HCMV NLV：C7 / C25','发现阶段','少量界面主链变化线索','仅10个配对位置；侧链较弱'],
['C21','A6：Tax变体 / Tel1p / HuD','发现阶段','可检验配体依赖的环境响应','当前预测弱；多变体对应待核'],
['C22','HBV–HLA-E三种肽','发现阶段','可检验非经典HLA环境响应','尚未形成完整区域机制证据'],
['C23','gp100锚定变体 / 两TCR','发现阶段','可检验锚点远端动态响应','未深复核；两受体非工程前后'],
['C24','HIV TCR55原肽 / Pep20','发现阶段','当前同类型区域证据不足','不包含受力与膜信号机制'],
['C25','MAGE-A3 / Titin','发现阶段','可检验交叉识别的环境响应','多处替换；不能预测临床毒性'],
['C26','HIV 868：SL9逃逸系列','发现阶段','可检验未替换界面动态差异','不能解释呈递量；未深复核'],
['C27','1E6自身 / 病原多配体','发现阶段','可检验识别热点周边约束','多处替换；配体与编号待核'],
['C28','公共 / 私人TCR','待生物关系核实','仅有发现性数值线索','寡聚和构建混杂；缺已核原文'],
['C29','P14病毒 / 自身肽','待生物关系核实','整肽变化有初筛排序线索','幅度收益小；缺已核原文'],
['C30','DMF5 / α24β17','发现阶段','少量界面双读数线索','非亲本关系；完整位点待核'],
['C31','KIX四种伙伴组成','专门MD对照','共同高活动环定位','状态差值方向未恢复；单轨迹'],
['C32','泛素–EDD UBA','留出结构实例','同类型侧链多样性区分','不直接给逐角分布或结合机制'],
]
assert len(rows)==32 and len({r[0] for r in rows})==32
for c in d['cases']:
    if c['id']=='KIX':c['catalog_case_id']='C31'
    if c['id']=='UBIQUITIN':c['catalog_case_id']='C32'
    c['paragraph_roles']=['生物背景与实验对照','模型实际证据']+(['机制联系与额外证据'] if len(c['paragraphs'])==4 else [])+['具体限制与未解决问题']
    assert len(c['paragraphs'])==len(c['paragraph_roles'])
    for i,t in enumerate(c['paragraphs']):
        t=t.replace('另一既有Full模型','另一训练对照模型').replace('另一既有Full模型','另一训练对照模型')
        c['paragraphs'][i]=t
p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
assessment=dict(date='2026-10-21',headers=['编号','对照','证据阶段','首要模型能力','首要缺口'],rows=rows,
               row_count=32,tcr_series=30,existing_non_tcr_cases=2,
               note='“可检验”表示生物问题与可用数据，不表示模型已证明机制；发现阶段的小区域结果需完整位点复核。负结果按实际证据写明，不自动归为成功能力。',
               source_narratives_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
               formatting={'orientation':'portrait','suggested_column_widths_cm':[1.0,3.8,2.7,4.5,4.7],'native_three_line_table':True,'repeat_header':True,'no_vertical_borders':True,'no_internal_horizontal_rules':True,'allow_page_split_between_rows_only':True,'font_size_pt':9,'no_cell_shading':True})
(OUT/'case_assessment_table.json').write_text(json.dumps(assessment,ensure_ascii=False,indent=2),encoding='utf-8')
reader='RMSF描述一个残基主链在一段时间里离开自身平均位置的幅度，单位为Å，数值较大表示这一位置活动较大。Hχ则描述同一侧链多个转角组合的取向有多丰富、多分散：较低意味着观测到的取向更集中，较高意味着有更多取向选择。两项一起看，可以区分“骨架较稳而侧链仍有多种朝向”和“骨架与侧链都受较强限制”等局部状态，但它们不会自动告诉我们哪条接触造成变化，更不是亲和力或结合自由能。模型每个残基只输出一个联合Hχ，图中的逐χ角度曲线、具体构象占据或接触比例若来自轨迹，需分别标明。文中变化排序相关只看模型能否把变化较大、较小的位置排对顺序，因此还需看变化幅度及相对“所有位置均不变”是否提供了信息；较高相关，尤其仅来自五六个配对位置时，不能单独作为机制成功的结论。'
(OUT/'reader_terms_paragraph.md').write_text(reader+'\n',encoding='utf-8')
audit=dict(status='passed_with_explicit_scope',date='2026-10-21',case_count=32,narrative_sha256=assessment['source_narratives_sha256'],
 checks=[
 '主TCR结果使用既有未含ITrace训练的Full第0折；KIX与泛素分别沿用原稿自身留出来源，不把32例宣称为同一测试队列。',
 'MEL5完整六环和同序列pMHC侧链分别讨论；未以TCR替换位置Hχ证明环境敏感性。',
 'SB27/SB47、PIK3CA、C18–20和C30均明示不同TCR；未写成单点工程因果。',
 'JM22与KRAS均将替换本位的化学类型影响与未替换环境响应分开。',
 'C18、C19、C20和C30的小n高相关配合位置数和范围限制；C29还保留幅度收益小。',
 'MAGE、AHIII、NYESO、NRAS、KRASG12V、PRAME、MART1长短肽和1G4成熟均明确弱结果、失败或训练模型间不一致。',
 'DMF4/DMF5与MEL5三抗原保留不利平行对照，不把单个较好对照泛化为整个系列。',
 'CatchBond未提供成功模型机制；RMSF/Hχ未当作力依赖寿命。',
 'Q72接触机制写作接触网络共同约束两读数，未写Hχ降低单向导致Gly54 RMSF降低。',
 '原文SPR、结构、原文MD、ITrace MD及本模型预测各自归属清楚；未用已结合态Hχ替代结合熵。',
 'C28/C29缺已核主源，明确生物关系/构建体仍待证，不补造功能数字。',
 'KIX共同活动区定位和组成状态差值分开，保留后者失败及每状态单轨迹条件。',
 '泛素为同结构位置的同类型差异案例，不推导结合前后功能。'],
 outstanding=['待根任务实际大图导出后，将图中显示的所有区域与本段结论逐例对齐。','待最终DOCX抽取及逐页视觉检查。'])
(OUT/'narrative_content_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 候选案例背景、模型证据与适用范围','',d['scope'],'',d['historical_alias_note'],'']
for c in sorted(d['cases'],key=lambda c:c['catalog_case_id']):
    lines += [f"## {c['catalog_case_id']}　{c['title']}",'',f"结构：{' / '.join(c['pdbs'])}；状态：{c['status']}。",'']
    for text in c['paragraphs']: lines +=[text,'']
    for s in c['sources']: lines +=[f"来源：[{s['title']}]({s['url']})。DOI：{s['doi'] or '当前未核实'}。{s.get('verification','原始内容或现稿已核')}；本地：`{s['local']}`。",'']
(OUT/'narratives.md').write_text('\n'.join(lines),encoding='utf-8')
shutil.copy2(__file__,OUT/Path(__file__).name)
print(json.dumps(dict(status='passed',rows=len(rows),max_cell_characters=max(len(s) for r in rows for s in r),out=str(OUT)),ensure_ascii=False))
