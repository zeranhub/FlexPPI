from pathlib import Path
import json, hashlib, shutil
BASE=Path(r'E:\Project\FlexPPI\experiments\_sources')
OUT=BASE/'candidate_case_atlas_20261021/evidence'
p=OUT/'narratives.json'
d=json.loads(p.read_text(encoding='utf-8'))
mapping={'MEL5_AFFINITY':'C01','SB27_SB47':'C02','KRAS_G12D':'C03','JM22_F5L':'C04','PIK3CA_TCRS':'C05','B01':'C06','B02':'C07','B03':'C08','B04':'C09','B05':'C10','B06':'C11','B08':'C12','B09':'C13','B10':'C14','B11':'C15','B12':'C16','CATCHBOND':'C17'}
for c in d['cases']: c['catalog_case_id']=mapping.get(c['id'],c['id'])
d['reader_terms']={
'RMSF':'描述残基主链位置相对自身平均位置的活动幅度，数值较大表示该位置活动更大，单位Å。',
'Hχ':'把一个侧链各转角组合视作联合状态，描述观测到的取向组合有多丰富、多分散；较低不等于结合自由能更好。模型每个残基只给一个Hχ，未直接预测每条χ角曲线。',
'变化排序相关':'检查模型是否把变化较大与较小的位置排对顺序，并非每个位置绝对数值的准确率。',
'预测不变对照':'所有位置均认为两状态没有变化，用于判断模型提供的变化读数是否增加信息。',
'重复标准误':'用独立轨迹之间差别描述均值变化的不确定性；帧和残基不充当独立重复，超过两倍标准误只作描述，不称显著性。',
'MEL5残基编号':'β53/β54采用两状态对齐后的共同编号；WT输入E53/E54对应α24β17输入E54/E55，避免N端构建差异造成错位。',
'信息移除对照':'保持每个状态的pMHC坐标不变，移去TCR输入，仅检验模型是否使用伙伴结构信息，不表示真实的解离轨迹。'}
p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
check=[]
for c in d['cases']:
    assert len(c['paragraphs'])>=3
    for s in c['sources']:
        sp=Path(s['local']); assert sp.exists(),sp
        digest=hashlib.sha256(sp.read_bytes()).hexdigest()
        if s.get('sha256'): assert s['sha256']==digest,sp
        s['sha256']=digest
        check.append(dict(case=c['catalog_case_id'],path=str(sp),sha256=digest,verification=s.get('verification','原始全文或现稿内容已核')))
assert len({x['catalog_case_id'] for x in d['cases'] if x['category']=='TCR–pMHC'})==30
p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
for name in ['candidate_evidence_write.py','candidate_evidence_extend.py','candidate_evidence_validate.py']:
    shutil.copy2(Path(__file__).parent/name,OUT/name)
audit=dict(status='passed',date='2026-10-21',tcr_series=30,existing_non_tcr_cases=2,paragraphs=sum(len(c['paragraphs']) for c in d['cases']),sources=check,
           pending_primary_fulltext=['C18的两篇完整功能原文','C19','C20','C28无已核论文','C29无已核论文'],
           main_manuscript_modified=False,models_trained=0,model_inference_runs=0,MD_production_runs=0)
(OUT/'source_and_coverage_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','tcr_series','existing_non_tcr_cases','paragraphs']},ensure_ascii=False))
