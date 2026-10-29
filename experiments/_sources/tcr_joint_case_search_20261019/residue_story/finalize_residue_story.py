from pathlib import Path
import hashlib,json
from datetime import datetime,timezone,timedelta
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
audit=read(ROOT/'trajectory_link/analysis_audit.json')
validation=read(ROOT/'trajectory_link/validation_and_inventory.json')
review=read(ROOT/'independent_review/final_evidence_chain_review.json')
assert audit['status']=='passed' and validation['status']=='passed'
assert review['status']=='passed_with_wording_constraints'
wording_updates=[]
for name,digest in review['source_SHA256'].items():
    current=sha(ROOT/'trajectory_link'/name)
    if current!=digest:
        assert name=='summarize_Q72_beta_link.py',name
        wording_updates.append({'file':name,'reviewed_sha256':digest,'current_sha256':current,'reason':'Report generator wording amended after review: higher contact occupancy rather than measured longer lifetime; all reviewed numerical data remain byte-identical.'})
for item in validation['files']:assert sha(item['path'])==item['sha256'],item['path']
pointer=read(ROOT.parent.parent/'completed_results_20261012/delivery_status.json')
assert sha(pointer['document'])==pointer['document_sha256']
result={'status':'completed_residue_level_case_analysis',
 'completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),
 'case':'MEL5 / alpha24beta17; HLA Q72 with TCR beta Phe53 / unchanged Gly54',
 'model_Q72_Hchi_nats':[1.792447805404663,1.3750499486923218],
 'model_betaGly54_RMSF_A':[1.1178646087646484,0.7468647956848145],
 'existing_MD_trajectories_reanalyzed':6,'saved_frames':6012,
 'new_MD_production':0,'new_training':0,'new_inference':0,
 'interpretation':'A remodeled contact patch provides a concrete local structural explanation for co-occurring side-chain state concentration and backbone-motion restriction; Hchi does not independently cause RMSF.',
 'wording':'Use higher contact occupancy, not measured longer contact lifetime. Three-replica Hchi means differ with overlapping ranges. Slower dissociation belongs to published SPR, not to the model or these bound-state trajectories.',
 'post_review_report_wording_updates':wording_updates,
 'document_modified':False,'current_document':pointer['document'],'document_sha256':pointer['document_sha256'],
 'prose_suggestion':str(ROOT/'Q72_Gly54具体叙述_建议稿.md'),
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='completion_status.json' and '__pycache__' not in p.parts}}
(ROOT/'completion_status.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:result[k] for k in ['status','case','existing_MD_trajectories_reanalyzed','saved_frames','document_modified','prose_suggestion']},ensure_ascii=False))
