from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
numeric=read(ROOT/'data_audit/candidate_numerical_audit.json')
assert numeric['status']=='passed' and all(r['passed'] for r in numeric['checks'])
assert (ROOT/'independent_review/case_argument_review.md').is_file()
assert (ROOT/'independent_review/case_argument_review.json').is_file()
report=ROOT/'TCR双输出案例_复核结论.md'
assert report.is_file()
pointer=read(ROOT.parent/'completed_results_20261012/delivery_status.json')
doc=Path(pointer['document']);assert sha(doc)==pointer['document_sha256']
manifest={'status':'completed_search_only','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),
 'recommended_case':'MEL5 parent / alpha24beta17, PDB 3HG1 / 4JFF',
 'joint_local_region':'TCR beta CDR2 backbone RMSF and neighboring HLA Q72 joint chi-state entropy; complete six-loop and pMHC-interface context retained',
 'alternative_same_peptide_case':'SB27 / SB47, PDB 2AK4 / 4JRY',
 'other_cases_rechecked':['JM22 WT / F5L, 1OGA / 5HHM','PIK3CA, 7L1D / 7RRG'],
 'new_model_training':0,'new_model_inference':0,'new_MD_production':0,'new_trajectory_label_computation':0,
 'document_modified':False,'delivery_pointer_modified':False,
 'current_document':str(doc),'current_document_sha256':sha(doc),
 'report':str(report),'report_sha256':sha(report),
 'checked_original_prediction_values':sum(r['checked_prediction_values'] for r in numeric['checks'] if 'checked_prediction_values' in r),
 'checked_existing_context_prediction_values':next(r['checked_raw_context_values'] for r in numeric['checks'] if 'checked_raw_context_values' in r),
 'source_hashes_verified':all(sha(s['path'])==s['sha256'] for s in numeric['sources']),
 'independent_argument_review':'independent_review/case_argument_review.md',
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='completion_status.json' and '__pycache__' not in p.parts}}
assert manifest['source_hashes_verified']
(ROOT/'completion_status.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:manifest[k] for k in ['status','recommended_case','report','document_modified','checked_original_prediction_values','checked_existing_context_prediction_values','source_hashes_verified']},ensure_ascii=False))
