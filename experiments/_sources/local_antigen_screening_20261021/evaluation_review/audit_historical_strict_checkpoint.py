"""Read-only audit of a locally available historical checkpoint and its manifest."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import os,json,hashlib
os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import torch
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
REF=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
CKPT=REF/'trained_checkpoints/strict_full/fold_0/epoch_040.pt'
STRICT=REF/'strict_external_payload/splits/strict_external_train.jsonl'
TRAIN2430=EXP/'family_folds/splits/fold_0_train.jsonl'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding='utf8'))
expected='b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f'
assert sha(CKPT)==expected
checkpoint=torch.load(str(CKPT),map_location='cpu',weights_only=False)
strict=[json.loads(l) for l in STRICT.read_text().splitlines() if l.strip()]
current=[json.loads(l) for l in TRAIN2430.read_text().splitlines() if l.strip()]
bystrict={r['sample_id']:r for r in strict};bycurrent={r['sample_id']:r for r in current}
assert len(bystrict)==2431 and len(bycurrent)==2430
onlystrict=sorted(set(bystrict)-set(bycurrent));onlycurrent=sorted(set(bycurrent)-set(bystrict))
assert onlystrict==['A00IL'] and onlycurrent==[]
common_diffs=[]
for sid in sorted(bycurrent):
    a,b=bystrict[sid],bycurrent[sid]
    core_a={k:v for k,v in a.items() if k!='metadata'}
    core_b={k:v for k,v in b.items() if k!='metadata'}
    assert core_a==core_b,(sid,'non-metadata changed')
    da,db=a['metadata'],b['metadata']
    # Split bookkeeping may change, but labels, input identity and provenance must not.
    for field in ['pdb_id','pdb_ids','label_sha256','original_record_sha256','original_label_sha256','training_metrics','residue_count']:
        assert da.get(field)==db.get(field),(sid,field)
    diffkeys=sorted(k for k in set(da)|set(db) if da.get(k)!=db.get(k))
    if diffkeys:common_diffs.append(dict(sample_id=sid,metadata_fields_only=diffkeys))
cfg=checkpoint['config']
assert checkpoint['epoch']==40 and cfg['seed']==7
assert cfg['data']['state_manifest']=='splits/strict_external_train.jsonl'
assert cfg['model']['information_ablation']=='full'
assert cfg['model']['metrics']==['rmsf','chi_entropy']
assert cfg['model']['sequence_backend']=='learned'
assert cfg['model']['enable_response_head']==False
assert cfg['model']['use_condition_features']==False
sources=dict(Counter(r['dataset'] for r in strict))
assert sources=={'atlas':1866,'dynarepo':565}
assert not any(r['dataset']=='itrace' for r in strict)
reported=read(EXP/'itrace_case_search_only_20261017/model_audit/remote_checkpoint_manifest_audit.json')
history=[c for c in reported['checkpoints'] if c['sha256']==expected]
assert len(history)==1
assert history[0]['config']==cfg
assert history[0]['manifests']['state_manifest']['sha256']==sha(STRICT)
inventory=read(REF/'sequence_audit/chain_inventory.json')
added=[c for c in inventory if c['sample_id']=='A00IL']
assert len(added)>=1
trainids=set(bystrict)
strictchains=[c for c in inventory if c['sample_id'] in trainids]
oldaudit=read(ROOT/'training_overlap_audit.json')
queryseqs={q['sequence'] for q in oldaudit['query_chain_audit']}
assert not any(c['sequence'] in queryseqs for c in added)
candidates=read(ROOT/'candidate_sequence_overlap.json')
assert not any(c['sequence']==r['sequence'] for c in added for r in candidates)

# Only the added record needs a new similarity calculation; all remaining
# chains were independently audited against the full 2430 base manifest.
from Bio.Align import PairwiseAligner
aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-1,open_gap_score=-3,extend_gap_score=-.5)
new_hits=[]
rawaudit=read(ROOT/'raw_template_overlap_audit.json')
for q in rawaudit['rows']:
    seq=q['sequence']
    for c in added:
        target=c['sequence']
        if seq==target:new_hits.append(dict(query=q,added_record=c,exact=True));continue
        if len(seq)<20 or min(len(seq),len(target))/max(len(seq),len(target))<.7:continue
        a=aligner.align(seq,target)[0];qb,tb=a.aligned
        if len(qb)==0:continue
        qc=int(qb[-1,1]-qb[0,0])/len(seq);tc=int(tb[-1,1]-tb[0,0])/len(target)
        if min(qc,tc)<.7:continue
        matches=sum(sum(aa==bb for aa,bb in zip(seq[int(qs):int(qe)],target[int(ts):int(te)])) for (qs,qe),(ts,te) in zip(qb,tb))
        x,y=a.coordinates
        columns=sum(max(int(x[i+1]-x[i]),int(y[i+1]-y[i])) for i in range(len(x)-1))
        identity=matches/columns
        if identity>=.3:new_hits.append(dict(pdb=q['pdb'],chain=q['chain'],scope=q['scope'],added_sample=c['sample_id'],added_chain=c['chain'],identity=identity,query_coverage=qc,target_coverage=tc))
assert new_hits==[],new_hits
out=dict(status='passed_read_only_historical_checkpoint_audit',audited_at_utc=datetime.now(timezone.utc).isoformat(),
    actual_checkpoint=dict(path=str(CKPT),sha256=sha(CKPT),bytes=CKPT.stat().st_size,epoch=checkpoint['epoch'],config=cfg),
    strict_train_manifest=dict(path=str(STRICT),sha256=sha(STRICT),records=len(strict),sources=sources,itrace_records=0,chains=len(strictchains),unique_chain_sequences=len({c['sequence'] for c in strictchains})),
    comparison_2430=dict(path=str(TRAIN2430),sha256=sha(TRAIN2430),records=len(current),common_records=2430,only_historical=onlystrict,only_2430=onlycurrent,
        common_input_label_protocol_and_paths_identical=True,common_label_hashes_identical=True,metadata_diff_keys=dict(Counter(k for d in common_diffs for k in d['metadata_fields_only'])),
        note='Family split bookkeeping differs; this does not make the two independently trained weight files identical or random replicates.'),
    added_record=bystrict['A00IL'],added_chain_inventory=added,
    added_record_exact_query_overlap=False,added_record_exact_candidate_peptide_overlap=False,added_record_raw_template_and_domain_30pct_70pct_both_hits=new_hits,
    inherited_2430_full_audit_sha256=sha(ROOT/'training_overlap_audit.json'),raw_template_audit_sha256=sha(ROOT/'raw_template_overlap_audit.json'),
    same_weights_as_current_manuscript_2430=False,
    justified_scope='Local pilot using the available 2431 historical strict no-ITrace joint model; an implementation and retrospective ranking test, not the exact current-manuscript checkpoint.',
    model_inference_run=False,training_run=False,remote_access=False)
(ROOT/'historical_strict_2431_checkpoint_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:out[k] for k in ['status','strict_train_manifest','comparison_2430','added_chain_inventory','same_weights_as_current_manuscript_2430']},ensure_ascii=False))
