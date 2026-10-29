"""Sequence-only overlap audit; reads no new model scores or response labels."""
from pathlib import Path
from collections import defaultdict
import csv,hashlib,importlib.util,json,sys,time
import Bio
from Bio.Align import PairwiseAligner,substitution_matrices
R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
B=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028\runtime')
F=Path(r'E:\Project\FlexPPI\experiments\_sources\family_folds')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
module_path=B/'family_input_audit.py'
spec=importlib.util.spec_from_file_location('alignment_functions',module_path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
inventory=json.loads((B/'family_input_audit_reference_chain_inventory.json').read_text())
assert len(inventory)==4973
family=list(csv.DictReader((F/'record_family_folds.csv').open()))
fold={(r['dataset'],r['sample_id']):int(r['fold']) for r in family}
refs=defaultdict(list)
for r in inventory:
 assert fold[r['dataset'],r['sample_id']]==r['fold']
 assert r['in_fold0_train']==(r['fold']!=0)
 refs[r['sequence']].append({k:r[k] for k in ['dataset','sample_id','chain','fold','component_id','in_fold0_train']})
prep=json.loads((R/'audits/preparation_audit.json').read_text())
states=json.loads((R/'states.json').read_text());queries=defaultdict(list)
for s in states:
 assert sha(s['path'])==s['sha256']
 for ch,seq0 in prep['raw_templates'][s['template']]['chain_sequences'].items():
  seq=s['sequence'] if ch=='C' else seq0
  queries[seq].append({'state_id':s['id'],'chain':ch})
matrix=substitution_matrices.load('BLOSUM62');aligners={}
for mode in ['local','global']:
 a=PairwiseAligner(mode=mode,open_gap_score=-11,extend_gap_score=-1);a.substitution_matrix=matrix;aligners[mode]=a
summaries=[];start=time.monotonic()
for i,(q,uses) in enumerate(sorted(queries.items(),key=lambda x:(-len(x[0]),x[0]))):
 hits=[]
 for t,rr in refs.items():
  if q==t:
   best={'method':'exact','identity':1.,'query_coverage':1.,'target_coverage':1.}
  elif len(q)<20 or len(t)<20:continue
  else:
   passing=[]
   for mode,a in aligners.items():
    z=module.align_stats(q,t,a,mode)
    if module.passes(z,q,t):passing.append(z)
   if not passing:continue
   best=max(passing,key=lambda x:x['identity']*min(x['query_coverage'],x['target_coverage']))
  hits.append({'target_sha256':hashlib.sha256(t.encode()).hexdigest(),'alignment':best,'references':rr})
 train_hits=[h for h in hits if any(r['in_fold0_train'] for r in h['references'])]
 summaries.append({'sequence':q,'length':len(q),'uses':uses,'passing_hits':hits,'potential_train_bridge':bool(train_hits),'train_hits':train_hits})
 if len(q)>=20:print('long chain',i+1,'length',len(q),'training_hits',len(train_hits),'seconds',round(time.monotonic()-start,1),flush=True)
report={'status':'complete','prediction_blind':True,'states':len(states),'unique_sequences':len(queries),'long_queries':sum(len(q)>=20 for q in queries),
 'all_training_pool_records':len(family),'all_training_pool_chains':len(inventory),'fold0_training_records':sum(r['fold']!='0' for r in family),
 'method':'Exact sequences at any length; exhaustive local/global BLOSUM62 alignment union for chains >=20 residues, >=30% identity and >=70% integer-span bidirectional coverage. Same transparent fallback as previous audit, not an MMseqs rerun.',
 'limitation':'Short novel peptide homology is not certified by whole-chain family rule. This fallback may conservatively flag bridges; absence does not prove identical original MMseqs result.',
 'queries_with_training_hits':[{'length':x['length'],'sequence':x['sequence'],'chains':sorted({u['chain'] for u in x['uses']}),'train_hits':x['train_hits']} for x in summaries if x['potential_train_bridge']],
 'queries':summaries,'biopython':Bio.__version__,'source_sha256':{str(p):sha(p) for p in [R/'states.json',R/'audits/preparation_audit.json',B/'family_input_audit_reference_chain_inventory.json',F/'record_family_folds.csv',module_path]},'seconds':time.monotonic()-start}
(R/'audits/current_model_training_overlap.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('COMPLETE','queries',len(queries),'train_bridge_queries',len(report['queries_with_training_hits']),flush=True)
