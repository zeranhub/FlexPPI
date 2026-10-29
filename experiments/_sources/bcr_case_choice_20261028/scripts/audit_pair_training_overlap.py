from pathlib import Path
from collections import defaultdict
import hashlib,importlib.util,json
from Bio.Align import PairwiseAligner,substitution_matrices
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
B=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028\runtime')
spec=importlib.util.spec_from_file_location('alignment_functions',B/'family_input_audit.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
refs=defaultdict(list)
for r in json.loads((B/'family_input_audit_reference_chain_inventory.json').read_text()):refs[r['sequence']].append({k:r[k] for k in ['sample_id','dataset','chain','fold','component_id','in_fold0_train']})
queries=defaultdict(list)
for s in json.loads((R/'pair_screen_states.json').read_text()):
 for ch,r in s['chain_records'].items():queries[r['sequence']].append({'state_id':s['id'],'chain':ch})
aligners={}
for mode in ['local','global']:
 a=PairwiseAligner(mode=mode,open_gap_score=-11,extend_gap_score=-1);a.substitution_matrix=substitution_matrices.load('BLOSUM62');aligners[mode]=a
out=[]
for i,(q,uses) in enumerate(queries.items()):
 hit=[]
 for t,rr in refs.items():
  if not any(r['in_fold0_train'] for r in rr):continue
  if q==t:stats={'method':'exact','identity':1.,'query_coverage':1.,'target_coverage':1.}
  else:
   passing=[z for mode,a in aligners.items() if m.passes(z:=m.align_stats(q,t,a,mode),q,t)]
   if not passing:continue
   stats=max(passing,key=lambda z:z['identity']*min(z['query_coverage'],z['target_coverage']))
  hit.append({'alignment':stats,'references':[r for r in rr if r['in_fold0_train']]})
 out.append({'sequence':q,'length':len(q),'uses':uses,'potential_train_hits':hit})
 print(i+1,len(queries),'length',len(q),'train_hits',len(hit),flush=True)
(R/'audits/paired_training_overlap.json').write_text(json.dumps({'status':'complete','unique_queries':len(queries),'queries_with_potential_train_hits':sum(bool(z['potential_train_hits']) for z in out),
 'method':'Exhaustive local/global BLOSUM62 union, >=30% identity and >=70% bidirectional coverage; not original MMseqs rerun.','independent_claim_limit':'Threshold coverage may change when cropping a domain. Do not certify strict novel-query independence from absence of these fallback hits alone.',
 'queries':out,'states_sha256':hashlib.sha256((R/'pair_screen_states.json').read_bytes()).hexdigest()},indent=2),encoding='utf-8')
