"""Additional numerical checks using an independent weighted-logistic objective."""
from pathlib import Path
from collections import OrderedDict
import csv,hashlib,json,math,sys
import numpy as np,pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
states=json.loads((R/'states.json').read_text());assert len(states)==564
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
def pdb_identity_and_backbone(path):
 ca=OrderedDict();bb=[]
 for l in Path(path).read_text().splitlines():
  if not l.startswith('ATOM  '):continue
  if l[12:16].strip()=='CA':ca[(l[21],l[22:27].strip())]=AA[l[17:20]]
  if l[12:16].strip() in ['N','CA','C','O']:bb.append(l[21:27]+l[12:16]+l[30:54])
 return ca,hashlib.sha256('\n'.join(bb).encode()).hexdigest()
reference_bb={s['id']:pdb_identity_and_backbone(s['path'])[1] for s in states if s['is_reference']}
counts={'full':0,'intrachain':0};max_intrachain_deltas=0.
for s in states:
 assert sha(s['path'])==s['sha256']
 ca,bb=pdb_identity_and_backbone(s['path']);assert bb==reference_bb[s['reference_state']]
 assert ''.join(aa for (c,a),aa in ca.items() if c=='C')==s['sequence']
 for v in ['full','intrachain']:
  z=np.load(R/'predictions'/v/(s['id']+'.npz'))
  assert str(z['input_sha256'])==s['sha256']
  assert list(zip(z['chain'],z['author']))==list(ca)
  assert ''.join(z['aa'])==''.join(ca.values())
  assert z['mean'].shape==(len(ca),2) and np.isfinite(z['mean'][z['model_mask']]).all()
  assert np.array_equal(z['hchi_valid'],z['model_mask']&z['geometry_mask']&~np.isin(z['aa'],['A','G']))
  counts[v]+=len(ca)
  if v=='intrachain' and not s['is_reference']:
   ref=np.load(R/'predictions'/v/(s['reference_state']+'.npz'))
   tcr=np.isin(z['chain'],['D','E'])
   max_intrachain_deltas=max(max_intrachain_deltas,float(np.max(np.abs(z['mean'][tcr]-ref['mean'][tcr]))))

fits=[]
for tv in ['primary','secondary']:
 O=R/'evaluation'/('full_'+tv);f=pd.read_csv(O/'features_and_labels.csv').set_index('candidate_id')
 memberships=json.loads((O/'fit_membership.json').read_text())
 reported=pd.read_csv(O/'all_heldout_predictions.csv')
 for record in memberships:
  if record['split']!='leave_system_out' or record['model'] not in ['B0_binding','B1_binding_static','B2_binding_static_dynamics']:continue
  train=f.loc[record['training_ids']];test=f.loc[record['test_ids']];cols=record['features']
  assert set(train.system).isdisjoint(set(test.system)) and len(train)==186 and len(test)==93
  x=train[cols].to_numpy(float);xt=test[cols].to_numpy(float);assert np.isfinite(x).all() and np.isfinite(xt).all()
  center=x.mean(axis=0);scale=x.std(axis=0);scale[scale==0]=1.
  x=(x-center)/scale;xt=(xt-center)/scale;y=train.functional_positive.to_numpy(int)
  weights=np.where(y==1,len(y)/(2*sum(y==1)),len(y)/(2*sum(y==0)))
  def objective(beta):
   lin=x@beta[:-1]+beta[-1];error=weights*(expit(lin)-y)
   loss=float(np.sum(weights*(np.logaddexp(0.,lin)-y*lin))+.5*np.sum(beta[:-1]**2))
   grad=np.r_[x.T@error+beta[:-1],error.sum()]
   return loss,grad
  opt=minimize(objective,np.zeros(x.shape[1]+1),jac=True,method='BFGS',options={'gtol':1e-7,'maxiter':2000})
  prob=expit(xt@opt.x[:-1]+opt.x[-1]);original=reported[(reported.split=='leave_system_out')&(reported.model==record['model'])&(reported.fold==record['fold'])].set_index('candidate_id').loc[test.index].score.to_numpy()
  ordered=test.assign(score=prob).sort_values(['score','candidate_id'],ascending=[False,True]) if 'candidate_id' in test.columns else test.assign(score=prob).reset_index().sort_values(['score','candidate_id'],ascending=[False,True])
  old=test.assign(score=original).reset_index().sort_values(['score','candidate_id'],ascending=[False,True])
  fits.append({'template_version':tv,'system':record['fold'],'model':record['model'],'max_probability_difference':float(np.max(np.abs(prob-original))),
    'independent_top10_hits':int(ordered.head(10).functional_positive.sum()),'original_top10_hits':int(old.head(10).functional_positive.sum()),
    'objective_gradient_max':float(np.max(np.abs(objective(opt.x)[1]))),'scipy_success':bool(opt.success),'scipy_message':str(opt.message)})
pd.DataFrame(fits).to_csv(R/'audits/independent_weighted_logistic_reproduction.csv',index=False)
report={'status':'checked','prediction_npz_files':len(states)*2,'input_PDBs_checked':len(states),'all_backbone_and_sequence_checks_passed':True,'residue_rows':counts,
 'intrachain_max_receptor_change_across_peptides':max_intrachain_deltas,'independent_logistic_fits':len(fits),
 'max_reproduced_probability_difference':max(r['max_probability_difference'] for r in fits),'top10_hit_counts_all_agree':all(r['independent_top10_hits']==r['original_top10_hits'] for r in fits),
 'numerical_note':'Independent BFGS minimization of weighted logistic sum loss plus 0.5*L2 penalty; sklearn uses a different stopping rule. Any probability or top10 differences retained, not hidden.',
 'source_sha256':{str(p):sha(p) for p in [R/'states.json',R/'protocol_frozen_before_new_predictions.json',Path(__file__)]}}
(R/'audits/final_numerical_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
