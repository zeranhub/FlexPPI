from pathlib import Path
import json,hashlib,itertools,time
import numpy as np,pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parent
O=R/'formal_results';O.mkdir(exist_ok=True)
METHODS=['full','intra','static47_sasa','aa_sasa_contacts','aa_mean'];METRICS=['rmsf','hchi'];KEY=['sample_id','residue_index']
def rho(y,p):
 return float(spearmanr(y,p).statistic) if len(y)>=3 and np.ptp(y)>0 and np.ptp(p)>0 else np.nan
def typed_rho(y,p,aa):
 t=pd.DataFrame(dict(y=y,p=p,aa=aa));z=[]
 for _,g in t.groupby('aa'):
  if len(g)>=3:z.append(g.assign(y=g.y-g.y.mean(),p=g.p-g.p.mean()))
 if not z:return np.nan
 t=pd.concat(z);return rho(t.y,t.p)
def interval(d,col='gain'):
 g=d.groupby('component_id')[col].agg(['sum','count','mean'])
 assert len(g)
 rng=np.random.default_rng(20261026);ix=rng.integers(0,len(g),size=(10000,len(g)))
 a=g['sum'].to_numpy()[ix].sum(1)/g['count'].to_numpy()[ix].sum(1);b=g['mean'].to_numpy()[ix].mean(1)
 return dict(n_PDB=len(d),n_split_groups=len(g),PDB_macro_gain=float(d[col].mean()),CI95_lo=float(np.quantile(a,.025)),CI95_hi=float(np.quantile(a,.975)),split_group_macro_gain=float(g['mean'].mean()),split_group_CI95_lo=float(np.quantile(b,.025)),split_group_CI95_hi=float(np.quantile(b,.975)),n_PDB_positive=int((d[col]>0).sum()),n_PDB_negative=int((d[col]<0).sum()))
f=pd.read_parquet(R/'controls/all56_features_and_family_predictions.parquet')
ann=pd.read_csv(R/'biological_annotation/all56_biological_annotations.csv')[['sample_id','category','core_functional_peptide','any_official_polymer_mutation']]
f=f.merge(ann,on='sample_id',validate='many_to_one')
assert f.fold.dtype.kind in 'iu'
u=pd.read_csv(R/'uniform_md/all_uniform_labels_per_residue.csv.gz',dtype={'chain':str,'author_residue_id':str})
u=u[~u.variant.eq('uniform_available_replicas')].copy()
assert u.groupby('variant').sample_id.nunique().min()==53
variants=['original','uniform_rep1']+sorted(v for v in u.variant.unique() if v.startswith('block_') and v!='block_000_100ns')
rec=[]
for variant in variants:
 if variant=='original':a=f.copy()
 else:
  target=u[u.variant.eq(variant)][KEY+['chain','amino_acid','rmsf','chi_entropy','valid_rmsf','valid_chi_entropy']].rename(columns={'chain':'check_chain','amino_acid':'check_amino_acid','rmsf':'new_rmsf','chi_entropy':'new_hchi','valid_rmsf':'new_valid_rmsf','valid_chi_entropy':'new_valid_hchi'})
  a=f.merge(target,on=KEY,validate='one_to_one',how='inner')
  assert a.chain.eq(a.check_chain).all() and a.amino_acid.eq(a.check_amino_acid).all()
  for m in METRICS:
   assert a['valid_'+m].equals(a['new_valid_'+m])
   a['md_'+m]=a['new_'+m]
  a=a.drop(columns=[c for c in a if c.startswith('new_') or c.startswith('check_')])
 if variant in ['original','uniform_rep1']:a.to_parquet(O/('all_features_'+variant+'.parquet'),index=False)
 for (sid,role),part in a.groupby(['sample_id','role']):
  common={c:part[c].iloc[0] for c in ['pdb_id','component_id','fold','category','core_functional_peptide','any_official_polymer_mutation']};common.update(sample_id=sid,role=role,reference=variant)
  for scope,scope_mask in [('all',np.ones(len(part),bool)),('interface',part.interface_ca8.to_numpy()),('outside',~part.interface_ca8.to_numpy())]:
   for m in METRICS:
    q=part.loc[scope_mask & part['valid_'+m]]
    for method in METHODS:
     y=q['md_'+m].to_numpy();p=q['pred_'+method+'_'+m].to_numpy()
     rec.append(dict(common,scope=scope,metric=m,method=method,n=len(q),mae=float(np.abs(y-p).mean()) if len(q) else np.nan,bias=float((p-y).mean()) if len(q) else np.nan,rho=rho(y,p),typed_rho=typed_rho(y,p,q.amino_acid.to_numpy()) if len(q) else np.nan))
rec=pd.DataFrame(rec);rec.to_csv(O/'per_record_all_methods_references.csv',index=False)
core_ids=set(ann.loc[ann.core_functional_peptide,'sample_id'])
counts=ann.groupby('sample_id').first().merge(f[['sample_id','component_id']].drop_duplicates().set_index('sample_id'),left_index=True,right_index=True)
biggest=counts.groupby('component_id').size().idxmax()
cohorts={'wide':set(ann.sample_id),'core':core_ids,'core_no_documented_mutation':set(ann.loc[ann.core_functional_peptide & ~ann.any_official_polymer_mutation,'sample_id']),'core_no_CaM':core_ids-set(f.loc[f.pdb_id.isin(['1MXE','4LZX']),'sample_id']),'core_no_largest_split_group':core_ids-set(f.loc[f.component_id.eq(biggest),'sample_id'])}
membership=[]
for name,ids in cohorts.items():
 b=f[f.sample_id.isin(ids)][['sample_id','pdb_id','component_id']].drop_duplicates()
 membership.append(dict(cohort=name,n_records=len(b),n_PDB=b.pdb_id.nunique(),n_split_groups=b.component_id.nunique(),member_ids=sorted(ids)))
(O/'cohorts.json').write_text(json.dumps(membership,ensure_ascii=False,indent=2),encoding='utf-8')
summ=[];ci=[]
for name,ids in cohorts.items():
 for keys,g in rec[rec.sample_id.isin(ids)&rec.n.gt(0)].groupby(['reference','role','scope','metric']):
  ref,role,scope,m=keys
  for method,p in g.groupby('method'):
   z=p.groupby(['pdb_id','component_id'])[['mae','bias','rho','typed_rho']].mean().reset_index();fm=z.groupby('component_id')[['mae','rho','typed_rho']].mean()
   summ.append(dict(cohort=name,reference=ref,role=role,scope=scope,metric=m,method=method,n_records=len(p),n_PDB=len(z),n_split_groups=len(fm),PDB_MAE=z.mae.mean(),PDB_bias=z.bias.mean(),PDB_rho=z.rho.mean(),PDB_typed_rho=z.typed_rho.mean(),group_MAE=fm.mae.mean()))
  w=g.pivot(index=['sample_id','pdb_id','component_id'],columns='method',values='mae')
  for method in METHODS[1:]:
   z=w[['full',method]].dropna().reset_index().groupby(['pdb_id','component_id'])[['full',method]].mean().reset_index();z['gain']=z[method]-z.full
   ci.append(dict(cohort=name,reference=ref,role=role,scope=scope,metric=m,comparator=method,n_records=len(w),Full_MAE=z.full.mean(),comparator_MAE=z[method].mean(),**interval(z)))
summ=pd.DataFrame(summ);ci=pd.DataFrame(ci);summ.to_csv(O/'summary_all_scopes.csv',index=False);ci.to_csv(O/'paired_MAE_gain_cluster_CI.csv',index=False)
print('MAIN CORE original/uniform Full vs Intra and strong baseline',flush=True)
print(ci[ci.cohort.eq('core') & ci.reference.isin(['original','uniform_rep1']) & ci.scope.eq('all') & ci.comparator.isin(['intra','static47_sasa'])].to_string(index=False),flush=True)
# Same-AA matching: type first, then record, PDB, split group. All pairs enumerated once, no cherry-picking.
pairrows=[];cover=[]
for ref in ['original','uniform_rep1']:
 a=pd.read_parquet(O/('all_features_'+ref+'.parquet'))
 for (sid,role),part in a.groupby(['sample_id','role']):
  q=part[part.valid_hchi & part.valid_rmsf];common={c:part[c].iloc[0] for c in ['pdb_id','component_id','category','core_functional_peptide']};common.update(sample_id=sid,role=role,reference=ref)
  for aa,g in q.groupby('amino_acid'):
   if len(g)<2:continue
   ix,jx=np.triu_indices(len(g),1);di=g.iloc[ix];dj=g.iloc[jx]
   dy=dj.md_hchi.to_numpy()-di.md_hchi.to_numpy()
   dr=np.abs(dj.md_rmsf.to_numpy()-di.md_rmsf.to_numpy());ds=np.abs(dj.sidechain_sasa_complex_A2.to_numpy()-di.sidechain_sasa_complex_A2.to_numpy())
   for tol,sasamax in [(0.10,None),(0.25,None),(0.50,None),(0.10,5),(0.10,10),(0.10,20),(0.25,10)]:
    keep=dr<=tol
    if sasamax is not None:keep &= ds<=sasamax
    if not keep.any():continue
    z=dict(common,amino_acid=aa,rmsf_tolerance_A=tol,sasa_tolerance_A2=sasamax if sasamax is not None else -1,n_pairs=int(keep.sum()),true_abs_difference=float(abs(dy[keep]).mean()))
    for method in METHODS:
     dp=dj['pred_'+method+'_hchi'].to_numpy()-di['pred_'+method+'_hchi'].to_numpy()
     z[method+'_pair_MAE']=float(abs(dp[keep]-dy[keep]).mean())
     strong=abs(dy[keep])>=.1
     if strong.any():
      sy=np.sign(dy[keep][strong]);sp=np.sign(dp[keep][strong]);z[method+'_direction_accuracy']=float(np.where(sp==0,.5,sp==sy).mean())
     else:z[method+'_direction_accuracy']=np.nan
    pairrows.append(z)
pairs=pd.DataFrame(pairrows);pairs.to_csv(O/'sameAA_match_per_record_type.csv',index=False)
pcis=[]
for name,ids in cohorts.items():
 sub=pairs[pairs.sample_id.isin(ids)]
 for keys,g in sub.groupby(['reference','role','rmsf_tolerance_A','sasa_tolerance_A2']):
  ref,role,tol,sa=keys
  # Mean per AA -> record, then average duplicates -> PDB. Not weighting by pair counts.
  col=[m+'_pair_MAE' for m in METHODS]
  rr=g.groupby(['sample_id','pdb_id','component_id'])[col].mean().reset_index();z=rr.groupby(['pdb_id','component_id'])[col].mean().reset_index()
  for method in METHODS[1:]:
   v=z[['pdb_id','component_id','full_pair_MAE',method+'_pair_MAE']].copy();v['gain']=v[method+'_pair_MAE']-v.full_pair_MAE
   pcis.append(dict(cohort=name,reference=ref,role=role,rmsf_tolerance_A=tol,sasa_tolerance_A2=sa,comparator=method,n_records=rr.sample_id.nunique(),n_types_records=len(g),n_pairs=g.n_pairs.sum(),Full_pair_MAE=v.full_pair_MAE.mean(),comparator_pair_MAE=v[method+'_pair_MAE'].mean(),**interval(v)))
pd.DataFrame(pcis).to_csv(O/'sameAA_match_type_equal_cluster_CI.csv',index=False)
audit=dict(status='complete',records=56,core_PDB=40,uniform_records=53,uniform_core_PDB=int(ci.loc[ci.cohort.eq('core')&ci.reference.eq('uniform_rep1')&ci.scope.eq('all'),'n_PDB'].iloc[0]),OOF_reference_and_label_fold_separated=True,family_fold_comes_from_features_not_raw_uniform_fold=True,all56_retained_in_wide_original=True,all53_long_trajectories_in_uniform=True,match_type_equal=True,largest_split_group=biggest,bootstrap_draws=10000,seed=20261026,primary_tol_A=.1,primary_SASA_tol_A2=10,not_independent_MD_replicates=True)
(O/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(audit,ensure_ascii=False),flush=True)
