import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np,pandas as pd,pyarrow.parquet as pq
from difflib import SequenceMatcher
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
R=Path('E:/Project/FlexPPI/paper/revision_v14_20261010/section33_expansion_20261025')
B=R/'feature_baseline';O=Path(__file__).resolve().parent
sid='A006W'
f=pq.read_table(B/'features_all.parquet',filters=[('sample_id','=',sid)]).to_pandas()
p=pq.read_table(B/'oof_predictions_all.parquet',filters=[('sample_id','=',sid)]).to_pandas()
pc=[c for c in p if c.startswith('pred_')]
f=f.drop(columns=[c for c in pc if c in f.columns])
f=f.merge(p[['sample_id','residue_index']+pc],on=['sample_id','residue_index'],validate='one_to_one')
path=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/dataset/structures/A006W.pdb')
atoms={};resnames={}
for l in path.read_text().splitlines():
 if l.startswith('ATOM') and l[12:16].strip()[0] not in 'HD':
  k=(l[21],l[22:27].strip());atoms.setdefault(k,[]).append([float(l[30:38]),float(l[38:46]),float(l[46:54])]);resnames[k]=l[17:20]
peptide=np.vstack([v for (c,r),v in atoms.items() if c in 'DE'])
dist=[]
for (c,r),a in atoms.items():
 if c in 'ABC':
  dist.append(dict(chain=c,author_residue_id=r,resname=resnames[(c,r)],p15_min_distance=float(cdist(a,peptide).min())))
d=pd.DataFrame(dist);f.author_residue_id=f.author_residue_id.astype(str);f=f.merge(d,on=['chain','author_residue_id'],how='left',validate='one_to_one')
f.to_csv(O/'pcna_all_residues_with_p15_distances.csv',index=False)
print(f.groupby('chain').agg(n=('residue_index','size'),p15contacts=('p15_min_distance',lambda x:(x<=5).sum())).to_string(),flush=True)
reference=f[f.chain=='A'].sort_values('residue_index')
maps={}
for ch in 'ABC':
 part=f[f.chain==ch].sort_values('residue_index');seq=''.join(part.amino_acid);refseq=''.join(reference.amino_acid)
 matcher=SequenceMatcher(None,refseq,seq,autojunk=False)
 mapping={}
 for a,b,n in matcher.get_matching_blocks():
  for offset in range(n):mapping[int(part.iloc[b+offset].residue_index)]=a+offset+1
 assert len(mapping)>=250
 maps[ch]=mapping
f['reference_position']=[maps.get(c,{}).get(int(i),np.nan) for c,i in zip(f.chain,f.residue_index)]
f.to_csv(O/'pcna_all_residues_with_p15_distances.csv',index=False)
q=f[f.chain.isin(list('ABC'))&f.valid_hchi&f.reference_position.notna()].copy()
out=[]
for ch1,ch2 in [('A','B'),('A','C'),('B','C')]:
 a=q[q.chain==ch1].set_index('reference_position');b=q[q.chain==ch2].set_index('reference_position');common=a.index.intersection(b.index);a=a.loc[common];b=b.loc[common]
 assert (a.amino_acid.values==b.amino_acid.values).all()
 for region,ix in [('all',np.ones(len(a),bool)),('p15_contact_union5',(a.p15_min_distance<=5)|(b.p15_min_distance<=5)),('p15_contact_union8',(a.p15_min_distance<=8)|(b.p15_min_distance<=8)),('asymmetric_contact',(a.p15_min_distance<=5)&(b.p15_min_distance>8)|((b.p15_min_distance<=5)&(a.p15_min_distance>8)))]:
  aa=a[ix];bb=b[ix];row=dict(pair=ch1+'-'+ch2,region=region,n=len(aa))
  if not len(aa):continue
  for metric in ['hchi','rmsf']:
   y=bb['md_'+metric].values-aa['md_'+metric].values;row[metric+'_md_delta_mean']=y.mean();row[metric+'_md_absdelta_mean']=abs(y).mean()
   for method in ['full','intra','aa_sasa_contacts','static47_sasa']:
    pred=bb['pred_'+method+'_'+metric].values-aa['pred_'+method+'_'+metric].values;row[metric+'_'+method+'_delta_mae']=abs(y-pred).mean();row[metric+'_'+method+'_delta_rho']=float(spearmanr(y,pred).statistic) if len(y)>=4 else np.nan
  out.append(row)
pd.DataFrame(out).to_csv(O/'pcna_equivalent_site_control_metrics.csv',index=False)
print(pd.DataFrame(out)[['pair','region','n','hchi_md_absdelta_mean','hchi_full_delta_mae','hchi_intra_delta_mae','hchi_aa_sasa_contacts_delta_mae','hchi_static47_sasa_delta_mae','hchi_full_delta_rho','rmsf_md_absdelta_mean']].to_string(index=False),flush=True)
