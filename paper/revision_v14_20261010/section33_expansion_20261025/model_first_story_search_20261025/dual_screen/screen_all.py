import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr

R=Path('E:/Project/FlexPPI/paper/revision_v14_20261010/section33_expansion_20261025')
B=R/'feature_baseline'
O=Path(__file__).resolve().parent
METHODS=['full','intra','aa_sasa_contacts','static47_sasa']
cols=['sample_id','chain','residue_index','author_residue_id','amino_acid','valid_hchi','valid_rmsf','sidechain_complete','md_hchi','md_rmsf','sidechain_sasa_complex_A2','sidechain_sasa_isolated_A2','interface_ca8','partner_ca_contacts_8A','distance_to_chain_end']
f=pq.read_table(B/'features_all.parquet',columns=cols).to_pandas()
pcols=['sample_id','residue_index']+['pred_'+m+'_'+v for m in METHODS for v in ['rmsf','hchi']]
p=pq.read_table(B/'oof_predictions_all.parquet',columns=pcols).to_pandas()
f=f.merge(p,on=['sample_id','residue_index'],validate='one_to_one')
meta=pd.read_csv(B/'cohort_records.csv').set_index('sample_id')
def rho(x,y):
    return float(spearmanr(x,y).statistic) if len(x)>=8 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
def centered(q,key):
    y=q[key].to_numpy(float).copy();aa=q.amino_acid.to_numpy()
    for a in set(aa):
        ix=aa==a;y[ix]-=y[ix].mean()
    return y
out=[]
for (sid,ch),q in f.groupby(['sample_id','chain'],sort=False):
    q=q[q.valid_hchi&q.valid_rmsf&q.sidechain_complete&(q.amino_acid!='P')&(q.distance_to_chain_end>3)].copy()
    if len(q)<60:continue
    counts=q.amino_acid.value_counts();q=q[q.amino_acid.isin(counts[counts>=4].index)]
    masks={'all':np.ones(len(q),bool),'core_sasa5':q.sidechain_sasa_complex_A2<=5,'low_rmsf_08':q.md_rmsf<=.8,'core_low_rmsf':(q.sidechain_sasa_complex_A2<=5)&(q.md_rmsf<=.8),'interface':q.interface_ca8,'interface_low_rmsf':q.interface_ca8&(q.md_rmsf<=.8),'surface_low_rmsf':(q.sidechain_sasa_complex_A2>=30)&(q.md_rmsf<=.8)}
    for mask,ix in masks.items():
        t=q[ix].copy();cs=t.amino_acid.value_counts();t=t[t.amino_acid.isin(cs[cs>=3].index)]
        if len(t)<25:continue
        y=t.md_hchi.to_numpy();r=t.md_rmsf.to_numpy();yc=centered(t,'md_hchi')
        d=dict(sample_id=sid,chain=ch,region=mask,n=len(t),n_aa=t.amino_acid.nunique(),hchi_sd=y.std(),hchi_iqr=np.quantile(y,.75)-np.quantile(y,.25),rmsf_median=np.median(r),rmsf_iqr=np.quantile(r,.75)-np.quantile(r,.25),rmsf_hchi_rho=rho(r,y),pdb_id=meta.loc[sid,'pdb_id'])
        for m in METHODS:
            yp=t['pred_'+m+'_hchi'].to_numpy();d[m+'_mae']=np.abs(y-yp).mean();d[m+'_rho']=rho(y,yp);d[m+'_typed_rho']=rho(yc,centered(t,'pred_'+m+'_hchi'));d[m+'_rmsf_mae']=np.abs(r-t['pred_'+m+'_rmsf']).mean()
        d['gain_strong']=min(d[m+'_mae'] for m in METHODS[2:])-d['full_mae']
        d['typed_gain']=d['full_typed_rho']-max(d[m+'_typed_rho'] for m in METHODS[2:])
        d['score']=d['gain_strong']*4+d['full_typed_rho']*.1+d['typed_gain']*.2+np.log(len(t))*.01
        out.append(d)
res=pd.DataFrame(out).sort_values('score',ascending=False)
desc=pd.read_csv('E:/Project/FlexPPI/experiments/non_tcr_case_screen_20261025/all_chain_metrics.csv')[['sample_id','chain','description']].drop_duplicates(['sample_id','chain'])
res=res.merge(desc,on=['sample_id','chain'],how='left')
res.to_csv(O/'all_chain_region_screen.csv',index=False)
short=res[(res.region!='all')&(res.n>=35)&(res.full_typed_rho>=.50)&(res.gain_strong>=.03)&(res.typed_gain>=.10)].copy()
short.to_csv(O/'shortlist_region_screen.csv',index=False)
top=short.drop_duplicates('pdb_id').head(30)
f[f.sample_id.isin(top.sample_id)].to_parquet(O/'top30_residues.parquet',index=False)
top.to_csv(O/'top30_regions.csv',index=False)
print(top[['sample_id','pdb_id','chain','region','n','full_mae','aa_sasa_contacts_mae','static47_sasa_mae','full_typed_rho','typed_gain','description']].to_string(index=False),flush=True)
