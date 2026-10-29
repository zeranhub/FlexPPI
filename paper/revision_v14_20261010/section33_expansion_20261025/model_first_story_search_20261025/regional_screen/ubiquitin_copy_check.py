from pathlib import Path
from difflib import SequenceMatcher
import numpy as np
import pandas as pd
from screen_regions import rho, METHODS
OUT=Path(__file__).resolve().parent
d=pd.read_csv(OUT/'seven_candidate_residue_predictions.csv')
a=d[(d.sample_id=='A006F')&(d.chain=='A')].reset_index(drop=True)
b=d[(d.sample_id=='A006F')&(d.chain=='C')].reset_index(drop=True)
blocks=SequenceMatcher(None,''.join(a.amino_acid),''.join(b.amino_acid),autojunk=False).get_matching_blocks()
ia=[i for m in blocks for i in range(m.a,m.a+m.size)]
ib=[i for m in blocks for i in range(m.b,m.b+m.size)]
q=pd.concat([a.iloc[ia].reset_index(drop=True).add_suffix('_A'),b.iloc[ib].reset_index(drop=True).add_suffix('_C')],axis=1)
assert (q.amino_acid_A==q.amino_acid_C).all()
out=[]
for met in ['rmsf','hchi']:
    valid=q['valid_'+met+'_A'] & q['valid_'+met+'_C']
    if met=='hchi':valid &= q.sidechain_complete_A & q.sidechain_complete_C & (q.amino_acid_A!='P')
    z=q[valid]
    y=z['md_'+met+'_C']-z['md_'+met+'_A']
    q.loc[valid,'md_C_minus_A_'+met]=y
    for m in METHODS:
        p=z['pred_'+m+'_'+met+'_C']-z['pred_'+m+'_'+met+'_A']
        q.loc[valid,m+'_C_minus_A_'+met]=p
        out.append(dict(metric=met,method=m,n=len(z),md_mean=y.mean(),pred_mean=p.mean(),delta_mae=abs(y-p).mean(),zero_mae=abs(y).mean(),rho=rho(y,p),directions=(np.sign(y)==np.sign(p)).mean()))
q.to_csv(OUT/'ubiquitin_3WWQ_copy_delta_residues.csv',index=False)
pd.DataFrame(out).to_csv(OUT/'ubiquitin_3WWQ_copy_delta_metrics.csv',index=False)
print(pd.DataFrame(out).to_string(index=False))
print(q[['author_residue_id_A','author_residue_id_C','amino_acid_A','md_C_minus_A_rmsf','full_C_minus_A_rmsf','md_C_minus_A_hchi','full_C_minus_A_hchi']].sort_values('md_C_minus_A_hchi').to_string(index=False))
