"""Small read-only check of the input structures behind the KRAS narrative."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
def read_ca(p):
    result={}
    for l in p.read_text().splitlines():
        if l.startswith('ATOM  ') and l[12:16].strip()=='CA':
            key=(l[21],l[22:27].strip())
            assert key not in result
            result[key]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return result
a=read_ca(ROOT/'replica_recheck/7ow5/7ow5_model_input.pdb')
b=read_ca(ROOT/'replica_recheck/7ow6/7ow6_model_input.pdb')
keys=list(a)
assert set(keys)==set(b)
fit=[k for k in keys if k[0]=='A' and int(k[1])<=180]
aa=np.array([a[k] for k in fit]);bb=np.array([b[k] for k in fit]);ac=aa.mean(0);bc=bb.mean(0)
u,s,vt=np.linalg.svd((bb-bc).T@(aa-ac));rotation=u@vt
if np.linalg.det(rotation)<0:u[:,-1]*=-1;rotation=u@vt
out={'alignment':'Fit input B=7OW6 to input A=7OW5 using first180 MHC C-alpha atoms, no scale change. All reported distances Angstrom.',
 'MHC_groove_n':len(fit),'static_RMSD_A':{}}
for label,select in [('MHC_groove',lambda k:k in fit),('peptide_all10',lambda k:k[0]=='C'),('TCR_variable_first125',lambda k:k[0] in 'DE' and int(k[1])<=125)]:
    kk=[k for k in keys if select(k)]
    delta=np.array([(b[k]-bc)@rotation+ac-a[k] for k in kk])
    out['static_RMSD_A'][label]={'n':len(kk),'RMSD':float(np.sqrt(np.mean(np.sum(delta**2,axis=1))))}
f=pd.read_csv(ROOT/'promising_cases_replica_validation/7ow5_7ow6_all_positions_three_replicas.csv')
p=f[f.role.eq('peptide')]
out['backbone_MD_and_predictions']=[{'position':int(r.position_a),'WT_AA':r.aa_a,'G12D_AA':r.aa_b,'MD_delta_A':float(r.rmsf_md_delta),'Full_delta_A':float(r.rmsf_full_delta),'old_strict_delta_A':float(r.rmsf_old_strict_delta)} for r in p.itertuples()]
(ROOT/'KRAS_story_static_and_dynamic_checks.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(out,ensure_ascii=False,indent=2))
