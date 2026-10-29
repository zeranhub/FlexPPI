from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
import numpy as np,pandas as pd
from scipy.stats import spearmanr

rows=[];details=[]
for sample,pepchain,pdbs in [('A0053','C',['3h9s']),('A0055','C',['2bnq','6q3s','6rp9','6rpa','6rpb']),('A005I','E',['2ak4','4jrx','4jry'])]:
    a=np.load(ROOT/f'{sample}.npz')
    ai=np.flatnonzero(a['chain']==pepchain)
    for pdb in pdbs:
        b=np.load(EXP/f'tcr_case_v40_20261017/family_predictions/full/itrace_{pdb}.npz')
        bi=np.flatnonzero(b['chain']=='C')
        assert ''.join(a['amino_acid'][ai])==''.join(b['amino_acid'][bi])
        for m,name in enumerate(['rmsf','Hchi']):
            valid=a['valid'][ai,m]&b['valid'][bi,m]
            idx=np.flatnonzero(valid)
            y=b['target'][bi[idx],m]-a['target'][ai[idx],m]
            p=b['mean'][bi[idx],m]-a['mean'][ai[idx],m]
            mae=np.mean(abs(y-p));zero=np.mean(abs(y))
            rows.append({'a':sample,'b':'itrace_'+pdb,'metric':name,'n':len(y),'change_rho':float(spearmanr(y,p).statistic),
                'change_MAE':float(mae),'zero_change_MAE':float(zero),'relative_MAE_reduction':float(1-mae/zero),
                'mean_MD_delta':float(np.mean(y)),'mean_model_delta':float(np.mean(p)),
                'sign_agreement':float(np.mean(np.sign(y)==np.sign(p))),
                'interpretation_scope':'cross_MD_protocol single local DynaRepo trajectory vs triplicate ITrace; target source effects must not be attributed exclusively to biological binding'})
            for k,t,pr in zip(idx,y,p):details.append({'a':sample,'b':'itrace_'+pdb,'position':int(k)+1,'aa':str(a['amino_acid'][ai[k]]),'metric':name,'MD_delta':float(t),'model_delta':float(pr)})
pd.DataFrame(rows).to_csv(ROOT/'cross_source_peptide_change_metrics.csv',index=False)
pd.DataFrame(details).to_csv(ROOT/'cross_source_peptide_all_positions.csv',index=False)
print(pd.DataFrame(rows).drop(columns='interpretation_scope').to_string(index=False))
