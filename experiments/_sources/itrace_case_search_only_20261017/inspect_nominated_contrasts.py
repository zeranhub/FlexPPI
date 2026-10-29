from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'label_sensitivity/runtime'))
import pandas as pd
x=pd.read_csv(ROOT/'numeric_rescreen/all_models_complete_region_comparison.csv')
pairs=[('1lp9','2j8u'),('1lp9','2jcc'),('1lp9','2uwe'),('7pdw','7pbc'),('8rym','8ryn'),('8rym','8ryo'),('2bnq','2bnr'),('8wte','8wul'),('8yiv','8yj2'),('7q99','7q9a'),('4qok','7q99'),('4qok','7q9a')]
parts=[]
for a,b in pairs:
    f=x[((x.a==a)&(x.b==b))|((x.a==b)&(x.b==a))]
    f=f[f.model.isin(['full','old_strict'])&f.scope.isin(['all_CDR','CDR3','whole_peptide','natural_pMHC_interface'])].copy()
    f['biological_contrast']=a+'→'+b
    parts.append(f)
out=pd.concat(parts,ignore_index=True)
out.to_csv(ROOT/'nominated_biological_contrast_model_comparison.csv',index=False,encoding='utf-8-sig')
columns=['biological_contrast','scope','model','n_rmsf','n_hchi','rho_delta_rmsf','rho_delta_hchi','relative_MAE_gain_rmsf','relative_MAE_gain_hchi','descriptor_above_2SEM_n_rmsf','descriptor_above_2SEM_n_hchi']
print(out[columns].round(3).to_string(index=False))
