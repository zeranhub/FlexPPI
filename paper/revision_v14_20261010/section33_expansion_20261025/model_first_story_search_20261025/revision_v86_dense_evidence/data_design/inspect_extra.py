from pathlib import Path
import pandas as pd,json
from scipy.stats import spearmanr
p=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v86_dense_evidence\data_design')
a=pd.read_csv(p/'all25_dynamics_and_static_lobe_contacts.csv');v=a[a.valid_chi_entropy]
print('Full below Intra rmsf',sum(a.pred_rmsf<a.intra_rmsf),'of',len(a),'Hchi',sum(v.pred_chi_entropy<v.intra_chi_entropy),'of',len(v))
print('exceptions',a.loc[a.pred_rmsf>=a.intra_rmsf,['CaMKI_native_position','amino_acid','pred_rmsf','intra_rmsf']].to_string(index=False))
print('Hchi exceptions',v.loc[v.pred_chi_entropy>=v.intra_chi_entropy,['CaMKI_native_position','amino_acid','pred_chi_entropy','intra_chi_entropy']].to_string(index=False))
for region,x in a.groupby('region'):
 print(region,len(x),'totalcontacts4',x.n_CaM_contact_residues4.sum(),'N',x.N_lobe_count_4A.sum(),'C',x.C_lobe_count_4A.sum(),'shieldfraction median',x.sidechain_shielding_fraction.median(),'sideSASA median',x.sidechain_SASA_complex_A2.median())
for aa,g in v.groupby('amino_acid'):
 if len(g)>1:
  print(aa,g[['CaMKI_native_position','pred_chi_entropy','intra_chi_entropy','sidechain_SASA_complex_A2','sidechain_shielding_fraction','target_chi_entropy']].to_string(index=False))
print('pair 303',a.loc[a.CaMKI_native_position==303, ['C_lobe_residues_4A','C_lobe_residues_4p5A','C_lobe_residues_5A']].to_dict('records'))
print('pair 316',a.loc[a.CaMKI_native_position==316, ['N_lobe_residues_4A','N_lobe_residues_4p5A','N_lobe_residues_5A']].to_dict('records'))