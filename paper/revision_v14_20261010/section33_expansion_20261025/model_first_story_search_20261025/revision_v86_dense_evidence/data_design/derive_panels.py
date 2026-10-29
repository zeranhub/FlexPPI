from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
src=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v85_1MXE\data_audit')
out=src.parent.parent/'revision_v86_dense_evidence'/'data_design'
out.mkdir(parents=True,exist_ok=True)
a=pd.read_csv(src/'peptide25_model_MD_structure_all_values.csv')
d=pd.read_csv(src/'all25_by_all144_CaM_residue_distances.csv')
# All 25 sites, every lobe and threshold retained; descriptive, static original coordinates.
rows=[]
for pos,g in d.groupby('CaMKI_native_position'):
 r=a[a.CaMKI_native_position==pos].iloc[0].to_dict()
 for lobe,z in g.groupby('CaM_lobe'):
  r[f'{lobe}_min_heavy_A']=z.minimum_heavy_atom_distance_A.min()
  for dist,c in [(4,'heavy_atom_contact4'),(4.5,'heavy_atom_contact4_5'),(5,'heavy_atom_contact5')]:
   v=z[z[c]]
   r[f'{lobe}_count_{str(dist).replace(".","p")}A']=len(v)
   r[f'{lobe}_residues_{str(dist).replace(".","p")}A']=';'.join(v.CaM_aa.astype(str)+v.CaM_native_id.astype(str))
 rows.append(r)
b=pd.DataFrame(rows)
b.to_csv(out/'all25_dynamics_and_static_lobe_contacts.csv',index=False)
d.to_csv(out/'all3600_CaMKI_CaM_distances.csv',index=False)
print('DISTANCE rows',len(d),'lobes',d.groupby('CaM_lobe').size().to_dict(),'dist range',d.minimum_heavy_atom_distance_A.min(),d.minimum_heavy_atom_distance_A.max())
for cutoff,col in [(4,'heavy_atom_contact4'),(4.5,'heavy_atom_contact4_5'),(5,'heavy_atom_contact5')]:
 x=d[d[col]]
 print('CONTACT',cutoff,'pairs',len(x),'CaMres',x.CaM_native_id.nunique(),'peptide',x.CaMKI_native_position.nunique(),'counts',x.groupby('CaM_lobe').size().to_dict())
 print(x.groupby('CaMKI_native_position').size().to_dict())
print(b[['CaMKI_native_position','amino_acid','region','N_lobe_count_4A','C_lobe_count_4A','linker_count_4A','sidechain_SASA_complex_A2','sidechain_shielding_fraction','pred_rmsf','intra_rmsf','pred_chi_entropy','intra_chi_entropy']].to_string(index=False))
print('valid Hchi',sum(a.valid_chi_entropy))
print('Full RMSF range',a.pred_rmsf.min(),a.pred_rmsf.max(),'Intra',a.intra_rmsf.min(),a.intra_rmsf.max())
print('Full Hchi range',a.loc[a.valid_chi_entropy,'pred_chi_entropy'].min(),a.loc[a.valid_chi_entropy,'pred_chi_entropy'].max(),'Intra',a.loc[a.valid_chi_entropy,'intra_chi_entropy'].min(),a.loc[a.valid_chi_entropy,'intra_chi_entropy'].max())
print('SASA',a.sidechain_SASA_complex_A2.min(),a.sidechain_SASA_complex_A2.max(),'shieldrange',a.sidechain_shielding_fraction.min(),a.sidechain_shielding_fraction.max())
# same-type position differences, all groups >=2; no cross-type inference
rep=a[a.valid_chi_entropy].groupby('amino_acid').filter(lambda x:len(x)>1).copy()
rep.to_csv(out/'all_same_type_sidechain_sites.csv',index=False)
summary=[]
for aa,g in rep.groupby('amino_acid'):
 r={'amino_acid':aa,'n':len(g),'positions':';'.join(g.CaMKI_native_position.astype(str)),'Full_Hchi_min':g.pred_chi_entropy.min(),'Full_Hchi_max':g.pred_chi_entropy.max(),'Intra_Hchi_min':g.intra_chi_entropy.min(),'Intra_Hchi_max':g.intra_chi_entropy.max(),'bound_sidechain_SASA_min':g.sidechain_SASA_complex_A2.min(),'bound_sidechain_SASA_max':g.sidechain_SASA_complex_A2.max(),'shield_fraction_min':g.sidechain_shielding_fraction.min(),'shield_fraction_max':g.sidechain_shielding_fraction.max()}
 summary.append(r)
print(pd.DataFrame(summary).to_string(index=False))
pd.DataFrame(summary).to_csv(out/'same_type_groups_summary.csv',index=False)