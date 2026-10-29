import pandas as pd
from pathlib import Path
p=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v88_two_panel_story')
df=pd.read_csv(p/'figure3/source_data/all25_dynamics_and_static_lobe_contacts.csv')
cols=['CaMKI_native_position','amino_acid','target_rmsf','pred_rmsf','intra_rmsf','target_chi_entropy','pred_chi_entropy','intra_chi_entropy','sidechain_shielding_fraction','minimum_heavy_atom_distance_A']
print(df[cols].to_string(index=False))
w=pd.read_csv(p/'data_audit/MD_rmsf_five_windows_25_residues.csv')
for a,b in [(294,302),(303,316),(317,318),(294,296),(299,315),(316,318)]:
 s=df[df.CaMKI_native_position.between(a,b)]
 print(a,b,'mean',s[['target_rmsf','pred_rmsf','intra_rmsf']].mean().to_dict(),'median',s[['target_rmsf','pred_rmsf','intra_rmsf']].median().to_dict())
 print('window',w[w.CaMKI_native_position.between(a,b)].groupby('block').rmsf.mean().to_dict())