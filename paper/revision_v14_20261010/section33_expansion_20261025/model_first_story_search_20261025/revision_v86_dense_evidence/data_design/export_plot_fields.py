from pathlib import Path
import pandas as pd
p=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v86_dense_evidence\data_design')
a=pd.read_csv(p/'all25_dynamics_and_static_lobe_contacts.csv')
a['site_label']=a.amino_acid+a.CaMKI_native_position.astype(str)
a['Full_RMSF']=a.pred_rmsf;a['Intra_RMSF']=a.intra_rmsf
cols=['site_label','CaMKI_native_position','amino_acid','region','anchor_span_303_316','Full_RMSF','Intra_RMSF','sidechain_SASA_complex_A2','sidechain_shielding_fraction','N_lobe_count_4A','C_lobe_count_4A','linker_count_4A']
a[cols].to_csv(p/'panel_d_all25_rmsf_paired.csv',index=False)
v=a[a.valid_chi_entropy].copy(); rep_order=['F','K','N','R','V']; singles=['I','S','W','Q','T','H','M']; order={aa:i for i,aa in enumerate(rep_order+singles)}
v['group_order']=v.amino_acid.map(order);v['Full_Hchi']=v.pred_chi_entropy;v['Intra_Hchi']=v.intra_chi_entropy;v['shielded_percent']=v.sidechain_shielding_fraction*100
v=v.sort_values(['group_order','CaMKI_native_position'])
assert len(v)==21
v[['site_label','CaMKI_native_position','amino_acid','group_order','Full_Hchi','Intra_Hchi','sidechain_SASA_complex_A2','shielded_percent','N_lobe_count_4A','C_lobe_count_4A','linker_count_4A']].to_csv(p/'panel_e_all21_Hchi_grouped_by_type.csv',index=False)
d=pd.read_csv(p/'all3600_CaMKI_CaM_distances.csv'); wide=d.pivot(index='CaMKI_native_position',columns='CaM_native_id',values='minimum_heavy_atom_distance_A')
wide.to_csv(p/'panel_f_all25_by144_distance_matrix_A.csv')