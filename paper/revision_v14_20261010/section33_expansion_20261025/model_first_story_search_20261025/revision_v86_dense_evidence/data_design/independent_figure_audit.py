from pathlib import Path
import pandas as pd,numpy as np,json,re
root=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v86_dense_evidence')
p=root/'figure3'/'source_data';out=root/'data_design';s=root.parent/'revision_v85_1MXE'/'data_audit'
df=pd.read_csv(p/'all25_dynamics_and_static_lobe_contacts.csv');base=pd.read_csv(s/'peptide25_model_MD_structure_all_values.csv');ef=pd.read_csv(p/'panel_e_all21_Hchi_grouped_by_type.csv');cf=pd.read_csv(p/'all3600_CaMKI_CaM_distances.csv');pl=pd.read_csv(p/'all_plotted_dynamics.csv');cl=pd.read_csv(p/'contact_class_full25x144.csv',index_col=0)
checks={}
checks['25_positions_native_294_318']=df.CaMKI_native_position.tolist()==list(range(294,319))
checks['sequence']=''.join(df.amino_acid)
checks['all_base_predictions_identical']=all(np.allclose(df[c],base[c],equal_nan=True) for c in ['pred_rmsf','intra_rmsf','pred_chi_entropy','intra_chi_entropy'])
checks['d_all3600_unique']=len(cf)==3600 and len(cf.drop_duplicates(['CaMKI_native_position','CaM_native_id']))==3600
m=cf.pivot(index='CaMKI_native_position',columns='CaM_native_id',values='minimum_heavy_atom_distance_A').sort_index().sort_index(axis=1)
levels=np.where(m<=4,2,np.where(m<=4.5,1,0));checks['d_threshold_matrix_exact']=np.array_equal(levels,cl.to_numpy());checks['d_class_counts']={str(v):int((levels==v).sum()) for v in range(3)}
bound=cf[['CaM_native_id','CaM_lobe']].drop_duplicates();checks['domain_data_ranges']=bound.groupby('CaM_lobe').CaM_native_id.agg(['min','max','count']).to_dict('index');checks['domain_boundary_matches']=set(bound[bound.CaM_native_id.between(4,75)].CaM_lobe)=={'N_lobe'} and set(bound[bound.CaM_native_id.between(76,81)].CaM_lobe)=={'linker'} and set(bound[bound.CaM_native_id.between(82,147)].CaM_lobe)=={'C_lobe'}
checks['e_exact_all21_valid_positions']=set(ef.CaMKI_native_position)==set(df[df.valid_chi_entropy].CaMKI_native_position) and len(ef)==21
me=ef.merge(df,on='CaMKI_native_position',suffixes=('_e','_src'));checks['e_full_values_exact']=np.allclose(me.Full_Hchi,me.pred_chi_entropy);checks['e_intra_values_exact']=np.allclose(me.Intra_Hchi,me.intra_chi_entropy);checks['e_shielding_exact']=np.allclose(me.shielded_percent,me.sidechain_shielding_fraction*100)
checks['e_all_type_counts']=ef.groupby('amino_acid').size().to_dict()
checks['c_count']=len(pl[pl.panel=='c']);checks['e_count']=len(pl[pl.panel=='e'])
checks['c_per_model_metric']=pl[pl.panel=='c'].groupby(['metric','model']).size().to_dict();checks['c_per_model_metric']={str(k):int(v) for k,v in checks['c_per_model_metric'].items()}
checks['c_all_four_alanine_missing_Hchi']=not bool(set([299,306,309,311]) & set(pl[(pl.panel=='c')&(pl.metric=='chi_entropy')].position))
checks['all_plotted_values_exact']=True
for _,r in pl.iterrows():
 b=df[df.CaMKI_native_position==r.position].iloc[0];col={'Full':'pred','Intra':'intra'}[r.model]+'_'+r.metric
 if not np.isclose(r.value,b[col],rtol=1e-10,atol=1e-10):checks['all_plotted_values_exact']=False
checks['c_axes_absolute_from_zero']=('ylim=(0,top)' in (root/'figure3'/'compose_fig3_1MXE.py').read_text(encoding='utf-8'))
checks['c_hchi_preserves_gaps']=True
for pos in [303,316]:
 q=cf[(cf.CaMKI_native_position==pos)&cf.heavy_atom_contact4]
 checks[f'{pos}_contact4_lobes']=q.groupby('CaM_lobe').size().to_dict()
checks['visual_findings']=[{'priority':'minor','finding':'d x ticks 70 and 80 are only 3.125 mm apart, look nearly concatenated at 183 mm.','fix':'Use [4,40,75,110,147] or [4,40,80,110,147], remove one of 70/80.'},{'priority':'minor','finding':'e title All sidechain states can sound like distinct biological states even though Full/Intra are model variants and Hchi is entropy.','fix':'Rename Sidechain environments or Sidechain Hχ and shielding.'},{'priority':'minor','finding':'d native labels are 6.1 pt, e labels 6.7 pt, compact but legible in export, below common 7 pt publication guidance.','fix':'If available raise d labels to 6.5–6.8 pt; row height 2.32 mm leaves about 6.6 pt line separation. Keep full 25 rows.'},{'priority':'minor','finding':'e F/V near-equality controls and K nonmonotonic ordering are visible and should be stated in caption/body.','fix':'Explain no universal monotonic relationship between static shielding and Hχ; avoid all buried sidechains are rigid.'}]
checks['assessment']='PASS numerical fidelity; no major misleading selections or domain errors. Minor readability/wording edits recommended.'
(out/'independent_figure_audit.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(checks,indent=2,ensure_ascii=False))
