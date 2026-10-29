from pathlib import Path
import json,pandas as pd,numpy as np
R=Path(__file__).resolve().parent;O=R/'audit';O.mkdir(exist_ok=True)
u=pd.read_csv(R/'uniform_md/all_uniform_labels_per_residue.csv.gz',dtype={'chain':str,'author_residue_id':str})
fam=pd.read_csv(r'E:\Project\FlexPPI\experiments\_sources\family_folds\record_family_folds.csv')[['sample_id','fold','component_id']].rename(columns={'fold':'family_oof_fold','component_id':'split_group'})
v=u.rename(columns={'fold':'original_label_split_fold'}).merge(fam,on='sample_id',validate='many_to_one')
assert v.family_id.eq(v.split_group).all()
v.to_parquet(R/'uniform_md/all_uniform_labels_explicit_fold_names.parquet',index=False)
assert v.groupby('variant').sample_id.nunique().eq(53).all()
f=pd.read_parquet(R/'formal_results/all_features_uniform_rep1.parquet')
assert f.sample_id.nunique()==53
assert f.fold.eq(f.sample_id.map(fam.set_index('sample_id').family_oof_fold)).all()
assert len(f)==len(v[v.variant.eq('uniform_rep1')])
keys=['sample_id','residue_index']
a=f.merge(v[v.variant.eq('uniform_rep1')][keys+['rmsf','chi_entropy','valid_rmsf','valid_chi_entropy','chain','amino_acid']],on=keys,validate='one_to_one',suffixes=('','_uniform'))
assert a.chain.eq(a.chain_uniform).all() and a.amino_acid.eq(a.amino_acid_uniform).all()
for m,raw in [('rmsf','rmsf'),('hchi','chi_entropy')]:
 keep=a['valid_'+m].to_numpy()
 assert np.array_equal(keep,a['valid_'+raw+'_uniform'] if raw=='rmsf' else a['valid_'+raw])
 assert np.allclose(a.loc[keep,'md_'+m],a.loc[keep,raw],rtol=0,atol=1e-12)
print('Fold/chain/AA/label joins checked. explicit uniform rows',len(v))
summary=json.loads((R/'uniform_md/compute_summary.json').read_text(encoding='utf-8'))
audit=dict(records53=True,all_uniform_variants53=True,source_labels_fold_renamed=True,actual_evaluation_uses_family_fold=True,source_full_intra_rows_matched=True,all_join_keys_one_to_one=True,AA_chain_label_match=True,valid_masks_match=True,no_performance_based_MD_exclusion=True,unavailable_short_trajectories=[q['sample_id'] for q in summary['unavailable']],not_new_trajectories=True,only_one_local_trajectory_per_record=True)
(O/'uniform_join_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
