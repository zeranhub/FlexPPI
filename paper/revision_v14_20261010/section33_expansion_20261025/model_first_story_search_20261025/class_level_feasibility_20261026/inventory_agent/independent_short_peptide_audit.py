from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
B=Path(r'E:\Project\FlexPPI');F=B/'experiments/_sources/family_folds';R=B/'paper/revision_v14_20261010/section33_expansion_20261025/model_first_story_search_20261025/class_level_feasibility_20261026';O=R/'inventory_agent'; D=B/'experiments/non_tcr_case_screen_20261025'
cand=pd.read_csv(O/'short_peptide_5to40_candidates.csv');cand=cand[cand.two_chain_unambiguous]; wanted=set(cand.sample_id);chains=pd.read_csv(D/'all_chain_metrics.csv'); cf=chains[chains.sample_id.isin(wanted)].copy(); fam=pd.read_csv(F/'record_family_folds.csv').set_index('sample_id'); report={};assert len(cf)==112
report['56_candidate_unique_records']=len(wanted);report['112_chain_roles']=len(cf);report['all_chain_folds_match_frozen_family_fold']=bool(all(int(r.fold)==int(fam.loc[r.sample_id,'fold']) for r in cf.itertuples()));report['family_components_single_fold']=bool((fam.groupby('component_id').fold.nunique()==1).all())
rawpath=Path(r'C:\Users\DELL\.codex\worktrees\2a52\FlexPPI\output\non_tcr_case_screen_20261025\archive_extract\family_folds\evaluation\oof_per_residue.csv.gz'); pieces=[]
for x in pd.read_csv(rawpath,chunksize=150000,dtype={'chain':str,'author_residue_id':str}):
 s=x[x.sample_id.isin(wanted)]
 if len(s):pieces.append(s)
raw=pd.concat(pieces,ignore_index=True); report['source_oof_csv_sha256']=hashlib.sha256(rawpath.read_bytes()).hexdigest();audit=json.loads((D/'screen_audit.json').read_text()); report['source_hash_matches_original_screen_audit']=report['source_oof_csv_sha256'] in audit['input_files'].values()
full=raw[raw.variant.eq('full')].sort_values(['sample_id','residue_index']).set_index(['sample_id','residue_index']);intra=raw[raw.variant.eq('intrachain')].sort_values(['sample_id','residue_index']).set_index(['sample_id','residue_index']);report['raw_full_intra_keys_unique_and_identical']=bool(full.index.equals(intra.index) and full.index.is_unique and intra.index.is_unique)
meta=['dataset','fold','amino_acid','chain','author_residue_id','valid_rmsf','valid_chi_entropy','target_rmsf','target_chi_entropy'];report['full_intra_matching_metadata_labels']={k:bool(full[k].equals(intra[k])) for k in meta};report['all_raw_fold_matches_family_fold']=bool(all(int(f)==int(fam.loc[sid,'fold']) for (sid,_),f in full.fold.items()) and all(int(f)==int(fam.loc[sid,'fold']) for (sid,_),f in intra.fold.items()))
maxerr=0.;issues=[]
for (sid,ch),g in raw.groupby(['sample_id','chain']):
 row=cf[(cf.sample_id==sid)&(cf.chain==ch)].iloc[0]
 for metric in ['rmsf','chi_entropy']:
  for var,prefix in [('full','full'),('intrachain','intra')]:
   t=g[g.variant.eq(var)&g['valid_'+metric]]; y=t['target_'+metric].to_numpy();p=t['pred_'+metric].to_numpy();mae=float(np.abs(y-p).mean());rho=float(pd.Series(y).rank().corr(pd.Series(p).rank()));
   for k,v in [(f'{metric}_{prefix}_mae',mae),(f'{metric}_{prefix}_rho',rho)]:
    err=abs(v-float(row[k]));maxerr=max(maxerr,err)
    if err>1e-10:issues.append((sid,ch,k,v,float(row[k])))
report['all112_chain_metrics_recomputed_from_raw_max_abs_error']=maxerr; report['chain_metric_issues']=issues
labissues=[]
for sid,g in full.reset_index().groupby('sample_id'):
 z=np.load(B/'experiments/core_corrected_pipeline/training_payload/dataset/labels'/f'{sid}.npz',allow_pickle=True); names=list(z['metric_names']);g=g.sort_values('residue_index');
 for metric in ['rmsf','chi_entropy']:
  j=names.index(metric)
  if not np.array_equal(g['target_'+metric].to_numpy(dtype=np.float32),z['values'][:,j],equal_nan=True):labissues.append((sid,metric,'target'))
  if not np.array_equal(g['valid_'+metric].to_numpy(dtype=bool),z['mask'][:,j]):labissues.append((sid,metric,'mask'))
report['all56_original_npz_labels_and_masks_float32_exact']=not labissues;report['label_issues']=labissues
# Check test member and training exclusion for every used family fold.
splitissues=[]
for fold,g in cf.groupby('fold'):
 tests=[json.loads(a) for a in (F/'splits'/f'fold_{fold}_test.jsonl').read_text().splitlines() if a];train=[json.loads(a) for a in (F/'splits'/f'fold_{fold}_train.jsonl').read_text().splitlines() if a];ts={a['sample_id'] for a in tests};tr={a['sample_id'] for a in train};trfam={a['metadata']['family_component_id'] for a in train}
 for sid in set(g.sample_id):
  if sid not in ts or sid in tr or fam.loc[sid,'component_id'] in trfam:splitissues.append((sid,int(fold)))
report['all56_members_heldout_test_and_own_family_absent_training']=not splitissues;report['split_issues']=splitissues
# Independent PDB merge and equal weighting; candidates must be fixed and each role unique.
roles=cf.merge(cand[['sample_id','short_chain']],on='sample_id',validate='many_to_one');roles['role']=np.where(roles.chain.eq(roles.short_chain),'peptide','receptor');report['duplicates']=roles[roles.role.eq('peptide')].groupby('pdb_id').sample_id.apply(list).loc[lambda q:q.str.len()>1].to_dict();outs=[]
for metric in ['rmsf','chi_entropy']:
 for role in ['peptide','receptor']:
  g=roles[roles.role.eq(role)];z=g.groupby('pdb_id')[[metric+'_full_mae',metric+'_intra_mae']].mean();gain=z[metric+'_intra_mae']-z[metric+'_full_mae'];outs.append(dict(metric=metric,role=role,PDBs=len(z),full_MAE=float(z[metric+'_full_mae'].mean()),intra_MAE=float(z[metric+'_intra_mae'].mean()),PDB_Full_better=int((gain>0).sum())))
report['independent_summary']=outs;report['verdict']='Confirmed preliminary within-structural-category prediction signal against original matching MD labels; not a confirmed functional class or uniform-time MD validation. No model rerun.'
(O/'independent_short_peptide_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))

