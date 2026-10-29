from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
source=Path(r'E:\Project\FlexPPI\experiments\non_tcr_case_screen_20261025')
profiles=Path(r'C:\Users\DELL\.codex\worktrees\2a52\FlexPPI\output\non_tcr_case_screen_20261025\all_non_tcr_profiles.csv.gz')
out=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\case_selection_audit_20261026');out.mkdir(exist_ok=True)
selected=[]
for chunk in pd.read_csv(profiles,chunksize=100000):
 g=chunk[chunk.sample_id.isin(['A001K','A004C','A007U'])]
 if len(g):selected.append(g)
d=pd.concat(selected,ignore_index=True);d.to_csv(out/'three_CaM_family_fold1_per_residue.csv',index=False)
assert set(d.fold)=={1}
records={}
for file in ['all_616_record_screen.csv','all_chain_metrics.csv','all_region_metrics.csv']:
 t=pd.read_csv(source/file);t=t[t.pdb_id.isin(['1MXE','2YGG','4LZX'])];t.to_csv(out/file,index=False);records[file]=t
q=d[(d.sample_id=='A001K')&(d.chain=='B')].copy();q['peptide_position']=np.arange(1,len(q)+1);q['CaMKI_native_position']=np.arange(294,319)
assert ''.join(q.amino_acid)=='IKKNFAKSKWKQAFNATAVVRHMRK'
assert q.loc[q.CaMKI_native_position==303,'amino_acid'].iloc[0]=='W'
assert q.loc[q.CaMKI_native_position==316,'amino_acid'].iloc[0]=='M'
q.to_csv(out/'CaMKI_294_318_family_fold1.csv',index=False)
checks=[]
for metric in ['rmsf','chi_entropy']:
 good=q['valid_'+metric].astype(bool)
 a=q.loc[good,'target_'+metric].to_numpy();full=q.loc[good,'pred_'+metric].to_numpy();intra=q.loc[good,'intra_'+metric].to_numpy()
 row=records['all_chain_metrics.csv'];row=row[(row.pdb_id=='1MXE')&(row.chain=='B')].iloc[0]
 for model,pred in [('full',full),('intra',intra)]:
  mae=float(np.mean(abs(a-pred)));assert abs(mae-row[metric+'_'+model+'_mae'])<1e-10
  checks.append(dict(metric=metric,model=model,n=int(good.sum()),target_mean=float(a.mean()),prediction_mean=float(pred.mean()),mae=mae))
old=[]
for variant in ['full','intrachain']:
 p=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\neural_exports\evaluation')/variant/'A001K.npz';z=np.load(p)
 pred=z['mean'][-25:,0];target=z['target'][-25:,0]
 old.append(dict(protocol='original_record_split_OOF_fold4_not_family_fold1',variant=variant,peptide_RMSF_mean=float(pred.mean()),peptide_RMSF_MAE=float(abs(pred-target).mean()),source=str(p)))
pd.DataFrame(checks).to_csv(out/'CaMKI_peptide_family_fold1_checked_metrics.csv',index=False)
audit=dict(decision='prefer_1MXE_for_primary_complex_context_case_not_proven_binding_induced_transition',family_fold1_checked=checks,screenshot_paragraph_matches_different_protocol=old,no_matching_isolated_1MXE_pair_in_existing_pair_delta_table=not pd.read_csv(source/'pair_delta_metrics.csv').pdb.eq('1MXE').any(),numbering='Native CaMKI 294-318 mapped by exact 25-residue sequence and the original structural paper; local preprocessed PDB numbering is different.',manuscript_modified=False,source_hashes={str(source/file):hashlib.sha256((source/file).read_bytes()).hexdigest() for file in ['all_616_record_screen.csv','all_chain_metrics.csv','all_region_metrics.csv']})
(out/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(pd.DataFrame(checks).to_string(index=False));print('OLD PARAGRAPH SOURCE',json.dumps(old));print(q[q.CaMKI_native_position.isin([296,298,303,316])][['CaMKI_native_position','amino_acid','target_rmsf','pred_rmsf','intra_rmsf','target_chi_entropy','pred_chi_entropy','intra_chi_entropy']].to_string(index=False));print('Saved',out)