from pathlib import Path
import json,csv
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_model_screen_20261028')
p=pd.read_csv(R/'all_states_per_residue_predictions.csv',dtype={'chain':str,'author_residue_id':str})
number=pd.read_csv(R/'reference/antibody_residue_numbering.csv',dtype={'chain':str,'author_residue_id':str})
p=p.merge(number,on=['state_id','chain','author_residue_id','aa'],how='left',validate='many_to_one')
p['author_number']=p.author_residue_id.str.extract(r'(-?\d+)',expand=False).astype(int)
cons=json.loads((R/'contrasts.json').read_text(encoding='utf-8'))
all_deltas=[];summaries=[]
for con in cons:
 for v in ['full','intrachain']:
  b=p[p.state_id.eq(con['before'])&p.variant.eq(v)];a=p[p.state_id.eq(con['after'])&p.variant.eq(v)]
  d=b.merge(a,on=['chain','author_residue_id'],suffixes=('_before','_after'),validate='one_to_one')
  d['same_aa']=d.aa_before.eq(d.aa_after);d['contrast']=con['id'];d['variant']=v;d['case']=con['case'];d['comparison_type']=con['type'];d['counterfactual']=con['counterfactual']
  d['delta_rmsf']=d.rmsf_after-d.rmsf_before;d['delta_hchi']=d.hchi_after-d.hchi_before
  d['valid_rmsf']=d.model_mask_before&d.model_mask_after&d.geometry_valid_before&d.geometry_valid_after
  d['valid_hchi']=d.valid_rmsf&d.same_aa&d.hchi_chemical_valid_before&d.hchi_chemical_valid_after
  all_deltas.append(d)
  scopes={'all_common_positions':np.ones(len(d),bool),'all_unchanged_positions':d.same_aa.to_numpy()}
  for ch in sorted(d.chain.unique()):scopes['chain_'+ch]=d.chain.eq(ch).to_numpy()
  fv=d.imgt_number_before.notna()&d.imgt_number_after.notna()
  scopes['complete_Fv']=fv.to_numpy()
  for role in ['H','L']:
   scopes['complete_V'+role]=(fv&d.antibody_role_before.eq(role)&d.antibody_role_after.eq(role)).to_numpy()
   for region in ['CDR1','CDR2','CDR3']:
    scopes[role+'_'+region]=(d.antibody_role_before.eq(role)&d.antibody_role_after.eq(role)&d.imgt_region_before.eq(region)&d.imgt_region_after.eq(region)).to_numpy()
  scopes['antigen_loop_73_78']=(d.chain.eq('A')&d.author_number_before.between(73,78)).to_numpy()
  scopes['H_outerloop_68_75']=(d.chain.eq('H')&d.author_number_before.between(68,75)).to_numpy()
  for scope,mask in scopes.items():
   for metric in ['rmsf','hchi']:
    z=d[mask&d['valid_'+metric]]
    if z.empty:continue
    summaries.append({'contrast':con['id'],'case':con['case'],'variant':v,'comparison_type':con['type'],'counterfactual':con['counterfactual'],'scope':scope,'metric':metric,'n':len(z),'mean_before':z[metric+'_before'].mean(),'mean_after':z[metric+'_after'].mean(),'mean_delta':z['delta_'+metric].mean(),'median_delta':z['delta_'+metric].median(),'mean_abs_delta':z['delta_'+metric].abs().mean(),'max_abs_delta':z['delta_'+metric].abs().max(),'n_decrease':int((z['delta_'+metric]<0).sum()),'n_increase':int((z['delta_'+metric]>0).sum())})
deltas=pd.concat(all_deltas,ignore_index=True);reg=pd.DataFrame(summaries)
deltas.to_csv(R/'all_contrasts_complete_per_position.csv',index=False)
reg.to_csv(R/'all_contrasts_predefined_regions.csv',index=False)
# All explicitly named NMR residues; unresolved identities remain auditable.
nmr=pd.read_csv(R/'reference/all_literature_named_NMR_sites.csv',dtype={'PDB_1C08_author_id':str,'paper_chain':str})
site=[]
for name,con in [('free_R71A_minus_free_WT','apo_HR71A'),('free_R45A_minus_free_WT','apo_LR45A')]:
 for r in nmr[nmr.contrast.eq(name)].to_dict('records'):
  z=deltas[deltas.contrast.eq(con)&deltas.chain.eq(r['paper_chain'])&deltas.author_residue_id.eq(str(r['PDB_1C08_author_id']))]
  for row in z.to_dict('records'):
   site.append({**r,'model_contrast':con,'variant':row['variant'],'actual_aa':row['aa_before'],'exact_reference_identity_eligible':bool(r['identity_matches']) and row['aa_before']==str(r['paper_aa']),'rmsf_before':row['rmsf_before'],'rmsf_after':row['rmsf_after'],'delta_rmsf':row['delta_rmsf'],'hchi_before':row['hchi_before'],'hchi_after':row['hchi_after'],'delta_hchi':row['delta_hchi'],'hchi_valid':row['valid_hchi']})
pd.DataFrame(site).to_csv(R/'all_named_NMR_sites_prediction_check.csv',index=False)
# Static side-chain network: values displayed, no expected entropy direction inferred.
network=[('H','34','W'),('H','29','I'),('H','53','Y'),('A','101','D'),('A','62','W')]
out=[]
for con in ['flex_joint_state_change_fv','flex_4X_apo_to_bound_fv']:
 for ch,author,aa in network:
  z=deltas[deltas.contrast.eq(con)&deltas.chain.eq(ch)&deltas.author_residue_id.eq(author)]
  for row in z.to_dict('records'):
   out.append({'contrast':con,'chain':ch,'author':author,'expected_aa':aa,'aa_before':row['aa_before'],'aa_after':row['aa_after'],'variant':row['variant'],'rmsf_before':row['rmsf_before'],'rmsf_after':row['rmsf_after'],'delta_rmsf':row['delta_rmsf'],'hchi_before':row['hchi_before'],'hchi_after':row['hchi_after'],'delta_hchi':row['delta_hchi'],'hchi_valid':row['valid_hchi'],'reference_is_static_geometry_not_entropy':True})
pd.DataFrame(out).to_csv(R/'complete_static_sidechain_network_predictions.csv',index=False)
# Real MD label evaluation: exact source structures, full positions and chain scopes.
l=pd.read_csv(R/'md_reference/existing_labels_per_residue.csv',dtype={'chain':str,'author_residue_id':str})
actual=[];mdrows=[]
for sample in ['A00EJ','A00HE','A00KW']:
 for v in ['full','intrachain']:
  a=p[p.state_id.eq('mdref_'+sample)&p.variant.eq(v)].merge(l[l.sample_id.eq(sample)],on=['model_index','chain','author_residue_id','aa'],suffixes=('_pred','_MD'),validate='one_to_one')
  assert len(a)==len(l[l.sample_id.eq(sample)])
  mdrows.append(a)
  scopes={'all':np.ones(len(a),bool)}
  for role in a.role.unique():scopes['role_'+role]=a.role.eq(role).to_numpy()
  for role in ['H','L']:
   for region in ['CDR1','CDR2','CDR3']:
    scopes[role+'_'+region]=(a.antibody_role.eq(role)&a.imgt_region.eq(region)).to_numpy()
  for scope,scope_mask in scopes.items():
   for metric,col,target,flag in [('rmsf','rmsf_pred','rmsf_MD','rmsf_valid'),('hchi','hchi','chi_entropy','chi_entropy_valid')]:
    m=scope_mask&a[flag]&a.model_mask&np.isfinite(a[col])&np.isfinite(a[target]);z=a[m]
    if len(z)<2:continue
    x=z[col].to_numpy();y=z[target].to_numpy()
    xx=x-z.groupby('aa')[col].transform('mean').to_numpy();yy=y-z.groupby('aa')[target].transform('mean').to_numpy()
    rho=float(spearmanr(x,y).statistic) if np.std(x)>0 and np.std(y)>0 else None
    centered=float(spearmanr(xx,yy).statistic) if np.std(xx)>0 and np.std(yy)>0 else None
    actual.append({'sample_id':sample,'pdb_id':a.pdb_id.iloc[0],'variant':v,'scope':scope,'metric':metric,'n':len(z),'mae':float(np.mean(abs(x-y))),'mean_prediction':float(np.mean(x)),'mean_MD':float(np.mean(y)),'spearman':rho,'same_aa_centered_spearman':centered,'independent_local_trajectories':1,'label_windows':int(z.label_windows.iloc[0])})
pd.concat(mdrows,ignore_index=True).to_csv(R/'exact_input_predictions_and_MD_labels.csv',index=False)
accuracy=pd.DataFrame(actual);accuracy.to_csv(R/'MD_absolute_prediction_metrics.csv',index=False)
show=reg[reg.variant.eq('full')&reg.scope.isin(['complete_Fv','H_CDR3','H_CDR1','antigen_loop_73_78'])&reg.metric.eq('rmsf')]
print(show[['contrast','scope','n','mean_before','mean_after','mean_delta']].to_string(index=False))
print('\nMD full metrics')
print(accuracy[accuracy.variant.eq('full')&accuracy.scope.eq('all')].to_string(index=False))
print('\nNMR named free-state qualitative hypothesis, not identical observable:')
sites=pd.DataFrame(site)
for (c,v),g in sites[sites.exact_reference_identity_eligible].groupby(['model_contrast','variant']):
 print(c,v,len(g),'RMSF decrease',int((g.delta_rmsf<0).sum()),'mean delta',g.delta_rmsf.mean())
