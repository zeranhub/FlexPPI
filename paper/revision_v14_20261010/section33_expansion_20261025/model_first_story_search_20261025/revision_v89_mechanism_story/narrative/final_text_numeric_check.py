from pathlib import Path
import json, hashlib
import pandas as pd, numpy as np
from scipy.stats import spearmanr
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
V=R/'revision_v89_mechanism_story'
d=pd.read_csv(V/'figure3/source_data/all25_dynamics_and_static_lobe_contacts.csv')
f=pd.read_csv(R/'case_selection_audit_20261026/three_CaM_family_fold1_per_residue.csv');f=f[f.sample_id.eq('A001K')]
assert len(f)==169
p={}
for name,g in [('core',d[d.CaMKI_native_position.between(303,316)]),('outside',d[~d.CaMKI_native_position.between(303,316)])]:
 for src,col in [('Full','pred_rmsf'),('Intra','intra_rmsf'),('MD','target_rmsf')]:p[name+'_'+src+'_RMSF']=float(g[col].mean())
for native in [297,308]:
 g=d[d.CaMKI_native_position.eq(native)].iloc[0]
 for key in ['pred_rmsf','pred_chi_entropy','target_chi_entropy','sidechain_shielding_fraction']:p[f'N{native}_{key}']=float(g[key])
for name,g in [('all169',f),('peptide25',d)]:
 for metric in ['rmsf','chi_entropy']:
  for model,col in [('Full','pred_'+metric),('Intra','intra_'+metric)]:
   q=g[g['valid_'+metric]];p[f'{name}_{metric}_{model}_n']=len(q);p[f'{name}_{metric}_{model}_MAE']=float((q[col]-q['target_'+metric]).abs().mean());p[f'{name}_{metric}_{model}_rho']=float(spearmanr(q[col],q['target_'+metric]).statistic)
contact=pd.read_csv(V/'pymol_mechanism/all_context_min_distances.csv')
for (structure,site,scope),g in contact.groupby(['structure','site','scope']):
 if scope not in ['sidechain','sidechain_heavy']:continue
 z=g[g.min_distance_A<=4.5].sort_values('min_distance_A')
 p[f'contacts_{structure}_{site}_{scope}']=z[['residue','resname','min_distance_A']].to_dict('records')
p['contact_scopes']=list(contact.scope.unique())
p['figure_png_SHA256']=hashlib.sha256((V/'figure3/fig3_mechanism_story.png').read_bytes()).hexdigest()
(V/'narrative/final_text_numeric_check.json').write_text(json.dumps(p,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(p,ensure_ascii=False,indent=2))