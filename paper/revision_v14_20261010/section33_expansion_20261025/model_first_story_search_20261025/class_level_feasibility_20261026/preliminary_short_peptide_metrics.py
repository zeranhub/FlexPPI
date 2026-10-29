from pathlib import Path
import pandas as pd, numpy as np, json,hashlib
ROOT=Path(__file__).resolve().parent
P=ROOT/'inventory_agent'
out=ROOT/'short_peptide_preliminary'
out.mkdir(exist_ok=True)
cand=pd.read_csv(P/'short_peptide_5to40_candidates.csv')
cand=cand[cand.two_chain_unambiguous].copy()
assert cand.sample_id.is_unique and len(cand)==56
chains=pd.read_csv(P/'all1535_chains_with_sequence_contacts_identity_family.csv')
assert not chains.duplicated(['sample_id','chain']).any()
selected=[]
for q in cand.itertuples():
    for chain,role in [(q.short_chain,'peptide'),(q.receptor_chains,'receptor')]:
        g=chains[(chains.sample_id==q.sample_id)&(chains.chain==chain)]
        assert len(g)==1,(q.sample_id,chain)
        g=g.copy();g['role']=role;selected.append(g)
x=pd.concat(selected,ignore_index=True)
assert len(x)==112
x.to_csv(out/'all56_records_two_roles.csv',index=False)
results=[]
for metric in ['rmsf','chi_entropy']:
    cols=[metric+'_full_mae',metric+'_intra_mae',metric+'_full_rho',metric+'_intra_rho']
    for role,g in x.groupby('role'):
        # Collapse duplicate PDB records first; do not let duplicate conditions count as independent examples.
        by=g.groupby(['pdb_id','component_id'],as_index=False)[cols].mean()
        by['abs_MAE_gain']=by[metric+'_intra_mae']-by[metric+'_full_mae']
        by['fractional_MAE_gain']=by.abs_MAE_gain/by[metric+'_intra_mae']
        by['rho_gain']=by[metric+'_full_rho']-by[metric+'_intra_rho']
        fm=by.groupby('component_id')[['abs_MAE_gain','fractional_MAE_gain','rho_gain']].mean()
        results.append(dict(metric=metric,role=role,n_records=len(g),n_unique_PDBs=by.pdb_id.nunique(),n_split_groups=by.component_id.nunique(),PDB_equal_Full_MAE=float(by[metric+'_full_mae'].mean()),PDB_equal_Intra_MAE=float(by[metric+'_intra_mae'].mean()),n_PDBs_Full_lower_MAE=int((by.abs_MAE_gain>0).sum()),PDB_fraction_Full_lower_MAE=float((by.abs_MAE_gain>0).mean()),PDB_equal_fractional_MAE_gain=float(by.fractional_MAE_gain.mean()),split_group_equal_fractional_MAE_gain=float(fm.fractional_MAE_gain.mean()),n_split_groups_mean_gain_positive=int((fm.abs_MAE_gain>0).sum()),PDB_equal_Full_rho=float(by[metric+'_full_rho'].mean()),PDB_equal_Intra_rho=float(by[metric+'_intra_rho'].mean()),PDB_equal_rho_gain=float(by.rho_gain.mean())))
        by.assign(metric=metric,role=role).to_csv(out/(metric+'_'+role+'_per_PDB.csv'),index=False)
summary=dict(selector='Membership fixed before looking at model metrics: short chain 5-40 aa, partner >=70 aa, distinct identity, both cross-chain CA8 contacts and valid RMSF/Hchi, two protein chains only.',data_scope='Existing family OOF predictions against corresponding ORIGINAL MD labels. Not uniform MD time windows; not an experiment and not a binding before/after comparison.',n_records=56,n_unique_PDBs=int(x.pdb_id.nunique()),n_split_groups=int(x.component_id.nunique()),no_new_model_run=True,results=results,caveats=['Biological eligibility of each short chain not manually confirmed. It may be a truncated construct.','Family connected components are split units, not biological class names.','No model-quality-based sample inclusion. This does not yet establish the biological pattern or benefit beyond burial/type/secondary-structure baselines.','Duplicates collapsed by PDB; remaining MD protocol, temperature and observation-time differences not controlled.','No inferential test or confidence intervals used in this feasibility pass.'])
(out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
