"""Bootstrap genuine all-chain components, retaining PDB-equal point estimates."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
systems=pd.read_csv(OUT/'test53_per_system_metrics.csv')
inventory=pd.read_csv(OUT/'paired_inventory.csv').query("status=='evaluated'")
# Both target sides of each PDB must share one conservative sequence component.
assert inventory.groupby('pdb_id').family_components.nunique().max()==1
mapping=inventory[['pdb_id','family_components']].drop_duplicates()
assert not mapping.family_components.str.contains(';').any()
systems=systems.merge(mapping,on='pdb_id',validate='many_to_one')
systems.to_csv(OUT/'test53_per_system_with_family_component.csv',index=False)
metrics=['mae_model','mae_zero','mae_gain_over_zero','spearman','direction_all','direction_strong','top20_abs_overlap']
rng=np.random.default_rng(20261025);rows=[]
for cohort,sub in [('all_same_fold_oof',systems),('family_clean_under_original_oof',systems[systems.family_clean])]:
    for metric,g in sub.groupby('metric'):
        clusters=sorted(g.family_components.unique())
        sampled=rng.integers(0,len(clusters),(10000,len(clusters)))
        result=dict(cohort=cohort,metric=metric,pdb_systems=g.pdb_id.nunique(),family_components=len(clusters),
                    bootstrap_unit='all-chain conservative sequence-family component',replicates=10000)
        for name in metrics:
            # Ratio estimate: every sampled cluster contributes all its systems;
            # each system remains equally weighted in the resampled mean.
            grouped=g.groupby('family_components')[name].agg(['sum','count']).reindex(clusters)
            sums=grouped['sum'].to_numpy()[sampled].sum(axis=1)
            counts=grouped['count'].to_numpy()[sampled].sum(axis=1)
            boot=np.divide(sums,counts,out=np.full(len(sums),np.nan),where=counts>0)
            result[name]=float(g[name].mean())
            result[name+'_ci_low'],result[name+'_ci_high']=np.nanquantile(boot,[.025,.975])
        rows.append(result)
pd.DataFrame(rows).to_csv(OUT/'test53_family_component_bootstrap_summary.csv',index=False)
source=pd.read_csv(OUT/'test53_corrected_full_oof_residue_deltas.csv.gz')
source=source.merge(mapping,on='pdb_id',validate='many_to_one')
source.to_csv(OUT/'test53_corrected_full_oof_residue_deltas_with_family.csv.gz',index=False,compression='gzip')
(OUT/'family_component_bootstrap_protocol.json').write_text(json.dumps(dict(
    clusters=int(mapping.family_components.nunique()),seed=20261025,replicates=10000,
    estimator='PDB-equal macro mean; chain-side metrics averaged within4ETP before bootstrap; resample full all-chain components, carry all constituent PDB systems together; divide accumulated metric sum by accumulated valid system count.',
    caveat='Family-clean is a conservatively screened subset of original OOF; no family-model training or inference was performed.',
    pdb_bootstrap_earlier_file='test53_macro_summary.csv contains explicitly PDB-level sensitivity intervals; use this component-level summary for primary intervals.'),indent=2),encoding='utf-8')
print(pd.DataFrame(rows)[['cohort','metric','pdb_systems','family_components','mae_gain_over_zero','mae_gain_over_zero_ci_low','mae_gain_over_zero_ci_high']].to_string(index=False))
