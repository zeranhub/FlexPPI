"""Keep directional accuracy interpretable under unequal increase/decrease classes."""
from pathlib import Path
import numpy as np
import pandas as pd
OUT=Path(__file__).resolve().parent
data=pd.read_csv(OUT/'test53_corrected_full_oof_residue_deltas_with_family.csv.gz')
rows=[]
for (pair,pdb,metric,family,clean),g in data.groupby(['pair_id','pdb_id','metric','family_components','family_clean']):
    strong=g[g.md_delta.abs()>=.1]
    increase=strong.md_delta>0;decrease=strong.md_delta<0
    model=float((np.sign(strong.md_delta)==np.sign(strong.model_delta)).mean())
    fixed=float(decrease.mean())
    rows.append(dict(pair_id=pair,pdb_id=pdb,metric=metric,family_components=family,family_clean=clean,
                     strong_residues=len(strong),always_decrease_accuracy=fixed,model_accuracy=model,
                     advantage_over_always_decrease=model-fixed,
                     balanced_accuracy=float(np.mean([(strong.loc[increase,'model_delta']>0).mean(),
                                                      (strong.loc[decrease,'model_delta']<0).mean()])) if increase.any() and decrease.any() else np.nan))
table=pd.DataFrame(rows);table.to_csv(OUT/'test53_direction_controls_per_pair.csv',index=False)
columns=['strong_residues','always_decrease_accuracy','model_accuracy','advantage_over_always_decrease','balanced_accuracy']
systems=table.groupby(['pdb_id','metric','family_components','family_clean'])[columns].mean().reset_index()
systems.to_csv(OUT/'test53_direction_controls_per_system.csv',index=False)
summary=[];rng=np.random.default_rng(20261025)
for cohort,sub in [('all_same_fold_oof',systems),('family_clean_under_original_oof',systems[systems.family_clean])]:
    for metric,g in sub.groupby('metric'):
        r=dict(cohort=cohort,metric=metric,pdb_systems=g.pdb_id.nunique())
        r.update(g[columns].mean().to_dict())
        grouped=g.groupby('family_components').advantage_over_always_decrease.agg(['sum','count'])
        draw=rng.integers(0,len(grouped),(10000,len(grouped)))
        boot=grouped['sum'].to_numpy()[draw].sum(axis=1)/grouped['count'].to_numpy()[draw].sum(axis=1)
        r['advantage_cluster_ci_low'],r['advantage_cluster_ci_high']=np.quantile(boot,[.025,.975])
        summary.append(r)
pd.DataFrame(summary).to_csv(OUT/'test53_direction_controls_summary.csv',index=False)
print(pd.DataFrame(summary).to_string(index=False))
