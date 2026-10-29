"""Compare frozen primary values with the previously used manuscript predictions."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser()
p.add_argument('--md-root',type=Path,required=True)
p.add_argument('--prior-csv',type=Path,required=True)
a=p.parse_args()
new=pd.read_csv(a.md_root/'frozen_predictions/residue_predictions.csv')
new=new[(new.method=='full') & (new.chain=='B')]
old=pd.read_csv(a.prior_csv)
assert set(old.model_fold)=={1}
rows=[]
for metric,column in [('rmsf','predicted_rmsf_angstrom'),('chi_entropy','predicted_hchi_nats')]:
    # Old prose may spell the metric hchi; use the observed CSV names only.
    if metric not in set(old.metric):
        alias={'chi_entropy':'hchi'}.get(metric,metric)
    else:alias=metric
    prior=old[old.metric==alias]
    if prior.empty:
        # The original KIX source table reports RMSF only. Do not invent an Hchi comparison.
        continue
    merged=prior.merge(new,on=['state','author_residue_number','amino_acid'],validate='one_to_one')
    assert len(merged)==len(prior)
    error=np.abs(merged.model_value-merged[column])
    rows.append(dict(metric=metric,positions=len(merged),maximum_absolute_difference=float(error.max())))
    if not np.allclose(merged.model_value,merged[column],rtol=2e-6,atol=2e-6):
        raise RuntimeError('Frozen input/checkpoint predictions differ from prior manuscript: '+str(rows[-1]))
assert any(row['metric']=='rmsf' and row['positions']==348 for row in rows)
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
output=a.md_root/'frozen_predictions/prior_comparison_audit.json'
output.write_text(json.dumps(dict(status='complete',comparison=rows,prior_sha256=sha(a.prior_csv),
   predictions_sha256=sha(a.md_root/'frozen_predictions/residue_predictions.csv')),indent=2),encoding='utf-8')
print(output)
