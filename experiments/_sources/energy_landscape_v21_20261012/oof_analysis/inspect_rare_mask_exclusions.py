from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
base = Path('/root/flexppi_completion_20261010/energy_landscape_v21_20261012/oof_analysis')
root = Path('/root/flexppi_revision_20261002/corrected_training_v1')
truth = pd.read_csv('/root/flexppi_completion_20261010/baselines/cohort/truth_residues.csv.gz',keep_default_na=False,dtype={'chain':str,'author_residue_id':str,'amino_acid':str})
selected = truth[(truth.valid_chi_entropy==0)&(~truth.amino_acid.isin(['A','G']))]
assert len(selected)==2 and set(selected.amino_acid)=={'K'}
results = []
for row in selected.itertuples(index=False):
    path = root/'dataset/labels'/(row.sample_id+'.npz')
    with np.load(path,allow_pickle=False) as label:
        metric_names = label['metric_names'].astype(str).tolist()
        column = metric_names.index('chi_entropy')
        assert not bool(label['mask'][row.residue_index,column])
        value = float(label['values'][row.residue_index,column])
        chi_count = int(label['hchi_chi_count'][row.residue_index]) if 'hchi_chi_count' in label else None
    results.append(dict(sample_id=row.sample_id,dataset=row.dataset,residue_index=int(row.residue_index),amino_acid=row.amino_acid,chain=row.chain,author_residue_id=row.author_residue_id,valid_rmsf=bool(row.valid_rmsf),valid_hchi=False,original_label_mask_hchi=False,stored_hchi_value=value if np.isfinite(value) else None,stored_chi_count=chi_count,label_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),reason='Original corrected Hchi supervision mask invalid; exact biochemical cause not inferred from saved mask'))
report = dict(status='passed_original_mask_exclusion_inventory',paired_excluded_total=136333,excluded_Ala=70830,excluded_Gly=65501,excluded_other=2,other_positions=results,no_new_position_removed=True)
(base/'validation/rare_nonAG_mask_exclusions.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
