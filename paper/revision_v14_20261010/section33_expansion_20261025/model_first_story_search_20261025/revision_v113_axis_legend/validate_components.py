"""Validate reusable components on an unsaved test canvas; no new figure release."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import native_components as native

ROOT = Path(__file__).resolve().parent
fig = plt.figure(figsize=(10,10))
canvas = (254.,254.)
rows = []
for case, panel in [('CaMKI','b'),('Ubiquitin','c')]:
    a, b, profiles = native.draw_case(fig, canvas, (20,20,100,30), case)
    for (metric, series), z in profiles.groupby(['metric','series'], sort=False):
        rows.extend(dict(panel=panel,kind='RMSF_A' if metric=='RMSF' else 'Hchi_nats',identity=str(int(i)),value=float(v),series=series)
                    for i,v in zip(z.residue_id,z.value))
    a, sasa = native.draw_sasa(fig, canvas, (20,60,100,12), case)
    rows.extend(dict(panel=panel,kind='static_sidechain_SASA_A2',identity=str(int(i)),value=float(v),series='starting_structure')
                for i,v in zip(sasa.residue_id,sasa.sidechain_SASA_A2))
_, scatter_values = native.draw_4chi_row(fig, canvas, [20,70,120,170], 130, 35.5, 25, 120)
rows.extend(scatter_values.to_dict(orient='records'))
old = pd.read_csv(ROOT/'frozen_all767_plotted_values.csv', float_precision='round_trip', dtype={'identity':str})
new = pd.DataFrame(rows)
keys = ['panel','kind','identity','series','chi']
new['chi'] = new['chi'].fillna(-1).astype(int)
old['chi'] = old['chi'].fillna(-1).astype(int)
new = new.sort_values(keys).reset_index(drop=True)
old = old.sort_values(keys).reset_index(drop=True)
assert len(new) == len(old) == 767
pd.testing.assert_frame_equal(new[keys], old[keys], check_dtype=False)
assert np.array_equal(new.value.to_numpy(), old.value.to_numpy())
native.molecule_pair(fig, canvas, (20,180,40,40), (70,180,40,40))
from plot_overlay import load_data
arrays, stats = load_data()
assert all(len(v)==824988 for v in arrays.values())
plt.close(fig)
report = {'all_767_native_values_exactly_equal_to_v110':True,'source_and_helper_imports_valid':True,
          'original_molecule_images_and_landmarks_load':True,'density_arrays_valid':True,
          'unsaved_test_canvas_only':True,'new_layout_or_artifact_rendered':False}
(ROOT/'native_component_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report))
