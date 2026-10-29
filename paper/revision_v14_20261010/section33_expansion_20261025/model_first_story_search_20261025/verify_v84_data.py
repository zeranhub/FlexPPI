from pathlib import Path
import pandas as pd,numpy as np,hashlib
root=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
out=root/'revision_v84_effect_audit'
f='cam_md_state_contrast_all_common_residues.csv'
assert (out/'source_data'/f).read_bytes()==(root/'revision_v83_model_story'/'figure3'/'source_data'/f).read_bytes()
df=pd.read_csv(out/'source_data'/f);pl=pd.read_csv(out/'source_data'/'plotted_values.csv');full=pl[pl.panel=='full']
assert len(full)==534
for state,suffix in [('NHE1','a'),('IQCG','b')]:
 for metric in ['rmsf','hchi']:
  source=df[df.metric==metric].set_index('mature_cam_number')['full0_'+suffix].sort_index()
  plotted=full[(full.state==state)&(full.metric==metric)].set_index('position').value.sort_index()
  assert np.array_equal(source.index.to_numpy(),plotted.index.to_numpy())
  assert np.allclose(source.to_numpy(),plotted.to_numpy(),rtol=0,atol=1e-12)
print('PASS: all 534 full-sequence plotted values match the unchanged v83 source; absolute axes retain the original zero baseline and range.')
print(pd.read_csv(out/'source_data'/'regional_effects_full_and_reference.csv')[['first','last','metric','n','full_NHE1','full_IQCG']].to_string(index=False))