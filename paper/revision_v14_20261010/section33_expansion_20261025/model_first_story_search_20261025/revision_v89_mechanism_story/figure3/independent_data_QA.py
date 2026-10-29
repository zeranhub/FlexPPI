from pathlib import Path
import pandas as pd, numpy as np, hashlib, json
ROOT=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
F=ROOT/'revision_v89_mechanism_story/figure3'
ref=ROOT/'revision_v88_two_panel_story/figure3/source_data/all25_dynamics_and_static_lobe_contacts.csv'
local=F/'source_data/all25_dynamics_and_static_lobe_contacts.csv'
points=F/'source_data/all_plotted_dynamics.csv'
d=pd.read_csv(ref).sort_values('CaMKI_native_position')
q=pd.read_csv(local).sort_values('CaMKI_native_position')
p=pd.read_csv(points)
assert len(d)==len(q)==25
assert ''.join(d.amino_acid)=='IKKNFAKSKWKQAFNATAVVRHMRK'
assert (d.CaMKI_native_position.to_numpy()==q.CaMKI_native_position.to_numpy()).all()
assert set(d.fold)=={1} and set(q.fold)=={1}
columns=['pred_rmsf','intra_rmsf','target_rmsf','pred_chi_entropy','intra_chi_entropy','target_chi_entropy','valid_rmsf','valid_chi_entropy']
for c in columns:
 assert np.array_equal(d[c].to_numpy(),q[c].to_numpy(),equal_nan=True), c
expected=[]
for metric in ['rmsf','chi_entropy']:
 for model,prefix in [('Full','pred'),('Intra','intra')]:
  g=d[d['valid_'+metric]]
  expected.extend({'panel':'b','metric':metric,'position':int(r.CaMKI_native_position),'model':model,'value':float(r[prefix+'_'+metric])} for _,r in g.iterrows())
e=pd.DataFrame(expected)
keys=['panel','metric','position','model']
assert len(e)==len(p)==92
assert not p.duplicated(keys).any()
m=p.merge(e,on=keys,how='outer',suffixes=('_plot','_frozen'),indicator=True,validate='one_to_one')
assert (m._merge=='both').all()
assert (m.value_plot==m.value_frozen).all()
assert len(p[p.metric=='rmsf'])==50 and len(p[p.metric=='chi_entropy'])==42
assert sorted(set(p.loc[p.metric=='chi_entropy','position']))==[x for x in range(294,319) if x not in [299,306,309,311]]
assert p.value.min()>=0
out={'reference':str(ref),'reference_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'local_source_sha256':hashlib.sha256(local.read_bytes()).hexdigest(),'local_source_byte_identical_to_frozen':ref.read_bytes()==local.read_bytes(),'all_plotted_points':92,'RMSF_values':50,'Hchi_values':42,'unique_RMSF_positions':25,'unique_Hchi_positions':21,'Hchi_omitted_Ala':[299,306,309,311],'all92_match_frozen_exactly':True,'max_absolute_numeric_difference':float((m.value_plot-m.value_frozen).abs().max()),'all_nonnegative':True,'fold':1,'PNG_sha256_reviewed':hashlib.sha256((F/'fig3_mechanism_story.png').read_bytes()).hexdigest()}
(F/'independent_data_QA.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(out,ensure_ascii=False,indent=2))