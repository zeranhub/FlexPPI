from pathlib import Path
import json, hashlib
import pandas as pd, numpy as np
from PIL import Image
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
F=R/'revision_v90_data_led_story/figure3'
ref=R/'revision_v89_mechanism_story/figure3/source_data/all25_dynamics_and_static_lobe_contacts.csv'
source=F/'source_data/all25_dynamics_and_static_lobe_contacts.csv'
d=pd.read_csv(ref).sort_values('CaMKI_native_position')
q=pd.read_csv(source).sort_values('CaMKI_native_position')
assert ref.read_bytes()==source.read_bytes()
assert len(d)==25 and ''.join(d.amino_acid)=='IKKNFAKSKWKQAFNATAVVRHMRK'
assert d.CaMKI_native_position.tolist()==list(range(294,319))
p=pd.read_csv(F/'source_data/all_plotted_dynamics.csv')
z=pd.read_csv(F/'source_data/repeated_Full_RMSF_expanded.csv')
expected=[]
for metric in ['rmsf','chi_entropy']:
 for model,prefix in [('Full','pred'),('Intra','intra')]:
  g=d[d['valid_'+metric]]
  expected.extend(dict(panel='b',metric=metric,position=int(r.CaMKI_native_position),model=model,value=float(r[prefix+'_'+metric])) for _,r in g.iterrows())
e=pd.DataFrame(expected);keys=['panel','metric','position','model']
assert len(e)==len(p)==92
assert not p.duplicated(keys).any()
m=p.merge(e,on=keys,how='outer',suffixes=('_plot','_frozen'),indicator=True,validate='one_to_one')
assert (m._merge=='both').all() and (m.value_plot==m.value_frozen).all()
assert len(z)==25 and z.position.tolist()==list(range(294,319))
assert set(z.model)=={'Full'} and set(z.metric)=={'rmsf'} and set(z.display)=={'expanded_y_scale'}
assert np.array_equal(z.value.to_numpy(),d.pred_rmsf.to_numpy())
assert np.array_equal(z.value.to_numpy(),p[(p.metric=='rmsf')&(p.model=='Full')].sort_values('position').value.to_numpy())
assert sorted(set(p.loc[p.metric=='chi_entropy','position']))==[x for x in range(294,319) if x not in [299,306,309,311]]
qa=json.loads((F/'figure_QA.json').read_text(encoding='utf-8'))
assert qa['panels']==['a','b'] and qa['no_mechanism_schematic']
assert qa['expanded_RMSF_ylim']==[0,2.3]
assert qa['original_axes'][0]['ylim']==[0,8.2] and qa['original_axes'][1]['ylim']==[0,4.3]
assert z.value.min()>=0 and z.value.max()<=2.3
assert qa['outside_text']==[]
sha=hashlib.sha256((F/'fig3_data_led_story.png').read_bytes()).hexdigest()
assert sha==qa['sha256']['png']
sc=(F/'compose_fig3_1MXE.py').read_text(encoding='utf-8')
assert "lo,hi=293.5,318.5" in sc
assert "a.set(xlim=(lo,hi),ylim=(0,2.3)" in sc
assert "df[prefix+'_'+metric].where(df['valid_'+metric])" in sc
out={'status':'pass','review_scope':'Actual final PNG viewed; independent frozen-source and all plotted-value audit; composer checked read-only. No root figure or data modified.','PNG_sha256':sha,'PNG_size_px':list(Image.open(F/'fig3_data_led_story.png').size),'canvas_mm':qa['canvas_mm'],'panels':['a','b'],'mechanism_schematic_removed':True,'frozen_reference':str(ref),'frozen_reference_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'source_byte_identical_to_V89':True,'unique_profile_values':92,'RMSF_unique_values':50,'Hchi_unique_values':42,'all92_exactly_equal_frozen':True,'max_numeric_difference':float((m.value_plot-m.value_frozen).abs().max()),'expanded_Full_RMSF_values':25,'all25_repeated_exactly_from_Full_profile':True,'total_curve_points_excluding_highlight_markers':117,'repeated_track_is_same_data_not_new_evidence':True,'native_xrange':[293.5,318.5],'native_positions':[294,318],'RMSF_overview_ylim':[0,8.2],'RMSF_expanded_ylim':[0,2.3],'Hchi_ylim':[0,4.3],'RMSF_Full_expanded_range':[float(z.value.min()),float(z.value.max())],'all_axes_start_at_zero':True,'no_axis_break_or_clipped_values':True,'omitted_invalid_Ala_Hchi':[299,306,309,311],'shade_span':[302.5,316.5],'shade_contains14_residues_303_316':True,'Full_expanded_scale_label_clear':True,'N297_N308_labels_correct':True,'K296_W303_M316_zoom_labels_correct':True,'palette_semantics_distinct':True,'no_need_for_separate_Lys_panel':True,'visual_findings':'The Full zoom makes the complete central low-fluctuation span readable without replacing the full-range Full/Intra comparison. Mechanistic inference belongs in the prose.','necessary_fixes':[],'min_font_from_renderer_pt':qa['min_font'],'outside_text_from_renderer':qa['outside_text']}
(F/'independent_final_QA.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(out,ensure_ascii=False,indent=2))