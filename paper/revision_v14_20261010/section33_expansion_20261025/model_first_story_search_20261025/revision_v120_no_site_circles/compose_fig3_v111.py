# Academic Figure Skill Asset Confirmation (existing native figure components)
# a -> genuine frozen PyMOL ribbon PNGs -> native source image placement
# b -> v107 transparent MD/Full hexbin + marginals -> parameter inheritance
# c/d -> v109 dual metric curves + static sidechain SASA -> original native values
# e -> v110 four-chi horizontal scatter row -> original native values
# PAGE DESIGN derives from actual PaperBanana candidate_id0, final two-critic run.
from pathlib import Path
from hashlib import sha256
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import native_components as n

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'figure3';OUT.mkdir(exist_ok=True)
PB=ROOT.parent/'paperbanana_layout_v111_20261027/generation_retry_20261027_155023'
candidate=json.loads((PB/'candidate_2_pipeline.json').read_text(encoding='utf-8'))
assert candidate['candidate_id']==0
assert (PB/'candidate_2_final.png').is_file()
n.apply_style()
W,H=183.,190.;canvas=(W,H)
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')

# PaperBanana A, row1: adjacent structures + compact global joint distribution.
n.panel(fig,canvas,'a','',3,2)
n.text(fig,canvas,30,3,'CaM–CaMKI\n(1MXE)',ha='center',linespacing=1.12)
n.text(fig,canvas,81,3,'Ubiquitin–UBA\n(2QHO)',ha='center',linespacing=1.12)
n.molecule_pair(fig,canvas,(7,13,47,38),(57,13,49,38))
n.component_legend(fig,canvas,[(15,'Analysed chain',n.PURPLE),(65,'Partner',n.GRAY)],56)
n.panel(fig,canvas,'b','RMSF–Hχ relationship',114,2)
density_audit=n.draw_density(fig,canvas,(125,23,45,29),xmax=10,gridsize=(50,35),marginal_mm=4.)

# PaperBanana A, row2: equal-sized case profiles with matched SASA tracks.
rows=[]
n.method_legend(fig,canvas,94,70)
for letter,old_panel,case,left,title in [('c','b','CaMKI',14,'CaMKI dynamics'),('d','c','Ubiquitin',108,'Ubiquitin dynamics')]:
    n.panel(fig,canvas,letter,title,left-11,64)
    a,b,profile_values=n.draw_case(fig,canvas,(left,76,63,25),case,xlabel=False)
    n.annotate_case_sites(a,case)
    for (metric,series),z in profile_values.groupby(['metric','series'],sort=False):
        rows.extend(dict(panel=old_panel,display_panel=letter,kind='RMSF_A' if metric=='RMSF' else 'Hchi_nats',
                    identity=str(int(i)),value=float(v),series=series)
                    for i,v in zip(z.residue_id,z.value))
    a,sasa=n.draw_sasa(fig,canvas,(left,105,63,10),case)
    rows.extend(dict(panel=old_panel,display_panel=letter,kind='static_sidechain_SASA_A2',identity=str(int(i)),
                value=float(v),series='starting_structure') for i,v in zip(sasa.residue_id,sasa.sidechain_SASA_A2))

# PaperBanana A, row3: all four measured angle-occupancy comparisons in one row.
n.panel(fig,canvas,'e','Per-χ occupancy association',3,129)
n.chi_legend(fig,canvas,145,132)
_,chi_values=n.draw_4chi_row(fig,canvas,[14,57.5,101,144.5],146,35.5,29,140.5)
chi_values['display_panel']='e';rows.extend(chi_values.to_dict(orient='records'))
n.chi_shared_labels(fig,canvas,3.1,160.5,97,182)

plotted=pd.DataFrame(rows)
assert len(plotted)==767
old=pd.read_csv(ROOT/'frozen_all767_plotted_values.csv',float_precision='round_trip',dtype={'identity':str})
keys=['panel','kind','identity','series','chi']
for table in [old,plotted]:table['chi']=table.chi.fillna(-1).astype(int)
check=plotted.sort_values(keys).reset_index(drop=True);old=old.sort_values(keys).reset_index(drop=True)
pd.testing.assert_frame_equal(check[keys],old[keys],check_dtype=False)
assert np.array_equal(check.value.to_numpy(),old.value.to_numpy())
plotted.to_csv(OUT/'all767_values_logical_and_display_panels.csv',index=False)

fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
for t in fig.findobj(mpl.text.Text):
    if t.get_visible() and t.get_text():
        b=t.get_window_extent(renderer)
        if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_PaperBanana_layout_native.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa={'canvas_mm':[W,H],'font_pt':{'axes_ticks_legend':9,'titles':10,'bold_letters':11},
 'actual_PaperBanana_run':str(PB),'selected_candidate_id':0,'selected_saved_image':'candidate_2_final.png',
 'candidate_image_sha256':sha256((PB/'candidate_2_final.png').read_bytes()).hexdigest(),
 'layout':'row1 structures and compact global density; row2 equal case trace+SASA blocks; row3 all4chi horizontal',
 'display_panel_mapping':{'a':'structures','b':'global RMSF-Hchi density','c':'CaMKI curves+SASA','d':'Ubiquitin curves+SASA','e':'4chi occupancy'},
 'all767_scalar_values_exactly_preserved':True,'profile_values':368,'SASA_values':97,'occupancy_coordinates':302,
 'case_curve_box_mm':{'CaMKI':[14,76,63,25],'Ubiquitin':[108,76,63,25]},
 'case_SASA_box_mm':{'CaMKI':[14,105,63,10],'Ubiquitin':[108,105,63,10]},
 'density':density_audit,'outside_text':outside,'generated_scientific_pixels_reused':False,
 'real_structure_renders_preserved':True,'Word_modified':False,'new_model_or_MD_run':False,
 'visual_review':'pending','sha256':{ext:sha256((OUT/f'fig3_PaperBanana_layout_native.{ext}').read_bytes()).hexdigest() for ext in ['png','pdf','svg']}}
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
plt.close(fig)
print(json.dumps(qa,ensure_ascii=False))
