# Native scientific figure, revised at the user's request from v111.
# Original predictions/MD/masks/structures are unchanged. Q62 is highlighted
# at its genuine data coordinates; Q2 remains in the unhighlighted population.
from pathlib import Path
from hashlib import sha256
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import native_components as n

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figure3';OUT.mkdir(exist_ok=True)
n.apply_style()
W,H=183.,168.;canvas=(W,H)
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')

n.panel(fig,canvas,'a','',3,2)
n.text(fig,canvas,26.5,3,'CaM–CaMKI\n(1MXE)',ha='center',linespacing=1.12)
n.text(fig,canvas,73,3,'Ubiquitin–UBA\n(2QHO)',ha='center',linespacing=1.12)
n.molecule_pair(fig,canvas,(4.5,12,44,37),(50,12,45,37))
n.component_legend(fig,canvas,[(11,'Analysed chain',n.PURPLE),(62,'Partner',n.GRAY)],53)
n.panel(fig,canvas,'b','RMSF–Hχ relationship',99,2)
density_audit=n.draw_density(fig,canvas,(111,15,59,36),xmax=10,gridsize=(50,35),marginal_mm=4.)

rows=[];case_axes={};sasa_axes={};curve_boxes={};sasa_boxes={}
n.method_legend(fig,canvas,146,63)
for letter,old_panel,case,left,width,title in [
    ('c','b','CaMKI',14,37.5,'CaMKI dynamics'),
    ('d','c','Ubiquitin',65,108,'Ubiquitin dynamics')]:
    n.panel(fig,canvas,letter,title,3 if case=='CaMKI' else 54,60)
    box=(left,70,width,26);sbox=(left,99,width,9)
    curve_boxes[case]=list(box);sasa_boxes[case]=list(sbox)
    a,b,profile_values=n.draw_case(fig,canvas,box,case,xlabel=False)
    case_axes[case]=(a,b)
    # Shared y-axis titles at the outer edges avoid a cramped central gutter.
    if case=='CaMKI':b.set_ylabel('')
    else:a.set_ylabel('')
    n.annotate_case_sites(a,case)
    for (metric,series),z in profile_values.groupby(['metric','series'],sort=False):
        rows.extend(dict(panel=old_panel,display_panel=letter,kind='RMSF_A' if metric=='RMSF' else 'Hchi_nats',
                    identity=str(int(i)),value=float(v),series=series)
                    for i,v in zip(z.residue_id,z.value))
    s,sasa=n.draw_sasa(fig,canvas,sbox,case);sasa_axes[case]=s
    s.sharex(a)
    if case=='Ubiquitin':s.set_ylabel('')
    rows.extend(dict(panel=old_panel,display_panel=letter,kind='static_sidechain_SASA_A2',identity=str(int(i)),
                value=float(v),series='starting_structure') for i,v in zip(sasa.residue_id,sasa.sidechain_SASA_A2))

n.panel(fig,canvas,'e','Per-χ occupancy association',3,119)
n.chi_legend(fig,canvas,146,122)
_,chi_values=n.draw_4chi_row(fig,canvas,[14,57.5,101,144.5],134.5,35.5,23,129)
chi_values['display_panel']='e';rows.extend(chi_values.to_dict(orient='records'))
n.chi_shared_labels(fig,canvas,3.1,146,97,163)

plotted=pd.DataFrame(rows)
assert len(plotted)==767
old=pd.read_csv(ROOT/'frozen_all767_plotted_values.csv',float_precision='round_trip',dtype={'identity':str})
keys=['panel','kind','identity','series','chi']
for table in [old,plotted]:table['chi']=table.chi.fillna(-1).astype(int)
check=plotted.sort_values(keys).reset_index(drop=True);old=old.sort_values(keys).reset_index(drop=True)
pd.testing.assert_frame_equal(check[keys],old[keys],check_dtype=False)
assert np.array_equal(check.value.to_numpy(),old.value.to_numpy())
plotted.to_csv(OUT/'all767_values_logical_and_display_panels.csv',index=False)

# Check actual rendered coordinate transforms, not just declared widths.
fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[];align={}
for case,(a,b) in case_axes.items():
    s=sasa_axes[case];lo,hi=a.get_xlim()
    assert np.array_equal(a.get_xlim(),s.get_xlim())
    samples=np.linspace(lo,hi,101)
    dx=np.max(np.abs(a.transData.transform(np.c_[samples,np.zeros(101)])[:,0]-
                     s.transData.transform(np.c_[samples,np.zeros(101)])[:,0]))
    assert dx<1e-8
    pixels_per_residue=(a.transData.transform((lo+1,0))-a.transData.transform((lo,0)))[0]
    align[case]={'maximum_curve_SASA_x_error_pixels':float(dx),'pixels_per_residue':float(pixels_per_residue),
                 'mm_per_residue':curve_boxes[case][2]/(hi-lo),'xlim':[lo,hi]}
assert abs(align['CaMKI']['mm_per_residue']-align['Ubiquitin']['mm_per_residue'])<1e-12
for t in fig.findobj(mpl.text.Text):
    if t.get_visible() and t.get_text():
        bb=t.get_window_extent(renderer)
        if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_compact_equal_residue_spacing.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa={'canvas_mm':[W,H],'previous_canvas_mm':[183,190],'height_reduction_percent':100*(190-H)/190,
 'font_pt':{'axes_ticks':9,'titles':10,'bold_letters':11,'density_legend':8},
 'layout':'larger density upper right; case width 25:72 with aligned shared x SASA; four chi horizontal',
 'display_panel_mapping':{'a':'structures','b':'global RMSF-Hchi density','c':'CaMKI curves+SASA','d':'Ubiquitin curves+SASA','e':'4chi occupancy'},
 'all767_scalar_values_exactly_preserved':True,'profile_values':368,'SASA_values':97,'occupancy_coordinates':302,
 'highlight_change':{'previous':[41,2],'current':[41,62],'Q2_retained_as_other':True,'Q62_actual_coordinates':True},
 'case_curve_box_mm':curve_boxes,'case_SASA_box_mm':sasa_boxes,'alignment':align,
 'density_area_change_percent':100*(59*36/(45*29)-1),'density':density_audit,'outside_text':outside,
 'real_structure_renders_preserved':True,'new_model_or_MD_run':False,'visual_review':'pending',
 'sha256':{ext:sha256((OUT/f'fig3_compact_equal_residue_spacing.{ext}').read_bytes()).hexdigest() for ext in ['png','pdf','svg']}}
(OUT/'figure_QA_v112.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
plt.close(fig)
print(json.dumps(qa,ensure_ascii=False))
