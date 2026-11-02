# Native scientific figure, revised at the user's request from v122.
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
W,H=183.,155.;canvas=(W,H)
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')

n.panel(fig,canvas,'a','',3,2)
n.text(fig,canvas,26.5,3,'CaM–CaMKI\n(1MXE)',ha='center',linespacing=1.12)
n.text(fig,canvas,73,3,'Ubiquitin–UBA\n(2QHO)',ha='center',linespacing=1.12)
n.molecule_pair(fig,canvas,(4.5,12,44,37),(50,12,45,37))
n.component_legend(fig,canvas,[(11,'Analysed chain',n.BLUE),(62,'Partner',n.GRAY)],53)
n.panel(fig,canvas,'b','RMSF–Hχ relationship',99,2)
density_audit=n.draw_density(fig,canvas,(111,15,59,36),xmax=10,gridsize=(50,35),marginal_mm=4.)

rows=[];case_axes={};sasa_axes={};curve_boxes={};sasa_boxes={}
n.method_legend(fig,canvas,146,63)
n.panel(fig,canvas,'c','',3,57)
for letter,old_panel,case,left,width,title in [
    ('c','b','CaMKI',14,46,'CaMKI dynamics'),
    ('c','c','Ubiquitin',69,104,'Ubiquitin dynamics')]:
    n.text(fig,canvas,9 if case=='CaMKI' else 69,57.25,title,10)
    box=(left,67,width,25);sbox=(left,94.5,width,8)
    curve_boxes[case]=list(box);sasa_boxes[case]=list(sbox)
    a,b,profile_values=n.draw_case(fig,canvas,box,case,xlabel=False)
    case_axes[case]=(a,b)
    # Shared y-axis titles at the outer edges avoid a cramped central gutter.
    if case=='CaMKI':
        b.set_ylabel('')
        a.yaxis.label.set(ha='center',va='center')
        a.yaxis.set_label_coords(3.6/W,1-(box[1]+box[3]/2)/H,transform=fig.transFigure)
    else:a.set_ylabel('')
    n.annotate_case_sites(a,case)
    for (metric,series),z in profile_values.groupby(['metric','series'],sort=False):
        rows.extend(dict(panel=old_panel,display_panel=letter,kind='RMSF_A' if metric=='RMSF' else 'Hchi_nats',
                    identity=str(int(i)),value=float(v),series=series)
                    for i,v in zip(z.residue_id,z.value))
    s,sasa=n.draw_sasa(fig,canvas,sbox,case);sasa_axes[case]=s
    s.sharex(a)
    if case=='Ubiquitin':s.set_ylabel('')
    else:
        s.yaxis.label.set(ha='center',va='center')
        s.yaxis.set_label_coords(3.6/W,1-(sbox[1]+sbox[3]/2)/H,transform=fig.transFigure)
    rows.extend(dict(panel=old_panel,display_panel=letter,kind='static_sidechain_SASA_A2',identity=str(int(i)),
                value=float(v),series='starting_structure') for i,v in zip(sasa.residue_id,sasa.sidechain_SASA_A2))

n.panel(fig,canvas,'d','Per-χ occupancy association',3,113)
n.chi_legend(fig,canvas,146,114.5)
_,chi_values=n.draw_4chi_row(fig,canvas,[14,57.5,101,144.5],122.5,35.5,21,118)
chi_values['display_panel']='d';rows.extend(chi_values.to_dict(orient='records'))
chi_labels=n.chi_shared_labels(fig,canvas,3.6,133,97,150.5)

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
# Axis titles were previously misplaced by automatic tick-label padding.
# Check actual visible text centres as well as the unchanged data transforms.
left_labels=[case_axes['CaMKI'][0].yaxis.label,sasa_axes['CaMKI'].yaxis.label,chi_labels[0]]
label_centres={t.get_text():float(t.get_window_extent(renderer).x0+t.get_window_extent(renderer).width/2)
               for t in left_labels}
assert max(label_centres.values())-min(label_centres.values())<1e-8,label_centres
for t in fig.findobj(mpl.text.Text):
    if t.get_visible() and t.get_text():
        bb=t.get_window_extent(renderer)
        if bb.x0<0 or bb.y0<0 or bb.x1>fig.bbox.width or bb.y1>fig.bbox.height:outside.append(t.get_text())
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_orange_sites_black_sasa.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa={'canvas_mm':[W,H],'previous_canvas_mm':[183,168],'height_reduction_percent':100*(168-H)/168,
 'font_pt':{'axes_ticks':9,'titles':10,'bold_letters':11,'density_legend':8},
 'layout':'compact four-panel composition; shorter row gaps and tighter four-chi row; matched left ylabel centres',
 'display_panel_mapping':{'a':'structures','b':'global RMSF-Hchi density','c':'CaMKI left and Ubiquitin right curves+SASA','d':'4chi occupancy'},
 'all767_scalar_values_exactly_preserved':True,'profile_values':368,'SASA_values':97,'occupancy_coordinates':302,
 'highlight_change':{'previous':[41,2],'current':[41,62],'Q2_retained_as_other':True,'Q62_actual_coordinates':True},
 'case_curve_box_mm':curve_boxes,'case_SASA_box_mm':sasa_boxes,'alignment':align,
 'left_ylabel_actual_bbox_centre_pixels':label_centres,'left_ylabel_figure_x_mm':3.6,
 'case_gutter_mm':9,'previous_case_gutter_mm':13.5,'previous_case_width_mm':{'CaMKI':37.5,'Ubiquitin':108},
 'SASA_bar_colour':n.BLUE,'SASA_bar_opacity':1.,'SASA_ylabel_colour':'#000000','SASA_ylabel_opacity':1.,
 'single_palette_source':'palette.py','Fig2_primary_BLUE':n.BLUE,'Fig2_primary_ORANGE':n.WARM,
 'non_density_blue_and_orange_elements_full_opacity':True,
 'exceptions':'Density opacity encodes cell frequency; PyMOL lighting provides 3D shading; point/line antialiasing affects edge pixels.',
 'method_legend':'paired blue/ochre lines, Full solid / MD dashed',
 'density_legend_labels':['MD','Full'],
 'density_legend_position':'outside whole coordinate plot, above marginal distribution, right aligned, horizontal, 8pt',
 'structure_label_colour':n.ACCENT,'structure_labels_near_sites':True,'structure_site_circles':False,'structure_site_leaders':['W303'],
 'highlights_colour':n.ACCENT,'green_removed_from_all_active_panels':True,
 'case_site_label_colour':n.ACCENT,'case_site_label_weight':'bold',
 'all_four_case_site_markers':{'colour':n.BLUE,'type':'vertical dotted line','width_pt':.6,'opacity':1.},
 'occupancy_block_mm':{'title_top':113,'axes_top':122.5,'axes_height':21,'shared_xlabel_top':150.5},
 'density_area_change_percent':100*(59*36/(45*29)-1),'density':density_audit,'outside_text':outside,
 'real_structure_renders':True,'structure_images_native_recoloured_only':True,'analysed_chain_colour':n.BLUE,
 'ubiquitin_camera_changed_in_this_revision':False,'structure_coordinates_unchanged':True,
 'new_model_or_MD_run':False,'visual_review':'pending',
 'sha256':{ext:sha256((OUT/f'fig3_orange_sites_black_sasa.{ext}').read_bytes()).hexdigest() for ext in ['png','pdf','svg']}}
(OUT/'figure_QA_v123.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
plt.close(fig)
print(json.dumps(qa,ensure_ascii=False))
