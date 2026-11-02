# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Actual molecular overview + one local view -> archived PyMOL render -> native raster placement
# (b,c) Residue profiles -> LineTrend/plot_trend.py -> param inherit
# (d upper/lower) Complete aligned delta maps -> heatmap/plot_composition.py -> param inherit
# Assets and companion PNGs inspected. Monthly counts, grouped fractions and clinical-stage
# counts have different data/axis semantics; preserve clean axes and grouping, use existing
# scientific source mappings. All user-requested values retained; native vectors for all charts.
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})
# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"
# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
import hashlib, json, shutil, sys
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,ConnectionPatch
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.text import Text
from matplotlib.transforms import Bbox

FIG=Path(__file__).resolve().parent
ROOT=FIG.parent
EXP=ROOT.parent
sys.path.insert(0,str(FIG))
from publication_style import BLUE, WARM, MD, MID, LIGHT, TEAL, SIGNED, signed_cmap, signed_norm
W,H=183,136
WT=WARM; ENG=BLUE; GRAY=MD
G4A=MD; I7A=MID
NA=LIGHT
STATES=['3hg1','4jff','4jfd','4jfe']
NAMES=['MEL5 · ELA','α24β17 · ELA','α24β17 · G4A','α24β17 · I7A']
COLORS=[WT,ENG,G4A,I7A]
MARKERS=['o','s','D','^']
STYLES=['-','-','--',':']
PAIRS=[('3hg1','4jff'),('4jff','4jfd'),('4jff','4jfe')]
PAIR_NAMES=['TCR maturation','Peptide G4A','Peptide I7A']
GROUPS=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
paths={
    'CDR_absolute':FIG/'source_data/source_CDR_absolute.csv',
    'interface_absolute':FIG/'source_data/source_interface_absolute.csv',
    'peptide_four_states':FIG/'source_data/source_peptide_four_states.csv',
    'overview':FIG/'结构原图/locator_trimmed.png',
    'local_structure':FIG/'结构原图/4jff_zoom_no_labels.png',
}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
cdr=pd.read_csv(paths['CDR_absolute'],float_precision='round_trip')
iface=pd.read_csv(paths['interface_absolute'],float_precision='round_trip')
pep=pd.read_csv(paths['peptide_four_states'],float_precision='round_trip')
assert [len(t) for t in [cdr,iface,pep]]==[45,23,40]

fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def axes(x,y,w,h): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(x,y,s,**kw): return fig.text(x/W,1-y/H,s,**kw)
def heading(letter,title,y,x=2,title_x=9):
    text(x,y,letter,fontsize=11,fontweight='bold',va='baseline')
    text(title_x,y,title,fontsize=10,va='baseline',fontweight='medium')
def crop(path,pad=12):
    # Freeze the exact original crop so colour changes cannot shift atom labels.
    frozen=json.loads((FIG/'original_figure_data_audit_v59.json').read_text(encoding='utf-8'))['crop_boxes']
    box=frozen['overview' if Path(path).name=='4jff_overview.png' else 'local']
    im=Image.open(path).convert('RGB')
    return im.crop(box),box

# (a) One complete-complex locator and two equally scaled local close-ups.
heading('a','Q72–CDR2β interface',5)
locator=Image.open(FIG/'结构原图/locator_trimmed.png').convert('RGB')
locator_meta=json.loads((FIG/'结构原图/locator_QA.json').read_text(encoding='utf-8'))
ovbox=locator_meta.get('white_margin_crop_box',[165,577,857,1294])
localbox=[0,0,1400,2100]
# Preserve the overview's physical aspect, and show only one intact complex.
oa=axes(1,20,39,39*locator.height/locator.width)
oa.imshow(locator);oa.axis('off')
roi=[280,213,376,374]
oa.add_patch(Rectangle((roi[0],roi[1]),roi[2]-roi[0],roi[3]-roi[1],
    fill=False,edgecolor='#3B3B3B',lw=.65))
text(27,18,'TCR',fontsize=8,color=BLUE,ha='center')
text(7,65.5,'HLA',fontsize=8,color=GRAY,ha='center')
labels=json.loads((FIG/'四位点无引线_布局.json').read_text(encoding='utf-8'))
local_axes=[]
for state,title,left in [('3hg1','MEL5',44),('4jff','Mature α24β17',80)]:
    text(left+16.5,12.5,title,fontsize=8.3,fontweight='bold',ha='center')
    ax=axes(left,16,33,49.5)
    ax.imshow(Image.open(FIG/f'结构原图/{state}_zoom_no_labels.png').convert('RGB'))
    ax.set_xlim(0,1400);ax.set_ylim(2100,0);ax.axis('off')
    ax.add_patch(Rectangle((0,0),1,1,transform=ax.transAxes,fill=False,
        edgecolor='#333333',lw=.65,clip_on=False))
    for r in labels['states'][state]['labels']:
        ax.text(*r['target_pixel'],r['text'],fontsize=8.0,fontweight='semibold',
            color=r['text_color'],ha='center',va='center',
            bbox=dict(facecolor='white',edgecolor='none',alpha=.94,pad=.25))
    local_axes.append(ax)
# These dashed lines communicate magnification, never residue identification.
for point,corner in [((roi[2],roi[1]),(0,1)),((roi[2],roi[3]),(0,0))]:
    fig.add_artist(ConnectionPatch(xyA=point,coordsA=oa.transData,xyB=corner,
        coordsB=local_axes[0].transAxes,color='#6A6F73',lw=.55,ls=(0,(3,2)),clip_on=False))
text(6,72.5,r'MEL5 → α24β17    $K_D$: 18 μM → 0.6 nM',fontsize=8.2)



# (b) All ten peptide positions for every existing observed bound structure.
heading('b','Peptide dynamics',5,x=119,title_x=126)
b=axes(130,18,48.5,34)
handles=[]
for state,name,col,marker,style in zip(STATES,NAMES,COLORS,MARKERS,STYLES):
    q=pep[pep.pdb.eq(state)].sort_values('id'); assert len(q)==10
    line,=b.plot(q.id,q.rmsf_full,color=col,ls=style,marker=marker,ms=2.2,lw=.9,
                mfc='white' if state=='3hg1' else col,mec=col,mew=.6,label=name)
    handles.append(line)
for pos in [4,10]:
    b.axvline(pos,color='#AAB2B8',lw=.5,ls='--',zorder=0)
    b.text(pos,1.012,f'P{pos}',ha='center',va='bottom',fontsize=8,color=GRAY)
b.set(xlim=(.7,10.3),ylim=(.57,1.02),xticks=range(1,11),yticks=[.6,.7,.8,.9,1],ylabel='Pred. RMSF (Å)',xlabel='Peptide position')
b.tick_params(labelsize=8,length=2,pad=2)
b.xaxis.label.set_size(9); b.yaxis.label.set_size(9)
b.xaxis.labelpad=2;b.yaxis.labelpad=3
state_legend=fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(119/W,1-63/H),ncol=2,
           fontsize=8,handlelength=1.6,handletextpad=.45,columnspacing=.7,labelspacing=.6,borderaxespad=0)

# (c) Paired backbone/side-chain response at the same aligned residues.
from compact_c import draw_compact_c
dqa=draw_compact_c(fig,(0,80,183,52),canvas_mm=(W,H))
d_arrays=dqa.pop('arrays')
dqa.pop('text_artists')
np.savez(FIG/'panel_c_arrays_in_group.npz',CDR=d_arrays['CDR'],pHLA=d_arrays['pHLA'])
fig.canvas.draw();renderer=fig.canvas.get_renderer()
legend_bb=state_legend.get_window_extent(renderer);xlabel_bb=b.xaxis.label.get_window_extent(renderer)
assert xlabel_bb.y0>legend_bb.y1,('Peptide xlabel overlaps legend',xlabel_bb,legend_bb)
texts=[]
for item in fig.findobj(match=Text):
    if not item.get_visible() or not item.get_text():continue
    bb=item.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    texts.append(dict(text=item.get_text(),font_pt=item.get_fontsize(),bounds=list(bb.extents),
        inside_canvas=bool(bb.x0>=-.001 and bb.y0>=-.001 and bb.x1<=1.001 and bb.y1<=1.001)))
assert all(t['inside_canvas'] for t in texts),[t for t in texts if not t['inside_canvas']]
assert min(t['font_pt'] for t in texts)>=7.5
stem=FIG/'figure4_去除原c_三面板'
for ext in ['png','pdf','svg']:
    with mpl.rc_context({'savefig.bbox':None}):fig.savefig(stem.with_suffix('.'+ext),dpi=600,bbox_inches=None)
for label,bounds in [('panel_a_整体与局部新比例',(0,0,116,75)),('panel_c_紧凑双指标热图',(0,80,183,52))]:
    x0,y0,ww,hh=bounds
    for ext in ['png','pdf','svg']:
        with mpl.rc_context({'savefig.bbox':None}):
            fig.savefig(FIG/(label+'.'+ext),dpi=600,bbox_inches=Bbox.from_bounds(x0/25.4,(H-y0-hh)/25.4,ww/25.4,hh/25.4))
plt.close(fig)
sources={key:{'path':str(paths[key]),'sha256':sha(paths[key])} for key in ['CDR_absolute','peptide_four_states']}
qa={'status':'pending_visual_QA','dimensions_mm':[W,H],'raster_dpi':600,
    'manuscript_modified':False,'native_structure_coordinates_changed':False,
    'structural_scale_mm':{'overview_previous':25,'overview_current':39,'local_previous':41,'local_current':33},
    'panel_b_data_unchanged':True,'original_panel_c_removed':True,'original_panel_d_relabelled_c':True,'panel_c':dqa,'sources':sources,
    'minimum_font_pt':min(t['font_pt'] for t in texts),'visible_text':texts}
(FIG/'整组核对.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'file':str(stem.with_suffix('.png')),'dimensions_mm':[W,H],'font_min':qa['minimum_font_pt']},ensure_ascii=False))

