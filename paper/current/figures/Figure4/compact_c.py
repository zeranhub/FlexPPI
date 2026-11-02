# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (regional prediction strips) heatmap/plot_comparison.py + project v124 -> param inherit
# (molecular panels) real MD/PyMOL renders in structure/ -> native raster placement
# (regional MD validation) same zero-centered residue strips -> param inherit
# The generic asset is a model/property matrix, not residue-aligned paired metrics.
# Preserve clean split blocks and direct labels; do not import its normalization or data.
# All 45 CDR and 30 pHLA interface positions in the stated scope are retained.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
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
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

"""Compact paired-metric maturation heatmap for Figure 4c.
All 45 CDR and 30 existing pHLA interface positions are retained.
Importing this module never exports files or modifies a manuscript.
"""
from pathlib import Path
from functools import lru_cache
from hashlib import sha256
import csv
import json
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

PAIR = Path(__file__).resolve().parents[4] / 'experiments/_sources/tcr_expanded_screen_v41_20261017/mel5_validation/corrected_pair_all_aligned_positions.csv'
OUT = Path(__file__).resolve().parent
GROUPS = ['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
BLUE='#2166AC'
WARM='#B87936'
BLACK='#222222'
MID='#62676B'
NA='#D9DDE0'
CMAP=LinearSegmentedColormap.from_list('compact_maturation_signed',[BLUE,'#F7F7F7',WARM])
NORMS={'rmsf':TwoSlopeNorm(vmin=-.85,vcenter=0,vmax=.85),
       'hchi':TwoSlopeNorm(vmin=-.45,vcenter=0,vmax=.45)}
BASE_W,BASE_H=183.,52.


@lru_cache(maxsize=1)
def source_data():
    """Return the exact original selection, its two metric arrays, and metadata."""
    with PAIR.open(encoding='utf-8-sig',newline='') as f:
        rows=[r for r in csv.DictReader(f) if (r['a'],r['b'])==('3hg1','4jff')]
    cdr=[r for g in GROUPS for r in rows if r['region']==g]
    iface=[r for r in rows if r['role'] in ['mhc','peptide'] and r['interface']=='True']
    assert len(cdr)==45 and len(iface)==30
    arrays=[]
    for data in [cdr,iface]:
        rmsf=np.array([float(r['rmsf_full_delta']) for r in data])
        hchi=np.array([float(r['hchi_full_delta']) if r['same_aa']=='True' and r['hchi_valid']=='True' else np.nan for r in data])
        assert all(r['rmsf_valid']=='True' for r in data)
        assert abs(rmsf).max()<=.85 and np.nanmax(abs(hchi))<=.45
        arrays.append(np.vstack([rmsf,hchi]))
    assert np.isfinite(arrays[0][1]).sum()==18
    assert np.isfinite(arrays[1][1]).sum()==23
    groups1=[(g,[r for r in cdr if r['region']==g]) for g in GROUPS]
    groups2=[('HLA α1',[r for r in iface if r['region']=='α1 helix']),
             ('HLA α2',[r for r in iface if r['region']=='α2 helix']),
             ('Peptide',[r for r in iface if r['role']=='peptide'])]
    assert sum(len(g) for _,g in groups2)==30
    return cdr,iface,arrays,groups1,groups2


def draw_compact_c(fig,bbox_mm=(0,118,183,52),canvas_mm=(183,174)):
    """Draw native vector cells/text in bbox (x,y,w,h), with y measured from top.

    Bbox height 52 mm is sufficient at 183 mm width. Font sizes remain in pt,
    so resizing the bbox changes spacing without silently shrinking text.
    """
    bx,by,bw,bh=map(float,bbox_mm)
    cw,ch=map(float,canvas_mm)
    sx,sy=bw/BASE_W,bh/BASE_H
    cdr,iface,arrays,cgroups,igroups=source_data()
    text_artists=[]
    def xy(x,y):return ((bx+x*sx)/cw,1-(by+y*sy)/ch)
    def txt(x,y,value,fs=8,**kwargs):
        a=fig.text(*xy(x,y),value,fontsize=fs,va=kwargs.pop('va','center'),**kwargs)
        text_artists.append(a)
        return a
    def rect(x,y,w,h,**kwargs):
        a=Rectangle(((bx+x*sx)/cw,1-(by+(y+h)*sy)/ch),w*sx/cw,h*sy/ch,
                    transform=fig.transFigure,**kwargs)
        fig.add_artist(a)
        return a
    def axmm(x,y,w,h):
        return fig.add_axes([(bx+x*sx)/cw,1-(by+(y+h)*sy)/ch,w*sx/cw,h*sy/ch])

    txt(2,3.4,'c',fs=11,fontweight='bold')
    txt(9,3.4,'TCR maturation: residue-resolved dynamics',fs=9.5,fontweight='medium')
    geometry=[]
    def band(data,split_groups,top,identity,number_all):
        x=27.0;width=138.0;gap=1.2
        cell=(width-gap*(len(split_groups)-1))/len(data)
        txt(2,top+1.65,identity,fs=8,fontweight='bold')
        for k,label in enumerate(['ΔRMSF (Å)',r'Δ$H_{\chi}$ (nats)']):
            txt(25.3,top+5.2+k*3.4,label,fs=7.8,ha='right')
        mapping={}
        for label,group in split_groups:
            start=x;span=cell*len(group)
            rect(start,top,span,3.3,facecolor='#F1F3F4',edgecolor='none')
            txt(start+span/2,top+1.65,label,fs=8,ha='center',fontweight='bold')
            for r in group:
                for k,metric in enumerate(['rmsf','hchi']):
                    valid=(metric=='rmsf') or (r['same_aa']=='True' and r['hchi_valid']=='True')
                    value=float(r[metric+'_full_delta']) if valid else np.nan
                    color=CMAP(NORMS[metric](value)) if valid else NA
                    rect(x,top+3.5+k*3.4,cell,3.4,facecolor=color,edgecolor='white',linewidth=.20)
                key=(r['role'],r['id_a'])
                mapping[key]=(x,cell)
                col=WARM if r['same_aa']!='True' or key==('mhc','72') else BLACK
                txt(x+cell/2,top+12.35,r['aa_a'],fs=7.5,ha='center',color=col,
                    fontweight='bold' if col==WARM else 'normal')
                geometry.append({'role':r['role'],'id':r['id_a'],'AA':r['aa_a'],
                    'left_mm':x,'width_mm':cell,'top_mm':top,
                    'same_aa':r['same_aa']=='True','hchi_valid':r['hchi_valid']=='True'})
                x+=cell
            if number_all:
                chosen=range(len(group))
            else:
                chosen=sorted(set([0,len(group)-1]+[j for j,r in enumerate(group)
                    if (r['role'],r['id_a'])==('tcr_beta','54')]))
            for j in chosen:
                r=group[j];key=(r['role'],r['id_a'])
                number_y=top+(18.3 if number_all and label=='HLA α2' and j%2 else 15.5)
                txt(start+(j+.5)*cell,number_y,r['id_a'],fs=7.5,ha='center',
                    color=WARM if key==('mhc','72') else MID,
                    fontweight='bold' if key==('mhc','72') else 'normal')
            x+=gap
        return mapping
    top1=9.4;top2=30.4
    cpos=band(cdr,cgroups,top1,'TCR CDRs',False)
    ipos=band(iface,igroups,top2,'pHLA',True)
    for mapping,key,top in [(cpos,('tcr_beta','54'),top1),(ipos,('mhc','72'),top2)]:
        left,cell=mapping[key]
        rect(left,top+3.5,cell,6.8,fill=False,edgecolor=WARM,linewidth=.8)

    # Two physical-unit scales are vertical and stacked in the right margin.
    bars=[]
    for metric,y,label in [('rmsf',12.4,'ΔRMSF\n(Å)'),('hchi',33.4,'ΔHχ\n(nats)')]:
        txt(174.0,y-4.2,label,fs=7.8,ha='center',linespacing=1.0)
        cbax=axmm(170.4,y,2.1,12.0)
        limit=.85 if metric=='rmsf' else .45
        cb=fig.colorbar(mpl.cm.ScalarMappable(norm=NORMS[metric],cmap=CMAP),
                       cax=cbax,orientation='vertical',ticks=[-limit,0,limit])
        cb.set_ticklabels([f'−{limit:.2f}','0',f'{limit:.2f}'])
        cb.outline.set_linewidth(.35)
        cb.ax.tick_params(labelsize=7.5,length=1.3,pad=1.0,width=.35)
        bars.append(cb)

    qa={
        'scope':'TCR maturation only: 4jff minus 3hg1',
        'residue_counts':{'CDR':45,'pHLA':30},
        'metric_rows_per_region':['Delta RMSF (angstrom)','Delta Hchi (nats)'],
        'hchi_comparable':{'CDR':18,'pHLA':23},
        'hchi_mask':'same_aa AND hchi_valid; all other cells retained as gray NaN, never zero',
        'gray_cells':{'CDR':27,'pHLA':7},
        'limits':{'rmsf':[-.85,.85],'hchi':[-.45,.45]},
        'colorbar_orientation':'vertical',
        'minimum_font_pt':7.5,
        'group_order':{'TCR':GROUPS,'pHLA':['HLA alpha1','HLA alpha2','Peptide']},
        'bbox_mm':list(bbox_mm),'canvas_mm':list(canvas_mm),
        'explanatory_sentences_in_figure':False,
        'sequence':'MEL5 reference AA; warm letters mark substitutions; Q72 also warm',
        'numbering':'MEL5 common TCR numbering; every pHLA author number retained. Seven HLA alpha2 numbers alternate two rows for spacing.',
        'source':str(PAIR),'source_sha256':sha256(PAIR.read_bytes()).hexdigest(),
        'coordinates_mm':geometry,
        'text_artists':text_artists,
        'arrays':{'CDR':arrays[0].copy(),'pHLA':arrays[1].copy()},
    }
    return qa


def export_standalone():
    """Produce inspectable native-vector previews and exact displayed source data."""
    import matplotlib.pyplot as plt
    from matplotlib.text import Text
    OUT.mkdir(exist_ok=True)
    fig=plt.figure(figsize=(183/25.4,56/25.4),facecolor='white')
    qa=draw_compact_c(fig,(0,2,183,52),(183,56))
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    text_checks=[]
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text():continue
        bb=t.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        text_checks.append({'text':t.get_text(),'font_pt':t.get_fontsize(),'bbox':bb.extents.tolist(),
                            'inside':bool(bb.x0>=0 and bb.x1<=1 and bb.y0>=0 and bb.y1<=1)})
    assert all(x['inside'] for x in text_checks),[x for x in text_checks if not x['inside']]
    assert min(x['font_pt'] for x in text_checks)>=7.5
    for ext in ['png','pdf','svg']:
        with mpl.rc_context({'savefig.bbox':None}):
            fig.savefig(OUT/f'panel_d_preview.{ext}',dpi=600,bbox_inches=None)
    plt.close(fig)
    cdr,iface,arrays,_,_=source_data()
    with (OUT/'panel_d_source_75_positions.csv').open('w',encoding='utf-8-sig',newline='') as f:
        fields=['role','id_a','region','aa_a','aa_b','same_aa','rmsf_valid','hchi_valid',
                'rmsf_full_a','rmsf_full_b','rmsf_full_delta','hchi_full_a','hchi_full_b','hchi_full_delta','hchi_displayed']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for r in cdr+iface:
            rr={k:r[k] for k in fields if k!='hchi_displayed'}
            rr['hchi_displayed']=float(r['hchi_full_delta']) if r['same_aa']=='True' and r['hchi_valid']=='True' else ''
            writer.writerow(rr)
    np.savez(OUT/'panel_d_arrays.npz',CDR=arrays[0],pHLA=arrays[1])
    qa.pop('text_artists');qa.pop('arrays')
    qa['text_checks']=text_checks;qa['status']='generated_pending_visual_QA';qa['manuscript_modified']=False
    (OUT/'panel_d_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:qa[k] for k in ['residue_counts','hchi_comparable','gray_cells','minimum_font_pt']},ensure_ascii=False))


if __name__=='__main__':
    export_standalone()

