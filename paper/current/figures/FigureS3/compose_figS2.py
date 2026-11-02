# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Existing verified PyMOL scene -> _case_validation_sources/ubiquitin -> native run
# (b) Existing verified real MD frame scenes -> _case_validation_sources/ubiquitin -> native run
# (c) Residue line profile -> existing compose_figure.py panel g -> param inherit
# (d) Cumulative state coverage -> existing compose_figure.py panel h -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
# The supplied scientific production source is authoritative. Skill LineTrend's
# publication timeline has incompatible semantics; only baseline styling applies.
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

from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.text import Text
from PIL import Image
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
SOURCE = Path(r'E:\Project\FlexPPI\experiments\_case_validation_sources')
W, H = 183., 130.
BLUE, WARM, TEAL, MD = '#2166AC', '#B87936', '#3A8D8C', '#666666'
# The explicitly requested 8 pt floor supersedes the baseline's 7 pt ticks.
mpl.rcParams.update({'xtick.labelsize':8,'ytick.labelsize':8,'mathtext.fontset':'dejavusans'})

def ax_mm(fig,x,y,w,h):
    return fig.add_axes([x/W,y/H,w/W,h/H])

def text(fig,x,y,label,**kwargs):
    kwargs.setdefault('fontsize',8); kwargs.setdefault('color',BLACK)
    return fig.text(x/W,y/H,label,**kwargs)

def title(fig,label,x,y,description):
    text(fig,x,y,label,fontsize=10,fontweight='bold')
    text(fig,x+6,y,description)

def image_panel(fig,name,x,y,w,h):
    path=HERE/'structure_renders'/name
    with Image.open(path) as im:
        pixels=np.asarray(im.convert('RGB'))
    ax=ax_mm(fig,x,y,w,h)
    ax.imshow(pixels); ax.set_axis_off()
    return ax

def main():
    data=HERE/'source_data';data.mkdir(exist_ok=True)
    for name in ['ubiquitin_fold0_all_residues.csv','ubiquitin_Q41_Q62_all27_joint_frequencies.csv','ubiquitin_Q41_Q62_validation_summary.csv']:
        if not (data/name).exists(): shutil.copy2(SOURCE/'data'/name,data/name)
    for name in ['representative_frames.csv','provenance.json','reference_environment_contacts.csv','reference_solvent_accessibility.csv','frame_geometry_audit.csv']:
        if not (data/name).exists(): shutil.copy2(SOURCE/'ubiquitin'/name,data/name)
    ub=pd.read_csv(data/'ubiquitin_fold0_all_residues.csv')
    ub=ub.loc[ub.chain.eq('A')].sort_values('residue_id').copy()
    rmsf=ub.loc[ub.valid_rmsf.eq(True)].sort_values('residue_id').copy()
    frequencies=pd.read_csv(data/'ubiquitin_Q41_Q62_all27_joint_frequencies.csv')
    summary=pd.read_csv(data/'ubiquitin_Q41_Q62_validation_summary.csv')
    assert len(rmsf)==72 and len(ub.loc[ub.valid_hchi.eq(True)])==66
    assert list(rmsf.residue_id)==list(range(1,73))
    assert np.isfinite(rmsf[['md_rmsf','pred_rmsf']]).all().all()
    assert len(frequencies)==54
    assert np.allclose(frequencies.groupby('residue').probability.sum(),1.)
    assert list(summary.states_covering_90_percent)==[7,12]
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    title(fig,'a',2,123,'2QHO ubiquitin–UBA complex')
    title(fig,'b',62,123,'MD conformational states')
    image_panel(fig,'ubiquitin_context.png',2,67,54,45)
    image_panel(fig,'q41_zoom.png',61,66,56,46.67)
    image_panel(fig,'q62_zoom.png',122,66,56,46.67)
    # Replaced the original tiny raster text with print-sized vector labels.
    text(fig,10,105,'UBA',ha='center',fontweight='bold',color=MD)
    text(fig,37,103.5,'Ub',ha='center',fontweight='bold',color=BLUE)
    text(fig,38,94,'Q41',ha='center',fontweight='bold',color=WARM)
    text(fig,46,75,'Q62',ha='center',fontweight='bold',color=TEAL)
    text(fig,89,115,'Q41',ha='center',color=WARM,fontweight='bold')
    text(fig,150,115,'Q62',ha='center',color=TEAL,fontweight='bold')
    title(fig,'c',2,58,'Full ubiquitin backbone profile')
    title(fig,'d',112,58,'Joint-state coverage')
    ax=ax_mm(fig,15,18,88,31)
    for source,color,label,marker in [('md',MD,'MD','o'),('pred',BLUE,'FlexPPI',None)]:
        ax.plot(rmsf.residue_id,rmsf[source+'_rmsf'],color=color,lw=1.1,label=label,
                marker=marker,markersize=2,markerfacecolor='white',markeredgewidth=.55,
                zorder=3 if source=='md' else 2)
    for number,color in [(41,WARM),(62,TEAL)]:
        ax.axvline(number,color=color,lw=.7,ls=':',zorder=1)
    ax.set_xlim(1,72);ax.set_ylim(0,1.5)
    ax.set_xticks([1,20,41,62,72]);ax.set_yticks([0,.5,1.,1.5])
    assert rmsf[['md_rmsf','pred_rmsf']].to_numpy().max() < ax.get_ylim()[1]
    ax.set_xlabel('Ubiquitin residue number',labelpad=3)
    ax.set_ylabel('RMSF (Å)',labelpad=3)
    ax.tick_params(length=2.4,pad=2)
    ax.legend(loc='lower left',bbox_to_anchor=(0,1.025),borderaxespad=0,ncol=2,handlelength=1.5,columnspacing=1.1)
    rho=float(spearmanr(rmsf.md_rmsf,rmsf.pred_rmsf).statistic)
    ax.text(.99,1.06,rf'$\rho={rho:.3f}$',transform=ax.transAxes,ha='right',va='bottom',fontsize=8)
    coverage=ax_mm(fig,125,18,52,31)
    coverage.axhline(90,color='#999999',ls=':',lw=.8)
    coverage.text(27,88,'90%',ha='right',va='top',fontsize=8,color=MD)
    coverage_records=[]
    for residue,color,marker,y in [('Q41',WARM,'o',63),('Q62',TEAL,'s',34)]:
        freq=frequencies.loc[frequencies.residue.eq(residue)].sort_values(['probability','joint_code'],ascending=[False,True])
        cumulative=freq.probability.cumsum().values*100
        n90=int(np.searchsorted(cumulative,90)+1)
        coverage.plot(range(1,28),cumulative,color=color,lw=1.2,marker=marker,ms=2.2,label=residue)
        coverage.plot([n90,n90],[0,cumulative[n90-1]],color=color,ls='--',lw=.7)
        coverage.text(n90+.5,y,f'{residue}: {n90}',fontsize=8,color=color)
        coverage_records.extend({'residue':residue,'rank':i+1,'joint_code':int(code),'probability':float(prob),'cumulative_percent':float(cum)} for i,(code,prob,cum) in enumerate(zip(freq.joint_code,freq.probability,cumulative)))
    coverage.set_xlim(1,27);coverage.set_ylim(0,104)
    coverage.set_xticks([1,7,12,20,27]);coverage.set_yticks([0,50,100])
    coverage.set_xlabel('Classes ranked by occupancy',labelpad=3)
    coverage.set_ylabel('MD frames covered (%)',labelpad=3)
    coverage.tick_params(length=2.4,pad=2)
    coverage.legend(loc='lower left',bbox_to_anchor=(0,1.025),borderaxespad=0,ncol=2,handlelength=1.2,columnspacing=1.)
    text(fig,2,2.8,'MD conformations; fold 0 held-out scalar predictions.',fontsize=8,color=MD)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    font_sizes=[item.get_fontsize() for item in fig.findobj(Text) if item.get_text()]
    assert min(font_sizes)>=8
    text_outside=[]
    figure_box=fig.bbox
    for item in fig.findobj(Text):
        if not item.get_text() or not item.get_visible():continue
        box=item.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > figure_box.x1+1 or box.y1 > figure_box.y1+1:
            text_outside.append(item.get_text())
    assert not text_outside,text_outside
    target=HERE/'figS2_ubiquitin_supplement'
    for extension in ('png','pdf','svg'):
        fig.savefig(str(target)+'.'+extension,bbox_inches=fig.bbox_inches,dpi=400 if extension=='png' else 600)
    pd.DataFrame(coverage_records).to_csv(data/'plotted_all_54_coverage_points.csv',index=False)
    rmsf[['residue_id','md_rmsf','pred_rmsf']].to_csv(data/'plotted_all_72_RMSF_positions.csv',index=False)
    report={'canvas_mm':[W,H],'PNG_dpi':400,'minimum_vector_text_pt':min(font_sizes),'RMSF_positions':len(rmsf),'RMSF_rho':rho,
            'coverage_points':len(coverage_records),'n90_Q41_Q62':[7,12],'Hchi_profile_repeated':False,
            'new_simulations_or_frames':False,'source_script':str(SOURCE/'compose_figure.py'),
            'representative_frame_indices':{'Q41':[2397,401,1171],'Q62':[924,1143,1467]},
            'native_vector_quantitative_panels':True,'text_outside_canvas':text_outside,
            'assets':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(data.iterdir()) if p.is_file()]}
    (HERE/'figure_QA_and_sources.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    plt.close(fig)
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':main()

