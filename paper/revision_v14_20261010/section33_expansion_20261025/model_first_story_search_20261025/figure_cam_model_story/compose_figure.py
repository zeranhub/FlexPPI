# Academic Figure Skill Asset Confirmation
# (a-c) image_plate + quant: real PyMOL simulation-input snapshots.
# (d,e) all-residue profiles: cross-type inherit from the existing CaM figure
# and LineTrend typography/direct labels. Timeline template data are incompatible.
# (f) scatter, (g,h) signed values: native scientific Matplotlib, param inherit.
# No template example data, point selection, fitted calibration, or invented atoms.
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
import json,hashlib,shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
OUT=Path(__file__).resolve().parent;P=OUT.parent/'partner_screen';W,H=183.,240.
MD='#424B50';FULL='#B2182B';NBLUE='#4C83A1';CTEAL='#52A49A';WARM='#B87936'
texts=[];axes=[]
def text(fig,x,y,s,size=8,**kw):
    t=fig.text(x/W,1-y/H,s,fontsize=size,color=kw.pop('color',BLACK),va=kw.pop('va','top'),**kw);texts.append(t);return t
def ax(fig,x,y,w,h):
    a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);axes.append(a);return a
def panel(fig,l,s,x,y):text(fig,x,y,l,9,fontweight='bold');text(fig,x+5,y,s,8)
def draw_image(fig,info,x,y,w,h):
    a=ax(fig,x,y,w,h);im=Image.open(OUT/'structure_renders'/info['file']);crop=info['crop'];a.imshow(im.crop(crop));a.axis('off');return a
def callout(a,info,key,xytext,fs=6.4):
    p=info['landmarks'][key];p=[p[0]-info['crop'][0],p[1]-info['crop'][1]]
    a.annotate(key,xy=p,xytext=xytext,textcoords='axes fraction',ha='center',va='center',fontsize=fs,
        bbox=dict(fc='white',ec='none',alpha=.88,pad=.45),arrowprops=dict(arrowstyle='-',lw=.45,color='#767D80',shrinkA=1.2,shrinkB=1.2))
def crop_pairs(man):
    for keys in [['whole_NHE1','whole_IQCG'],['M72_NHE1','M72_IQCG'],['M109_NHE1','M109_IQCG']]:
        boxes=[Image.open(OUT/'structure_renders'/man['assets'][k]['file']).getchannel('A').getbbox() for k in keys]
        ww,hh=Image.open(OUT/'structure_renders'/man['assets'][keys[0]]['file']).size
        box=[max(0,min(b[0] for b in boxes)-25),max(0,min(b[1] for b in boxes)-25),min(ww,max(b[2] for b in boxes)+25),min(hh,max(b[3] for b in boxes)+25)]
        for k in keys:man['assets'][k]['crop']=box
def main():
    data=pd.read_csv(P/'cam_md_state_contrast_all_common_residues.csv');mets=pd.read_csv(P/'final_cam_all_nine_methionines_with_stability.csv')
    groups=pd.read_csv(P/'final_cam_fixed_structural_group_metrics.csv');blocks=pd.read_csv(P/'cam_md_state_contrast_all_residues_five_blocks.csv')
    assert len(data)==267 and sum(data.metric=='hchi')==122 and sum(data.metric=='rmsf')==145
    man=json.loads((OUT/'structure_renders/render_manifest.json').read_text());crop_pairs(man)
    sd=OUT/'source_data';sd.mkdir(exist_ok=True)
    for src in ['cam_md_state_contrast_all_common_residues.csv','final_cam_all_nine_methionines_with_stability.csv','final_cam_fixed_structural_group_metrics.csv','cam_md_state_contrast_all_residues_five_blocks.csv','cam_md_state_contrast_region_metrics.csv','cam_md_mature_number_mapping.csv']:
        shutil.copy2(P/src,sd/src)
    fig=plt.figure(figsize=(W/25.4,H/25.4));fig.patch.set_facecolor('white')
    text(fig,1,2,'CaM | Backbone motion and side-chain entropy change differently',9.2,fontweight='bold')
    text(fig,1,8,'NHE1-based state → IQCG-based state · MD and frozen FlexPPI predictions',7.2,color='#5C656A')
    panel(fig,'a','Two actual simulation inputs',1,16)
    for i,label in enumerate(['NHE1','IQCG']):
        x=2+i*46;text(fig,x+21,21,label+'–CaM',7.5,ha='center',fontweight='bold')
        info=man['assets']['whole_'+label];a=draw_image(fig,info,x,25,44,45)
        for k,loc in [('M72',(.12,.75)),('M109',(.86,.48)),('M144',(.76,.12))]:callout(a,info,k,loc,6.0)
        text(fig,x+21,70,'Based on '+('2YGG' if label=='NHE1' else '4LZX'),6.3,ha='center',color='#626B70')
    for x,c,s in [(3,NBLUE,'N-lobe'),(26,CTEAL,'C-lobe'),(49,WARM,'Partner')]:
        text(fig,x,75,'■',7,color=c);text(fig,x+3.5,75,s,6.2,color='#626B70')
    for x,l,site in [(99,'b',72),(143,'c',109)]:
        panel(fig,l,'M'+str(site)+' environment',x,16)
        for i,label in enumerate(['NHE1','IQCG']):
            info=man['assets'][f'M{site}_{label}'];a=draw_image(fig,info,x+i*20.3,25,19.3,43)
            callout(a,info,'M'+str(site),(.55,.13),6.0)
            text(fig,x+9.6+i*20.3,70,label,6.4,ha='center',color='#626B70')
    text(fig,181,76,'Snapshots locate sites; they do not depict motion.',5.9,ha='right',color='#626B70')
    for x,l,metric,title in [(1,'d','rmsf','Backbone: all 145 common positions'),(96,'e','hchi','Side chains: all 122 valid residues')]:
        panel(fig,l,title,x,83);a=ax(fig,x+10,91,74,32)
        g=data[data.metric==metric].set_index('mature_cam_number').reindex(range(4,149))
        for lo,hi in [(20,31),(56,67),(93,104),(129,140)]:a.axvspan(lo,hi,color='#D9E8E5',alpha=.65,lw=0,zorder=0)
        a.axhline(0,color='#B7BDC1',lw=.5)
        a.plot(g.index,g.md_delta,color=MD,lw=1,label='MD')
        a.plot(g.index,g.full0_delta,color=FULL,lw=.9,label='Full')
        a.set_xlim(1,150);a.set_xticks([20,60,100,140]);a.set_xlabel('Mature CaM position',fontsize=7)
        a.set_ylabel('ΔRMSF (Å)' if metric=='rmsf' else 'ΔHχ (nats)',fontsize=7,labelpad=3)
        a.text(.98,.98,'ρ = '+('0.827' if metric=='rmsf' else '0.630'),transform=a.transAxes,ha='right',va='top',fontsize=6.7)
        if metric=='rmsf':a.legend(fontsize=6.1,loc='lower right',ncol=2,handlelength=1.5,columnspacing=.8)
        else:a.text(.01,.02,'Green shading: four EF-hand loops',transform=a.transAxes,fontsize=5.8,color='#5A746E',va='bottom')
    panel(fig,'f','Changes across the whole side-chain set',1,136)
    a=ax(fig,12,144,70,35);g=data[data.metric=='hchi'];colors=np.where(g.contact_b>0,CTEAL,'#9CA7AE')
    lo=min(g.md_delta.min(),g.full0_delta.min());hi=max(g.md_delta.max(),g.full0_delta.max());pad=.10*(hi-lo)
    limits=[float(lo-pad),float(hi+pad)]
    a.axhline(0,color='#B7BDC1',lw=.5);a.axvline(0,color='#B7BDC1',lw=.5);a.plot(limits,limits,color='#C4C9CC',lw=.6,ls='--')
    a.scatter(g.md_delta,g.full0_delta,c=colors,s=12,edgecolors='white',linewidths=.25,alpha=.9)
    for site,offset in [(72,(-12,8)),(109,(8,-12)),(144,(8,4))]:
        r=g[g.mature_cam_number==site].iloc[0];a.annotate('M'+str(site),(r.md_delta,r.full0_delta),xytext=offset,textcoords='offset points',fontsize=6,arrowprops=dict(arrowstyle='-',lw=.4,color='#777777'))
    a.set_xlim(limits);a.set_ylim(limits);a.set_xlabel('MD ΔHχ (nats)',fontsize=7);a.set_ylabel('Predicted ΔHχ (nats)',fontsize=7)
    a.text(.02,.98,'58/70 larger changes: correct direction\n|MD ΔHχ| ≥ 0.1 nats',transform=a.transAxes,ha='left',va='top',fontsize=6.1)
    panel(fig,'g','Contact regions change in opposite directions',96,136)
    a=ax(fig,109,144,69,35);names=['NHE1_contact_only','IQCG_contact_only','contact_in_both_states','outside_both_CA8_interfaces']
    gg=groups[(groups.metric=='hchi')].set_index('region').loc[names]
    xx=np.arange(4);a.axhline(0,color='#B7BDC1',lw=.5)
    a.bar(xx-.17,gg.MD_delta,width=.30,color=MD,label='MD');a.bar(xx+.17,gg.Full_delta,width=.30,color=FULL,label='Full')
    a.set_xticks(xx,['NHE1 only\nn = 8','IQCG only\nn = 21','Both\nn = 3','Neither\nn = 90'],fontsize=6.1)
    a.set_ylabel('Mean ΔHχ (nats)',fontsize=7,labelpad=3);a.set_ylim(-.28,.21)
    a.legend(loc='upper left',fontsize=6,ncol=2,handlelength=1,columnspacing=.8)
    text(fig,178,190,'Contact: another chain Cα within 8 Å. Teal scatter points: IQCG-contact positions.',5.8,ha='right',color='#626B70')
    panel(fig,'h','All nine methionines: local side-chain changes despite a common sequence',1,198)
    a=ax(fig,12,206,166,20);m=mets[mets.metric=='hchi'].sort_values('mature_cam_number');assert len(m)==9
    xx=np.arange(9);a.axhline(0,color='#B7BDC1',lw=.5)
    a.bar(xx-.17,m.md_delta,width=.30,color=MD,label='MD');a.bar(xx+.17,m.full0_delta,width=.30,color=FULL,label='Full')
    a.set_xticks(xx,['M'+str(p) for p in m.mature_cam_number],fontsize=6.7);a.set_ylabel('ΔHχ (nats)',fontsize=7,labelpad=3)
    a.set_ylim(-.58,.64);a.set_yticks([-.4,0,.4]);a.legend(loc='upper left',fontsize=6,ncol=2,handlelength=1.2)
    text(fig,181,234,'Training-set illustration: fixed fold 0; both states seen during training. MD: one 500 ns trajectory per state.',5.8,ha='right',color='#626B70')
    text(fig,181,237,'All changes are IQCG minus NHE1. Hχ is torsional occupancy entropy; these are MD results, not experimental S².',5.8,ha='right',color='#626B70')
    fig.canvas.draw();renderer=fig.canvas.get_renderer();issues=[]
    for t in texts:
        b=t.get_window_extent(renderer)
        if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:issues.append(t.get_text())
    for ext,dpi in [('png',400),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig_cam_model_supported_story.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
    (OUT/'figure_QA.json').write_text(json.dumps(dict(out_of_canvas_text=issues,physical_size_mm=[W,H],PNG_dpi=400,
        data_rows=len(data),valid_sidechains=122,whole_common_backbone=145,methionines=9,figure_kind='Training-set interpretive case',
        visual_QA='Pending image inspection',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
    assert not issues,issues
    plt.close(fig)
if __name__=='__main__':main()
