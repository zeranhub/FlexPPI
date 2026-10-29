# Asset Confirmation Table (see output/tcr_real_states_20260915/figure_contract.md)
# 3D panels | deposited-coordinate PyMOL scenes | param inherit; native asset mismatch
# Profiles | LineTrend/plot_trend.py, plot_sweep.py + previews | param inherit
# Fold points | GroupedBarChart/plot_GroupedBarChartv1.py + preview | param inherit
# Matrices | heatmap/plot_comparison.py + preview | param inherit; no column normalization
# All companion previews and source scripts were inspected before implementation.

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

from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/tcr_real_states_20260915'
ASSETS=OUT/'structure_assets'
COLORS=['#2166AC','#B35806','#762A83','#666666','#A3A3A3']
MARKERS=['o','s','^','D','v']
STYLES=['-','--',':','-.','--']
APO='#526581';HOLO='#69A69B';SITE='#762A83'
BIO_B4=[f'b4_{p}_apo_holo_full' for p in ['ab','cd','ef']]
BIO_1G4=[f'1g4_{p}_holo_full' for p in ['8tw4','8tw6','9bbc']]
EXPORT_AUDIT=[]


def export(fig, stem, claim):
    fig.canvas.draw()
    save_cns_figure(fig, OUT/stem)
    fig.savefig(OUT/(stem+'.svg'),bbox_inches='tight',dpi=300)
    texts=[x for x in fig.findobj(mpl.text.Text) if x.get_visible() and x.get_text()]
    EXPORT_AUDIT.append(dict(file=stem,claim=claim,width_mm=float(fig.get_figwidth()*25.4),height_mm=float(fig.get_figheight()*25.4),
        minimum_text_pt=min(t.get_fontsize() for t in texts),n_text_objects=len(texts),
        raster_content='Only physically rendered 3D structures; all plotted data and text are vector.',
        no_data_downsampling=True))
    plt.close(fig)


def heading(fig, letter, title, x, y):
    fig.text(x,y,letter,fontweight='bold',fontsize=10,va='bottom')
    fig.text(x+.025,y,title,fontsize=8,va='bottom')


def title(fig, name, detail):
    fig.text(.025,.98,name,fontsize=11,fontweight='bold',va='top')
    subtitle_y=.98-5.2/(fig.get_figheight()*25.4)
    fig.text(.025,subtitle_y,detail,fontsize=7.5,color='#555555',va='top')


def structure(fig, box, name, local=False):
    ax=fig.add_axes(box)
    im=Image.open(ASSETS/(name+'.png'))
    ax.imshow(im);ax.set_axis_off()
    # Tighten only blank camera margins; every nonwhite rendered pixel remains.
    pixels=np.asarray(im.convert('RGB'))
    yy,xx=np.where(np.any(pixels<245,axis=2))
    ax.set_xlim(max(0,xx.min()-55),min(im.width,xx.max()+55))
    ax.set_ylim(min(im.height,yy.max()+55),max(0,yy.min()-55))
    system=name.split('_')[0]
    meta=json.loads((ASSETS/(system+'_render_provenance.json')).read_text())
    asset=next(a for a in meta['assets'] if a['name']==name)
    landmarks=[l for l in asset['landmarks'] if l['object']=='holo']
    if local:
        # Positions are fixed in the camera frame; leader lines point to verified
        # deposited CA atoms. No hand-drawn protein geometry is added.
        anchors=([(.04,.14),(.04,.51),(.04,.91),(.72,.66)] if system=='b423' else [(.10,.90),(.78,.18)])
        for point, xytext in zip(landmarks,anchors):
            ax.annotate(point['label'],xy=point['pixel_xy_top_left'],xycoords='data',xytext=xytext,textcoords='axes fraction',
                fontsize=7,color=SITE,ha='left',va='center',arrowprops=dict(arrowstyle='-',lw=.6,color=SITE),
                bbox=dict(boxstyle='round,pad=.12',fc='white',ec='none',alpha=.9))
    else:
        xy=np.array([p['pixel_xy_top_left'] for p in landmarks])
        if len(xy):
            bounds=np.r_[xy.min(0)-70,xy.max(0)+70]
            ax.add_patch(Rectangle(bounds[:2],*(bounds[2:]-bounds[:2]),fill=False,ec=SITE,lw=.8,ls='--'))
            ax.annotate('Local view →',xy=xy.mean(0),xytext=(.03,.1),textcoords='axes fraction',fontsize=6.5,color=SITE,
                arrowprops=dict(arrowstyle='-',color=SITE,lw=.6),bbox=dict(fc='white',ec='none',pad=.2))
    return ax


def overlay_legend(fig,y,x=.03):
    handles=[Line2D([0],[0],color=APO,lw=2),Line2D([0],[0],color=HOLO,lw=2),Line2D([0],[0],color='#C3B69F',lw=2)]
    fig.legend(handles,['Apo TCR','Holo TCR','pMHC'],loc='lower left',bbox_to_anchor=(x,y),ncol=3,
        fontsize=6.5,handlelength=1.2,columnspacing=1.0,borderaxespad=0)


def forest(ax,rd,contrasts,regions,metric,show_y=True,mark_fold3=False):
    offsets=np.linspace(-.28,.28,len(contrasts))
    for j,cid in enumerate(contrasts):
        sub=rd[(rd.contrast==cid)&(rd.metric==metric)]
        for i,region in enumerate(regions):
            q=sub[sub.region==region].sort_values('fold')
            assert len(q)==5,(cid,region)
            vals=q.mean_delta.to_numpy(); y=i+offsets[j]
            ax.scatter(vals,y+np.linspace(-.016,.016,5),s=7,c=COLORS[j],alpha=.35,linewidths=0,zorder=2)
            ax.errorbar(vals.mean(),y,xerr=vals.std(ddof=1),fmt=MARKERS[j],color=COLORS[j],ms=3.4,lw=.8,capsize=1.6,zorder=3)
            if mark_fold3:
                v=q.loc[q.fold==3,'mean_delta'].iloc[0]
                ax.plot(v,y,'D',mfc='white',mec=COLORS[j],mew=.7,ms=3.1,zorder=4)
    ax.axvline(0,color='#888888',lw=.65,ls='--',zorder=0)
    ax.set_yticks(range(len(regions)))
    ax.set_yticklabels([r.replace('H3-associated 127-138','H3 window').replace('Other distal beta','Other distal β').replace('Lock-site windows','Lock-site windows').replace('Other distal noncontact','Distal noncontact') for r in regions] if show_y else [])
    ax.set_ylim(len(regions)-.5,-.5)
    ax.set_xlabel('Predicted ΔRMSF (Å)' if metric=='rmsf' else 'Predicted ΔHχ (nats)')
    ax.tick_params(axis='y',length=0)
    ax.spines['left'].set_visible(False)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4))


def contrast_legend(fig,labels,y,x=.16,ncol=None):
    handles=[Line2D([0],[0],marker=MARKERS[j],color=COLORS[j],lw=.9,ms=3) for j in range(len(labels))]
    fig.legend(handles,labels,loc='lower left',bbox_to_anchor=(x,y),ncol=ncol or len(labels),
        fontsize=6.5,columnspacing=1.2,handlelength=1.4,borderaxespad=0)


def profile(ax,df,contrasts,chain,metric,band=True,labels=None):
    sub=df[df.common_all & (df.chain==chain)]
    if metric=='hchi':
        sub=sub[sub.hchi_eligible]
    for j,cid in enumerate(contrasts):
        g=sub[sub.contrast==cid].groupby(['author_number','author_id'])['delta_'+metric].agg(['mean','std']).reset_index().sort_values('author_number')
        x=g.author_number.to_numpy();ym=g['mean'].to_numpy();sd=g['std'].to_numpy()
        # Break curves at missing author numbers or omitted A/G/P, never join
        # across an unresolved stretch as though it were predicted.
        splits=np.r_[0,np.where(np.diff(x)>1)[0]+1,len(x)]
        for a,b in zip(splits[:-1],splits[1:]):
            ax.plot(x[a:b],ym[a:b],color=COLORS[j],ls=STYLES[j],lw=.85,
                marker='.' if b-a==1 else None,ms=1.5)
            if band:
                ax.fill_between(x[a:b],ym[a:b]-sd[a:b],ym[a:b]+sd[a:b],color=COLORS[j],alpha=.09,lw=0)
    ax.axhline(0,color='#777777',lw=.6,ls='--')
    if chain=='F':
        ax.axvspan(126.5,138.5,color=SITE,alpha=.08,lw=0)
    if chain=='A':
        for a,b in [(98.5,109.5),(176.5,187.5)]:
            ax.axvspan(a,b,color=SITE,alpha=.08,lw=0)
    ax.set_ylabel('ΔRMSF (Å)' if metric=='rmsf' else 'ΔHχ (nats)')
    ax.margins(x=.005)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(3))


def b4_main(df,rd,geometry):
    fig=plt.figure(figsize=(183/25.4,210/25.4))
    title(fig,'B4.2.3 | testing a distal signaling site','Real apo/holo coordinates · frozen FlexPPI checkpoints · literature-defined region')
    heading(fig,'a','Experimental states',.02,.924)
    structure(fig,[.01,.665,.285,.25],'b423_experimental_overlay')
    heading(fig,'b','Local Cα/Cβ interface',.32,.924)
    structure(fig,[.30,.69,.34,.225],'b423_distal_local',True)
    overlay_legend(fig,.645)
    fig.text(.325,.676,'5IW1 (AB) → 5IVX; one whole-TCR fit',fontsize=6.2,color='#555555')
    heading(fig,'c','Deposited contacts',.70,.924)
    ax=fig.add_axes([.75,.72,.22,.172])
    pairs=['E121alpha-K134beta charged atoms','Y125alpha-N133beta polar atoms','Y125alpha-K134beta polar atoms']
    order=['b4_5iw1_ab_full','b4_5iw1_cd_full','b4_5iw1_ef_full','b4_5ivx_full_full']
    for i,pair in enumerate(pairs):
        q=geometry[geometry.pair==pair].set_index('state').loc[order]
        ax.plot(range(4),q.distance_A,color=COLORS[i],marker=MARKERS[i],ms=3,lw=.9,ls=STYLES[i])
    ax.set_xticks(range(4),['AB','CD','EF','Holo']);ax.set_ylabel('Atom distance (Å)');ax.set_ylim(2.7,5.5)
    ax.legend(['E121α–K134β','Y125α–N133β','Y125α–K134β'],fontsize=6.5,loc='upper left',bbox_to_anchor=(-.23,-.20),handlelength=1.8)
    heading(fig,'d','Functional sites versus crystal-copy controls',.02,.624)
    contrasts=BIO_B4+['b4_apo_ab_vs_cd_full','b4_apo_ab_vs_ef_full']
    contrast_legend(fig,['AB→holo','CD→holo','EF→holo','AB→CD','AB→EF'],.600,x=.13)
    regions=['N133','K134','H3-associated 127-138','Other distal beta']
    forest(fig.add_axes([.16,.410,.32,.170]),rd,contrasts,regions,'rmsf',mark_fold3=True)
    forest(fig.add_axes([.65,.410,.32,.170]),rd,contrasts,regions,'hchi',show_y=False,mark_fold3=True)
    fig.text(.13,.343,'Small points: five checkpoints; symbols ± SD: mean; open diamonds: fold 3 (5IVX held out).',fontsize=6.6,color='#555555')
    heading(fig,'e','Whole β-chain predictions; H3 window shaded',.02,.312)
    ax1=fig.add_axes([.105,.203,.865,.091]);ax2=fig.add_axes([.105,.075,.865,.091])
    profile(ax1,df,BIO_B4,'F','rmsf');profile(ax2,df,BIO_B4,'F','hchi')
    ax1.tick_params(labelbottom=False);ax2.set_xlabel('β-chain author residue number')
    contrast_legend(fig,['AB→holo','CD→holo','EF→holo'],.304,x=.57)
    fig.text(.025,.018,'Predictions, not measured dynamics. 5IVX occurs in ITrace; this is not independent external validation.',fontsize=6.8,color='#555555')
    export(fig,'fig01_b423_real_states','Local contact reorganization is resolved; model changes at literature sites require crystal-copy and training-overlap controls.')


def assay_panel(ax):
    p=ROOT/'artifacts/tcr_pmhc_reassessment_20260915/literature/notti2025_functional_source_audit.json'
    audit=json.loads(p.read_text())
    raw=[]
    for i,item in enumerate(audit['functional_checks']):
        for j,key in enumerate(['WT','mutant_S104C_V182C']):
            vals=np.array(item[key],dtype=float)/item['WT_mean']*100
            x=i+(j-.5)*.22
            ax.scatter(x+np.linspace(-.025,.025,3),vals,s=11,color=[APO,SITE][j],marker=['o','s'][j],zorder=3)
            ax.errorbar(x,vals.mean(),yerr=vals.std(ddof=1),color=[APO,SITE][j],lw=.8,capsize=2,fmt='none',zorder=2)
            for k,v in enumerate(vals):
                raw.append(dict(assay=item['assay'],genotype=key,observation=k+1,raw_value=item[key][k],percent_assay_WT_mean=float(v)))
        ax.plot([i-.11,i+.11],[100,item['mutant_mean']/item['WT_mean']*100],color='#BFBFBF',lw=.8,zorder=1)
    ax.set_xticks(range(3),['Tetramer','Peptide / APC','PMA / iono.'])
    ax.set_ylabel('Response (% of assay WT mean)');ax.set_ylim(0,125)
    ax.axhline(100,color='#C5C5C5',lw=.6,ls='--')
    ax.legend([Line2D([0],[0],marker='o',color=APO,lw=0,ms=3),Line2D([0],[0],marker='s',color=SITE,lw=0,ms=3)],
        ['WT','S104C/V182C'],fontsize=6.5,loc='lower left',bbox_to_anchor=(0,1.02),ncol=2,borderaxespad=0)
    pd.DataFrame(raw).to_csv(OUT/'plotted_literature_functional_source.csv',index=False)


def coverage_panel(ax,df,cid,metric='rmsf'):
    sub=df[df.common_all & df.contrast.isin([cid,cid[:-5]+'_common'])]
    table=sub.pivot_table(index=['chain','author_id','regions'],columns='contrast',values='delta_'+metric,aggfunc='mean')
    x=table[cid];y=table[cid[:-5]+'_common']
    chosen=np.array(['Lock-site windows' in r for r in table.index.get_level_values('regions')])
    ax.scatter(x[~chosen],y[~chosen],s=6,c='#B9C3CC',alpha=.55,linewidths=0)
    ax.scatter(x[chosen],y[chosen],s=13,c=SITE,alpha=.8,linewidths=.2,edgecolors='white')
    bounds=[min(x.min(),y.min(),0)-.08,max(x.max(),y.max(),0)+.08]
    ax.plot(bounds,bounds,color='#777777',ls='--',lw=.7)
    ax.set_xlim(bounds);ax.set_ylim(bounds)
    ax.set_xlabel('Full-context ΔRMSF (Å)');ax.set_ylabel('Common-context ΔRMSF (Å)')
    ax.set_aspect('equal',adjustable='box')
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4));ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(4))


def ig4_main(df,rd,geometry):
    fig=plt.figure(figsize=(183/25.4,205/25.4))
    title(fig,'1G4–CD3 | separating conformation from antigen specificity','Real nanodisc states · complete observed protein context · explicit environment controls')
    heading(fig,'a','Real receptor states',.02,.924)
    structure(fig,[.006,.659,.30,.257],'1g4_experimental_overlay')
    heading(fig,'b','WT lock-site neighborhood',.33,.924)
    structure(fig,[.315,.70,.325,.21],'1g4_distal_local',True)
    overlay_legend(fig,.642)
    fig.text(.328,.684,'S104 / V182 shown in WT; no mutant built.',fontsize=6.4,color='#555555')
    heading(fig,'c','Lock-site separation',.71,.924)
    ax=fig.add_axes([.77,.73,.20,.16])
    q=geometry[geometry.pair=='S104-V182 CA distance'].set_index('state').loc[[f'1g4_{p}_full_full' for p in ['8tw4','8tw6','9bbc','9c3e']]]
    ax.plot(range(4),q.distance_A,color=APO,lw=.9)
    for i,v in enumerate(q.distance_A):
        ax.scatter(i,v,color=[COLORS[0],COLORS[1],GREY,HOLO][i],marker=MARKERS[i],s=20,zorder=3)
    ax.set_xticks(range(4),['ND-I','ND-II','GDN','Holo'],rotation=25,ha='right');ax.set_ylabel('WT Cα–Cα distance (Å)');ax.set_ylim(5.5,12)
    fig.text(.72,.641,'ND: nanodisc apo\nGDN: detergent apo\nCD3 is shown in light grey in a.',fontsize=6.6,color='#555555',linespacing=1.5)
    heading(fig,'d','Predicted regional changes and unbound-state controls',.02,.615)
    cs=[BIO_1G4[0],BIO_1G4[1],'1g4_8tw4_vs_detergent_full',BIO_1G4[2]]
    contrast_legend(fig,['ND-I→holo','ND-II→holo','ND-I→GDN','GDN→holo'],.590,x=.16)
    regions=['pMHC interface','Lock-site windows','Distal CD3 contact','Other distal noncontact']
    forest(fig.add_axes([.19,.407,.30,.162]),rd,cs,regions,'rmsf')
    forest(fig.add_axes([.67,.407,.30,.162]),rd,cs,regions,'hchi',show_y=False)
    fig.text(.16,.340,'Identical residue masks; points: five checkpoints; summary symbols and error bars: mean ± SD.',fontsize=6.6,color='#555555')
    heading(fig,'e','Coverage sensitivity',.02,.305)
    ax=fig.add_axes([.12,.115,.27,.161]);coverage_panel(ax,df,BIO_1G4[0])
    fig.text(.10,.048,'ND-I→holo; all 366 matched sites\nPurple: lock-site windows; line: identity.',fontsize=6.4,color='#555555')
    heading(fig,'f','Published functional data (not prediction)',.48,.305)
    ax=fig.add_axes([.585,.098,.385,.165]);assay_panel(ax)
    fig.text(.575,.047,'Source Data: Notti et al. (2025); n = 3 per group.\nNormalization uses each assay’s WT mean; raw data supplied.',fontsize=6.2,color='#555555')
    fig.text(.025,.012,'State association is not proof of antigen-specific dynamics. Missing membrane/coverage and training homology remain limitations.',fontsize=6.7,color='#555555')
    export(fig,'fig02_1g4_real_states','Real structural opening and functional constraint coexist, but model signal is not specific to antigen and is coverage sensitive.')


def full_profiles(df,case):
    b4=case=='B4.2.3';bio=BIO_B4 if b4 else BIO_1G4
    chains=['E','F'] if b4 else ['A','B']
    fig,axs=plt.subplots(4,2,figsize=(183/25.4,227/25.4),gridspec_kw=dict(wspace=.3,hspace=.62))
    fig.subplots_adjust(left=.095,right=.975,top=.885,bottom=.086)
    title(fig,case+' | complete matched-residue profiles','Left: full observed context · right: common-coverage sensitivity · same residue set in both')
    for col,suffix in enumerate(['_full','_common']):
        cs=[c[:-5]+suffix for c in bio]
        for row,(ch,metric) in enumerate([(chains[0],'rmsf'),(chains[0],'hchi'),(chains[1],'rmsf'),(chains[1],'hchi')]):
            profile(axs[row,col],df,cs,ch,metric)
            axs[row,col].set_title(('α' if ch==chains[0] else 'β')+' chain · '+('full context' if col==0 else 'common context'),loc='left')
            axs[row,col].set_xlabel('Author residue number')
            axs[row,col].text(-.12,1.08,chr(97+row*2+col),transform=axs[row,col].transAxes,fontweight='bold',fontsize=9)
    # Corresponding physical metrics share y limits across full/common inputs.
    for row in range(4):
        lo=min(axs[row,j].get_ylim()[0] for j in range(2));hi=max(axs[row,j].get_ylim()[1] for j in range(2))
        for j in range(2):axs[row,j].set_ylim(lo,hi)
    contrast_legend(fig,['AB→holo','CD→holo','EF→holo'] if b4 else ['ND-I→holo','ND-II→holo','GDN→holo'],.913,x=.24)
    fig.text(.095,.012,'Lines: five-checkpoint means; bands: ± SD, not confidence intervals. Gaps are not interpolated. Hχ omits A/G/P.\nAll matched residues are retained; purple bands mark predeclared literature-site windows.',fontsize=6.6,color='#555555')
    export(fig,'figS01_b423_complete_profiles' if b4 else 'figS02_1g4_complete_profiles','Complete profiles including missing-residue gaps and the same coverage controls.')


def structural_maps():
    # Slightly narrower canvas keeps the tight PDF bounding box below 183 mm.
    fig=plt.figure(figsize=(181/25.4,119/25.4))
    title(fig,'Where the model predicts a change','Mean holo-minus-apo RMSF; common TCR positions only; unmatched positions and partners in grey')
    for i,(name,label) in enumerate([('b423','B4.2.3: AB→holo'),('1g4','1G4: ND-I→holo')]):
        heading(fig,chr(97+i),label,.025+i*.5,.88)
        ax=fig.add_axes([.01+i*.50,.13,.48,.73]);ax.imshow(Image.open(ASSETS/(name+'_delta_rmsf.png')));ax.set_axis_off()
    cmap=mpl.colors.LinearSegmentedColormap.from_list('physical_delta',DIVERGING)
    ax=fig.add_axes([.36,.065,.28,.023]);fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1,1),cmap=cmap),cax=ax,orientation='horizontal',extend='both',ticks=[-1,0,1])
    ax.set_xlabel('Predicted ΔRMSF (Å); color saturated beyond ±1 Å',fontsize=7)
    export(fig,'figS03_predicted_structure_maps','Spatial model-output maps; not observed motion or an independent validation.')


def existing_md_figure():
    d=pd.read_csv(OUT/'b423_existing_holo_md_comparison.csv',dtype={'author_id':str})
    d=d[d.valid & d.chain.isin(['E','F'])].copy()
    d=d[(d.metric!='hchi')|(~d.amino_acid.isin(list('AGP')))]
    d['author_number']=d.author_id.str.extract(r'(-?\d+)').astype(int)
    fig,axs=plt.subplots(3,2,figsize=(183/25.4,205/25.4),gridspec_kw=dict(hspace=.68,wspace=.34))
    fig.subplots_adjust(left=.11,right=.975,top=.87,bottom=.095)
    title(fig,'B4.2.3 | held-out holo prediction versus existing MD','Frozen fold 3 · deposited PDB input · absolute flexibility, not the apo→holo change')
    for row,metric in enumerate(['rmsf','hchi']):
        for col,ch in enumerate(['E','F']):
            ax=axs[row,col];q=d[d.metric.eq(metric)&d.chain.eq(ch)].sort_values('author_number')
            x=q.author_number.to_numpy()
            cuts=np.r_[0,np.where(np.diff(x)>1)[0]+1,len(x)]
            for field,color,style in [('md_reference','#333333','-'),('predicted_deposited_fold3',COLORS[0],'--')]:
                vals=q[field].to_numpy()
                for a,b in zip(cuts[:-1],cuts[1:]):
                    ax.plot(x[a:b],vals[a:b],color=color,ls=style,lw=.8,marker='.' if b-a==1 else None,ms=1.5)
            if ch=='F':ax.axvspan(126.5,138.5,color=SITE,alpha=.08,lw=0)
            ax.set_title(('α' if ch=='E' else 'β')+' chain',loc='left')
            ax.set_xlabel('Author residue number');ax.set_ylabel('RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)')
            ax.margins(x=.01);ax.set_ylim(bottom=0)
    metrics=pd.read_csv(OUT/'b423_holo_md_metrics.csv')
    for col,metric in enumerate(['rmsf','hchi']):
        ax=axs[2,col];q=d[d.metric.eq(metric)]
        selected=q.chain.eq('F')&q.author_id.isin(['133','134'])
        ax.scatter(q.loc[~selected,'md_reference'],q.loc[~selected,'predicted_deposited_fold3'],s=7,color=COLORS[0],alpha=.30,linewidths=0)
        ax.scatter(q.loc[selected,'md_reference'],q.loc[selected,'predicted_deposited_fold3'],s=22,color=SITE,edgecolors='white',linewidths=.5,zorder=4)
        for j,(_,r) in enumerate(q[selected].iterrows()):
            ax.annotate(r.amino_acid+r.author_id,xy=(r.md_reference,r.predicted_deposited_fold3),xytext=(12,-12 if j==0 else 12),
                textcoords='offset points',fontsize=6.5,color=SITE,arrowprops=dict(arrowstyle='-',lw=.5,color=SITE))
        maxv=max(q.md_reference.max(),q.predicted_deposited_fold3.max())*1.05
        ax.plot([0,maxv],[0,maxv],color='#888888',lw=.7,ls='--');ax.set_xlim(0,maxv);ax.set_ylim(0,maxv)
        ax.set_xlabel('MD reference '+('RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)'))
        ax.set_ylabel('Predicted '+('RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)'))
        m=metrics[(metrics.scope=='TCR alpha beta')&(metrics.metric==metric)&(metrics.predictor=='predicted_deposited_fold3')].iloc[0]
        ax.text(.05,.96,f'Spearman ρ = {m.spearman:.3f}\nMAE = {m.mae:.3f} '+('Å' if metric=='rmsf' else 'nats')+f'\nn = {m.n_residues} residues',
            transform=ax.transAxes,va='top',fontsize=7,linespacing=1.4)
        ax.set_aspect('equal',adjustable='box')
        ax.set_anchor('W')
    for row in range(2):
        high=max(axs[row,c].get_ylim()[1] for c in range(2))
        for c in range(2):axs[row,c].set_ylim(0,high)
    for i,ax in enumerate(axs.flat):ax.text(-.15,1.08,chr(97+i),transform=ax.transAxes,fontweight='bold',fontsize=9)
    fig.legend([Line2D([0],[0],color='#333333',lw=1),Line2D([0],[0],color=COLORS[0],lw=1,ls='--')],
        ['Existing MD aggregate','FlexPPI: held-out fold 3'],loc='lower left',bbox_to_anchor=(.24,.905),ncol=2,fontsize=7)
    fig.text(.035,.02,'ITrace 5IVX: three-replica label protocol; no new MD was run. Same-AA-centered Hχ correlation = 0.844 (descriptive).\nCombination-level holdout, not unseen-family validation. No matched apo MD reference is available in this comparison.',fontsize=6.6,color='#555555')
    export(fig,'figS04_b423_heldout_holo_MD','Fold-3 absolute holo-state flexibility agrees with existing MD aggregates; this is not validation of the binding-induced delta.')


def main():
    df=pd.read_csv(OUT/'residue_predictions.csv',dtype={'author_id':str})
    rd=pd.read_csv(OUT/'region_fold_statistics.csv')
    geometry=pd.read_csv(OUT/'local_experimental_geometry.csv')
    b4_main(df,rd,geometry);ig4_main(df,rd,geometry)
    full_profiles(df,'B4.2.3');full_profiles(df,'1G4-CD3');structural_maps();existing_md_figure()
    (OUT/'figure_export_audit.json').write_text(json.dumps(dict(figures=EXPORT_AUDIT,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_sha256={f:hashlib.sha256((OUT/f).read_bytes()).hexdigest() for f in ['residue_predictions.csv','region_fold_statistics.csv','local_experimental_geometry.csv']},
        asset_route='Parameter inheritance; transformations incompatible with scientific variables were not reused.'),indent=2),encoding='utf-8')
    print('Exported six figures in PDF, PNG and SVG.',flush=True)


if __name__=='__main__':
    main()
