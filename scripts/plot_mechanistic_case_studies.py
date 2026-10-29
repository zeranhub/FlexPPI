# Asset confirmation: LineTrend/plot_trend.py + trend_by_month.png inspected;
# heatmap/plot_comparison.py + results_comparison_optimization.png inspected.
# Cross-type inheritance: thin dark lines/light SD ribbons, aligned cells.
# NOT inherited: cumulative y variables, per-column normalization, method ranking.
# Scientific overrides: shared zero-centered signed scales; genuine PyMOL raster
# structures; numerical panels/text remain vector. No correlations are invented.
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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/mechanistic_case_studies_20260913'
OUT=ROOT/'output/flexppi_mechanistic_case_studies_20260913'
DATA=pd.read_csv(BASE/'residue_predictions.csv',dtype={'author_id':str,'before_author':str})
REG=pd.read_csv(BASE/'region_predictions.csv')
CMAP=LinearSegmentedColormap.from_list('flexppi_delta',DIVERGING)
BLUE,RED,GREEN,ORANGE,PURPLE,GRAY=CATEGORICAL
HCHI=r'$\Delta H_\chi$ (nats)'
RMSF=r'$\Delta$RMSF ($\AA$)'
MANIFEST=[]

def panel(ax,letter,title):
    ax.text(-.09,1.07,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
    ax.set_title(title,loc='left',pad=8,fontweight='bold')

def image_panel(fig,rect,name,kind,letter,title):
    rect=[rect[0],.605,rect[2],.265]
    ax=fig.add_axes(rect)
    img=plt.imread(OUT/'structure_assets'/f'{name}_{kind}.png')
    # Crop only external white canvas, never protein content. Keeps structure
    # raster native; this is a layout operation, not biological image retouching.
    mask=np.any(img[:,:,:3]<.975,axis=2)
    yy,xx=np.where(mask)
    pad=35
    img=img[max(0,yy.min()-pad):min(img.shape[0],yy.max()+pad+1),max(0,xx.min()-pad):min(img.shape[1],xx.max()+pad+1)]
    ax.imshow(img);ax.axis('off');panel(ax,letter,title)
    if name=='kix' and kind=='whole':
        ax.text(.05,.01,'cMyb',color=PURPLE,transform=ax.transAxes,fontsize=8,weight='bold')
        ax.text(.78,.01,'MLL',color='#9A7B43',transform=ax.transAxes,fontsize=8,weight='bold')
    if name=='pdz3' and kind=='whole':
        ax.text(.48,.01,'Orange: deleted 396–402',color='#AA6B20',transform=ax.transAxes,fontsize=7)
    if name.startswith('tcr') and kind=='whole':
        ax.text(.05,.01,'Purple: selected NMR-perturbed sites',color=PURPLE,transform=ax.transAxes,fontsize=7)
    return ax

def canvas(title,subtitle):
    fig=plt.figure(figsize=(183/25.4,185/25.4))
    fig.text(.07,.975,title,fontsize=12,fontweight='bold',va='top',color=BLACK)
    fig.text(.07,.938,subtitle,fontsize=8,color='#555555',va='top')
    return fig

def footer(fig,text):
    fig.text(.07,.055,text,fontsize=7,va='top',linespacing=1.35,color='#555555')

def colorbar(fig,metric,rect=(.14,.562,.68,.012)):
    lim=1 if metric=='rmsf' else .2
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-lim,lim),cmap=CMAP),cax=fig.add_axes(rect),orientation='horizontal',ticks=[-lim,0,lim])
    cb.outline.set_visible(False)
    cb.set_label((RMSF if metric=='rmsf' else HCHI)+'   |   blue: lower flexibility; red: higher',fontsize=7,labelpad=2)
    cb.ax.tick_params(length=2,pad=1)

def profile(ax,contrast,metric,chain=None,label=None,color=BLUE,style='-',shade=True):
    df=DATA[DATA.contrast==contrast].copy()
    if chain:
        df=df[df.chain==chain]
    if metric=='chi_entropy':
        df.loc[df.aa.isin(['A','G','P']),'delta_'+metric]=np.nan
    agg=df.groupby(['chain','author_id','author_number'],sort=False)['delta_'+metric].agg(['mean','std']).reset_index()
    # Explicit observed-residue index preserves insertion-coded thrombin60A-I.
    x=np.arange(len(agg));y=agg['mean'].values;s=agg['std'].values
    ax.plot(x,y,color=color,lw=1.1,ls=style,label=label)
    if shade:
        ax.fill_between(x,y-s,y+s,color=color,alpha=.14,lw=0)
    ax.axhline(0,color='#999999',lw=.6,ls='--',zorder=0)
    step=max(1,len(agg)//5)
    ticks=list(range(0,len(agg),step))
    ax.set_xticks(ticks,[agg.author_id.iloc[i] for i in ticks])
    ax.set_xlim(-1,len(agg))
    ax.set_ylabel(RMSF if metric=='rmsf' else HCHI)
    ax.set_xlabel('PDB residue (observed positions)')
    return agg

def shade_region(ax,agg,low,high,color=ORANGE,label=None):
    ix=np.flatnonzero((agg.author_number>=low)&(agg.author_number<=high))
    if len(ix):
        ax.axvspan(ix.min()-.4,ix.max()+.4,color=color,alpha=.14,lw=0,zorder=-1)
        if label:
            ax.text(np.mean(ix),.99,label,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=7)

def dots(ax,groups,labels,colors=None,ylabel=RMSF):
    for i,values in enumerate(groups):
        v=np.asarray(values,dtype=float);v=v[np.isfinite(v)]
        c=(colors or [BLUE]*len(groups))[i]
        ax.scatter(i+np.linspace(-.10,.10,len(v)),v,s=15,color=c,alpha=.7,zorder=3)
        ax.errorbar(i,v.mean(),yerr=v.std(ddof=1),fmt='D',ms=4,color=c,elinewidth=1,capsize=3,zorder=5,mec='white',mew=.4)
    ax.axhline(0,color=GREY,lw=.6,ls='--')
    ax.set_xticks(range(len(labels)),labels);ax.set_ylabel(ylabel)
    ax.set_xlim(-.5,len(labels)-.5)

def region_values(contrast,region,metric):
    return REG[(REG.contrast==contrast)&(REG.region==region)&(REG.metric==metric)].sort_values('fold').delta.to_numpy()

def save(fig,name,caption):
    save_cns_figure(fig,OUT/name)
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight',dpi=300)
    MANIFEST.append(dict(name=name,caption=caption,pdf=str(OUT/(name+'.pdf')),png=str(OUT/(name+'.png')),svg=str(OUT/(name+'.svg'))))
    plt.close(fig)

def pdz():
    fig=canvas('PDZ3 | deletion response without full dynamic rescue','Four-state perturbation screen • 5HEB coordinates • five frozen checkpoints')
    image_panel(fig,[.055,.62,.46,.28],'pdz3','whole','a','Where the helix was removed')
    image_panel(fig,[.58,.62,.37,.28],'pdz3','local','b','Local packing around the deletion')
    colorbar(fig,'chi_entropy')
    ax=fig.add_axes([.10,.30,.43,.18]);panel(ax,'c','Side-chain response along the core')
    ag=profile(ax,'pdz_deletion_free','chi_entropy',label='Free',color=RED)
    profile(ax,'pdz_deletion_bound','chi_entropy',label='Bound',color=BLUE,style='--',shade=False)
    shade_region(ax,ag,389,395,label='near deletion');ax.legend(loc='upper left',fontsize=7,ncol=2)
    ax=fig.add_axes([.66,.30,.29,.18]);panel(ax,'d','Deletion effect across states')
    vals=[region_values(c,'shared core','chi_entropy') for c in ['pdz_deletion_free','pdz_deletion_bound','pdz_deletion_apo']]
    dots(ax,vals,['Free','Bound','Apo\n1BFE'],[RED,BLUE,GRAY],HCHI)
    ax=fig.add_axes([.10,.115,.85,.095]);panel(ax,'e','Fold-resolved deletion response: localized, not a uniform shift')
    raw=DATA[DATA.contrast=='pdz_deletion_free'].copy()
    raw.loc[raw.aa.isin(['A','G','P']),'delta_chi_entropy']=np.nan
    columns=raw[raw.fold==0].author_id.tolist()
    mat=raw.pivot(index='fold',columns='author_id',values='delta_chi_entropy')[columns]
    ax.imshow(mat,cmap=CMAP,vmin=-.2,vmax=.2,aspect='auto',interpolation='none')
    ax.set_yticks(range(5),range(1,6));ax.set_ylabel('Checkpoint')
    ticks=list(range(0,len(columns),18));ax.set_xticks(ticks,[columns[i] for i in ticks]);ax.set_xlabel('PDB residue; white gaps: no chi summary (A/G/P)')
    footer(fig,'Dots: 5 checkpoints; diamonds ± SD. Common residues only; Hχ excludes A/G/P. No biological replicates.\nObserved CRIPT7/9; deletion and stripped states are modeled controls. Homology overlap flagged; not OOF validation.')
    save(fig,'fig01_pdz3_deletion','PDZ3 fixed-coordinate deletion screen. Orange helix396–402 is shown as a positional reference, not retained in deletion inference. Shared-core Hchi rises with deletion in free and bound states; the experimental rescue criterion is not reproduced. 1BFE I328V sequence-matched apo sensitivity has missing303–305. NMR summary is qualitative; no S2 values were invented.')

def tcr():
    fig=canvas('A6 TCR | separate partner effects from crystal differences','3QH3 → sequence-matched 4FTV • exact matched residues • not native wild-type coordinates')
    image_panel(fig,[.055,.62,.43,.28],'tcr','whole','a','Independent apo-to-holo difference')
    image_panel(fig,[.56,.62,.40,.28],'tcr_context','whole','b','Partner-only effect, same backbone')
    colorbar(fig,'rmsf')
    ax=fig.add_axes([.10,.30,.43,.18]);panel(ax,'c','Beta-chain response decomposition')
    for c,l,col,sty in [('tcr_apo_holo','Total',RED,'-'),('tcr_geometry','Backbone/construct',GRAY,'--'),('tcr_context','Partner',BLUE,':')]:
        ag=profile(ax,c,'rmsf','E',l,col,sty,shade=False)
    shade_region(ax,ag,218,231,label='218–231')
    sites=json.loads((BASE/'tcr/selected_nmr_site_mapping.json').read_text())
    for site in sites:
        ix=np.flatnonzero(ag.author_id==site.get('pdb_author_id'))
        if len(ix):
            ax.plot(ix,[.96]*len(ix),'|',color=PURPLE,ms=7,transform=ax.get_xaxis_transform())
    ax.legend(loc='lower left',fontsize=7)
    ax=fig.add_axes([.66,.30,.29,.18]);panel(ax,'d','Constant-domain response')
    labels=['Cα\ntotal','Cβ\ntotal','Cα\npartner','Cβ\npartner']
    groups=[region_values(c,r,'rmsf') for c,r in [('tcr_apo_holo','C alpha'),('tcr_apo_holo','C beta'),('tcr_context','C alpha'),('tcr_context','C beta')]]
    dots(ax,groups,labels,[RED,RED,BLUE,BLUE])
    ax=fig.add_axes([.10,.115,.85,.095]);panel(ax,'e','Partner-only response is concentrated near the binding interface')
    rr=DATA[DATA.contrast=='tcr_context'].groupby(['chain','author_id'],sort=False).agg(distance=('partner_distance','mean'),delta=('delta_rmsf','mean')).reset_index()
    for ch,col,lab,marker in [('D',BLUE,'Alpha','o'),('E',PURPLE,'Beta','^')]:
        r=rr[rr.chain==ch];ax.scatter(r.distance,r.delta,c=col,s=8,alpha=.65,label=lab,marker=marker)
    ax.axhline(0,color=GREY,lw=.6,ls='--');ax.set_ylabel(RMSF);ax.set_xlabel('Minimum heavy-atom distance to pMHC ($\AA$)')
    ax.legend(loc='lower right',ncol=2,fontsize=7)
    footer(fig,'Dots n=5 checkpoints; diamonds ± SD. Purple markers: six sequence-verified NMR sites, not a complete mask.\nBeta 218–231 remains an exploratory window. NMR perturbation is not signed RMSF. Exact training overlap detected.')
    save(fig,'fig02_tcr_context_decomposition','TCR total change decomposes into independent stripped-crystal difference plus partner context on the same exactly matched residues. Constant-domain changes are essentially absent for the partner-only intervention; no claim of successful reproduction of distal allostery. Engineered cysteine background and holo CDR backbone remain caveats; beta218–231 is not an experimentally validated mask.')

def thrombin():
    fig=canvas('Thrombin | adding TM4 does not recover the proposed loop response','TM56 versus TM456 • same 1DX5 backbone • sequence reversions recorded')
    image_panel(fig,[.055,.62,.46,.28],'thrombin','whole','a','TM4 increment mapped onto thrombin')
    image_panel(fig,[.59,.62,.37,.28],'thrombin','local','b','60s-loop neighborhood')
    colorbar(fig,'rmsf')
    ax=fig.add_axes([.10,.30,.43,.18]);panel(ax,'c','Heavy-chain response')
    ag=profile(ax,'th_tm456_effect','rmsf','M','TM456',BLUE,shade=True)
    profile(ax,'th_tm56_effect','rmsf','M','TM56',PURPLE,'--',shade=False)
    profile(ax,'th_tm4_increment','rmsf','M','TM4 increment',ORANGE,':',shade=False)
    shade_region(ax,ag,60,64,label='60s');ax.legend(loc='lower right',fontsize=7)
    ax=fig.add_axes([.66,.30,.29,.18]);panel(ax,'d','60s loop: no positive increment')
    dots(ax,[region_values(c,'60s loop','rmsf') for c in ['th_tm56_effect','th_tm456_effect','th_tm4_increment']],['TM56','TM456','TM4\nincrement'],[PURPLE,BLUE,ORANGE])
    ax=fig.add_axes([.10,.115,.85,.095]);panel(ax,'e','The additional TM4 domain has little effect across functional loops')
    labs=['60s loop','Exosite I','170s loop','220s loop']
    dots(ax,[region_values('th_tm4_increment',r,'rmsf') for r in labs],labs,[ORANGE]*4)
    footer(fig,'Dots n=5 checkpoints; diamonds ± SD. Loop60–64 includes insertion-coded residues60A–60I.\nInhibitor-conditioned crystal backbone; inhibitor/Ca chemistry is not encoded. Homology overlap flagged.')
    save(fig,'fig03_thrombin_tm4_control','Thrombin partner context from1DX5 A/M/I; TM56 author388–462 and TM456345–462. Reversion I60I to T follows SEQADV, rather than incorrectly changing author91. The additional TM4 domain gives little60s response and does not reproduce the selected literature direction. Protein-only screen, not a protein C activity prediction.')

def kix():
    fig=canvas('KIX | MLL preferentially stabilizes the L12 loop','All 20 deposited 2AGH conformers × 5 checkpoints • MLL 840–858 • paired context changes')
    image_panel(fig,[.055,.62,.45,.28],'kix','whole','a','Two partners, one flexible loop')
    image_panel(fig,[.58,.62,.38,.28],'kix','local','b','L12 loop 614–621: MLL response')
    colorbar(fig,'rmsf')
    ax=fig.add_axes([.10,.30,.43,.18]);panel(ax,'c','L12 response for all 20 conformers')
    k=REG[(REG.case=='kix')&(REG.region=='L12 loop')&(REG.metric=='rmsf')&~REG.contrast.str.contains('full31')].copy()
    k['partner']=k.contrast.str.extract(r'm\d+_(\w+)_effect')[0].values
    # Connect paired conformer means; overlay checkpoint means across conformers.
    p=k.groupby(['source_model','partner']).delta.mean().unstack()[['myb','mll','both']]
    for _,r in p.iterrows():
        ax.plot(range(3),r.to_numpy(),color=GRAY,alpha=.36,lw=.6,marker='o',ms=2)
    groups=[k[k.partner==p].groupby('fold').delta.mean().to_numpy() for p in ['myb','mll','both']]
    dots(ax,groups,['cMyb','MLL','Both'],[PURPLE,BLUE,BLACK]);ax.set_ylabel('L12 '+RMSF)
    ax=fig.add_axes([.66,.30,.29,.18]);panel(ax,'d','Side-chain response in the loop')
    kh=REG[(REG.case=='kix')&(REG.region=='L12 loop')&(REG.metric=='chi_entropy')&~REG.contrast.str.contains('full31')].copy()
    kh['partner']=kh.contrast.str.extract(r'm\d+_(\w+)_effect')[0].values
    dots(ax,[kh[kh.partner==p].groupby('fold').delta.mean().to_numpy() for p in ['myb','mll','both']],['cMyb','MLL','Both'],[PURPLE,BLUE,BLACK],HCHI)
    ax=fig.add_axes([.10,.115,.43,.10]);panel(ax,'e','Loop selectivity across the sequence')
    ag=profile(ax,'kix_m01_mll_effect','rmsf','B','MLL',BLUE,shade=False)
    profile(ax,'kix_m01_myb_effect','rmsf','B','cMyb',PURPLE,'--',shade=False)
    shade_region(ax,ag,614,621);ax.set_xlabel('PDB residue (conformer1)')
    ax.legend(loc='lower left',ncol=2,fontsize=7)
    ax=fig.add_axes([.66,.115,.29,.10]);panel(ax,'f','Non-additivity, not affinity')
    wide=k.pivot(index=['source_model','fold'],columns='partner',values='delta')
    interaction=wide['both']-wide['mll']-wide['myb']
    dots(ax,[interaction.groupby('fold').mean().values],['Both − MLL − cMyb'],[BLACK])
    footer(fig,'Thin gray lines: 20 paired conformer means. Colored points: 5 checkpoint means over all 20; diamonds ± SD.\nQualitative MD-pattern agreement only. Near-identical KIX training sequence detected; not independent validation.')
    save(fig,'fig04_kix_partner_selectivity','KIX all20 NMR conformer fixed-coordinate contexts. L12-loop614–621 responds more strongly to MLL than cMyb, consistent with the selected qualitative literature pattern. Model non-additivity is not thermodynamic cooperativity, nor proof of a dynamic communication pathway. Near-identical training overlap precludes calling the fivefold average OOF validation.')

def reserve(case,num,title,subtitle,contrast,metric,chain,low,high,localtitle):
    fig=canvas(title,subtitle)
    image_panel(fig,[.055,.62,.46,.28],case,'whole','a','Spatial response map')
    image_panel(fig,[.59,.62,.37,.28],case,'local','b',localtitle)
    colorbar(fig,metric)
    ax=fig.add_axes([.10,.30,.43,.18]);panel(ax,'c','Backbone response')
    ag=profile(ax,contrast,'rmsf',chain);shade_region(ax,ag,low,high)
    ax=fig.add_axes([.66,.30,.29,.18]);panel(ax,'d','Side-chain response')
    ag=profile(ax,contrast,'chi_entropy',chain);shade_region(ax,ag,low,high)
    ax=fig.add_axes([.10,.115,.85,.095]);ax.axis('off');panel(ax,'e','Evidence boundary')
    if case=='giv':
        line1='Switch II rigidification is detected, but it is also close to the GIV contact surface.'
        line2='This is not a matched W211A HDX test and does not predict GDP release.'
        foot='6MHF has 18 observed GIV residues; GDP is not encoded. Close G-protein homologs occur in training.'
    else:
        line1='L49V is a sequence-only intervention; backbone coordinates are held fixed.'
        line2='A marginal Hχ change cannot establish the correlated-motion mechanism in the paper.'
        foot='UCH 3–238 extracted from 8H1T; C91S reverted by sequence; missing loops and chemical context remain unmodeled.'
    ax.text(0,.70,line1,fontsize=8);ax.text(0,.26,line2,fontsize=8,color=RED)
    footer(fig,'Lines: mean ± SD across 5 checkpoints. Hχ excludes A/G/P. No biological replication or matched experimental error.\n'+foot)
    save(fig,f'fig{num:02d}_{case}_exploratory',line1+' '+line2+' '+foot)

def overview():
    fig=plt.figure(figsize=(183/25.4,115/25.4))
    fig.text(.05,.95,'Mechanistic screens | agreement is not independent validation',fontsize=12,fontweight='bold')
    ax=fig.add_axes([.05,.14,.90,.72]);ax.axis('off')
    rows=[['PDZ3','Deletion loosens shared core','Partial; bound rescue fails','Low-identity flags'],
          ['A6 TCR','Opposite distal responses','Not reproduced by context','Exact sequence'],
          ['Thrombin','TM4-specific 60s response','Not reproduced','Near-identical'],
          ['KIX','MLL-selective loop ordering','Qualitative agreement','Near-identical'],
          ['Gαi / GIV','Switch II stabilization','Detected; contact-proximal','Close homologs'],
          ['BAP1','L49V disrupts coupling','Not testable with Hχ alone','Low-identity flag']]
    table=ax.table(cellText=rows,colLabels=['Protein','Literature criterion','Model screen','Training overlap'],cellLoc='left',colLoc='left',colWidths=[.15,.31,.34,.20],bbox=[0,0,1,1])
    table.auto_set_font_size(False);table.set_fontsize(8)
    for (r,c),cell in table.get_celld().items():
        cell.set_edgecolor('white');cell.set_linewidth(2)
        cell.set_facecolor('#EDF2F5' if r==0 else ('#F5F7F8' if r%2 else 'white'))
        if r==0:cell.set_text_props(weight='bold')
        if c==2 and r==4:cell.set_text_props(color=BLUE,weight='bold')
        if c==2 and r in [2,3]:cell.set_text_props(color=RED)
    fig.text(.05,.055,'102 input states; 510 checkpoint evaluations; 20 KIX conformers. No new MD, fitting or synthetic experimental data.',fontsize=7)
    save(fig,'fig00_evidence_overview','Summary of prespecified literature screens and training-overlap limitations. Descriptive screen outcomes are not binary biological validation and are not selected on best checkpoint.')

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    pdz();tcr();thrombin();kix()
    reserve('giv',5,'Gαi / GIV | a detected response near Switch II','Observed6MHF protein context • exploratory, not a matched HDX replication','giv_binding','rmsf','A',201,216,'Switch II around W211')
    reserve('bap1',6,'BAP1 | mutation sensitivity is not dynamic coupling','L49V sequence-only screen • shared UCH backbone • no modeled mutant relaxation','bap_mutation_free','chi_entropy','L',45,53,'L49 and aromatic neighbors')
    overview()
    (OUT/'figure_manifest.json').write_text(json.dumps(MANIFEST,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Exported',len(MANIFEST),'figures as PNG/PDF/SVG')

if __name__=='__main__':
    main()
