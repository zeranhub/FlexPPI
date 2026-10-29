# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Ordered residue heatmap -> heatmap/plot_comparison.py -> param inherit
# (b) Protein plate -> authentic PyMOL renders from render_pilot.py -> native image
# (c) State/trajectory lines -> LineTrend/plot_trend.py -> param inherit
# (d) Model agreement -> heatmap/plot_comparison.py -> param inherit
# Semantic mismatch prevents native template execution: no method rankings,
# cumulative publication counts, or column-specific min/max normalization here.
# Retained: aligned heatmap cells, no frame, dark lines, light band hierarchy.
# Scientific overrides: common units/scales, all residues/folds, no fabricated CI.
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
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# Preserve the requested final canvas width; no implicit 0.1 inch export pad.
mpl.rcParams['savefig.pad_inches']=0
mpl.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'artifacts/kix_md_20260914/pilot_results_20260916'
OUT=ROOT/'output/flexppi_kix_md_pilot_20260916'
STATES=['apo','myb','mll','both']
LABELS=['Apo','cMyb','MLL','Both']
CONTRASTS=['myb_minus_apo','mll_minus_apo','both_minus_apo','mll_minus_myb']
CLABELS=['cMyb - apo','MLL - apo','Both - apo','MLL - cMyb']
MD='#222222'; MODEL='#2166AC'; L12='#F1A340'
STATE_COLORS=['#666666','#762A83','#1B7837','#B2182B']
DIV=LinearSegmentedColormap.from_list('signed',DIVERGING)
SEQ=LinearSegmentedColormap.from_list('rank',SEQUENTIAL)
SEQ.set_bad('#D4D4D4');DIV.set_bad('#D4D4D4')
R=pd.read_csv(BASE/'relative_residues.csv',dtype={'model_fold':str})
D=pd.read_csv(BASE/'relative_deltas.csv',dtype={'model_fold':str})
REG=pd.read_csv(BASE/'relative_regions.csv',dtype={'model_fold':str})
MET=pd.read_csv(BASE/'relative_metrics.csv',dtype={'model_fold':str})
DM=pd.read_csv(BASE/'relative_delta_metrics.csv',dtype={'model_fold':str})
W=pd.read_csv(BASE/'quality/residue_window_metrics.csv')
TRACE=pd.read_csv(BASE/'quality/trajectory_trace_binned.csv')
THERMO=pd.read_csv(BASE/'quality/thermodynamics_binned.csv')
META=json.loads((BASE/'quality/summary.json').read_text())
FIGS=[]

def subset(frame, **kwargs):
    mask=np.ones(len(frame),bool)
    for k,v in kwargs.items():mask &= frame[k].eq(v)
    return frame.loc[mask].sort_values('author_residue_number') if 'author_residue_number' in frame else frame.loc[mask]

def panel(fig, x, y, letter, title):
    fig.text(x,y,letter,fontsize=11,fontweight='bold',va='bottom')
    fig.text(x+.029,y,title,fontsize=8.5,fontweight='bold',va='bottom')

def footer(fig, text):
    fig.text(.045,.017,text,fontsize=6.6,color='#555555',va='bottom')

def save(fig,name):
    # An invisible fixed canvas prevents tight export from changing the 183 mm width.
    fig.add_artist(Rectangle((0,0),1,1,transform=fig.transFigure,fill=False,edgecolor='none'))
    save_cns_figure(fig,str(OUT/name))
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight',dpi=300,pad_inches=0)
    FIGS.append(name)
    plt.close(fig)

def image_ax(fig,rect,name):
    ax=fig.add_axes(rect)
    img=plt.imread(OUT/'structure_assets'/name)
    yy,xx=np.where(np.any(img[:,:,:3]<.975,axis=2));p=20
    img=img[max(0,yy.min()-p):min(img.shape[0],yy.max()+p),max(0,xx.min()-p):min(img.shape[1],xx.max()+p)]
    ax.imshow(img);ax.axis('off')
    return ax

def state_plot(ax, metric, normalization, ylabel):
    data=subset(REG,metric=metric,normalization=normalization)
    for f in map(str,range(5)):
        ys=[subset(data,state=s,model_fold=f).model_mean.iloc[0] for s in STATES]
        ax.plot(range(4),ys,color=MODEL,lw=.55,alpha=.22,marker='o',ms=2)
    e=subset(data,model_fold='ensemble_mean')
    ax.plot(range(4),[subset(e,state=s).md_mean.iloc[0] for s in STATES],color=MD,lw=1.6,marker='o',ms=4,label='MD: one trajectory')
    ax.plot(range(4),[subset(e,state=s).model_mean.iloc[0] for s in STATES],color=MODEL,lw=1.6,marker='s',ms=4,label='FlexPPI: five-model mean')
    ax.set_xticks(range(4),LABELS);ax.set_ylabel(ylabel);ax.set_xlim(-.2,3.2)

def rank_heatmaps(fig, bottom, height, left=.16, width=.75):
    for j,metric in enumerate(['rmsf','hchi']):
        ax=fig.add_axes([left+j*(width/2+.035),bottom,width/2-.035,height])
        arrays=[];labels=[]
        for s,l in zip(STATES,LABELS):
            group=subset(R,state=s,metric=metric,model_fold='ensemble_mean',normalization='percentile').set_index('author_residue_number')
            for col,label in [('md_value','MD'),('model_value','FlexPPI')]:
                arrays.append(group[col].reindex(range(586,673)).to_numpy());labels.append(l+' / '+label)
        ax.imshow(np.array(arrays),aspect='auto',cmap=SEQ,vmin=0,vmax=100,extent=[585.5,672.5,7.5,-.5],interpolation='nearest')
        ax.add_patch(Rectangle((613.5,-.5),8,8,fill=False,edgecolor=L12,lw=1))
        ax.set_xticks([586,610,630,650,672]);ax.set_yticks(range(8),labels if j==0 else [])
        ax.tick_params(length=0,labelsize=6.5);ax.set_frame_on(False)
        ax.set_title('Main chain: RMSF' if j==0 else r'Side chain: $H_\chi$',fontsize=8,pad=5)
        ax.set_xlabel('KIX residue number',fontsize=7)
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0,100),cmap=SEQ),cax=fig.add_axes([.925,bottom,.011,height]))
    cb.set_ticks([0,50,100]);cb.ax.tick_params(labelsize=6,length=2);cb.set_label('Within-state percentile',fontsize=7)

def figure_absolute():
    fig=plt.figure(figsize=(183/25.4,175/25.4))
    fig.text(.045,.967,'KIX | State-wise flexibility maps',fontsize=11,fontweight='bold')
    panel(fig,.045,.928,'a','The flexible and rigid regions are recovered across four states')
    rank_heatmaps(fig,.715,.172)
    panel(fig,.045,.64,'b','Spatial localization in the ternary complex')
    for x,source,label in [(.04,'md','MD'),(.285,'model','FlexPPI')]:
        fig.text(x+.11,.611,label,ha='center',fontsize=8,fontweight='bold')
        image_ax(fig,[x,.402,.235,.198],f'both_rank_{source}_whole.png')
        image_ax(fig,[x,.192,.235,.19],f'both_rank_{source}_local.png')
    fig.text(.045,.398,'Whole KIX + partners',fontsize=6.8,color='#555555')
    fig.text(.045,.174,'L12 close-up (614–621); identical reference coordinates',fontsize=6.7,color='#555555')
    fig.text(.045,.147,'cMyb: muted purple    MLL: tan    Sticks: L12',fontsize=6.7,color='#555555')
    panel(fig,.565,.64,'c','L12 remains a relatively flexible region')
    ax=fig.add_axes([.638,.438,.305,.154]);state_plot(ax,'rmsf','percentile','Mean RMSF percentile')
    ax.set_ylim(72,98)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,1.32),fontsize=6.4,handlelength=1.4)
    panel(fig,.565,.348,'d','State-wise ranking agreement')
    ax=fig.add_axes([.638,.145,.305,.163])
    for metric,color,marker in [('rmsf',MODEL,'o'),('hchi','#762A83','s')]:
        for f in map(str,range(5)):
            vals=[subset(MET,state=s,metric=metric,model_fold=f,normalization='raw').spearman.iloc[0] for s in STATES]
            ax.scatter(np.arange(4),vals,s=7,c=color,alpha=.22,marker=marker)
        vals=[subset(MET,state=s,metric=metric,model_fold='ensemble_mean',normalization='raw').spearman.iloc[0] for s in STATES]
        ax.plot(range(4),vals,color=color,marker=marker,lw=1.2,ms=4,label='RMSF' if metric=='rmsf' else r'$H_\chi$')
    ax.set_xticks(range(4),LABELS);ax.set_ylim(.65,1.01);ax.set_ylabel('Spearman correlation')
    ax.legend(loc='upper center',bbox_to_anchor=(.5,1.25),ncol=2,fontsize=7)
    footer(fig,'One 100 ns MD trajectory per state; faint marks = model checkpoints, not MD replicates.\nRMSF: 87 residues; Hχ: 74 (A/G/P excluded; grey). Near-identical KIX sequence occurs in training.\nState-wise ranking is not evidence of correctly predicting binding-induced changes.')
    save(fig,'fig01_statewise_flexibility')

def figure_response():
    fig=plt.figure(figsize=(183/25.4,179/25.4))
    fig.text(.045,.967,'KIX | Binding responses: amplitude versus relative ranking',fontsize=10.5,fontweight='bold')
    panel(fig,.045,.92,'a','L12 amplitude does not follow the MD trend')
    panel(fig,.545,.92,'b','L12 ranking is largely preserved')
    ax=fig.add_axes([.115,.708,.337,.17]);state_plot(ax,'rmsf','raw',r'Mean RMSF ($\AA$)')
    ax.legend(loc='lower left',bbox_to_anchor=(-.02,1.01),fontsize=6.1,handlelength=1.3)
    ax=fig.add_axes([.628,.708,.325,.17]);state_plot(ax,'rmsf','percentile','Mean RMSF percentile');ax.set_ylim(72,98)
    panel(fig,.045,.627,'c','Scale-free comparison recovers some, but not all, response patterns')
    for j,metric in enumerate(['rmsf','hchi']):
        ax=fig.add_axes([.245+j*.37,.397,.305,.164])
        array=[]
        for con in CONTRASTS:
            array.append([subset(DM,contrast=con,metric=metric,model_fold='ensemble_mean',region='KIX',normalization=mode).spearman.iloc[0] for mode in ['raw','percentile']])
        array=np.asarray(array)
        ax.imshow(array,cmap=DIV,vmin=-1,vmax=1,aspect='auto')
        for (i,k),v in np.ndenumerate(array):ax.text(k,i,f'{v:.2f}',ha='center',va='center',fontsize=8,color='white' if abs(v)>.65 else BLACK)
        ax.set_xticks([0,1],['Raw difference','Rank shift']);ax.set_yticks(range(4),CLABELS if j==0 else [])
        ax.tick_params(length=0,labelsize=7);ax.set_frame_on(False)
        ax.set_title('RMSF' if metric=='rmsf' else r'$H_\chi$',pad=5)
    fig.text(.05,.361,'Values: Spearman correlation across all matched residues; rank shift = change in within-state percentile.',fontsize=6.7)
    panel(fig,.045,.297,'d','Ternary - apo: opposite L12 amplitude responses')
    for x,source,label in [(.075,'md','MD'),(.46,'model','FlexPPI')]:
        image_ax(fig,[x,.058,.30,.205],f'both_delta_{source}_local.png')
        fig.text(x+.15,.263,label,ha='center',fontsize=8,fontweight='bold')
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-4.2,4.2),cmap=DIV),cax=fig.add_axes([.855,.092,.016,.14]))
    cb.set_label(r'$\Delta$RMSF ($\AA$)',fontsize=7);cb.set_ticks([-4,0,4])
    fig.text(.81,.059,'Blue: lower\nRed: higher',fontsize=6.5)
    footer(fig,'Same input geometry; colors map MD-derived/predicted changes, not coordinate displacements.\nRank normalization is a post-readout sensitivity analysis; it removes global shifts. MD n = 1 per state.')
    save(fig,'fig02_binding_response')

def figure_profiles():
    fig,axes=plt.subplots(4,2,figsize=(183/25.4,205/25.4))
    fig.subplots_adjust(left=.1,right=.97,top=.89,bottom=.09,hspace=.49,wspace=.28)
    for i,(state,label) in enumerate(zip(STATES,LABELS)):
        for j,metric in enumerate(['rmsf','hchi']):
            ax=axes[i,j]
            data=subset(R,state=state,metric=metric,normalization='raw')
            ax.axvspan(613.5,621.5,color=L12,alpha=.15,zorder=0)
            for f in map(str,range(5)):
                g=subset(data,model_fold=f).set_index('author_residue_number').reindex(range(586,673))
                ax.plot(range(586,673),g.model_value,color=MODEL,lw=.5,alpha=.22)
            e=subset(data,model_fold='ensemble_mean').set_index('author_residue_number').reindex(range(586,673))
            ax.plot(range(586,673),e.md_value,color=MD,lw=1,label='MD')
            ax.plot(range(586,673),e.model_value,color=MODEL,lw=1,label='FlexPPI mean')
            ax.set_xlim(586,672);ax.set_ylim(0,11 if metric=='rmsf' else 4.5)
            ax.set_title(label+' | '+('RMSF' if metric=='rmsf' else r'$H_\chi$'),loc='left',fontweight='bold')
            ax.set_ylabel(r'RMSF ($\AA$)' if metric=='rmsf' else r'$H_\chi$ (nats)')
            ax.set_xticks([586,610,630,650,672])
            if i==3:ax.set_xlabel('KIX residue number')
            rho=subset(MET,state=state,metric=metric,model_fold='ensemble_mean',normalization='raw').spearman.iloc[0]
            ax.text(.98,.94,f'ρ = {rho:.3f}',transform=ax.transAxes,ha='right',va='top',fontsize=7)
    fig.text(.06,.965,'Full residue profiles | Four states × two observables',fontsize=11,fontweight='bold')
    fig.legend(handles=[Line2D([],[],color=MD,label='MD: one 100 ns trajectory'),Line2D([],[],color=MODEL,label='FlexPPI: mean of five checkpoints')],loc='upper center',bbox_to_anchor=(.53,.94),ncol=2,fontsize=7)
    footer(fig,'All 87 KIX positions are retained; Hχ curves break at excluded A/G/P (74 evaluated residues). Orange: L12.\nThin blue curves: individual checkpoints. Common axes across states; no smoothing or fitted rescaling.')
    save(fig,'figS01_full_residue_profiles')

def figure_deltas():
    fig,axes=plt.subplots(4,2,figsize=(183/25.4,210/25.4))
    fig.subplots_adjust(left=.12,right=.97,top=.88,bottom=.09,hspace=.49,wspace=.3)
    for i,(contrast,label) in enumerate(zip(CONTRASTS,CLABELS)):
        for j,metric in enumerate(['rmsf','hchi']):
            ax=axes[i,j]
            data=subset(D,contrast=contrast,metric=metric,normalization='raw')
            ax.axvspan(613.5,621.5,color=L12,alpha=.15)
            ax.axhline(0,color=GREY,lw=.6,ls='--')
            for f in map(str,range(5)):
                g=subset(data,model_fold=f).set_index('author_residue_number').reindex(range(586,673))
                ax.plot(range(586,673),g.model_delta,color=MODEL,lw=.5,alpha=.2)
            e=subset(data,model_fold='ensemble_mean').set_index('author_residue_number').reindex(range(586,673))
            ax.plot(range(586,673),e.md_delta,color=MD,lw=1)
            ax.plot(range(586,673),e.model_delta,color=MODEL,lw=1)
            ax.set_xlim(586,672);ax.set_xticks([586,610,630,650,672]);ax.set_ylim(-4.2,4.2) if metric=='rmsf' else ax.set_ylim(-2.2,2.2)
            ax.set_ylabel(r'$\Delta$RMSF ($\AA$)' if metric=='rmsf' else r'$\Delta H_\chi$ (nats)')
            ax.set_title(label,loc='left',fontweight='bold')
            if i==3:ax.set_xlabel('KIX residue number')
    fig.text(.06,.965,'Full response profiles | No selection of favorable residues',fontsize=10.7,fontweight='bold')
    fig.legend(handles=[Line2D([],[],color=MD,label='MD difference'),Line2D([],[],color=MODEL,label='FlexPPI mean difference')],loc='upper center',bbox_to_anchor=(.52,.94),ncol=2,fontsize=7)
    footer(fig,'Raw differences use the frozen before/after state predictions. MD n = 1 per state; thin lines = checkpoints.\nPositive: higher flexibility in the first-named state. Negative: lower. Orange: literature-defined L12.')
    save(fig,'figS02_full_binding_differences')

def figure_quality():
    fig,axes=plt.subplots(3,2,figsize=(183/25.4,189/25.4))
    fig.subplots_adjust(left=.115,right=.97,top=.87,bottom=.13,hspace=.57,wspace=.33)
    for i,(state,label,color) in enumerate(zip(STATES,LABELS,STATE_COLORS)):
        for ax,metric in [(axes[0,0],'KIX_CA_RMSD_angstrom'),(axes[0,1],'L12_CA_RMSD_after_KIX_fit_angstrom')]:
            g=subset(TRACE,state=state,metric=metric)
            ax.plot(g.end_ns_inclusive,g['mean'],color=color,lw=.7,alpha=.8,label=label)
        g=W[(W.state==state)&(W.chain=='B')&W.author_residue_number.between(614,621)&(W.window_kind=='cumulative')].groupby('end_ns_inclusive').rmsf_angstrom.mean()
        axes[1,0].plot(g.index,g.values,color=color,lw=1.2,marker='o',ms=2.5)
        hh=W[(W.state==state)&(W.chain=='B')&W.author_residue_number.between(614,621)&(W.window_kind=='disjoint_half')].groupby('start_ns_exclusive').rmsf_angstrom.mean()
        axes[1,1].plot([i-.14,i+.14],hh.values,color=color,lw=1)
        axes[1,1].scatter(i-.14,hh.iloc[0],c=color,s=24,marker='o')
        axes[1,1].scatter(i+.14,hh.iloc[1],edgecolors=color,facecolors='white',s=24,marker='s')
        t=subset(THERMO,state=state,metric='temperature_kelvin')
        axes[2,1].plot(t.end_ns_inclusive,t['mean'],color=color,lw=.55,alpha=.7)
    for state,chain,color,ls in [('myb','A','#762A83','-'),('mll','C','#1B7837','-'),('both','A','#B2182B','-'),('both','C','#B2182B','--')]:
        g=subset(TRACE,state=state,metric=f'partner_{chain}_native_residue_contact_fraction')
        axes[2,0].plot(g.end_ns_inclusive,g['mean'],color=color,ls=ls,lw=.65,label=('Both / ' if state=='both' else '')+('cMyb' if chain=='A' else 'MLL'))
    titles=['KIX RMSD from first production frame','L12 RMSD after whole-KIX alignment','Cumulative L12 RMSF','L12 sensitivity to the analysis window','Initial interface-contact retention','Production temperature']
    for letter,ax,title in zip('abcdef',axes.flat,titles):
        ax.set_title(title,loc='left',fontweight='bold',fontsize=7.5)
        ax.text(-.16,1.07,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
        ax.set_xlabel('Production time (ns)');ax.set_xlim(0,100)
    for ax in axes[0,:]:ax.set_ylabel(r'RMSD ($\AA$)')
    axes[1,0].set_ylabel(r'Mean RMSF ($\AA$)')
    axes[1,1].set_ylabel(r'Mean RMSF ($\AA$)');axes[1,1].set_xlim(-.5,3.5);axes[1,1].set_xticks(range(4),LABELS);axes[1,1].set_xlabel('')
    axes[1,1].legend(handles=[Line2D([],[],marker='o',color=GREY,ls='',label='0–50 ns'),Line2D([],[],marker='s',markerfacecolor='white',color=GREY,ls='',label='50–100 ns')],loc='upper center',bbox_to_anchor=(.5,1.25),ncol=2,fontsize=6)
    axes[2,0].set_ylabel('Retained fraction');axes[2,0].set_ylim(0,1)
    axes[2,0].legend(loc='upper center',bbox_to_anchor=(.5,1.31),ncol=2,fontsize=6)
    axes[2,1].set_ylabel('Temperature (K)');axes[2,1].axhline(298,color=GREY,lw=.5,ls='--')
    fig.text(.06,.965,'Trajectory quality | Stable thermodynamics ≠ converged regional dynamics',fontsize=10,fontweight='bold')
    fig.legend(handles=[Line2D([],[],color=c,label=l) for c,l in zip(STATE_COLORS,LABELS)],loc='upper center',bbox_to_anchor=(.52,.935),ncol=4,fontsize=7)
    footer(fig,'Every production frame included: traces show 0.1 ns bin means (100 frames/bin), not independent replicates.\nNative contacts: fixed initial ≤4.5 Å heavy-atom contacts; new contacts excluded. Retention is not a binding assay.\nCumulative windows overlap. Disjoint halves are a sensitivity check, not independent-repeat uncertainty.')
    save(fig,'figS03_trajectory_quality')

def figure_relative():
    fig,axes=plt.subplots(2,2,figsize=(183/25.4,157/25.4))
    fig.subplots_adjust(left=.12,right=.97,top=.88,bottom=.13,hspace=.44,wspace=.35)
    for j,metric in enumerate(['rmsf','hchi']):
        ax=axes[0,j]
        for mode,color,marker,offset in [('raw',MD,'o',-.13),('percentile',MODEL,'s',.13)]:
            g=subset(DM,metric=metric,normalization=mode,region='KIX')
            for f in map(str,range(5)):
                ys=[subset(g,contrast=c,model_fold=f).spearman.iloc[0] for c in CONTRASTS]
                ax.scatter(np.arange(4)+offset,ys,color=color,s=9,alpha=.22,marker=marker)
            ys=[subset(g,contrast=c,model_fold='ensemble_mean').spearman.iloc[0] for c in CONTRASTS]
            ax.scatter(np.arange(4)+offset,ys,color=color,s=29,marker=marker,label='Raw difference' if mode=='raw' else 'Percentile shift')
        ax.axhline(0,color=GREY,lw=.6,ls='--');ax.set_ylim(-1,1);ax.set_ylabel('Response Spearman correlation')
        ax.set_title('RMSF' if metric=='rmsf' else r'$H_\chi$',loc='left',fontweight='bold')
        ax=axes[1,j]
        g=subset(DM,metric=metric,normalization='percentile',region='KIX',model_fold='ensemble_mean')
        for offset,col,color,marker,label in [(-.13,'zero_response_mae',GREY,'o','No-change baseline'),(.13,'mae',MODEL,'s','FlexPPI rank shift')]:
            ys=[subset(g,contrast=c)[col].iloc[0] for c in CONTRASTS]
            ax.scatter(np.arange(4)+offset,ys,color=color,s=28,marker=marker,label=label)
        ax.set_ylabel('Rank-shift MAE (percentile points)');ax.set_ylim(0,17)
        ax.set_title('Does prediction improve on no change?',loc='left',fontsize=8)
    for ax in axes.flat:
        ax.set_xticks(range(4),CLABELS,rotation=22,ha='right',fontsize=6.8);ax.set_xlim(-.5,3.5)
    for ax in axes[:,1]:ax.legend(loc='upper center',bbox_to_anchor=(.45,1.29),ncol=2,fontsize=6.5)
    for letter,ax in zip('abcd',axes.flat):ax.text(-.17,1.06,letter,transform=ax.transAxes,fontsize=11,fontweight='bold')
    fig.text(.055,.965,'Scale-free response audit | Better ranking is not universal',fontsize=11,fontweight='bold')
    footer(fig,'Rank shift = after percentile − before percentile; no MD-fitted calibration. All states and checkpoints retained.\nNo-change baseline predicts zero rank shift; lower MAE is better. Faint dots are model checkpoints, not MD error.\nPost-readout sensitivity analysis; MD n = 1/state. RMSF n = 87; Hχ n = 74.')
    save(fig,'figS04_scale_free_response')

def validate():
    assert len(R)==4*6*(87+74)*4 and len(D)==4*6*(87+74)*4
    assert R[['md_value','model_value']].notna().all().all()
    for name,info in META['outputs'].items():assert hashlib.sha256((BASE/'quality'/name).read_bytes()).hexdigest()==info['sha256']
    # Warn-only weak signal is a result, never a reason to alter observed data.
    for frame,cols in [(R,['md_value','model_value']),(D,['md_delta','model_delta'])]:assert np.isfinite(frame[cols].to_numpy()).all()
    # Figure profiles use explicit limits: no silent clipping.
    raw=subset(R,normalization='raw')
    assert raw.loc[raw.metric=='rmsf',['md_value','model_value']].max().max()<11
    assert raw.loc[raw.metric=='hchi',['md_value','model_value']].max().max()<4.5
    dd=subset(D,normalization='raw')
    assert dd.loc[dd.metric=='rmsf',['md_delta','model_delta']].abs().max().max()<4.2
    assert dd.loc[dd.metric=='hchi',['md_delta','model_delta']].abs().max().max()<2.2

if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    validate()
    figure_absolute();figure_response();figure_profiles();figure_deltas();figure_quality();figure_relative()
    (OUT/'figure_manifest.json').write_text(json.dumps(dict(figures=FIGS,source_directory=str(BASE),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rmsf_residues=87,hchi_residues=74,md_repeats_per_state=1,model_checkpoints=5),indent=2)+'\n')
    print('\n'.join(FIGS))
