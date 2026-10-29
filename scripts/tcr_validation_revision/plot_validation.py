# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) residue profile -> LineTrend/plot_trend.py -> param inherit (cumulative
#     publication/event input is semantically incompatible with residue data)
# (b) paired residue strip -> cross-type inherit -> param inherit
# (c,d) calibration scatter -> MarginalDensity/plot_VariableCorrelation.py ->
#     param inherit (log-genome regression/inset tests are not validation axes)
# (e) baseline errors -> BarComparison/plot_comparison_GeneRegulatory.py ->
#     param inherit (single-complex descriptive errors, no replicate bar SD)
# structure panels -> existing project PyMOL layers -> native rendered image
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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
import json, re, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize,TwoSlopeNorm
from matplotlib.lines import Line2D
from PIL import Image
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/tcr_prediction_reference_revision_20260915'
ASSETS=OUT/'structure_assets'
REF=CATEGORICAL[0];PRED='#D67A3A';AA='#999999';SILVER='#D8DCDF'
CMAP=LinearSegmentedColormap.from_list('reference_blue',SEQUENTIAL)
DCMAP=LinearSegmentedColormap.from_list('signed_error',DIVERGING)
EXPORTS=[]

def finish(fig,name):
    # Export fixed journal dimensions: avoid tight-bbox changing final width.
    with mpl.rc_context({'savefig.bbox':None}):
        for ext in ['pdf','svg','png']:fig.savefig(OUT/(name+'.'+ext),dpi=300,facecolor='white')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    outside=[]
    for obj in fig.findobj(mpl.text.Text):
        if not obj.get_visible() or not obj.get_text().strip():continue
        bb=obj.get_window_extent(renderer)
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1>fig.bbox.width+1 or bb.y1>fig.bbox.height+1:outside.append(obj.get_text())
    EXPORTS.append(dict(name=name,width_mm=fig.get_figwidth()*25.4,height_mm=fig.get_figheight()*25.4,outside_text=outside))
    plt.close(fig);print('exported',name,flush=True)

def title(fig,text,sub):
    fig.text(.045,.971,text,fontsize=10,fontweight='bold',va='top')
    fig.text(.045,.943,sub,fontsize=7,color='#555555',va='top')

def head(fig,label,text,x,y):
    fig.text(x,y,label,fontsize=10,fontweight='bold',va='bottom')
    fig.text(x+.035,y,text,fontsize=8,va='bottom')

def number(s):return int(re.match(r'-?\d+',str(s))[0])

def combined_x(df,offset=225):return np.array([number(a)+(offset if c in ['F','beta'] else 0) for c,a in zip(df.chain,df.author_id)])

def profile(ax,df,metric,md='md_reference',pred='prediction',offset=225,sd=None):
    q=df[df.metric==metric].copy();q['x']=combined_x(q,offset)
    for chain,g in q.groupby('chain',sort=False):
        g=g.sort_values('x');x=g.x.to_numpy();y=g[md].to_numpy();p=g[pred].to_numpy()
        if sd is not None:ax.fill_between(x,y-g[sd].to_numpy(),y+g[sd].to_numpy(),color=REF,alpha=.13,lw=0)
        ax.plot(x,y,color=REF,lw=.95,label='MD reference')
        ax.plot(x,p,color=PRED,lw=.95,ls='--',label='FlexPPI')
    ax.set_xlim(0,offset+245);ax.set_xticks([1,50,100,150,200,offset+1,offset+50,offset+100,offset+150,offset+200])
    ax.set_xticklabels(['1','50','100','150','200','1','50','100','150','200'])
    ax.axvline(offset-12,color=SILVER,lw=.6)
    ax.text(.19,1.03,'TCR α',transform=ax.transAxes,color='#555555',fontsize=7)
    ax.text(.73,1.03,'TCR β',transform=ax.transAxes,color='#555555',fontsize=7)
    ax.set_ylabel('RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)')

def agreement(ax,q,metric,x='md_reference',y='prediction',labels=False):
    upper=float(max(q[x].max(),q[y].max())*1.07);lower=min(0.,float(min(q[x].min(),q[y].min())*1.07))
    ax.plot([lower,upper],[lower,upper],color=GREY,ls='--',lw=.6,zorder=0)
    ax.scatter(q[x],q[y],s=10,alpha=.5,c=PRED,edgecolor='none')
    ax.set(xlim=(lower,upper),ylim=(lower,upper),xlabel='MD reference',ylabel='FlexPPI prediction')
    ax.set_aspect('equal',adjustable='box')
    err=q[y]-q[x];rho=spearmanr(q[x],q[y]).statistic
    ax.text(.04,.96,f'ρ = {rho:.3f}\nMAE = {abs(err).mean():.3f}\nn = {len(q)} residues',transform=ax.transAxes,va='top',fontsize=6.5)
    if labels:
        row=q.loc[abs(err).idxmax()]
        ax.annotate('α'+str(row.author_id),(row[x],row[y]),xytext=(-25,12),textcoords='offset points',fontsize=6.5,
                    arrowprops=dict(arrowstyle='-',lw=.5,color=GREY),color=BLACK)

def b423_main(df,metrics):
    valid=df[df.evaluation_mask].copy();fig=plt.figure(figsize=(183/25.4,180/25.4))
    title(fig,'B4.2.3 | direct prediction–MD validation','5IVX holo · frozen fold 3 · absolute flexibility, not the apo→holo change')
    fig.legend([Line2D([],[],color=REF,lw=1.2),Line2D([],[],color=PRED,lw=1.2,ls='--')],['MD reference','FlexPPI'],
               loc='upper center',bbox_to_anchor=(.53,.916),ncol=2,fontsize=7)
    head(fig,'a','Backbone flexibility at every TCR residue',.045,.858)
    ax=fig.add_axes([.095,.663,.855,.173]);profile(ax,valid,'rmsf');ax.set_ylim(0,7.4)
    q=valid[(valid.chain=='E')&(valid.author_id=='196')&(valid.metric=='rmsf')].iloc[0]
    ax.annotate('α196: large underprediction',(196,q.md_reference),xytext=(240,5.4),fontsize=6.5,
                arrowprops=dict(arrowstyle='-',color=BLACK,lw=.5))
    head(fig,'b','Side-chain flexibility: paired residue-by-residue values',.045,.589)
    ax=fig.add_axes([.095,.480,.855,.088])
    h=valid[valid.metric=='hchi'];arr=np.full((2,471),np.nan)
    for r in h.itertuples():
        ix=number(r.author_id)+(225 if r.chain=='F' else 0)
        arr[:,ix]=[r.md_reference,r.prediction]
    cmap=CMAP.copy();cmap.set_bad('#E7E9EB')
    ax.imshow(arr,aspect='auto',cmap=cmap,vmin=0,vmax=4,interpolation='none',extent=(-.5,470.5,1.5,-.5))
    ax.set_yticks([0,1],['MD','FlexPPI']);ax.tick_params(axis='y',length=0);ax.set_xticks([1,100,200,226,325,425],['α1','100','200','β1','100','200'])
    ax.spines[['bottom','left']].set_visible(False)
    ca=fig.add_axes([.77,.594,.18,.009]);cb=fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0,4),cmap=CMAP),cax=ca,orientation='horizontal',ticks=[0,2,4]);cb.ax.tick_params(labelsize=6,length=2);cb.outline.set_visible(False)
    fig.text(.72,.591,'nats',fontsize=6,ha='right')
    fig.text(.095,.433,'Author residue number · grey cells: unscored / absent; no interpolation',fontsize=6.5,color='#666666')
    head(fig,'c','RMSF agreement (Å)',.045,.385);head(fig,'d','Hχ agreement (nats)',.37,.385);head(fig,'e','Error versus baseline',.715,.385)
    agreement(fig.add_axes([.095,.126,.23,.235]),valid[valid.metric=='rmsf'],'rmsf',labels=True)
    agreement(fig.add_axes([.42,.126,.23,.235]),valid[valid.metric=='hchi'],'hchi')
    for i,metric in enumerate(['rmsf','hchi']):
        ax=fig.add_axes([.795,.283-i*.153,.155,.065])
        g=metrics[(metrics.scope=='TCR alpha beta')&(metrics.metric==metric)].set_index('predictor')
        vals=[g.loc['FlexPPI','mae'],g.loc['AA-only','mae']]
        ax.barh([1,0],vals,color=[PRED,AA],height=.48)
        ax.set_yticks([1,0],['FlexPPI','AA-only'],fontsize=6.5);ax.set_xlim(0,max(vals)*1.37)
        ax.set_xticks([]);ax.spines[['bottom','left']].set_visible(False);ax.tick_params(length=0)
        for y,v in zip([1,0],vals):ax.text(v+max(vals)*.055,y,f'{v:.3f}',va='center',fontsize=6.5)
        ax.text(0,1.23,'RMSF MAE (Å)' if metric=='rmsf' else 'Hχ MAE (nats)',transform=ax.transAxes,fontsize=7)
        fig.text(.745,.257-i*.153,f"{g.loc['FlexPPI','mae_reduction_vs_AA_percent']:.1f}% lower MAE",fontsize=7,color=PRED)
    fig.text(.045,.057,'AA-only: means fitted on fold-3 training labels only. One complex; n denotes residues, not biological repeats.',fontsize=6.3,color='#555555')
    fig.text(.045,.035,'Combination-level holdout; homologous training proteins remain. Reference: existing ITrace MD, not experimental dynamics.',fontsize=6.3,color='#555555')
    finish(fig,'fig01_b423_prediction_vs_MD')

def image_layer(fig,bounds,name,meta,labels=False):
    ax=fig.add_axes(bounds);im=Image.open(ASSETS/(name+'.png'));ax.imshow(im);ax.axis('off')
    if labels:
        for r in meta['landmarks']:
            x,y=r['pixel'];label=('N133' if r['author_id']=='133' else 'K134')
            if r['author_id']=='196':continue
            ax.annotate(label,(x,y),xytext=(20,-25 if label=='N133' else 25),textcoords='offset points',fontsize=6.5,
                        arrowprops=dict(arrowstyle='-',lw=.6,color=BLACK),color=BLACK)
    if name.startswith('b423_rmsf_'):
        for r in meta['landmarks']:
            if r['author_id'] not in ['196','213']:continue
            label='α196' if r['chain']=='E' else 'β213'
            ax.annotate(label,r['pixel'],xytext=(-37,5) if label=='α196' else (22,-10),textcoords='offset points',fontsize=6.4,
                        arrowprops=dict(arrowstyle='-',lw=.5,color=BLACK),color=BLACK)
    return ax

def structure_plate(df):
    path=ASSETS/'render_provenance.json'
    if not path.exists():return
    audit=json.loads(path.read_text());meta={r['name']:r for r in audit['assets']}
    fig=plt.figure(figsize=(183/25.4,208/25.4))
    title(fig,'B4.2.3 | where predictions agree, and where they fail','Identical 5IVX coordinates and camera · scalars from MD reference, model prediction and their difference')
    cols=[.04,.355,.67]
    for x,lab in zip(cols,['MD reference','FlexPPI','Prediction − MD']):fig.text(x+.14,.897,lab,ha='center',fontsize=9,fontweight='bold')
    for row,metric in enumerate(['rmsf','hchi']):
        y=.644-row*.294
        for col,(x,kind) in enumerate(zip(cols,['md','pred','error'])):
            name=f'b423_{metric}_{kind}';image_layer(fig,[x,y,.285,.231],name,meta[name])
            fig.text(x+.002,y+.22,chr(97+row*3+col),fontsize=10,fontweight='bold')
        fig.text(.014,y+.13,'RMSF' if metric=='rmsf' else 'Hχ',rotation=90,va='center',fontsize=9)
        r=meta[f'b423_{metric}_md']['range'];er=meta[f'b423_{metric}_error']['range']
        for bounds,ran,cm in [([.115,y-.015,.43,.008],r,LinearSegmentedColormap.from_list('struct',['#DBEBF5','#6BAED6','#08306B'])),([.70,y-.015,.235,.008],er,DCMAP)]:
            cb=fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(*ran),cmap=cm),cax=fig.add_axes(bounds),orientation='horizontal');cb.outline.set_visible(False);cb.ax.tick_params(labelsize=6,length=2)
        fig.text(.56,y-.015,'Å' if metric=='rmsf' else 'nats',fontsize=6)
    head(fig,'g','Local side-chain comparison: β127–138',.045,.270)
    for x,kind in zip(cols,['md','pred','error']):
        name='b423_local_hchi_'+kind;image_layer(fig,[x,.097,.285,.155],name,meta[name],labels=True)
    q=df[(df.metric=='hchi')&(df.chain=='F')&df.author_id.isin(['133','134'])].set_index('author_id')
    for x,kind,col in zip(cols,['md','pred','error'],['md_reference','prediction','model_error']):
        fig.text(x+.14,.086,f"N133: {q.loc['133',col]:+.3f}   K134: {q.loc['134',col]:+.3f} nats" if kind=='error' else f"N133: {q.loc['133',col]:.3f}   K134: {q.loc['134',col]:.3f} nats",ha='center',fontsize=6.3)
    fig.text(.045,.047,'Blue error: underprediction; red error: overprediction. Each metric uses one shared full-range scale; no value clipping.',fontsize=6.3,color='#555555')
    fig.text(.045,.027,'Grey: partner / unscored residue. Local Hχ panels inherit row-2 scales. These are scalar maps, not three different conformations.',fontsize=6.3,color='#555555')
    finish(fig,'fig02_b423_structure_reference_prediction_error')

def b423_supplement(df):
    v=df[df.evaluation_mask];fig=plt.figure(figsize=(183/25.4,160/25.4))
    title(fig,'B4.2.3 | validation beyond amino-acid type','Training-only priors · local successes and failures retained · no test-fitted calibration')
    fig.legend([Line2D([],[],color=REF,lw=1.2),Line2D([],[],color=PRED,lw=1.2,ls='--')],['MD reference','FlexPPI'],loc='upper center',bbox_to_anchor=(.54,.916),ncol=2,fontsize=7)
    head(fig,'a','Hχ after subtracting the training amino-acid prior',.045,.86)
    q=v[v.metric=='hchi'];ax=fig.add_axes([.105,.49,.33,.32])
    agreement(ax,q,'hchi',x='aa_adjusted_md_reference',y='aa_adjusted_prediction')
    ax.set_xlabel('MD Hχ − training AA mean (nats)');ax.set_ylabel('Predicted Hχ − training AA mean (nats)')
    head(fig,'b','Within-type prediction errors',.50,.86)
    ax=fig.add_axes([.56,.49,.39,.31]);types=sorted(q.amino_acid.unique())
    for j,aa in enumerate(types):
        g=q[q.amino_acid==aa];xx=np.full(len(g),j)
        ax.scatter(xx,g.model_error,s=5,alpha=.45,color=PRED,edgecolors='none')
        ax.plot([j-.22,j+.22],[g.model_error.mean()]*2,color=BLACK,lw=.9)
    ax.axhline(0,color=GREY,lw=.6,ls='--');ax.set_xticks(range(len(types)),types);ax.set_ylabel('Prediction − MD (nats)');ax.set_xlabel('Amino-acid type')
    for k,(chain,start,end,title_,metric) in enumerate([('F',127,138,'Predefined β127–138 window','hchi'),('E',185,196,'α-chain terminus: largest RMSF error','rmsf')]):
        head(fig,chr(99+k),title_,.045+.455*k,.363)
        ax=fig.add_axes([.10+.455*k,.123,.36,.196]);g=v[(v.chain==chain)&(v.metric==metric)].copy();g['number']=g.author_id.map(number);g=g[g.number.between(start,end)].sort_values('number')
        ax.plot(g.number,g.md_reference,color=REF,lw=1,marker='o',ms=2.8,label='MD reference');ax.plot(g.number,g.prediction,color=PRED,lw=1,marker='s',ms=2.8,ls='--',label='FlexPPI')
        ax.set_ylabel('Hχ (nats)' if metric=='hchi' else 'RMSF (Å)');ax.set_xlabel('Author residue number')
    fig.text(.045,.036,'Dots: all scored residues within each amino-acid type; black ticks: mean signed error. Local windows supplement full-residue main panels.',fontsize=6.3,color='#555555')
    finish(fig,'figS01_b423_AA_controls_and_local_errors')

def ig4_figures():
    if not (OUT/'1g4_validation_audit.json').exists():return
    states=pd.read_csv(OUT/'1g4_state_summary.csv',dtype={'author_id':str});delta=pd.read_csv(OUT/'1g4_delta_paired_residues.csv',dtype={'author_id':str})
    priors=pd.read_csv(ROOT/'artifacts/evaluation_atlas_hchi_v1/analysis/training_fold_aa_priors.csv')
    priors=priors[priors.fold==0].copy();priors.metric=priors.metric.replace({'chi_entropy':'hchi'})
    states=states.merge(priors[['amino_acid','metric','training_mean']],on=['amino_acid','metric'],validate='many_to_one')
    for metric in ['rmsf','hchi']:
        unit='Å' if metric=='rmsf' else 'nats';pretty='RMSF' if metric=='rmsf' else 'Hχ';q=delta[delta.metric==metric].copy()
        fig=plt.figure(figsize=(183/25.4,193/25.4))
        title(fig,f'1G4 | {pretty} predictions tested against matched MD',
              'Free TCR versus TCR–pMHC (no CD3) · fixed 10 ns snapshot inputs · 20–100 ns MD reference')
        fig.legend([Line2D([],[],color=REF,lw=1.2),Line2D([],[],color=PRED,lw=1.2,ls='--')],['MD reference','FlexPPI'],loc='upper center',bbox_to_anchor=(.54,.915),ncol=2,fontsize=7)
        head(fig,'a','Binding-induced change: prediction versus reference',.045,.858)
        ax=fig.add_axes([.10,.64,.85,.185]);q['x']=combined_x(q,225)
        for chain,g in q.groupby('chain',sort=False):
            g=g.sort_values('x');x=g.x.to_numpy();y=g.md_delta.to_numpy();e=g.md_delta_se.to_numpy()
            ax.fill_between(x,y-e,y+e,color=REF,alpha=.16,lw=0)
            ax.plot(x,y,color=REF,lw=.9);ax.plot(x,g.prediction_delta,color=PRED,lw=.9,ls='--')
        ax.axhline(0,color=GREY,ls=':',lw=.6);ax.axvline(214,color=SILVER,lw=.6)
        ax.set_xlim(0,470);ax.set_xticks([1,100,200,226,325,425],['α1','100','200','β1','100','200']);ax.set_ylabel(f'Δ{pretty} ({unit})');ax.set_xlabel('MD author residue number')
        fig.text(.10,.572,'Δ = holo − apo; shading: MD mean difference ± SE across independent 10-run state ensembles.',fontsize=6.4,color='#555555')
        head(fig,'b','Are changes predicted correctly?',.045,.536);head(fig,'c','Absolute apo prediction',.375,.536);head(fig,'d','Absolute holo prediction',.695,.536)
        for i in range(3):
            ax=fig.add_axes([.10+.31*i,.262,.225,.24])
            if i==0:agreement(ax,q,metric,x='md_delta',y='prediction_delta');ax.axhline(0,color=SILVER,lw=.5);ax.axvline(0,color=SILVER,lw=.5)
            else:
                g=states[(states.metric==metric)&(states.state==('apo' if i==1 else 'holo'))]
                agreement(ax,g,metric,x='md_mean',y='prediction_mean')
                ax.text(.04,.69,f'AA-only MAE = {abs(g.training_mean-g.md_mean).mean():.3f}',transform=ax.transAxes,fontsize=6.2,color='#666666')
            ax.set_xlabel(f'MD {"Δ" if i==0 else ""}{pretty} ({unit})');ax.set_ylabel(f'Predicted {"Δ" if i==0 else ""}{pretty} ({unit})')
        head(fig,'e','Does the model beat a no-change prediction?',.045,.178)
        ax=fig.add_axes([.37,.095,.43,.070]);errors=[float(abs(q.prediction_delta-q.md_delta).mean()),float(abs(q.md_delta).mean())]
        ax.barh([1,0],errors,color=[PRED,AA],height=.5);ax.set_yticks([1,0],['FlexPPI','Predict no change']);ax.set_xlim(0,max(errors)*1.25);ax.tick_params(axis='y',length=0)
        ax.spines[['bottom','left']].set_visible(False);ax.set_xticks([])
        for y,value in zip([1,0],errors):ax.text(value+max(errors)*.025,y,f'{value:.3f} {unit}',fontsize=7,va='center')
        fig.text(.37,.080,'Mean absolute error of the binding-induced change (lower is better)',fontsize=6.5)
        fig.text(.045,.042,'Both predictions and references are 10-run state means. MD-derived inputs; not an MD-free PDB deployment.',fontsize=6.3,color='#555555')
        fig.text(.045,.022,'Fold 0 excludes 2BNQ/2BNR; homologs remain. AA-only: training priors. Published MD: Knapp et al., PLoS Comput Biol (2019).',fontsize=6.3,color='#555555')
        finish(fig,f'fig03_1g4_{metric}_prediction_vs_MD' if metric=='rmsf' else 'figS02_1g4_hchi_prediction_vs_MD')

def main():
    p=argparse.ArgumentParser();p.add_argument('--only',choices=['all','b423','1g4'],default='all');args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.only in ['all','b423']:
        df=pd.read_csv(OUT/'validation_paired_residues.csv',dtype={'author_id':str});metrics=pd.read_csv(OUT/'validation_metrics.csv')
        b423_main(df,metrics);structure_plate(df);b423_supplement(df)
    if args.only in ['all','1g4']:ig4_figures()
    path=OUT/'figure_export_audit.json';old=json.loads(path.read_text()) if path.exists() else []
    d={r['name']:r for r in old};d.update({r['name']:r for r in EXPORTS});path.write_text(json.dumps(list(d.values()),indent=2),encoding='utf-8')

if __name__=='__main__':main()
