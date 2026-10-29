# Asset confirmation: quantitative-grid revision of the existing v13 scripts.
# Figure 2 ECDF / box summaries: original fig2/plot.py, visual adaptation.
# Figure 3 scatter: prodesign/plot_reference_amplitude.py, parameter inheritance.
# Figure 4 paired data: external prepare_and_plot_external_layout_v2.py, visual adaptation.
# Figure S1 composition: BarComposition/plot_rewriting.py, cross-type inheritance.
# Existing companion PNGs inspected; frozen CSVs replace original data mappings.
# User's explicit red-bracket revision instructions authorize the present layout.
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
import hashlib, json, shutil, xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, NullFormatter, NullLocator
from matplotlib import font_manager
from PIL import Image

OUT = Path(__file__).resolve().parent
BASE = Path(r"C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002")
MAIN = BASE / 'figures/corrected_main_results_preparation/completed/20261003T212155Z_0caa2efe'
EXT = BASE / 'figures/corrected_external_results_preparation/completed/20261003T163520Z_03a3b308'
COHORT = BASE / 'figures/corrected_cohort_lineage'
LABEL = BASE / 'figures/corrected_label_variation'
M = ['full', 'intrachain', 'static_sequence_complex', 'aa_only']
ML = dict(zip(M, ['Full', 'Intra', 'Static', 'AA']))
MC = dict(zip(M, [CATEGORICAL[0], CATEGORICAL[3], CATEGORICAL[2], CATEGORICAL[5]]))
LS = dict(zip(M, ['-', '--', '-.', ':']))
TASKS = ['rmsf', 'chi_entropy']
TN = {'rmsf':'RMSF', 'chi_entropy':'Hχ'}
UN = {'rmsf':'Å', 'chi_entropy':'nats'}
SRC = ['atlas', 'dynarepo', 'itrace']
SC = dict(zip(SRC+["thermopcd"], [CATEGORICAL[0], CATEGORICAL[3], CATEGORICAL[4], BLACK]))
SN = dict(zip(SRC+["thermopcd"], ['ATLAS', 'DynaRepo', 'ITrace', 'ThermoPCD']))
SM = dict(zip(SRC, ['o', '^', 's']))
AA = list('AGCSTVNDHILFPWYQEMRK')
EVIDENCE = {}
MANIFEST = {'figures':{}, 'confidence_intervals':[], 'statistical_tests':[],
            'new_training':0, 'new_inference':0, 'random_subsampling':0}
font_manager.findfont('Arial', fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Arial', 'Helvetica', 'Liberation Sans'],
    'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':8,
    'axes.titlesize':8.5, 'axes.labelsize':8, 'axes.unicode_minus':False})

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def read(path, audit=None):
    path=Path(path); actual=sha(path)
    if audit:
        a=json.loads(Path(audit).read_text(encoding='utf-8'))
        assert actual==a['output_sha256'][path.name], f'Upstream hash mismatch: {path}'
        EVIDENCE[str(audit)]=sha(audit)
    EVIDENCE[str(path)]=actual
    return pd.read_csv(path)

def write(d,name):d.to_csv(OUT/name,index=False,float_format='%.17g',na_rep='')

def label(ax,letter,title):
    ax.set_title(title,loc='left',pad=9)
    ax.text(-.15,1.055,letter,transform=ax.transAxes,fontsize=11,weight='bold',ha='left',va='bottom')

def grid(height,rows=3,bottom=.105):
    fig=plt.figure(figsize=(183/25.4,height/25.4))
    gs=fig.add_gridspec(rows,2,left=.105,right=.94,bottom=bottom,top=.875 if rows==1 else .93,
        wspace=.40,hspace=.65)
    return fig,[fig.add_subplot(gs[i//2,i%2]) for i in range(rows*2)]

def numeric_ticks(ax,which='x'):
    axis=ax.xaxis if which=='x' else ax.yaxis
    axis.set_major_formatter(FuncFormatter(lambda v,p:format(v,'g')))
    axis.set_minor_locator(NullLocator())

def export(fig,stem,caption,panels,stats,sources):
    fig.canvas.draw(); renderer=fig.canvas.get_renderer()
    outside_view_ticks=set()
    for ax in fig.axes:
        for axis in [ax.xaxis,ax.yaxis]:
            lo,hi=sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if tick.get_loc()<lo-1e-10 or tick.get_loc()>hi+1e-10:
                    outside_view_ticks.update([tick.label1,tick.label2])
    texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text() and t not in outside_view_ticks]
    outside=[]
    for t in texts:
        bb=t.get_window_extent(renderer)
        if bb.x0 < -.5 or bb.y0 < -.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:outside.append(t.get_text())
    assert not outside, (stem,'text outside canvas',outside)
    assert min(t.get_fontsize() for t in texts)>=7.5
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/(stem+'.'+ext),bbox_inches=fig.bbox_inches,dpi=300)
    svg=ET.parse(OUT/(stem+'.svg')); ns={'s':'http://www.w3.org/2000/svg'}
    assert len(svg.findall('.//s:image',ns))==0, 'Rasterized SVG layer'
    im=Image.open(OUT/(stem+'.png')); assert min(im.info.get('dpi',(0,0)))>=299
    dpi=im.info['dpi']; dims=[float(v*25.4) for v in fig.get_size_inches()]
    pixels=np.asarray(im.convert('RGB')); density=[]
    for ax in fig.axes:
        p=ax.get_position(); x0=int(p.x0*im.width);x1=int(p.x1*im.width)
        y0=int((1-p.y1)*im.height);y1=int((1-p.y0)*im.height)
        crop=pixels[y0:y1,x0:x1];density.append(float(np.mean(np.any(crop<245,axis=2))))
    qa={'width_mm':dims[0],'height_mm':dims[1],'dpi':list(dpi),'minimum_text_pt':min(t.get_fontsize() for t in texts),
        'text_outside_canvas':outside,'svg_raster_images':0,'svg_editable_text':len(svg.findall('.//s:text',ns)),
        'pixel_content_density':density,'visual_review':'pending full-image inspection',
        'font':'Microsoft YaHei / Arial; TrueType PDF fonts, editable SVG text',
        'all_data_preserved':True, 'no_ci_or_p_values':True}
    write(pd.DataFrame(stats),stem+'_statistics.csv')
    (OUT/(stem+'_caption_zh.md')).write_text(caption,encoding='utf-8')
    MANIFEST['figures'][stem]={'caption_zh':caption,'panels':panels,'statistics_csv':stem+'_statistics.csv',
        'inputpaths':sources,'input_sha256':{s:EVIDENCE[s] for s in sources},'qa':qa,
        'outputs':{ext:{'path':str(OUT/(stem+'.'+ext)),'sha256':sha(OUT/(stem+'.'+ext))} for ext in ['png','pdf','svg']}}
    plt.close(fig)

def prepare():
    per=read(MAIN/'per_record_metrics.csv',MAIN/'render_audit.json')
    assert sha(MAIN/'per_record_metrics.csv')=='22cc931321f0962e7f89d736b9036dea5e0802fce89fd906ea191f56ea1e83be'
    per=per[per.variant.isin(M)].copy()
    assert len(per)==3063*8 and not per.duplicated(['sample_id','variant','metric']).any()
    expected={'full':[.577104,.165866],'intrachain':[.635527,.169617],'static_sequence_complex':[.677328,.186988]}
    for method,vals in expected.items():
        for task,v in zip(TASKS,vals):
            d=per[(per.variant==method)&(per.metric==task)]
            assert len(d)==3063 and np.isclose(d.mae.mean(),v,atol=.0000006,rtol=0)
            assert d['count'].sum()==(961321 if task=='rmsf' else 824988)
    chains=read(COHORT/'all_chain_lengths_and_overlap.csv')
    assert len(chains)==5068 and not chains.duplicated(['sample_id','chain']).any()
    counts=chains.groupby('sample_id').size().rename('input_chains')
    per=per.merge(counts,on='sample_id',validate='many_to_one')
    assert per.input_chains.notna().all()
    per['chain_group']=np.where(per.input_chains==1,'1',np.where(per.input_chains==2,'2','3+'))
    labels=read(LABEL/'per_record_labels.csv',LABEL/'label_analysis_audit.json')
    aa=read(MAIN/'amino_acid_support_and_macro_metrics.csv',MAIN/'render_audit.json')
    pairs=read(MAIN/'amino_acid_paired_mae_support_and_macros.csv',MAIN/'render_audit.json')
    external=read(EXT/'all_32_per_record_metrics.csv',EXT/'render_audit.json')
    write(per,'fig2_all_3063_four_method_metrics.csv')
    write(chains,'figS1_all_5068_chain_observations.csv')
    return per,chains,labels,aa,pairs,external

def fig2(per):
    fig,axs=grid(188);stats=[];plot=[]
    handles=[Line2D([],[],color=MC[m],ls=LS[m],lw=1.4,label=ML[m]) for m in M]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.995),ncol=4,columnspacing=2,handlelength=2.5)
    for j,task in enumerate(TASKS):
        ax=axs[j]
        for method in M:
            d=per[(per.variant==method)&(per.metric==task)].sort_values('mae',kind='mergesort')
            fr=np.arange(1,len(d)+1)/len(d)
            ax.step(d.mae,fr,where='post',color=MC[method],ls=LS[method],lw=1.4)
            stats.append({'panel':'ab'[j],'metric':task,'variant':method,'n_records':len(d),
                'valid_residues':int(d['count'].sum()),'mean_mae':d.mae.mean(),'median_mae':d.mae.median(),
                'mean_spearman':d.spearman.mean(),'correlation_defined_records':int(d.spearman.notna().sum())})
            for r,f in zip(d.itertuples(),fr):plot.append({'panel':'ab'[j],'sample_id':r.sample_id,'variant':method,'metric':task,'x':r.mae,'y':f})
        ax.set(ylim=(0,1.02),ylabel='Fraction of records',xlabel=f'MAE ({UN[task]})')
        if task=='rmsf':
            ax.set_xscale('symlog',linthresh=.1);ax.set_xticks([0,.1,1,10]);numeric_ticks(ax)
            ax.set_xlabel('MAE (Å)')
        ax.set_xlim(0,per[(per.metric==task)].mae.max()*1.04)
        label(ax,'ab'[j],TN[task]+': error distribution')
        ax=axs[2+j]
        values=[per[(per.variant==m)&(per.metric==task)].spearman.dropna().to_numpy() for m in M]
        # Descriptive quantiles: boxes Q1-Q3, middle line median, whiskers full range.
        bp=ax.boxplot(values,orientation='horizontal',positions=np.arange(4),widths=.50,patch_artist=True,
            whis=(0,100),showfliers=True,medianprops={'color':BLACK,'lw':1.0},
            whiskerprops={'lw':.8},capprops={'lw':.8})
        for i,m in enumerate(M):
            bp['boxes'][i].set(facecolor=MC[m],alpha=.25,edgecolor=MC[m],linewidth=.9)
            mu=float(np.mean(values[i]));ax.scatter(mu,i,s=21,color=MC[m],edgecolors='white',linewidths=.6,zorder=4)
            ax.text(1.012,i,f'{mu:.3f}',transform=ax.get_yaxis_transform(),ha='left',va='center',fontsize=8)
            q=np.quantile(values[i],[0,.25,.5,.75,1])
            stats.append({'panel':'cd'[j],'metric':task,'variant':m,'n_records':len(values[i]),
                'minimum':q[0],'q1':q[1],'median':q[2],'q3':q[3],'maximum':q[4],'mean':mu})
        ax.set_yticks(range(4),[ML[m] for m in M]);ax.invert_yaxis();ax.set_xlabel('Within-record Spearman ρ')
        ax.set_xlim(min(np.min(v) for v in values)-.03,1.01);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
        label(ax,'cd'[j],TN[task]+': residue ranking')
        ax=axs[4+j]
        d=per[(per.metric==task)&per.variant.isin(['full','intrachain'])]
        wide=d.pivot(index=['sample_id','dataset','input_chains','chain_group'],columns='variant',values='mae').reset_index()
        wide['gain']=wide.intrachain-wide.full
        write(wide,f'fig2_{task}_all_3063_chain_gains.csv')
        bins=[]
        for i,g in enumerate(['1','2','3+']):
            v=wide[wide.chain_group==g].sort_values('sample_id');bins.append(v.gain.to_numpy())
            x=i+((np.arange(len(v))%67)/66-.5)*.40
            ax.scatter(x,v.gain,s=2.2,color=CATEGORICAL[0],alpha=.32,linewidths=0)
            ax.scatter(i,v.gain.mean(),s=30,color=BLACK,marker='D',edgecolor='white',linewidths=.5,zorder=5)
            stats.append({'panel':'ef'[j],'metric':task,'input_chain_group':g,'n_records':len(v),
                'mean_gain':v.gain.mean(),'median_gain':v.gain.median(),'positive_records':int((v.gain>0).sum()),
                'zero_records':int((v.gain==0).sum()),'negative_records':int((v.gain<0).sum()),
                'minimum_gain':v.gain.min(),'maximum_gain':v.gain.max()})
        bp=ax.boxplot(bins,positions=[0,1,2],widths=.48,showfliers=False,patch_artist=True,
            medianprops={'color':BLACK,'lw':1},whiskerprops={'color':GREY,'lw':.6},capprops={'color':GREY,'lw':.6})
        for b in bp['boxes']:b.set(facecolor='none',edgecolor=BLACK,linewidth=.8)
        ax.axhline(0,color=GREY,lw=.8,ls='--');ax.set_yscale('symlog',linthresh=.04 if j==0 else .003)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:format(v,'g')))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_xticks([0,1,2],[g+'\nn = '+str(len(v)) for g,v in zip(['1','2','3+'],bins)])
        ax.set_xlabel('Number of input chains');ax.set_ylabel('Intra − Full MAE ('+UN[task]+')')
        label(ax,'ef'[j],TN[task]+': effect of chain count')
    write(pd.DataFrame(plot),'fig2_ecdf_all_24504_points.csv')
    caption='图2：五折折外预测的整体准确率及输入链数分层。Full为完整结构模型，Intra为分别训练的链内模型，Static为静态特征回归对照，AA为氨基酸类型均值对照。（a,b）两项目标逐记录MAE的经验累积分布；每条曲线包含全部3,063条记录，RMSF横轴在近零区域线性、远端为对数。（c,d）逐记录Spearman相关的分布，箱体为四分位区间，线为中位数，须线覆盖完整范围，圆点及右侧数字为记录等权均值。（e,f）按静态输入中的实际链数分为1、2、3条及以上，展示同一记录的Intra减Full MAE；正值表示Full误差较低，每点为一条记录，黑色菱形为组内记录等权均值。差值轴在零附近线性、两端为对数。1链组的差值反映两种分别训练模型的预测差异，不能解释为加入伙伴链带来的改善。所有分布和分层均为描述性统计，未绘制置信区间或进行显著性检验。'
    export(fig,'fig2_accuracy_and_chain_context',caption,['a RMSF MAE ECDF','b Hχ MAE ECDF','c RMSF Spearman','d Hχ Spearman','e RMSF gain by actual chain count','f Hχ gain by actual chain count'],stats,[str(MAIN/'per_record_metrics.csv'),str(COHORT/'all_chain_lengths_and_overlap.csv')])

def fig3(per,labels,aa,pairs):
    fig,axs=grid(202);stats=[];plotted=[]
    # Bottom left logical panel e is split into two aligned axes to keep units separate.
    rect=axs[4].get_position();fig.delaxes(axs[4]);h=(rect.height-.028)/2
    e1=fig.add_axes([rect.x0,rect.y0+h+.028,rect.width,h]);e2=fig.add_axes([rect.x0,rect.y0,rect.width,h])
    sources=[Line2D([],[],ls='',marker=SM[s],color=SC[s],markersize=4,label=SN[s]) for s in SRC]
    fig.legend(handles=sources,loc='upper center',bbox_to_anchor=(.53,.995),ncol=3,columnspacing=2.3)
    full=per[per.variant=='full'].merge(labels,on=['sample_id','dataset','biological_system'],validate='many_to_one')
    assert len(full)==6126
    for j,task in enumerate(TASKS):
        d=full[full.metric==task].copy();prefix='rmsf' if j==0 else 'hchi'
        d['reference_mean']=d[prefix+'_label_mean'];assert (d.reference_mean>0).all()
        d['relative_mae_percent']=100*d.mae/d.reference_mean
        assert len(d)==3063 and np.isfinite(d.relative_mae_percent).all()
        for rowoff,ycol in [(0,'mae'),(2,'relative_mae_percent')]:
            ax=axs[rowoff+j]
            for s in SRC:
                v=d[d.dataset==s]
                ax.scatter(v.reference_mean,v[ycol],s=4,marker=SM[s],color=SC[s],alpha=.46,linewidths=0)
                rho=float(v.reference_mean.rank().corr(v[ycol].rank()))
                stats.append({'panel':'abcd'[rowoff+j],'metric':task,'dataset':s,'n_records':len(v),
                    'spearman_reference_error':rho,'mean_y':v[ycol].mean(),'minimum_y':v[ycol].min(),'maximum_y':v[ycol].max()})
                for r in v.itertuples():plotted.append({'panel':'abcd'[rowoff+j],'sample_id':r.sample_id,'dataset':s,'metric':task,'x':r.reference_mean,'y':getattr(r,ycol)})
            ax.set_xlabel('Mean reference '+TN[task]+' ('+UN[task]+')')
            ax.set_ylabel('Full MAE ('+UN[task]+')' if rowoff==0 else 'Full MAE / reference mean (%)')
            if j==0:
                ax.set_xscale('log');ax.set_yscale('log');ax.set_xticks([.5,1,2,5,10,15]);numeric_ticks(ax)
                ax.set_yticks([.1,.2,.5,1,2,5,10] if rowoff==0 else [10,20,50,100,200]);numeric_ticks(ax,'y')
                ax.set_xlabel('Mean reference RMSF (Å)')
                ax.set_ylabel(ax.get_ylabel())
            label(ax,'abcd'[rowoff+j],TN[task]+(': amplitude and error' if rowoff==0 else ': relative error'))
    for j,ax in enumerate([e1,e2]):
        task=TASKS[j];v=aa[(aa.variant=='full')&(aa.metric==task)&(aa.dataset=='all')].set_index('amino_acid').reindex(AA)
        assert v.valid_residues.sum()==(961321 if j==0 else 824988)
        ax.plot(range(20),v.record_macro_mae,color=MC['full'],lw=1.1,marker='o',ms=2.5)
        ax.set_ylabel(TN[task]+' MAE\n('+UN[task]+')');ax.set_xlim(-.6,19.6)
        ax.set_xticks(range(20),AA if j==1 else ['']*20);ax.tick_params(axis='x',length=2)
        if j==1:ax.set_xlabel('Amino acid (ordered by χ count)')
        for a,r in v.iterrows():stats.append({'panel':'e','metric':task,'amino_acid':a,'n_records_with_type':int(r.records_valid),'valid_residues':int(r.valid_residues),'full_record_macro_mae':r.record_macro_mae})
    label(e1,'e','Full: within-type error')
    ax=axs[5]
    handles=[]
    for m in ['intrachain','static_sequence_complex']:
        for j,task in enumerate(TASKS):
            v=pairs[(pairs.comparator==m)&(pairs.metric==task)&(pairs.dataset=='all')].set_index('amino_acid').reindex(AA)
            valid=v.paired_records_defined>0
            assert (v.loc[valid,'comparator_record_mean_common_support']>0).all()
            v['relative_gain_percent']=100*v.record_macro_improvement/v.comparator_record_mean_common_support
            ax.plot(range(20),v.relative_gain_percent,color=MC[m],ls='-' if j==0 else '--',marker='o' if j==0 else '^',ms=2.6,lw=1.0)
            handles.append(Line2D([],[],color=MC[m],ls='-' if j==0 else '--',marker='o' if j==0 else '^',ms=3,lw=1,label=TN[task]+' / '+ML[m]))
            for a,r in v.iterrows():stats.append({'panel':'f','metric':task,'comparator':m,'amino_acid':a,'n_paired_records':int(r.paired_records_defined),'valid_residues':int(r.valid_residues),'relative_gain_percent':r.relative_gain_percent,'full_mae_common_support':r.full_record_mean_common_support,'comparator_mae_common_support':r.comparator_record_mean_common_support})
    ax.axhline(0,color=GREY,ls=':',lw=.8);ax.set_ylabel('MAE reduction (%)');ax.set_xticks(range(20),AA);ax.set_xlim(-.6,19.6)
    ax.set_yscale('symlog',linthresh=10)
    ax.set_yticks([-5,0,5,10,20,50,100]);numeric_ticks(ax,'y')
    ax.set_ylabel('MAE reduction (%)')
    ax.set_xlabel('Amino acid (ordered by χ count)');label(ax,'f','Full: within-type gain')
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(0,-.27),ncol=2,fontsize=7.5,columnspacing=.8,handlelength=1.4)
    fig.subplots_adjust(bottom=.095)  # individual axes set above retain their positions
    write(full,'fig3_all_6126_full_metrics_and_reference.csv');write(pd.DataFrame(plotted),'fig3_all_12252_scatter_points.csv')
    write(aa[(aa.variant=='full')&(aa.dataset=='all')],'fig3_full_amino_acid_source.csv')
    write(pairs[(pairs.comparator.isin(['intrachain','static_sequence_complex']))&(pairs.dataset=='all')],'fig3_amino_acid_paired_source.csv')
    caption='图3：主链波动幅度、侧链状态多样性与预测误差的关系，以及同类残基内的结构信息收益。（a,b）全部3,063条主评价记录的参考标签均值与Full MAE，每点为一条记录；（c,d）相同记录的MAE除以各自参考标签均值并乘100%，用于区分绝对误差与相对动态幅度的误差，RMSF相关坐标采用对数轴。三个来源仅用颜色与形状标识，不包含普通外部集成预测。（e）Full在同一氨基酸类型内的记录等权MAE，RMSF与Hχ用上下两条独立坐标显示，分别使用Å和nats，不混合两种单位；Ala和Gly无侧链χ角，其Hχ不适用。（f）在相同有效记录支持上，相对Intra及Static的MAE降幅，定义为100%×（对照MAE−Full MAE）/对照MAE；颜色区分对照，圆点实线表示RMSF，三角虚线表示Hχ，百分比坐标在±10%内线性、远端为对数。类型内分数先在每条记录中计算，再对存在该类型有效标签的记录等权平均。所有记录和有效氨基酸类型均纳入，未抽样、绘制置信区间或进行显著性检验。'
    export(fig,'fig3_amplitude_and_residue_types',caption,['a RMSF amplitude vs full MAE','b Hχ amplitude vs full MAE','c RMSF amplitude vs relative MAE','d Hχ amplitude vs relative MAE','e full within-type MAE, two separate unit axes','f paired within-type relative percent gain'],stats,[str(MAIN/'per_record_metrics.csv'),str(LABEL/'per_record_labels.csv'),str(MAIN/'amino_acid_support_and_macro_metrics.csv'),str(MAIN/'amino_acid_paired_mae_support_and_macros.csv')])

def fig4(external):
    d=external[external.protocol.eq('strict_single_seed7')&external.variant.ne('backflip2_official')].copy()
    d['variant']=d.variant.replace({'matched_static':'static_sequence_complex','matched_aa':'aa_only'})
    assert len(d)==32*8 and d.family_components.eq(1).all()
    write(d,'fig4_all_32_strict_single_model_metrics.csv')
    fig,axs=grid(100,rows=1);stats=[]
    fig.subplots_adjust(top=.865)
    for j,task in enumerate(TASKS):
        ax=axs[j];wide=d[d.metric==task].pivot(index='sample_id',columns='variant',values='mae').reindex(columns=M).sort_index()
        assert len(wide)==32 and np.isfinite(wide).all().all()
        for _,r in wide.iterrows():ax.plot(range(4),r,color='#C5C5C5',lw=.65,alpha=.7,zorder=1)
        for i,m in enumerate(M):
            vals=wide[m];ax.scatter(np.full(32,i),vals,color=MC[m],s=15,alpha=.70,edgecolors='white',linewidths=.35,zorder=3)
            mu=vals.mean();ax.scatter(i,mu,s=50,marker='D',color=MC[m],edgecolors=BLACK,linewidths=.6,zorder=4)
            stats.append({'panel':'ab'[j],'metric':task,'variant':m,'n_records':32,'family_components':1,'mean_mae':mu,'median_mae':vals.median(),'min_mae':vals.min(),'max_mae':vals.max(),'valid_residues':int(d[(d.metric==task)&(d.variant==m)]['count'].sum())})
        ax.set_xticks(range(4),[ML[m]+'\n'+f'{wide[m].mean():.3f}' for m in M]);ax.set_xlim(-.3,3.3)
        ax.set_ylim(wide.to_numpy().min()*.91,wide.to_numpy().max()*1.09)
        if j==0:
            ax.set_yscale('log');ax.set_yticks([.2,.5,1,2]);numeric_ticks(ax,'y')
        ax.set_ylabel('MAE ('+UN[task]+')');label(ax,'ab'[j],TN[task]+'：严格外部单模型验证')
        if j==0:ax.set_ylabel('MAE (Å；对数)')
    fig.text(.54,.99,'32 条记录；1 个序列关联组；菱形 = 记录等权均值',ha='center',va='top',fontsize=8)
    caption='图4：移除相关训练记录后，在32条ThermoPCD记录上的严格外部验证。（a,b）RMSF与Hχ的逐记录MAE；Full、Intra、Static和AA均使用严格保留的2,431条训练记录建立对应模型或对照。每条灰线连接同一记录的四种方法，小圆点为该记录误差，大菱形及方法名下方数字为32条记录等权均值；RMSF采用对数误差轴以同时保留低误差记录和完整的高误差尾部。Full与Intra各为种子7的单模型，未与普通五模型集成方案混合。32条外部记录归属同一序列关联组，因此这里只进行描述性比较，不把它们视为32个独立家族重复，不绘制置信区间或报告显著性检验。'
    export(fig,'fig4_strict_external_paired',caption,['a RMSF paired MAE, strict single seed7','b Hχ paired MAE, strict single seed7'],stats,[str(EXT/'all_32_per_record_metrics.csv')])

def figS1(chains):
    records=read(BASE/'experiments/corrected_pipeline/sequence_audit/record_overlap_audit.csv')
    comp=read(COHORT/'all_component_sizes.csv');strict=read(COHORT/'strict_membership_by_source.csv')
    assert len(records)==3095 and len(comp)==1664 and comp.records.sum()==3095
    assert records[records.role=='external'].component_id.nunique()==1
    rec=records.merge(chains.groupby('sample_id').size().rename('input_chains'),on='sample_id',validate='one_to_one')
    rec['chain_group']=np.where(rec.input_chains==1,'1',np.where(rec.input_chains==2,'2','3+'))
    write(rec,'figS1_all_3095_record_composition.csv');write(comp,'figS1_all_1664_components.csv');write(strict,'figS1_strict_membership_by_source.csv')
    fig,axs=grid(143,rows=2,bottom=.155);stats=[];colors=[CATEGORICAL[0],CATEGORICAL[3],CATEGORICAL[4]]
    ax=axs[0]
    for i,s in enumerate(SRC+['thermopcd']):
        v=rec[rec.dataset==s];left=0
        for k,g in enumerate(['1','2','3+']):
            n=int(v.chain_group.eq(g).sum());frac=100*n/len(v)
            ax.barh(i,frac,left=left,color=colors[k],height=.58,edgecolor='white',linewidth=.5)
            if frac>=10:ax.text(left+frac/2,i,str(n),ha='center',va='center',color='white' if k!=1 else BLACK,fontsize=8)
            stats.append({'panel':'a','dataset':s,'input_chain_group':g,'n_records':n,'fraction_percent':frac})
            left+=frac
    ax.set_yticks(range(4),[SN[s]+'\nn = '+str(len(rec[rec.dataset==s])) for s in SRC+['thermopcd']]);ax.invert_yaxis();ax.set_xlim(0,100)
    ax.set_xlabel('记录比例 (%)');ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    label(ax,'a','输入结构的链组成')
    ax.legend(handles=[Patch(color=colors[i],label=g+' 条链') for i,g in enumerate(['1','2','3+'])],loc='upper left',bbox_to_anchor=(0,-.28),ncol=3,fontsize=7.5,columnspacing=.8,handlelength=1)
    ax=axs[1]
    for i,s in enumerate(SRC+['thermopcd']):
        v=chains[chains.dataset==s].sort_values('length');y=np.arange(1,len(v)+1)/len(v)
        ax.step(v.length,y,where='post',color=SC[s],ls=['-','--','-.',':'][i],lw=1.2)
        stats.append({'panel':'b','dataset':s,'n_chain_observations':len(v),'min_length':int(v.length.min()),'median_length':v.length.median(),'max_length':int(v.length.max())})
    ax.set_xscale('log');ax.set_xticks([1,10,100,1000]);numeric_ticks(ax);ax.set_ylim(0,1.02)
    ax.set_xlabel('链长（残基数；对数）');ax.set_ylabel('累计链比例');label(ax,'b','全部 5,068 条链的长度分布')
    ax.legend(handles=[Line2D([],[],color=SC[s],ls=['-','--','-.',':'][i],lw=1.1,label=SN[s]) for i,s in enumerate(SRC+['thermopcd'])],loc='upper left',bbox_to_anchor=(0,-.28),ncol=2,fontsize=7.5,columnspacing=1)
    ax=axs[2];ax.plot(comp['rank'],comp.records,color=GREY,lw=.7);ax.scatter(comp['rank'],comp.records,s=3,color=GREY,linewidths=0)
    e=comp[comp.external_records>0];assert len(e)==1 and e.records.iloc[0]==664
    ax.scatter(e['rank'],e.records,s=30,color=BLACK,marker='D',zorder=5)
    ax.annotate('632 训练 + 32 外部',(1,664),xytext=(5,230),fontsize=8,arrowprops={'arrowstyle':'-','lw':.65,'color':BLACK})
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xticks([1,10,100,1000]);numeric_ticks(ax);ax.set_yticks([1,10,100,1000]);numeric_ticks(ax,'y')
    ax.set_ylim(.8,1100);ax.set_xlabel('关联组规模排序（对数）');ax.set_ylabel('每组记录数（对数）');label(ax,'c','序列与体系关联的 1,664 个组')
    stats.append({'panel':'c','n_components':1664,'n_records':3095,'external_component_records':664,'external_records':32,'related_training_records':632})
    ax=axs[3]
    for i,s in enumerate(SRC):
        own=strict[strict.dataset==s];keep=int(own.loc[~own.purged_training_member,'records'].sum());remove=int(own.loc[own.purged_training_member,'records'].sum())
        ax.barh(i,keep,color=CATEGORICAL[0],height=.56);ax.barh(i,remove,left=keep,color='#BBBBBB',height=.56,hatch='///',edgecolor='white',linewidth=.4)
        if keep:ax.text(keep/2,i,str(keep),ha='center',va='center',color='white',fontsize=8)
        if remove:ax.text(keep+remove/2,i,str(remove),ha='center',va='center',color=BLACK,fontsize=8)
        stats.append({'panel':'d','dataset':s,'retained_training_records':keep,'removed_training_records':remove})
    assert sum(r.get('retained_training_records',0) for r in stats)==2431
    ax.set_yticks(range(3),[SN[s] for s in SRC]);ax.invert_yaxis();ax.set_xlabel('训练记录数');ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    label(ax,'d','严格外部训练：移除相关记录')
    ax.legend(handles=[Patch(color=CATEGORICAL[0],label='保留 2,431'),Patch(facecolor='#BBBBBB',hatch='///',label='移除 632')],loc='upper left',bbox_to_anchor=(0,-.27),ncol=2,fontsize=7.5)
    caption='图S1：数据组成及严格外部验证的序列关系。（a）三个主评价来源和ThermoPCD外部集的实际输入链数组成，柱内数字为记录数；主评价共3,063条，外部32条。（b）全部5,068条输入链的链长经验累积分布，链长为结构中保留的氨基酸数。（c）按序列和体系关系连接的1,664个关联组及其规模，完整显示全部3,095条记录；黑色菱形所在组包含632条训练记录及全部32条外部记录。（d）按预设关联规则移除632条相关训练记录后，各来源保留及移除的数量，严格训练集合共2,431条。关联组用于说明当前数据关系，不将普通五折划分写成家族完全独立的测试，也不把32条外部记录视为独立家族。此图已删除原标签修复和覆盖率审计面板，所有组成数据均完整纳入。'
    export(fig,'figS1_composition_and_sequence_relations',caption,['a actual input chain composition','b all chain-length ECDF','c all sequence/system component sizes','d strict retained/removed records'],stats,[str(COHORT/'all_chain_lengths_and_overlap.csv'),str(BASE/'experiments/corrected_pipeline/sequence_audit/record_overlap_audit.csv'),str(COHORT/'all_component_sizes.csv'),str(COHORT/'strict_membership_by_source.csv')])

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    per,chains,labels,aa,pairs,external=prepare()
    fig2(per);fig3(per,labels,aa,pairs);fig4(external);figS1(chains)
    # Retain copies of the original plotting scripts only in the new output folder.
    refs=[BASE/'figures/figure_revision_20261005/fig2/plot.py',BASE/'figures/figure_revision_20261005/fig3/plot.py',BASE/'figures/prodesign_inspired_analysis_20261006/plot_reference_amplitude.py',BASE/'figures/corrected_cohort_lineage/make_corrected_cohort_lineage.py']
    provenance=OUT/'original_source_scripts';provenance.mkdir(exist_ok=True)
    for i,p in enumerate(refs):shutil.copyfile(p,provenance/(str(i+1)+'_'+p.name));EVIDENCE[str(p)]=sha(p)
    for p,h in EVIDENCE.items():assert sha(p)==h, 'Input changed during render'
    MANIFEST['source_sha256']=EVIDENCE;MANIFEST['script_sha256']=sha(__file__)
    MANIFEST['data_scope']={'primary_records':3063,'primary_valid_rmsf':961321,'primary_valid_hchi':824988,'external_records':32,'external_family_components':1}
    (OUT/'figure_manifest.json').write_text(json.dumps(MANIFEST,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({'generated':list(MANIFEST['figures']),'main_records':3063,'strict_external':32,'source_files':len(EVIDENCE)},ensure_ascii=False))

if __name__=='__main__' and False:main()
