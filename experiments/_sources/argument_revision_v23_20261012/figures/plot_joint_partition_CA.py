# Copy-first native Matplotlib visual adaptation of redraw_quantitative.py; original SHA e2a2710ee16499d8f577a8d4c48a7779714f11fc9c3d6471ce5c4d285789cc2d.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'completed_results_20261012' / 'plot_runtime'))
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
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter,MaxNLocator,NullLocator
from matplotlib import font_manager
from PIL import Image

OUT=Path(__file__).resolve().parent
TASK=OUT.parents[1]
TRAJ=TASK/'energy_landscape_v21_20261012'/'trajectory_analysis'
V22=TASK/'energy_landscape_v22_maintext_20261012'
SOURCE_SCRIPT=TASK/'energy_landscape_v21_20261012'/'figures'/'plot_partition_sampling.py'
STEM='fig8_joint_partition_CA'
SRC=['atlas','dynarepo','itrace']
SC=dict(zip(SRC,[CATEGORICAL[0],CATEGORICAL[3],CATEGORICAL[4]]))
SN=dict(zip(SRC,['ATLAS','DynaRepo','ITrace']))
SM=dict(zip(SRC,['o','^','s']))
SL=dict(zip(SRC,['-','--',':']))
EVIDENCE={}
POINTS=[]
LINES=[]
STATS=[]
font_manager.findfont('Microsoft YaHei',fallback_to_default=False)
mpl.rcParams.update({'font.sans-serif':['Microsoft YaHei','Arial'],'xtick.labelsize':7.5,
    'ytick.labelsize':7.5,'legend.fontsize':8,'axes.titlesize':8.5,'axes.labelsize':8,'axes.unicode_minus':False})


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def read(path,expected=None):
    actual=sha(path)
    if expected:
        assert actual==expected,(path,'frozen input hash changed')
    EVIDENCE[str(path)]=actual
    return pd.read_csv(path)


def write(frame,name):
    frame.to_csv(OUT/name,index=False,float_format='%.17g',na_rep='')


def label(ax,letter,title):
    ax.set_title(title,loc='left',pad=9)
    ax.text(-.18,1.055,letter,transform=ax.transAxes,fontsize=11,weight='bold',ha='left',va='bottom')


def prepare():
    audit=json.loads((TRAJ/'completion_audit.json').read_text(encoding='utf-8'))
    independent=json.loads((TRAJ/'independent_verification_audit.json').read_text(encoding='utf-8'))
    ca_audit=json.loads((V22/'representative_Fstar_and_CA_maintext_audit.json').read_text(encoding='utf-8'))
    assert audit['status']=='complete' and audit['records']==83 and audit['trajectories']==203
    assert independent['status']==ca_audit['status']=='passed'
    assert independent['analysis_plan_sha256']==audit['plan_sha256']
    for path in [TRAJ/'completion_audit.json',TRAJ/'independent_verification_audit.json',
                 TRAJ/'analysis_plan.json',V22/'representative_Fstar_and_CA_maintext_audit.json',SOURCE_SCRIPT]:
        EVIDENCE[str(path)]=sha(path)
    source_path=TRAJ/'distribution_record_macro.csv'
    full=read(source_path,audit['output_tables']['distribution_record_macro']['sha256'])
    selected=full[(full.replica_mode=='replica1')&full.shift_degrees.isin([0,60])].copy()
    assert len(selected)==166 and not selected.duplicated(['sample_id','shift_degrees']).any()
    assert selected.groupby('shift_degrees').size().to_dict()=={0:83,60:83}
    assert selected[selected.shift_degrees==0].dataset.value_counts().to_dict()=={'atlas':30,'itrace':30,'dynarepo':23}
    needed=['h_joint','ca_between_joint_fraction','residues']
    assert np.isfinite(selected[needed]).all().all()
    denominator_path=V22/'CA_joint_state_decomposition_per_record.csv'
    den=read(denominator_path,ca_audit['new_csv_outputs'][denominator_path.name])
    joined=selected.merge(den,on=['sample_id','dataset','shift_degrees'],validate='one_to_one')
    assert len(joined)==166
    np.testing.assert_allclose(joined.ca_between_joint_fraction,joined.between_fraction_record_mean,rtol=0,atol=1e-14)
    np.testing.assert_array_equal(joined.residues.to_numpy(),joined.valid_residues.to_numpy())
    assert joined.groupby('shift_degrees').valid_residues.sum().to_dict()=={0:33493,60:33493}
    ids0=joined[joined.shift_degrees==0][['sample_id','valid_residues']].sort_values('sample_id')
    ids60=joined[joined.shift_degrees==60][['sample_id','valid_residues']].sort_values('sample_id')
    np.testing.assert_array_equal(ids0.to_numpy(),ids60.to_numpy())
    for shift,expected in [(0,[1.484190947286289,.062420048040040955]),(60,[.9618438253602354,.07769458266934855])]:
        part=joined[joined.shift_degrees==shift]
        # Hχ rounded published means independently checked to six decimals;
        # CA fraction has its exact independently audited denominator and mean.
        assert np.isclose(part.h_joint.mean(),expected[0],atol=1e-6,rtol=0)
        assert np.isclose(part.ca_between_joint_fraction.mean(),expected[1],atol=1e-14,rtol=0)
    joined['ca_between_joint_percent']=100*joined.ca_between_joint_fraction
    write(joined,STEM+'.all83_source_data.csv')
    return joined


def line(ax,panel,sample,dataset,x,y,scope='record',color='#C5C5C5',lw=.6,alpha=.7,ls='-'):
    x=np.asarray(x,dtype=float)
    y=np.asarray(y,dtype=float)
    rendered,=ax.plot(x,y,color=color,lw=lw,alpha=alpha,ls=ls,zorder=1 if scope=='record' else 5)
    np.testing.assert_array_equal(rendered.get_xdata(),x)
    np.testing.assert_array_equal(rendered.get_ydata(),y)
    uid=len(LINES)
    for vertex,(px,py) in enumerate(zip(x,y)):
        LINES.append(dict(line_id=uid,panel=panel,scope=scope,sample_id=sample,dataset=dataset,vertex=vertex,x=px,y=py))


def scatter(ax,panel,sample_ids,sources,shift,x,y,scope='record',color=GREY,marker='o',size=8,alpha=.56):
    x=np.asarray(x,dtype=float)
    y=np.asarray(y,dtype=float)
    collection=ax.scatter(x,y,color=color,marker=marker,s=size,alpha=alpha,edgecolors='white' if scope=='mean' else 'none',
                          linewidths=.4 if scope=='mean' else 0,zorder=6 if scope=='mean' else 3)
    np.testing.assert_allclose(np.asarray(collection.get_offsets()),np.column_stack([x,y]),rtol=0,atol=1e-14)
    for sample,source,px,py in zip(sample_ids,sources,x,y):
        POINTS.append(dict(panel=panel,scope=scope,sample_id=sample,dataset=source,shift_degrees=shift,x=float(px),y=float(py)))


def main():
    source=prepare()
    fig=plt.figure(figsize=(183/25.4,101/25.4))
    gs=fig.add_gridspec(1,2,left=.12,right=.96,bottom=.235,top=.785,wspace=.43)
    axs=[fig.add_subplot(gs[0,i]) for i in range(2)]
    handles=[Line2D([],[],color=SC[s],marker=SM[s],ls=SL[s],lw=1.1,markersize=4,label=SN[s]+' 均值') for s in SRC]
    handles.append(Line2D([],[],color=BLACK,marker='D',lw=1.4,markersize=4,label='总体均值'))
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.985),ncol=4,columnspacing=1.1,handlelength=2)
    titles=['完整联合 χ 状态占据熵','侧链状态条件下的 Cα 波动']
    labels=['完整联合 Hχ (nats)','状态间 / 总 Cα 均方波动 (%)']
    columns=['h_joint','ca_between_joint_percent']
    for panel,ax,column,title,ylabel in zip('ab',axs,columns,titles,labels):
        wide=source.pivot(index=['sample_id','dataset'],columns='shift_degrees',values=column).reindex(columns=[0,60]).sort_index()
        assert len(wide)==83 and np.isfinite(wide).all().all()
        jitter=(np.arange(83)%17/16-.5)*.12
        for index,((sample,dataset),row) in enumerate(wide.iterrows()):
            line(ax,panel,sample,dataset,[jitter[index],1+jitter[index]],row.to_numpy())
        for shift,xpos in [(0,0),(60,1)]:
            scatter(ax,panel,[x[0] for x in wide.index],[x[1] for x in wide.index],shift,xpos+jitter,wide[shift].to_numpy())
        offsets={'atlas':-.105,'dynarepo':-.045,'itrace':.105,'all':0}
        for dataset in SRC+['all']:
            part=wide if dataset=='all' else wide.xs(dataset,level='dataset')
            means=part.mean(axis=0).to_numpy()
            x=np.asarray([0,1])+offsets[dataset]
            color=BLACK if dataset=='all' else SC[dataset]
            marker='D' if dataset=='all' else SM[dataset]
            linestyle='-' if dataset=='all' else SL[dataset]
            line(ax,panel,'record_mean',dataset,x,means,'mean',color=color,lw=1.3 if dataset=='all' else 1.0,alpha=1,ls=linestyle)
            for shift,px,py in zip([0,60],x,means):
                scatter(ax,panel,['record_mean'],[dataset],shift,[px],[py],'mean',color,marker,size=30 if dataset=='all' else 24,alpha=1)
                v=part[shift].to_numpy()
                STATS.append(dict(panel=panel,metric=column,dataset=dataset,shift_degrees=shift,n_records=len(part),
                    effective_residues=int(source[(source.shift_degrees==shift)&(source.dataset.isin(SRC if dataset=='all' else [dataset]))].valid_residues.sum()),
                    mean=float(v.mean()),record_sample_sd=float(v.std(ddof=1)),median=float(np.median(v)),minimum=float(v.min()),maximum=float(v.max())))
        ax.set_xticks([0,1],['原分区','整体平移 60°'])
        ax.set_xlim(-.27,1.27)
        limits=wide.to_numpy()
        spread=limits.max()-limits.min()
        ax.set_ylim(max(0,limits.min()-.1*spread),limits.max()+.12*spread)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('同 83 条记录；第一条轨迹')
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value,position:format(value,'g')))
        ax.yaxis.set_minor_locator(NullLocator())
        label(ax,panel,title)
    fig.text(.12,.078,'灰点和连线：全部记录；黑菱形：记录等权均值；完整联合 χ 状态保留 Pro。',ha='left',va='center',fontsize=7.5)
    fig.text(.12,.035,'统一窗口 0–100 ns；200 ps 间隔；每条轨迹 501 帧；有效位置共 33,493 个。',ha='left',va='center',fontsize=7.5)
    points=pd.DataFrame(POINTS)
    lines=pd.DataFrame(LINES)
    stats=pd.DataFrame(STATS)
    assert points.groupby(['panel','scope']).size().to_dict()=={('a','record'):166,('a','mean'):8,('b','record'):166,('b','mean'):8}
    assert points[points.scope=='record'].groupby('panel').sample_id.nunique().to_dict()=={'a':83,'b':83}
    for index,panel in enumerate('ab'):
        data=points[points.panel==panel]
        assert data.x.between(*axs[index].get_xlim()).all() and data.y.between(*axs[index].get_ylim()).all()
    write(points,STEM+'.points.csv')
    write(lines,STEM+'.paired_line_vertices.csv')
    write(stats,STEM+'.statistics.csv')
    caption='图8：同一轨迹中侧链占据熵与Cα波动分解的角区间敏感性。83条冻结记录各取第一条轨迹，统一为0–100 ns、200 ps间隔的501帧。（a）原分区与整体平移60°等宽分区的完整联合Hχ，保留Pro。（b）以相同联合χ状态对链内对齐Cα的总均方波动（RMSF²）分解，显示状态均值间项占总量的比例。两面板使用相同33,493个有效位置，先在记录内平均，再对记录等权平均。灰线配对同一记录，彩色符号为来源均值，黑菱形为总体均值。'
    (OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf-8')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    outside_ticks=set()
    for ax in axs:
        for axis in [ax.xaxis,ax.yaxis]:
            lo,hi=sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not lo-1e-10<=tick.get_loc()<=hi+1e-10:
                    outside_ticks.update([tick.label1,tick.label2])
    texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text() and t not in outside_ticks]
    outside=[]
    for text in texts:
        box=text.get_window_extent(renderer)
        if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.width+.5 or box.y1>fig.bbox.height+.5:
            outside.append(text.get_text())
    assert not outside and min(t.get_fontsize() for t in texts)>=7.5
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/(STEM+'.'+ext),bbox_inches=fig.bbox_inches,dpi=300)
    preview=Image.open(OUT/(STEM+'.png'))
    assert min(preview.info['dpi'])>=299
    svg=ET.parse(OUT/(STEM+'.svg'))
    ns={'s':'http://www.w3.org/2000/svg'}
    assert not svg.findall('.//s:image',ns)
    pixels=np.asarray(preview.convert('RGB'))
    density=[]
    for ax in axs:
        p=ax.get_position()
        crop=pixels[int((1-p.y1)*preview.height):int((1-p.y0)*preview.height),int(p.x0*preview.width):int(p.x1*preview.width)]
        density.append(float(np.any(crop<245,axis=2).mean()))
    result=dict(status='rendered_pending_visual_QA',figure='8',scientific_claim='The same angular partition defines jointHchi occupancy and the descriptive sidechain-state-conditioned fraction ofCAmean-square fluctuation',
        archetype='quantitative_grid',copy_first_source=str(SOURCE_SCRIPT),copy_first_source_sha256=sha(SOURCE_SCRIPT),
        plot_code_sha256=sha(__file__),input_sha256=EVIDENCE,source_rows=166,records=83,
        primary='replica1',all83records_preserved=True,random_subsampling=0,all_model_outputs_and_primary_labels_unchanged=True,
        shared_effective_residues=33493,source_records={'atlas':30,'dynarepo':23,'itrace':30},
        sources_valid_positions={'atlas':5660,'dynarepo':6384,'itrace':21449},
        center='Per-residue occupancy H or between/total variance ratio; arithmetic residual-position mean withinrecord;83recordequal mean',
        spread='Individual83recordpoints; sampleSD,median,min,max documentedinstatistics; noCIorerrorbars',
        tests='Descriptive paired comparison; no significance test orpvalue',
        units={'a':'nats','b':'percent ofalignedCAmean-squarefluctuation'},
        CA_denominator='Same33493χ+CAvalidpositions; scalar RMSF² inchain-aligned coordinate time series; includesPro jointcodes; noenergetic barrier orcausal contribution asserted',
        dimensions_mm=[float(v*25.4) for v in fig.get_size_inches()],minimum_text_pt=min(t.get_fontsize() for t in texts),
        png_dpi=list(preview.info['dpi']),png_pixels=list(preview.size),svg_editable_text=len(svg.findall('.//s:text',ns)),svg_raster_images=0,
        text_outside_canvas=outside,all_points_inside_axes=True,panel_content_density=density,
        actual_record_point_counts=points[points.scope=='record'].groupby('panel').size().to_dict(),
        all_rendered_scatter_offsets_and_line_vertices_asserted_against_source=True,caption_zh=caption,
        outputs={})
    for ext in ['png','pdf','svg','all83_source_data.csv','points.csv','paired_line_vertices.csv','statistics.csv','caption_zh.md']:
        file=OUT/(STEM+'.'+ext)
        result['outputs'][ext]=dict(path=str(file),sha256=sha(file))
    (OUT/(STEM+'.audit.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    plt.close(fig)
    print(json.dumps({k:v for k,v in result.items() if k not in ['outputs','input_sha256','caption_zh']},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
