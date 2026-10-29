from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullLocator, MaxNLocator
from matplotlib import font_manager
from PIL import Image

OUT=Path(__file__).resolve().parent
INPUT=OUT.parent/'trajectory_analysis'
STEM='figS12_partition_sampling'
BASE_SCRIPT=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\figures\redraw_quantitative.py')
BASE_PREVIEW=BASE_SCRIPT.parent/'fig4_strict_external_paired.png'
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


def read(path,expected=None,**kwargs):
    actual=sha(path)
    if expected is not None:
        assert actual==expected,(path,'input hash changed')
    EVIDENCE[str(path)]=actual
    return pd.read_csv(path,**kwargs)


def write(frame,name):
    frame.to_csv(OUT/name,index=False,float_format='%.17g',na_rep='')


def label(ax,letter,title):
    ax.set_title(title,loc='left',pad=9)
    ax.text(-.18,1.055,letter,transform=ax.transAxes,fontsize=11,weight='bold',ha='left',va='bottom')


def numeric_ticks(ax):
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value,position:format(value,'g')))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))


def scatter(ax,panel,identity,scope,category,x,y,color,marker='o',size=8,alpha=.6):
    x=np.asarray(x,dtype=float)
    y=np.asarray(y,dtype=float)
    assert x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all()
    collection=ax.scatter(x,y,color=color,marker=marker,s=size,alpha=alpha,edgecolors='white' if scope=='mean' else 'none',
                          linewidths=.4 if scope=='mean' else 0,zorder=6 if scope=='mean' else 3)
    np.testing.assert_allclose(np.asarray(collection.get_offsets()),np.column_stack([x,y]),rtol=0,atol=1e-14)
    for row,px,py in zip(identity,x,y):
        POINTS.append(dict(panel=panel,scope=scope,category=category,x=float(px),y=float(py),**row))


def paired_line(ax,panel,sample,dataset,x,y,pair_scope='record',color='#C5C5C5',lw=.55,alpha=.64,ls='-'):
    x=np.asarray(x,dtype=float)
    y=np.asarray(y,dtype=float)
    assert x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all()
    line,=ax.plot(x,y,color=color,lw=lw,alpha=alpha,ls=ls,zorder=2 if pair_scope=='record' else 5)
    np.testing.assert_array_equal(line.get_xdata(),x)
    np.testing.assert_array_equal(line.get_ydata(),y)
    for vertex,(px,py) in enumerate(zip(x,y)):
        LINES.append(dict(panel=panel,pair_scope=pair_scope,sample_id=sample,dataset=dataset,vertex=vertex,x=float(px),y=float(py)))


def stat(panel,category,dataset,values):
    values=np.asarray(values,dtype=float)
    assert np.isfinite(values).all() and len(values)>0
    STATS.append(dict(panel=panel,category=category,dataset=dataset,n_records=len(values),mean=float(values.mean()),
        record_sample_sd=float(values.std(ddof=1)),median=float(np.median(values)),minimum=float(values.min()),maximum=float(values.max())))


def tight_y(ax,values,zero=False):
    values=np.asarray(values,dtype=float)
    lo,hi=float(values.min()),float(values.max())
    spread=hi-lo
    ax.set_ylim(0 if zero else max(0,lo-spread*.1),hi+spread*.12)
    numeric_ticks(ax)


def prepare():
    audit_path=INPUT/'completion_audit.json'
    audit=json.loads(audit_path.read_text(encoding='utf-8'))
    assert audit['status']=='complete' and audit['records']==83 and audit['trajectories']==203
    verify_path=INPUT/'independent_verification_audit.json'
    verify=json.loads(verify_path.read_text(encoding='utf-8'))
    assert verify['status']=='passed' and verify['analysis_plan_sha256']==audit['plan_sha256']
    derivative_path=INPUT/'delivery_summary_audit.json'
    derivative=json.loads(derivative_path.read_text(encoding='utf-8'))
    assert derivative['status']=='complete' and derivative['frozen_plan_sha256']==audit['plan_sha256']
    for path in [audit_path,verify_path,derivative_path,INPUT/'analysis_plan.json',BASE_SCRIPT,BASE_PREVIEW]:
        EVIDENCE[str(path)]=sha(path)
    full=read(INPUT/'distribution_record_macro.csv',audit['output_tables']['distribution_record_macro']['sha256'])
    chi1=read(INPUT/'nonpro_chi1_record_macro.csv',derivative['outputs']['nonpro_chi1_record_macro.csv'])
    full=full[full.replica_mode=='replica1'].copy()
    chi1=chi1[chi1.replica_mode=='replica1'].copy()
    assert len(full)==len(chi1)==83*4
    assert not full.duplicated(['sample_id','shift_degrees']).any()
    assert not chi1.duplicated(['sample_id','shift_degrees']).any()
    for frame in [full,chi1]:
        assert set(frame.shift_degrees)=={0,30,60,90}
        assert frame.groupby('shift_degrees').size().to_dict()=={0:83,30:83,60:83,90:83}
        assert frame[frame.shift_degrees==0].dataset.value_counts().to_dict()=={'atlas':30,'itrace':30,'dynarepo':23}
    assert (chi1.nonpro_chi1_residues>0).all()
    repeat=read(INPUT/'repeat_jsd_per_residue.csv.gz',audit['output_tables']['repeat_jsd_per_residue']['sha256'],
                dtype={'sample_id':str,'chain':str,'author_residue_id':str})
    annotations=read(INPUT/'per_residue_trajectory.csv.gz',audit['output_tables']['per_residue_trajectory']['sha256'],
        usecols=['sample_id','residue_index','chi1_conventional_eligible'],dtype={'sample_id':str})
    annotations=annotations.drop_duplicates()
    assert not annotations.duplicated(['sample_id','residue_index']).any()
    repeat=repeat.merge(annotations,on=['sample_id','residue_index'],validate='many_to_one')
    assert repeat.chi1_conventional_eligible.notna().all()
    repeat=repeat[repeat.chi1_conventional_eligible==1].copy()
    per_pair=repeat.groupby(['sample_id','dataset','fold','shift_degrees','replica_first','replica_second'])['chi1_repeat_jsd'].mean().reset_index()
    assert per_pair.groupby(['sample_id','shift_degrees']).size().eq(3).all()
    repeats=per_pair.groupby(['sample_id','dataset','fold','shift_degrees'])['chi1_repeat_jsd'].mean().reset_index()
    assert len(repeats)==60*4 and set(repeats.dataset)=={'atlas','itrace'}
    matched=chi1[chi1.sample_id.isin(repeats.sample_id)].merge(repeats,on=['sample_id','dataset','fold','shift_degrees'],validate='one_to_one')
    assert len(matched)==60*4
    assert matched.groupby('shift_degrees').size().eq(60).all()
    write(full,STEM+'.all83_jointH_source.csv')
    write(chi1,STEM+'.all83_nonpro_chi1_source.csv')
    write(per_pair,STEM+'.all60_nonpro_chi1_repeat_pair_source.csv')
    write(matched,STEM+'.all60_nonpro_chi1_time_repeat_source.csv')
    return full,chi1,matched


def wide(frame,column,offsets):
    result=frame[frame.shift_degrees.isin(offsets)].pivot(index=['sample_id','dataset'],columns='shift_degrees',values=column)
    result=result.reindex(columns=offsets).sort_index()
    assert result.notna().all().all()
    return result


def summary_overlay(ax,panel,values,categories,positions,ordinal=True):
    sources=[s for s in SRC if s in values.index.get_level_values('dataset')]
    mean_offsets=dict(zip(sources,np.linspace(-.105,.105,len(sources))))
    for source in sources+['all']:
        subset=values if source=='all' else values.xs(source,level='dataset')
        y=np.asarray(subset.mean(axis=0),dtype=float)
        x=np.asarray(positions,dtype=float)+(0 if source=='all' or not ordinal else mean_offsets[source])
        color=BLACK if source=='all' else SC[source]
        marker='D' if source=='all' else SM[source]
        line_style='-' if source=='all' else SL[source]
        paired_line(ax,panel,'record_mean',source,x,y,'mean',color=color,lw=1.4 if source=='all' else 1.0,alpha=1,ls=line_style)
        for category,px,py,raw in zip(categories,x,y,subset.to_numpy().T):
            scatter(ax,panel,[dict(sample_id='record_mean',dataset=source,n_records=len(subset))],'mean',category,[px],[py],color,marker,size=32 if source=='all' else 23,alpha=1)
            stat(panel,category,source,raw)


def paired_panel(ax,panel,frame,column,title,ylabel,zero=False):
    values=wide(frame,column,[0,60])
    assert len(values)==83
    jitter=(np.arange(len(values))%17/16-.5)*.12
    for index,((sample,source),row) in enumerate(values.iterrows()):
        paired_line(ax,panel,sample,source,[jitter[index],1+jitter[index]],row.to_numpy())
    for category,x in [(0,0),(60,1)]:
        ids=[dict(sample_id=sample,dataset=source,shift_degrees=category) for sample,source in values.index]
        scatter(ax,panel,ids,'record',str(category),x+jitter,values[category].to_numpy(),GREY,size=7.5,alpha=.56)
    summary_overlay(ax,panel,values,['0','60'],[0,1])
    ax.set_xticks([0,1],['原分区','常规 χ1 中心'])
    ax.set_xlim(-.27,1.27)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('非 Pro χ1；83 条记录')
    tight_y(ax,values.to_numpy(),zero=zero)
    label(ax,panel,title)


def draw(full,chi1,matched):
    fig=plt.figure(figsize=(183/25.4,166/25.4))
    gs=fig.add_gridspec(2,2,left=.12,right=.96,bottom=.15,top=.875,wspace=.43,hspace=.67)
    axs=[fig.add_subplot(gs[i//2,i%2]) for i in range(4)]
    handles=[Line2D([],[],color=SC[s],marker=SM[s],ls=SL[s],lw=1.1,markersize=4,label=SN[s]+' 均值') for s in SRC]
    handles.append(Line2D([],[],color=BLACK,marker='D',lw=1.4,markersize=4,label='总体均值'))
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.985),ncol=4,columnspacing=1.1,handlelength=2)
    values=wide(full,'h_joint',[0,30,60,90])
    assert len(values)==83
    ax=axs[0]
    for (sample,source),row in values.iterrows():
        paired_line(ax,'a',sample,source,[0,30,60,90],row.to_numpy())
    for shift in [0,30,60,90]:
        ids=[dict(sample_id=sample,dataset=source,shift_degrees=shift) for sample,source in values.index]
        scatter(ax,'a',ids,'record',str(shift),np.full(83,shift),values[shift].to_numpy(),GREY,size=7.5,alpha=.5)
    summary_overlay(ax,'a',values,['0','30','60','90'],[0,30,60,90],ordinal=False)
    ax.set_xticks([0,30,60,90],['0°','30°','60°','90°'])
    ax.set_xlim(-7,97)
    ax.set_xlabel('等宽角区间整体平移；83 条记录')
    ax.set_ylabel('完整联合 Hχ (nats)')
    tight_y(ax,values.to_numpy())
    label(ax,'a','联合 χ 状态的角区间敏感性')
    paired_panel(axs[1],'b',chi1,'h_chi1','χ1 占据熵与区间中心','χ1 边缘熵 (nats)')
    paired_panel(axs[2],'c',chi1,'apparent_chi1_changes','χ1 保存帧类别变化','变化次数 / 501 个保存帧',zero=True)
    ax=axs[3]
    rows=[]
    for (sample,source),part in matched.groupby(['sample_id','dataset'],sort=True):
        part=part.set_index('shift_degrees')
        rows.append(dict(sample_id=sample,dataset=source,half0=part.loc[0,'half_chi1_jsd'],half60=part.loc[60,'half_chi1_jsd'],
                         repeat0=part.loc[0,'chi1_repeat_jsd'],repeat60=part.loc[60,'chi1_repeat_jsd']))
    values=pd.DataFrame(rows).set_index(['sample_id','dataset'])[['half0','half60','repeat0','repeat60']]
    assert len(values)==60 and np.isfinite(values).all().all()
    jitter=(np.arange(60)%17/16-.5)*.12
    x_positions=[0,1,3,4]
    for index,((sample,source),row) in enumerate(values.iterrows()):
        paired_line(ax,'d',sample,source,[jitter[index],1+jitter[index]],row[['half0','half60']].to_numpy())
        paired_line(ax,'d',sample,source,[3+jitter[index],4+jitter[index]],row[['repeat0','repeat60']].to_numpy())
    for category,x in zip(values.columns,x_positions):
        ids=[dict(sample_id=sample,dataset=source,comparison=category) for sample,source in values.index]
        scatter(ax,'d',ids,'record',category,x+jitter,values[category].to_numpy(),GREY,size=7.5,alpha=.56)
    # The two comparison types are not connected: each paired line only joins
    # original and shifted60 definitions within the same sampling comparison.
    summary_overlay(ax,'d',values[['half0','half60']],['half0','half60'],[0,1])
    summary_overlay(ax,'d',values[['repeat0','repeat60']],['repeat0','repeat60'],[3,4])
    ax.set_xticks(x_positions,['半程\n原分区','半程\n常规中心','重复间\n原分区','重复间\n常规中心'])
    ax.set_xlim(-.4,4.4)
    ax.set_ylabel('χ1 占据分布 JSD (nats)')
    ax.set_xlabel('同 60 条多重复记录；非 Pro χ1')
    tight_y(ax,values.to_numpy(),zero=True)
    label(ax,'d','χ1 占据分布：时间与重复')
    fig.text(.12,.055,'灰点与连线：全部记录；彩色符号：来源均值；黑菱形：记录等权均值。',ha='left',va='center',fontsize=7.5)
    fig.text(.12,.027,'统一轨迹窗口 0–100 ns；每 200 ps 读取一帧；主分析每记录使用第一条轨迹。',ha='left',va='center',fontsize=7.5)
    return fig,axs


def export(fig,axs):
    points=pd.DataFrame(POINTS)
    lines=pd.DataFrame(LINES)
    stats=pd.DataFrame(STATS)
    assert points.groupby(['panel','scope']).size().to_dict()=={('a','record'):332,('a','mean'):16,('b','record'):166,('b','mean'):8,('c','record'):166,('c','mean'):8,('d','record'):240,('d','mean'):12}
    assert points[points.scope=='record'].groupby('panel').sample_id.nunique().to_dict()=={'a':83,'b':83,'c':83,'d':60}
    record=points[points.scope=='record']
    for index,panel in enumerate('abcd'):
        selected=record[record.panel==panel]
        xlo,xhi=axs[index].get_xlim()
        ylo,yhi=axs[index].get_ylim()
        assert selected.x.between(xlo,xhi).all() and selected.y.between(ylo,yhi).all()
    write(points,STEM+'.points.csv')
    write(lines,STEM+'.paired_line_vertices.csv')
    write(stats,STEM+'.statistics.csv')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    outside_ticks=set()
    for ax in fig.axes:
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
    assert not outside,('text outside canvas',outside)
    assert min(t.get_fontsize() for t in texts)>=7.5
    for extension in ['png','pdf','svg']:
        fig.savefig(OUT/(STEM+'.'+extension),bbox_inches=fig.bbox_inches,dpi=300)
    svg=ET.parse(OUT/(STEM+'.svg'))
    namespace={'s':'http://www.w3.org/2000/svg'}
    assert not svg.findall('.//s:image',namespace)
    preview=Image.open(OUT/(STEM+'.png'))
    assert min(preview.info['dpi'])>=299
    pixels=np.asarray(preview.convert('RGB'))
    density=[]
    for ax in axs:
        p=ax.get_position()
        crop=pixels[int((1-p.y1)*preview.height):int((1-p.y0)*preview.height),int(p.x0*preview.width):int(p.x1*preview.width)]
        density.append(float(np.any(crop<245,axis=2).mean()))
    caption='图 S12：侧链占据统计对角区间定义及轨迹采样的依赖。（a）83条冻结记录在0°、30°、60°和90°整体平移的等宽角区间下，完整联合χ状态熵的逐记录均值；Pro及所有原有定义χ角保留在联合状态中，各种平移均作为角区间敏感性比较。（b）非Pro残基的χ1边缘占据熵，比较原分区与以−60°、+60°和180°为中心的常规χ1分区；圆周接缝两侧的trans区间合并为同一态。（c）同一χ1序列在两种分区下的相邻保存帧类别变化次数，每条轨迹包含501帧；该计数按保存帧定义。（d）在同一批60条ATLAS/ITrace多重复记录中，比较χ1前后两个不重叠半程的占据分布Jensen–Shannon散度，以及三条重复轨迹全部三对占据分布的平均散度；均排除Pro，单位为nats。半程比较使用第一条轨迹，时间段分别为0–49.8 ns和50–100 ns。所有轨迹统一取0–100 ns、200 ps间隔的501个已保存帧；（a–c）主分析每条记录均使用第一条轨迹，来源为ATLAS 30条、DynaRepo 23条、ITrace 30条；全部203条可用轨迹另保留用于重复分析，其中ATLAS/ITrace各30条记录均有三条轨迹，DynaRepo各记录有一条轨迹。数值先在轨迹内对相应有效残基取平均，再对记录等权取平均。灰点与配对连线覆盖全部对应记录，彩色符号为来源均值，黑菱形为总体记录均值；（b–c）的83条记录均存在有效非Pro χ1残基。'
    (OUT/(STEM+'.caption_zh.md')).write_text(caption+'\n',encoding='utf-8')
    metadata=dict(status='rendered_pending_visual_QA',figure='S12',scientific_claim='Defined angular partitions change observed occupancy entropy and saved-frame category changes; sampling comparisons characterize occupancy distributions',
        archetype='quantitative_grid',copy_first_source=str(BASE_SCRIPT),copy_first_source_sha256=sha(BASE_SCRIPT),copy_first_preview_sha256=sha(BASE_PREVIEW),
        plot_code_sha256=sha(__file__),input_sha256=EVIDENCE,caption_zh=caption,
        traceability='Every individual scatter point and paired-line vertex is recorded with sample/source/category; all saved scatter offsets and line vertices asserted against source values',
        dimensions_mm=[float(v*25.4) for v in fig.get_size_inches()],minimum_text_pt=min(t.get_fontsize() for t in texts),
        png_dpi=list(preview.info['dpi']),png_pixels=list(preview.size),editable_svg_text=len(svg.findall('.//s:text',namespace)),svg_raster_images=0,
        text_outside_canvas=outside,all_points_inside_axes=True,panel_content_density=density,
        record_points=record.groupby('panel').size().to_dict(),record_counts=record.groupby('panel').sample_id.nunique().to_dict(),
        n_definition='Record means, not residues treated as replicates;83records primary and60matched multi-replica records forpaneld',
        center='Arithmetic record mean; per-record average over applicable residues; pairwise-repeat comparisons average allthree repeat pairs',
        spread='Individual record values visible; sampleSD,median,range in statistics.csv; no CI/errorbars drawn',
        tests='Descriptive comparisons; no significance test, p-value or multiple-testing correction',
        source_subsampling=0,new_record_selection=0,Pro_conventional_chi1_exclusion='Inherited explicit non-Pro interpretation rule; complete83 retained forpanelsb/c',
        outputs={})
    for extension in ['png','pdf','svg','points.csv','paired_line_vertices.csv','statistics.csv','caption_zh.md',
                      'all83_jointH_source.csv','all83_nonpro_chi1_source.csv','all60_nonpro_chi1_repeat_pair_source.csv','all60_nonpro_chi1_time_repeat_source.csv']:
        file=OUT/(STEM+'.'+extension)
        metadata['outputs'][extension]=dict(path=str(file),sha256=sha(file))
    (OUT/(STEM+'.audit.json')).write_text(json.dumps(metadata,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    plt.close(fig)
    print(json.dumps({k:v for k,v in metadata.items() if k not in ['input_sha256','outputs','caption_zh']},ensure_ascii=False,indent=2))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    full,chi1,matched=prepare()
    fig,axs=draw(full,chi1,matched)
    export(fig,axs)


if __name__=='__main__':
    main()
