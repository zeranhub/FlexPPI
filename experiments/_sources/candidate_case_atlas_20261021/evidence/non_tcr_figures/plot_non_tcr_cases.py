# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) Existing verified PyMOL coordinates and stored camera scenes -> molecular native render
# (c,d,e) LineTrend project production composition -> param inherit; all KIX observations retained
# (b,c,d) Ubiquitin original_fig7 + figure_merge_v43 production panels -> visual adapt; all66/151/432 observations retained
# RULE: native molecular renders loaded as images; quantitative lines/points are vector data artists.
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
import sys,json,hashlib,re,shutil
BASE=Path(r'E:\Project\FlexPPI\experiments\_sources')
sys.path.insert(0,str(BASE/'completed_results_20261012/plot_runtime'))
sys.path.insert(1,str(BASE/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
OUT=BASE/'candidate_case_atlas_20261021/evidence/non_tcr_figures'
OLD=BASE/'reader_revision_v39_20261016/english_figures/fig6'
UBOLD=BASE/'figure_merge_v43_20261018/figure5/source/original_fig7'
W,H=183.,205.
BLUE,ORANGE,PURPLE=CATEGORICAL[0],CATEGORICAL[3],CATEGORICAL[4]
MD='#666666'
mpl.rcParams.update({'font.sans-serif':['Arial','Liberation Sans'],'font.size':9,
 'axes.titlesize':9,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,
 'legend.fontsize':9,'savefig.bbox':None,'axes.unicode_minus':False})
read=lambda p:pd.read_csv(p,float_precision='round_trip')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def scaffold(caseid):
    dest=OUT/caseid;dest.mkdir(exist_ok=True,parents=True)
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    axes=[];points=[]
    def ax(x,y,w,h):
        a=fig.add_axes([x/W,y/H,w/W,h/H]);axes.append(a);return a
    def text(x,y,s,**kw):
        kw.setdefault('fontsize',9);kw.setdefault('color',BLACK)
        return fig.text(x/W,y/H,s,**kw)
    def panel(letter,x,y,title):
        text(x,y,letter,fontsize=12,fontweight='bold');text(x+8,y,title,fontweight='bold')
    def picture(path,x,y,w,h):
        a=ax(x,y,w,h);a.imshow(plt.imread(path));a.axis('off');return a
    return dest,fig,axes,points,ax,text,panel,picture

def save(dest,fig,axes,points,audit,caption):
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    texts=[x for x in fig.findobj(mpl.text.Text) if x.get_visible() and x.get_text()]
    outside=[]
    for t in texts:
        b=t.get_window_extent(renderer)
        if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(t.get_text())
    assert not outside,outside
    assert min(t.get_fontsize() for t in texts)>=9
    assert not [t.get_text() for t in texts if re.search('[\u3400-\u9fff]',t.get_text())]
    stem='fig_'+dest.name.lower()+'_case'
    for ext in ['png','pdf','svg']:fig.savefig(dest/(stem+'.'+ext),dpi=450,bbox_inches=None)
    pd.DataFrame(points).to_csv(dest/'all_plotted_points.csv',index=False,float_format='%.17g')
    (dest/'caption_zh.md').write_text(caption+'\n',encoding='utf-8')
    (dest/'visible_text.json').write_text(json.dumps([t.get_text() for t in texts],ensure_ascii=False,indent=2),encoding='utf-8')
    audit.update(status='numerical_passed_visual_pending',canvas_mm=[W,H],dpi=450,minimum_font_pt=9,
                 visible_language='English',axes=[list(a.get_position().bounds) for a in axes],
                 point_rows=len(points),text_outside_canvas=outside,
                 source_paths_unmodified=True,statistical_tests='None; descriptive complete-region evaluation',
                 output_sha256={ext:sha(dest/(stem+'.'+ext)) for ext in ['png','pdf','svg']})
    (dest/'statistics_and_data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    plt.close(fig)
    return dict(caseid=dest.name,png=str(dest/(stem+'.png')),pdf=str(dest/(stem+'.pdf')),svg=str(dest/(stem+'.svg')),caption=str(dest/'caption_zh.md'),audit=str(dest/'statistics_and_data_audit.json'))

def kix():
    dest,fig,axes,points,ax,text,panel,picture=scaffold('C31')
    data=read(OLD/'data/kix_fold1_all_residues.csv');data=data[data.metric=='rmsf'].copy()
    metrics=read(OLD/'data/kix_fold1_case_metrics.csv');metrics=metrics[metrics.metric=='rmsf']
    states=['apo','myb','mll','both'];names=['Apo','cMyb','MLL','Both']
    assert len(data)==348 and data.groupby('state').size().to_dict()=={s:87 for s in states}
    text(2,200,'C31  CBP KIX: locating a mobile connecting loop',fontsize=11,fontweight='bold')
    panel('a',2,190,'Binding partners')
    picture(OUT/'structure_renders/kix_context.png',2,145,56,42)
    for x,label,col in [(3,'KIX',MD),(20,'cMyb',CATEGORICAL[2]),(41,'MLL',PURPLE)]:text(x,142,label,color=col)
    panel('b',63,190,'L12 (614–621), both-partner state')
    picture(OUT/'structure_renders/kix_l12_md.png',63,145,54,42)
    picture(OUT/'structure_renders/kix_l12_model.png',123,145,54,42)
    text(89,142,'MD reference',ha='center');text(149,142,'Model prediction',ha='center',color=BLUE)
    panel('c',2,134,'All 87 KIX positions in each state')
    matrix=[];labels=[]
    for s,n in zip(states,names):
        p=data[data.state==s].sort_values('author_residue_number')
        assert list(p.author_residue_number)==list(range(586,673))
        for key,name in [('md_percentile','MD'),('model_percentile','model')]:
            matrix.append(p[key].to_numpy());labels.append(n+' '+name)
            for _,r in p.iterrows():points.append(dict(panel='c',state=s,residue=int(r.author_residue_number),kind=name,x=int(r.author_residue_number),y=float(r[key])))
    a=ax(24,100,93,28)
    im=a.imshow(matrix,origin='upper',aspect='auto',cmap=LinearSegmentedColormap.from_list('mobility',SEQUENTIAL),vmin=0,vmax=100,extent=(585.5,672.5,7.5,-.5),interpolation='nearest')
    a.set_yticks(range(8),labels);a.set_xticks([586,614,646,672]);a.tick_params(length=2,pad=2)
    a.set_xlabel('KIX residue number',labelpad=2)
    for y in [1.5,3.5,5.5]:a.axhline(y,color='white',lw=1.2)
    a.add_patch(Rectangle((613.5,-.5),8,8,fill=False,lw=1.0,ec=ORANGE))
    for spine in a.spines.values():spine.set_visible(False)
    cbax=ax(120,103,2.5,23);cb=fig.colorbar(im,cax=cbax,ticks=[0,50,100]);cb.ax.tick_params(length=2,pad=2)
    cb.set_label('RMSF percentile',labelpad=2,fontsize=9)
    panel('d',139,134,'L12 ranking')
    dm=metrics.set_index('state').loc[states]
    a=ax(146,102,32,26)
    a.axhspan(80,100,color=ORANGE,alpha=.12,lw=0)
    a.plot(range(4),dm.l12_md_mean_percentile,'o-',color=MD,lw=1.0,ms=3)
    a.plot(range(4),dm.l12_model_mean_percentile,'s--',color=BLUE,lw=1.,ms=3)
    a.set_xlim(-.3,3.3);a.set_ylim(75,101);a.set_yticks([80,90,100]);a.set_xticks(range(4),names,rotation=45,ha='right')
    a.set_ylabel('Mean percentile',labelpad=2);a.tick_params(length=2,pad=2)
    text(139,89,'Top-18 overlap:',fontsize=9)
    text(139,84,'14, 16, 11, 14 / 18',fontsize=9,color=BLUE)
    for s in states:
        r=dm.loc[s];points.append(dict(panel='d',state=s,kind='L12_mean_MD_percentile',y=float(r.l12_md_mean_percentile)))
        points.append(dict(panel='d',state=s,kind='L12_mean_model_percentile',y=float(r.l12_model_mean_percentile)))
    panel('e',2,81,'Absolute model and MD RMSF, complete KIX')
    handles=[Line2D([0],[0],color=MD,lw=1),Line2D([0],[0],color=BLUE,lw=1.1)]
    fig.legend(handles,['MD reference','Model prediction'],loc='upper right',bbox_to_anchor=(.98,.405),ncol=2,fontsize=9,handlelength=1.5,columnspacing=1.)
    ymax=float(max(data.md_value.max(),data.model_value.max()))*1.04
    for i,(s,n) in enumerate(zip(states,names)):
        left=17 if i%2==0 else 109;bottom=46 if i<2 else 11
        a=ax(left,bottom,69,23)
        p=data[data.state==s].sort_values('author_residue_number')
        a.axvspan(613.5,621.5,color=ORANGE,alpha=.15,lw=0)
        a.plot(p.author_residue_number,p.md_value,color=MD,lw=.9,marker='o',ms=1.8,mfc='white',mew=.4)
        a.plot(p.author_residue_number,p.model_value,color=BLUE,lw=1.0)
        a.set_title(n+'  ρ='+format(float(dm.loc[s,'spearman']),'.3f'),loc='left',pad=3)
        a.set_xlim(586,672);a.set_ylim(0,ymax);a.set_xticks([586,610,640,672]);a.tick_params(length=2,pad=2)
        if i%2==0:a.set_ylabel('RMSF (Å)',labelpad=2)
        if i>=2:a.set_xlabel('KIX residue number',labelpad=2)
        for _,r in p.iterrows():
            for kind,key in [('MD','md_value'),('model','model_value')]:points.append(dict(panel='e',state=s,residue=int(r.author_residue_number),kind=kind,x=int(r.author_residue_number),y=float(r[key])))
    caption='图C31：KIX共同活动区域的预测。（a）真实2AGH来源结构，cMyb和MLL结合KIX的不同表面；（b）双伙伴状态的L12（614–621），分别按MD与模型RMSF百分位着色，使用与c相同的0–100蓝色色标；（c）四种伙伴组成中KIX全部87位的MD/模型百分位，橙框为L12；（d）L12八位平均百分位及全部最高20%集合的交集数，四种状态依次为Apo、cMyb、MLL与Both；（e）全部87位绝对RMSF曲线。MD每状态一条100 ns轨迹，模型为排除KIX相关记录的第1折预测。该图检验活动区定位；四状态L12差值方向没有被模型恢复，不解释结合协同。'
    audit=dict(case='C31',model='chemical_identity_corrected_v1 full fold1; KIX A00AQ excluded from training',n_states=4,n_residues_per_state=87,n_L12=8,MD_independent_trajectories_per_state=1,MD_duration_ns=100,metric='Spearman within each complete KIX chain; exact percentages across87 positions',spread='None; no replicate uncertainty claimed',top20_counts=dm.top_overlap.astype(int).tolist(),top20_n=18,spearman=dm.spearman.to_dict(),all_source_points_retained=True,subsampling=False,source_sha256={str(OLD/'data'/name):sha(OLD/'data'/name) for name in ['kix_fold1_all_residues.csv','kix_fold1_case_metrics.csv']},boundary='Common flexible-loop localization, not cooperative state differences or binding thermodynamics')
    return save(dest,fig,axes,points,audit,caption)

def ub():
    dest,fig,axes,points,ax,text,panel,picture=scaffold('C32')
    rows=read(UBOLD/'source/all66_residue_summary.csv');defined=read(UBOLD/'experiment/defined_angles_long.csv')
    hist=pd.read_csv(UBOLD/'source/histograms_501frames.csv')
    corr=read(UBOLD/'experiment/per_chi_correlations.csv');corr=corr[(corr.predictor=='predicted_joint_Hchi')&(corr.metric=='fine_entropy')]
    assert len(rows)==66 and len(defined)==151 and len(hist)==432 and len(corr)==4
    text(2,200,'C32  Ubiquitin: side-chain diversity in one complex',fontsize=11,fontweight='bold')
    panel('a',2,190,'2QHO-derived input')
    picture(OUT/'structure_renders/ubiquitin_context.png',2,145,57,42)
    text(3,142,'Ubiquitin',color=MD);text(32,142,'EDD UBA',color=CATEGORICAL[2])
    text(6,136,'Q2',color=PURPLE);text(24,136,'Q41',color=ORANGE)
    panel('b',64,190,'All 66 valid side-chain positions')
    a=ax(77,146,100,34)
    p=rows.sort_values('author_residue_id')
    a.plot(p.author_residue_id,p.frozen_MD_joint_Hchi_nats,color=MD,lw=.8,ls=':',marker='o',ms=2.3,mfc='white',mew=.5,label='MD reference')
    a.plot(p.author_residue_id,p.predicted_Hchi_nats,color=BLUE,lw=1.,marker='s',ms=2.0,mew=0,label='Model prediction')
    for site,col in [(2,PURPLE),(41,ORANGE)]:
        a.axvspan(site-.5,site+.5,color=col,alpha=.14,lw=0)
        a.text(site,-.22,'Q'+str(site),transform=a.get_xaxis_transform(),ha='center',color=col,fontsize=9)
    missing=sorted(set(range(1,73))-set(p.author_residue_id))
    a.plot(missing,[-.03]*len(missing),'|',color='#AAAAAA',ms=4,transform=a.get_xaxis_transform(),clip_on=False)
    a.set_xlim(.5,72.5);a.set_ylim(-.1,4.3);a.set_yticks([0,2,4]);a.set_xticks([1,20,40,60,72]);a.set_ylabel('Joint Hχ (nats)',labelpad=2)
    a.tick_params(length=2,pad=2)
    a.legend(loc='lower right',bbox_to_anchor=(1,1.02),ncol=2,fontsize=9,handlelength=1.3,columnspacing=.8,borderaxespad=0)
    rho=float(spearmanr(p.frozen_MD_joint_Hchi_nats,p.predicted_Hchi_nats).statistic);mae=float(np.mean(np.abs(p.frozen_MD_joint_Hchi_nats-p.predicted_Hchi_nats)))
    text(78,135,'MAE '+format(mae,'.3f')+' nats   ρ '+format(rho,'.3f'))
    for _,r in p.iterrows():
        for name,key in [('MD','frozen_MD_joint_Hchi_nats'),('model','predicted_Hchi_nats')]:points.append(dict(panel='b',residue_id=int(r.author_residue_id),kind=name,x=int(r.author_residue_id),y=float(r[key])))
    panel('c',2,126,'Joint Hχ versus individual-angle MD occupancy')
    expected={1:66,2:52,3:23,4:10}
    for chi,left,bottom in [(1,23,96),(2,113,96),(3,23,62),(4,113,62)]:
        a=ax(left,bottom,64,23)
        q=defined[defined.chi=='chi'+str(chi)];assert len(q)==expected[chi]
        st=corr[corr.chi=='chi'+str(chi)].iloc[0]
        a.set_title('χ'+str(chi)+'   n='+str(len(q))+'   ρ='+format(float(st.rho_AA_centered),'.2f'),loc='left',pad=3)
        for label,g in [('other',q[~q.author_residue_id.isin([2,41])]),('Q2',q[q.author_residue_id==2]),('Q41',q[q.author_residue_id==41])]:
            color=PURPLE if label=='Q2' else ORANGE if label=='Q41' else BLUE
            a.scatter(g.predicted_Hchi_nats_AA_centered,g.fine_entropy_four_origins_mean_nats_AA_centered,s=12 if label=='other' else 24,color=color,edgecolors='white',linewidths=.25,alpha=.7 if label=='other' else 1,zorder=2 if label=='other' else 3)
        a.axhline(0,color=GREY,lw=.5,ls=':');a.axvline(0,color=GREY,lw=.5,ls=':')
        a.set_xlim(-.7,.5);a.set_ylim(-1,.75);a.set_xticks([-.5,0,.5]);a.set_yticks([-.5,0,.5]);a.tick_params(length=2,pad=2)
        if chi in [1,2]:a.tick_params(labelbottom=False)
        for _,r in q.iterrows():points.append(dict(panel='c',residue_id=int(r.author_residue_id),chi=chi,kind='occupancy_pair',x=float(r.predicted_Hchi_nats_AA_centered),y=float(r.fine_entropy_four_origins_mean_nats_AA_centered)))
    text(6,93,'MD entropy deviation (nats)',rotation=90,ha='center',va='center')
    text(101,53,'Predicted joint Hχ deviation (nats)',ha='center')
    panel('d',2,46,'MD free energy at prediction-selected glutamines')
    text(19,40,'Q41: lowest Hχ 2.000',color=ORANGE)
    text(103,40,'Q2: highest Hχ 2.990',color=PURPLE)
    maxf=float(hist.F_star.max());widths=[]
    for chi,left in enumerate([19,76,133],1):
        a=ax(left,9,43,25)
        a.set_title('χ'+str(chi),loc='left',pad=3)
        for j,(site,col,ls) in enumerate([('Q41',ORANGE,'-'),('Q2',PURPLE,'--')]):
            p=hist[(hist.residue==site)&(hist.chi=='chi'+str(chi))].sort_values('angle_center_degrees')
            assert len(p)==72
            a.plot(p.angle_center_degrees,p.F_star,color=col,ls=ls,lw=1.1)
            good=np.isfinite(p.F_star)&(p.F_star<=1)
            a.broken_barh([(float(v),5.) for v in p.loc[good,'angle_lower_degrees']],(-.65+j*.3,.16),facecolors=col,edgecolors='none')
            width=5*int(good.sum());widths.append(dict(residue=site,chi=chi,width_degrees=width))
            a.text(.03,.96-j*.16,site+' '+str(width)+'°',transform=a.transAxes,color=col,fontsize=9,va='top')
            for _,r in p.iterrows():points.append(dict(panel='d',residue_id=int(site[1:]),chi=chi,kind='MD_Fstar',x=float(r.angle_center_degrees),y=float(r.F_star),low_Fstar_le1=bool(np.isfinite(r.F_star) and r.F_star<=1)))
        a.axhline(1,color=GREY,ls=':',lw=.7)
        a.set_xlim(-180,180);a.set_ylim(-.8,maxf+1.8);a.set_xticks([-180,-60,60,180]);a.tick_params(length=2,pad=2)
        a.set_xlabel('Angle (°)',labelpad=2)
        if chi==1:a.set_ylabel('MD F*',labelpad=2)
        else:a.tick_params(labelleft=False)
    assert len(points)==132+151+432
    caption='图C32：泛素同类型侧链取向多样性的预测。（a）真实2QHO来源泛素–EDD UBA模型输入，显示Q2与Q41的位置；（b）泛素全部66个Hχ有效位置的MD五窗口均值参考与留出预测，线连接有效观测，轴下短灰刻度标示没有可定义Hχ的六位；（c）全部151个残基–角度观测，χ1–χ4分别66、52、23和10点，横轴仍是同一个联合Hχ模型读数，纵轴是MD单角细箱熵；两轴在各角层内按氨基酸类型去均值。Q2/41为紫/橙色；（d）由六个谷氨酰胺中模型最低Q41与最高Q2选定的MD经验投影自由能F*=−ln(p̂/p̂max)，每角72个5°箱、共432箱，采用0–100 ns的501帧；零计数箱断开，F*≤1的角区以底部色条和宽度数字表示。模型没有直接预测逐角曲线、具体构象或结合自由能。'
    audit=dict(case='C32',model='chemical_identity_corrected_v1 full OOF fold0; DynaRepo A003R held out',input_structure=str(OLD/'ubiquitin/inputs/A003R_input_structure.pdb'),input_ubiquitin_CA=72,n_all_valid_Hchi=66,n_per_chi=expected,n_all_angle_pairs=151,n_Fstar_bins=432,n_zero_count_Fstar_bins=int(hist.F_star.isna().sum()),MD_frames_for_angles=501,MD_time_window_ns=[0,100],MD_frame_spacing_ps=200,MD_Hchi_reference='Mean of five consecutive1000-frame windows in original5000 label frames',spread='No independent-repeat uncertainty claimed; one DynaRepo trajectory',mae_nats=mae,spearman=rho,lowF_widths=widths,all_source_points_retained=True,subsampling=False,source_sha256={str(p):sha(p) for p in [UBOLD/'source/all66_residue_summary.csv',UBOLD/'source/histograms_501frames.csv',UBOLD/'experiment/defined_angles_long.csv',UBOLD/'experiment/per_chi_correlations.csv']},boundary='Within-complex same-type side-chain diversity; joint scalar is not per-angle prediction or binding entropy')
    return save(dest,fig,axes,points,audit,caption)

if __name__=='__main__':
    results=[kix(),ub()]
    (OUT/'non_tcr_figure_manifest.json').write_text(json.dumps(dict(status='numerical_passed_visual_pending',figures=results),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=False))
