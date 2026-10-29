from pathlib import Path
import json,hashlib,sys,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image
OUT=ROOT/'figures';EXP=ROOT/'experiment';SOURCE=ROOT/'source'
OLD=ROOT.parent/'figure7_all_chi_v35_20261015'
STEM='fig7_rich_chi_scatter_and_MD_curves'
OLDSTEM='fig7_all_chi_associations_and_MD_curves'
read=lambda p:pd.read_csv(p,float_precision='round_trip')
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks=[]
def check(name,value):
    assert bool(value),name
    checks.append({'check':name,'pass':True})
new=read(OUT/(STEM+'.points.csv'))
old=read(OLD/'figures'/(OLDSTEM+'.points.csv'))
audit=readj(OUT/(STEM+'.audit.json'))
check('723_total_rows',len(new)==723)
check('all_counts',new.groupby('panel').size().to_dict()=={'a_Hchi':132,'b_scatter':151,'b_summary':8,'c':432})
a=new[new.panel=='a_Hchi'].reset_index(drop=True)
b=old[old.panel=='a_Hchi'].reset_index(drop=True)
pd.testing.assert_frame_equal(a[['residue_id','kind','x','y']],b[['residue_id','kind','x','y']],check_dtype=False,check_exact=True)
check('a_all132_values_exactly_v35',True)
before=old[old.kind=='MD_Fstar'].sort_values(['chi','residue_id','x']).reset_index(drop=True)
after=new[new.panel=='c'].sort_values(['chi','residue_id','x']).reset_index(drop=True)
pd.testing.assert_frame_equal(after[['residue_id','chi','x','y','low_Fstar_le1']],before[['residue_id','chi','x','y','low_Fstar_le1']],check_dtype=False,check_exact=True)
check('c_all432_values_exactly_v35',True)
defined=read(EXP/'defined_angles_long.csv')
scatter=new[new.panel=='b_scatter']
for _,r in defined.iterrows():
    chi=int(str(r.chi).replace('chi',''))
    match=scatter[(scatter.chi==chi)&(scatter.residue_id==r.author_residue_id)]
    check('b_one_'+str(r.residue)+'_chi'+str(chi),len(match)==1)
    p=match.iloc[0]
    check('b_coordinates_'+str(r.residue)+'_chi'+str(chi),float(p.x)==float(r.predicted_Hchi_nats_AA_centered) and float(p.y)==float(r.fine_entropy_four_origins_mean_nats_AA_centered))
check('all_Pro_observations_preserved',int(scatter.is_proline.sum())==int(defined.is_proline.sum()))
check('no_undefined_angles_zero_imputed',len(scatter)==151 and len(read(EXP/'all_angle_slots_long.csv'))==264)
n_by_chi={1:66,2:52,3:23,4:10}
zero_overlaps={}
corr=read(EXP/'per_chi_correlations.csv')
corr=corr[(corr.predictor=='predicted_joint_Hchi')&corr.metric.isin(['fine_entropy','low_Fstar_width'])]
for chi,n in n_by_chi.items():
    p=scatter[scatter.chi==chi]
    d=defined[defined.chi=='chi'+str(chi)]
    check('chi'+str(chi)+'_all_n',len(p)==len(d)==n)
    zero_overlaps[chi]=int(((p.x==0)&(p.y==0)).sum())
    for metric,field in [('fine_entropy','fine_entropy_four_origins_mean_nats_AA_centered'),('low_Fstar_width','width_four_origins_mean_degrees_AA_centered')]:
        rho=float(np.corrcoef(d.predicted_Hchi_nats_AA_centered.rank(method='average'),d[field].rank(method='average'))[0,1])
        source=corr[(corr.chi=='chi'+str(chi))&(corr.metric==metric)].iloc[0]
        shown=new[(new.panel=='b_summary')&(new.chi==chi)&(new.kind==metric)].iloc[0]
        check('summary_original_'+str(chi)+'_'+metric,float(shown.rho)==float(source.rho_AA_centered))
        check('summary_independent_rank_'+str(chi)+'_'+metric,abs(float(shown.rho)-rho)<1e-14)
check('all_zero_overlap_observations_retained',zero_overlaps=={1:3,2:3,3:1,4:0})
check('chi4_negative_width_kept',float(new[(new.panel=='b_summary')&(new.chi==4)&(new.kind=='low_Fstar_width')].rho.iloc[0])<0)
check('all_scatter_points_within_common_xlim',scatter.x.between(-.7,.5).all())
check('all_scatter_points_within_common_ylim',scatter.y.between(-1,.75).all())
check('common_linear_axes',audit['scatter_xlim']==[-.7,.5] and audit['scatter_ylim']==[-1,.75])
check('full_fixed_correlation_colorscale',audit['correlation_color_scale']==[-1,1])
check('only_abc_labels',audit['panel_letters']==['a','b','c'])
check('only_Q2_Q41_sequence_labels',[r['residue_id'] for r in audit['a_labels']]==[2,41])
check('original_c_limits',audit['c_curve_xlim']==[-180,180] and audit['c_curve_ylim']==readj(OLD/'figures'/(OLDSTEM+'.audit.json'))['c_curve_ylim'] if 'c_curve_ylim' in readj(OLD/'figures'/(OLDSTEM+'.audit.json')) else audit['c_curve_xlim']==[-180,180] and audit['c_curve_ylim']==[-.8,6.628313737302301])
check('no_outside_text',audit['text_outside_canvas']==[])
check('font_at_least7point5',audit['minimum_font_pt']>=7.5)
check('unchanged_66_metrics',audit['frozen_overall_metrics']==readj(OLD/'source/fig7_all_residues_Gln_and_MD_free_energy.statistics.json'))
copied={}
for folder in ['source','experiment']:
    for p in sorted((ROOT/folder).iterdir()):
        if p.is_file():
            check('copy_hash_'+folder+'/'+p.name,sha(p)==sha(OLD/folder/p.name))
            copied[str(p.relative_to(ROOT))]=sha(p)
check('SVG_vector_only',not ET.parse(OUT/(STEM+'.svg')).findall('.//{http://www.w3.org/2000/svg}image'))
im=Image.open(OUT/(STEM+'.png'))
info=im.info;actualmm=[im.width/info['dpi'][0]*25.4,im.height/info['dpi'][1]*25.4]
check('300dpi',min(info['dpi'])>=299)
audit.update(status='complete',visual_review='passed after direct full-PNG inspection of two rendered iterations: abc numbering, all four scatter panels, eight-value overview, legend, Q2/Q41 colors, and all three c curves clear; no overlap or clipping.',independent_check_count=len(checks),independent_checks=checks,actual_export_pixels=list(im.size),actual_export_dpi=list(info['dpi']),actual_export_mm=actualmm,export_crop_policy='183x160mm configured canvas with inherited tight outer crop.',scatter_axis_size_mm=[30,23],scatter_origin_overlap_counts=zero_overlaps,source_copy_sha256=copied,all723_original_or_derived_values_verified=True,no_single_angle_prediction_accuracy_claim=True,no_smoothing_jitter_or_filtering=True)
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
statistics={'a':{'n':66,'unit':'residue, not independent protein','prediction':'Frozen full-model fold0 per-residue joint Hchi, original test output.','reference':'Mean of original five consecutive 1000-frame windows; same MD trajectory, not independent MD replicates.','metrics':audit['frozen_overall_metrics'],'MAE_definition':'mean absolute prediction-reference difference, nats','rho_definition':'Spearman rank correlation; type-corrected rho frozen from earlier analysis','spread_CI_test':'No interval or new significance test.'},'b':{'n_defined_by_chi':n_by_chi,'all_defined_observations':151,'undefined_slots':113,'x':'Frozen model joint Hchi minus mean prediction of same amino-acid type within chi layer, nats.','y':'MD per-angle fine-bin entropy minus mean of same amino-acid type within chi layer, nats.','model_output_scope':'One joint Hchi per residue; no single-angle model prediction or per-angle prediction error.','MD_metrics':'Shannon entropy in 5-degree bins and total width of Fstar<=1 boxes; both averaged across four bin origins 0/1.25/2.5/3.75 degrees.','time_scope':'Original 0-100ns, 501 frames at 200ps intervals.','rho':'Average-rank Spearman of all within-type-centered paired values per chi layer; all eight entropy/width associations retained.','zero_ties':'Singleton types and constant groups retained; zero-overlap counts 3/3/1/0. No jitter or deduplication.','scatter_range':{'x':[-.7,.5],'y':[-1,.75]},'summary_color_range':[-1,1],'Pro':'All three valid Pro residues retained in all defined chi layers.','pvalues_intervals':'None. No new training or MD run; observations within one protein.'},'c':{'n':432,'definition':'2 residues x 3 angles x 72 five-degree bins; retains original NaN zero-count bins','selection':'Minimum Q41 and maximum Q2 frozen predicted joint Hchi among all six ubiquitin glutamines.','time_scope':'Same 501 frames, origin0 bins','Fstar':'-ln(p/pmax), dimensionless; zero-count bins remain disconnected','low_F_width':'5deg times number of finite boxes with Fstar<=1, origin0','axes':{'x':[-180,180],'y':audit['c_curve_ylim']}},'source_traceability':{'a':'source/all66_residue_summary.csv and source/ubiquitin_all66_chi_501frames.npz','b_scatter':'experiment/defined_angles_long.csv, centered fields','b_summary':'experiment/per_chi_correlations.csv, predicted_joint_Hchi rows, rho_AA_centered','c':'source/histograms_501frames.csv'},'export':{'configured_mm':[183,160],'actual_mm':actualmm,'actual_pixels':list(im.size),'dpi':list(info['dpi']),'vector':['PDF embedded TrueType','SVG editable text, no raster image nodes']}}
(OUT/(STEM+'.statistics.json')).write_text(json.dumps(statistics,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
qa=f'''# 图7 v36 核验与统计说明

最终图采用a/b/c三个标号：a为全部66残基的132个参考/预测点；b为χ1–χ4全部151个有效角的散点及八项相关性概览，四个散点轴均为30×23mm、共享线性范围；c为三条逐角投影的小轴组，仅一个c标号。b包含n=66/52/23/10，八个ρ全部保留，χ4宽度负值−0.183原样显示；Q2紫、Q41橙跨面板对应。散点轴名明确同类型偏差，不称作单角预测误差。

独立核验共{len(checks)}项通过。723行点表=132+151+8+432；a的132原值、c的432原值均与v35严格相同，151角逐行对应保存来源；八个ρ从保存的同类型中心化值独立以平均秩重算一致。原点重叠数χ1/χ2/χ3/χ4=3/3/1/0，全部保留，没有抖动、筛点、去重或平滑；Pro与χ4全部保留。复制的全部source/experiment文件SHA256一致，未修改v35。

已直接查看两次完整PNG，最终版本中a、b、c标题、参考/预测图例、散点标题、共享轴标签、相关性数字、c三轴角标和端点说明均清晰，无碰撞、裁切或覆盖。最小字体7.5pt；无y=x线、回归线或新p值。PDF与SVG维持矢量；SVG无位图节点。

配置画布183×160mm，沿用tight外沿裁边。PNG为{im.width}×{im.height}，约300dpi，实际外沿约{actualmm[0]:.2f}×{actualmm[1]:.2f}mm。统计定义、样本单位、模型折号、MD时段、四箱原点和逐文件来源见statistics.json；source/experiment完整保留151角和113缺失槽的原始处理说明。没有新MD生产、训练或显著性检验。
'''
(OUT/(STEM+'.QA_statistics_zh.md')).write_text(qa,encoding='utf8')
print(json.dumps({'status':'complete','checks':len(checks),'points':len(new),'point_counts':audit['point_counts'],'pixels':list(im.size),'actual_mm':actualmm,'zero_overlaps':zero_overlaps},ensure_ascii=False))
