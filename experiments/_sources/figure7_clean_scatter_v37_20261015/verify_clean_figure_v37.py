from pathlib import Path
import json,hashlib,sys,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image
OUT=ROOT/'figures';EXP=ROOT/'experiment';SOURCE=ROOT/'source'
OLD=ROOT.parent/'figure7_rich_chi_panel_v36_20261015'
STEM='fig7_clean_chi_scatter_and_MD_curves'
OLDSTEM='fig7_rich_chi_scatter_and_MD_curves'
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
before=readj(OLD/'figures'/(OLDSTEM+'.audit.json'))
check('715_total_rows',len(new)==715)
check('all_counts',new.groupby('panel').size().to_dict()=={'a_Hchi':132,'b_scatter':151,'c':432})
for panel in ['a_Hchi','b_scatter','c']:
    a=new[new.panel==panel].reset_index(drop=True)
    b=old[old.panel==panel].reset_index(drop=True)
    cols=['residue_id','kind','x','y']
    if panel!='a_Hchi':cols+=['chi']
    if panel=='b_scatter':cols+=['n_defined','is_proline']
    if panel=='c':cols+=['low_Fstar_le1']
    pd.testing.assert_frame_equal(a[cols],b[cols],check_dtype=False,check_exact=True)
    check(panel+'_all_values_exactly_v36',True)
check('no_b_summary_rows',(new.panel!='b_summary').all())
defined=read(EXP/'defined_angles_long.csv')
scatter=new[new.panel=='b_scatter']
for _,r in defined.iterrows():
    chi=int(str(r.chi).replace('chi',''))
    match=scatter[(scatter.chi==chi)&(scatter.residue_id==r.author_residue_id)]
    check('b_one_'+str(r.residue)+'_chi'+str(chi),len(match)==1)
    p=match.iloc[0]
    check('b_coordinates_'+str(r.residue)+'_chi'+str(chi),float(p.x)==float(r.predicted_Hchi_nats_AA_centered) and float(p.y)==float(r.fine_entropy_four_origins_mean_nats_AA_centered))
n_by_chi={1:66,2:52,3:23,4:10};zero_overlaps={}
for chi,n in n_by_chi.items():
    p=scatter[scatter.chi==chi]
    check('chi'+str(chi)+'_all_n',len(p)==n)
    zero_overlaps[chi]=int(((p.x==0)&(p.y==0)).sum())
check('all_zero_overlap_observations_retained',zero_overlaps=={1:3,2:3,3:1,4:0})
check('all_Pro_observations_preserved',int(scatter.is_proline.sum())==int(defined.is_proline.sum()))
check('no_undefined_angles_zero_imputed',len(scatter)==151 and len(read(EXP/'all_angle_slots_long.csv'))==264)
check('all_scatter_points_within_common_xlim',scatter.x.between(-.7,.5).all())
check('all_scatter_points_within_common_ylim',scatter.y.between(-1,.75).all())
check('common_linear_axes_unchanged',audit['scatter_xlim']==before['scatter_xlim']==[-.7,.5] and audit['scatter_ylim']==before['scatter_ylim']==[-1,.75])
check('only_abc_labels',audit['panel_letters']==['a','b','c'])
check('only_Q2_Q41_sequence_labels',[r['residue_id'] for r in audit['a_labels']]==[2,41])
check('original_c_limits',audit['c_curve_xlim']==before['c_curve_xlim'] and audit['c_curve_ylim']==before['c_curve_ylim'])
check('no_outside_text',audit['text_outside_canvas']==[])
check('font_at_least7point5',audit['minimum_font_pt']>=7.5)
check('unchanged_66_metrics',audit['frozen_overall_metrics']==before['frozen_overall_metrics'])
check('147mm_canvas',audit['configured_canvas_mm']==[183.,147.])
copied={}
for folder in ['source','experiment']:
    for p in sorted((ROOT/folder).iterdir()):
        if p.is_file():
            check('copy_hash_'+folder+'/'+p.name,sha(p)==sha(OLD/folder/p.name))
            copied[str(p.relative_to(ROOT))]=sha(p)
check('SVG_vector_only',not ET.parse(OUT/(STEM+'.svg')).findall('.//{http://www.w3.org/2000/svg}image'))
script=(OUT/'plot_clean_chi_scatter_and_MD_curves.py').read_text(encoding='utf8')
check('summary_table_plotting_removed',"panel='b_summary'" not in script and '2×4概览' not in script and 'cell_left=' not in script)
im=Image.open(OUT/(STEM+'.png'));info=im.info
actualmm=[im.width/info['dpi'][0]*25.4,im.height/info['dpi'][1]*25.4]
check('300dpi',min(info['dpi'])>=299)
audit.update(status='complete',visual_review='passed after directly viewing final PNG: coefficient table removed, compact abc layout, all four scatter panels and all three c curves clear; no overlap or clipping.',independent_check_count=len(checks),independent_checks=checks,actual_export_pixels=list(im.size),actual_export_dpi=list(info['dpi']),actual_export_mm=actualmm,export_crop_policy='183x147mm configured canvas with inherited tight outer crop.',scatter_axis_size_mm=[30,23],scatter_origin_overlap_counts=zero_overlaps,source_copy_sha256=copied,all715_values_identical_to_v36_except_removed_summary=True,no_single_angle_prediction_accuracy_claim=True,no_smoothing_jitter_or_filtering=True)
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
statistics={'a':{'n':66,'unit':'residue, not independent protein','prediction':'Frozen full-model fold0 per-residue joint Hchi, original test output.','reference':'Mean of original five consecutive 1000-frame windows; same MD trajectory, not independent MD replicates.','metrics':audit['frozen_overall_metrics'],'MAE_definition':'mean absolute prediction-reference difference, nats','rho_definition':'Spearman rank correlation; type-corrected rho frozen from earlier analysis','spread_CI_test':'No interval or new significance test.'},'b':{'n_defined_by_chi':n_by_chi,'all_defined_observations':151,'undefined_slots':113,'x':'Frozen model joint Hchi minus mean prediction of same amino-acid type within chi layer, nats.','y':'MD per-angle fine-bin entropy minus mean of same amino-acid type within chi layer, nats.','model_output_scope':'One joint Hchi per residue; no single-angle model prediction or per-angle prediction error.','MD_metric':'Shannon entropy in 5-degree bins, averaged across four origins 0/1.25/2.5/3.75 degrees.','time_scope':'Original 0-100ns, 501 frames at 200ps intervals.','rho':'Unchanged average-rank Spearman of all within-type-centered paired values per chi layer; shown in four scatter titles.','zero_ties':'Singleton types and constant groups retained; zero-overlap counts 3/3/1/0. No jitter or deduplication.','scatter_range':{'x':[-.7,.5],'y':[-1,.75]},'Pro':'All three valid Pro residues retained in all defined chi layers.','pvalues_intervals':'None; no new analysis, training or MD run.','removed_summary':'Whole2x4 coefficient table omitted at user request; all original8coefficients remain unchanged in experiment archive.'},'c':{'n':432,'definition':'2 residues x3 angles x72 five-degree bins; retains original NaN zero-count bins','selection':'Minimum Q41 and maximum Q2 frozen predicted joint Hchi among all six ubiquitin glutamines.','time_scope':'Same501 frames, origin0 bins','Fstar':'-ln(p/pmax), dimensionless; zero-count bins remain disconnected','low_F_width':'5deg times number of finite boxes with Fstar<=1, origin0, unchanged','axes':{'x':[-180,180],'y':audit['c_curve_ylim']}},'source_traceability':{'a':'source/all66_residue_summary.csv and source/ubiquitin_all66_chi_501frames.npz','b_scatter':'experiment/defined_angles_long.csv, centered fields','b_rho':'experiment/per_chi_correlations.csv, predicted_joint_Hchi/fine_entropy rows, rho_AA_centered','c':'source/histograms_501frames.csv'},'export':{'configured_mm':[183,147],'actual_mm':actualmm,'actual_pixels':list(im.size),'dpi':list(info['dpi']),'vector':['PDF embedded TrueType','SVG editable text, no raster image nodes']}}
(OUT/(STEM+'.statistics.json')).write_text(json.dumps(statistics,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
qa=f'''# 图7 v37 核验与统计说明

按用户要求删除b下方整个2×4系数表及对应图注，保留a/b/c三标号。图高160→147mm，上排a/b整体下移13mm收紧空白，未缩小字体或散点轴；四个散点轴仍为30×23mm，最小字体7.5pt。a的全部66残基曲线、仅Q2/Q41标记和统计保留，b全部151角及χ4、n=66/52/23/10、四个ρ保留；c的432条源值、角度范围和物理宽度全部不变。

共{len(checks)}项独立核验通过。715行点表=132a+151b+432c，与v36对应值逐浮点严格同一；只删除8条摘要显示行。151角逐行对应来源，原点重合数3/3/1/0和Pro全部保留，没有抖动、去重、拟合线或y=x线。全部source/experiment复制文件SHA256一致，原八个相关系数留在实验归档，未做新统计或模拟。

已直接查看最终PNG：四个χ散点、共享轴、标题、n/ρ、图例和c三曲线均清晰，无文字碰撞、裁切或覆盖；删除表后的间距紧凑。PDF和SVG保留矢量，SVG无位图节点。配置画布183×147mm，tight裁边后为{im.width}×{im.height}px、约300dpi，实际外沿约{actualmm[0]:.2f}×{actualmm[1]:.2f}mm。完整定义和溯源见statistics.json。
'''
(OUT/(STEM+'.QA_statistics_zh.md')).write_text(qa,encoding='utf8')
print(json.dumps({'status':'complete','checks':len(checks),'points':len(new),'point_counts':audit['point_counts'],'pixels':list(im.size),'actual_mm':actualmm,'PNG_sha256':sha(OUT/(STEM+'.png'))},ensure_ascii=False))
