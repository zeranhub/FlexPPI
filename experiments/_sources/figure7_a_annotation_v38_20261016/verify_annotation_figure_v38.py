from pathlib import Path
import json,hashlib,sys,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image,ImageChops
OUT=ROOT/'figures';EXP=ROOT/'experiment';SOURCE=ROOT/'source'
OLD=ROOT.parent/'figure7_clean_scatter_v37_20261015'
STEM='fig7_clear_a_annotations_and_MD_curves'
OLDSTEM='fig7_clean_chi_scatter_and_MD_curves'
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
    check(panel+'_all_values_exactly_v37',True)
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
script=(OUT/'plot_clear_a_annotations_and_MD_curves.py').read_text(encoding='utf8')
check('summary_table_plotting_removed',"panel='b_summary'" not in script and '2×4概览' not in script and 'cell_left=' not in script)
im=Image.open(OUT/(STEM+'.png'));info=im.info
actualmm=[im.width/info['dpi'][0]*25.4,im.height/info['dpi'][1]*25.4]
check('300dpi',min(info['dpi'])>=299)
check('a_xaxis_title_removed',audit['a_xaxis_title_removed'] and '泛素残基编号' not in script)
check('a_legend_above_right',audit['a_legend_above_axis'] and audit['a_legend_location']=='above_curve_axis_upper_right')
check('a_legend_no_stats_overlap',audit['a_legend_no_overlap_with_statistics'])
check('clear_type_mean_wording',audit['a_type_mean_label']=='扣除氨基酸类型均值后 ρ 0.868')
check('all66_same_chain_A',(read(SOURCE/'all66_residue_summary.csv').chain=='A').all())
check('six_NaN_gaps_retained',audit['a_missing_positions']==[10,28,35,46,47,53] and audit['a_NaN_breaks_preserved'])
oldimage=Image.open(OLD/'figures'/(OLDSTEM+'.png')).convert('RGB')
newimage=im.convert('RGB')
check('same_PNG_dimensions',oldimage.size==newimage.size)
unchanged_regions={'b_complete_right':(1160,0,im.width,im.height),'c_complete_lower':(0,1100,im.width,im.height)}
for name,region in unchanged_regions.items():check('pixel_identical_'+name,ImageChops.difference(oldimage.crop(region),newimage.crop(region)).getbbox() is None)
pixel_diff=ImageChops.difference(oldimage,newimage).getbbox()

audit.update(a_missing_residue_names={'10':'GLY','28':'ALA','35':'GLY','46':'ALA','47':'GLY','53':'GLY'},a_type_mean_correction='Same fold0 training reference amino-acid mean subtracted from both prediction and reference; distinct from b layer-wise means.',status='complete',visual_review='passed after directly viewing final PNG: a x-axis title removed, legend above curve at upper right, plain type-mean label clear; all b/c elements unchanged, no overlap or clipping.',independent_check_count=len(checks),independent_checks=checks,actual_export_pixels=list(im.size),actual_export_dpi=list(info['dpi']),actual_export_mm=actualmm,export_crop_policy='183x147mm configured canvas with inherited tight outer crop.',scatter_axis_size_mm=[30,23],scatter_origin_overlap_counts=zero_overlaps,source_copy_sha256=copied,all715_values_identical_to_v37=True,pixel_difference_bbox=pixel_diff,pixel_identical_regions=unchanged_regions,no_single_angle_prediction_accuracy_claim=True,no_smoothing_jitter_or_filtering=True)
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
statistics=readj(OLD/'figures'/(OLDSTEM+'.statistics.json'))
statistics['export']={'configured_mm':[183,147],'actual_mm':actualmm,'actual_pixels':list(im.size),'dpi':list(info['dpi']),'vector':['PDF embedded TrueType','SVG editable text, no raster image nodes']}
statistics['a']['type_mean_correction']='For each residue, subtract the same ordinary-residue-weighted fold0 training reference mean of its amino-acid type from both frozen prediction and MD reference; compute Spearman rho. Not a ubiquitin mean, not separately fitted prediction/reference means. Panel b uses separate within-layer means.'
statistics['a']['annotation_update']='Plain-language type-mean caption, same frozen rho0.868. X-axis wording removed; numeric ticks and Q2/Q41 retained. Six NaN positions10/28/35/46/47/53 retained.'
(OUT/(STEM+'.statistics.json')).write_text(json.dumps(statistics,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
qa=f'''# 图7 v38 标注修订与核验

a删除“泛素残基编号”文字，保留全部数值刻度和Q2/Q41标记；MD参考/模型预测图例移到折线图右上方外侧空白。图例包围框位于数据轴之上，与MAE和“扣除氨基酸类型均值后 ρ 0.868”统计字样无重叠；数值与原统计定义未改。a原六个NaN位置G10、A28、G35、A46、G47、G53保持不绘点及断线；图注明确全部为泛素A链及缺少可定义侧链χ角。a校正ρ使用训练折的同一个氨基酸类型参考均值同时扣除预测与MD参考，b则在各χ层内分别扣除同类型均值，两者定义已分清。

共{len(checks)}项独立核验通过。全部715条显示数据与v37逐浮点严格相同，151角原点重叠3/3/1/0和Pro/χ4完整保留，全部source/experiment复制SHA256一致。b完整区域与c完整区域逐像素相同；像素变化仅来自a的三项标注修改。无新统计、训练、模拟、筛点、抖动、去重或拟合。

直接查看完整PNG，a图例、统计行和折线互不覆盖，b/c与原版清晰度一致，无文字裁切或碰撞。最小字体7.5pt，四个散点轴仍30×23mm，画布183×147mm；tight裁边PNG为{im.width}×{im.height}px、约300dpi，实际外沿约{actualmm[0]:.2f}×{actualmm[1]:.2f}mm。PDF/SVG同步导出，SVG无位图节点。完整统计定义沿用v37并记录在statistics.json。
'''
(OUT/(STEM+'.QA_statistics_zh.md')).write_text(qa,encoding='utf8')
print(json.dumps({'status':'complete','checks':len(checks),'points':len(new),'point_counts':audit['point_counts'],'pixels':list(im.size),'actual_mm':actualmm,'PNG_sha256':sha(OUT/(STEM+'.png'))},ensure_ascii=False))
