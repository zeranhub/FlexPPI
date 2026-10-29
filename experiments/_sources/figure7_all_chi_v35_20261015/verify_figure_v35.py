from pathlib import Path
import sys, json, hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image,ImageChops
OUT=ROOT/'figures';EXP=ROOT/'experiment'
STEM='fig7_all_chi_associations_and_MD_curves'
OLD=ROOT.parent/'figure7_sequence_order_v34_20261015'
OLDSTEM='fig7_sequence_order_predictions_and_MD_curves'
read=lambda p:pd.read_csv(p,float_precision='round_trip')
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks=[]
def check(name,value):
    assert bool(value),name
    checks.append({'check':name,'pass':True})
points=read(OUT/(STEM+'.points.csv'));oldpoints=read(OLD/'figures'/(OLDSTEM+'.points.csv'))
check('572_plot_rows',len(points)==572)
check('all_panel_counts',points.groupby('panel').size().to_dict()=={'a_Hchi':132,'a_chi_association':8,'b':144,'c':144,'d':144})
for panel in ['a_Hchi','b','c','d']:
    a=points[points.panel==panel].reset_index(drop=True)
    b=oldpoints[oldpoints.panel==panel].reset_index(drop=True)
    cols=['panel','residue_id','kind','x','y']
    if panel!='a_Hchi':cols+=['chi','low_Fstar_le1']
    pd.testing.assert_frame_equal(a[cols],b[cols],check_dtype=False,check_exact=True)
    check('unchanged_'+panel+'_source_values',True)
source=read(EXP/'per_chi_correlations.csv')
chosen=source[(source.predictor=='predicted_joint_Hchi')&source.metric.isin(['fine_entropy','low_Fstar_width'])]
check('exactly8_source_associations',len(chosen)==8)
defined=read(EXP/'defined_angles_long.csv')
slots=read(EXP/'all_angle_slots_long.csv')
check('all151_defined_angles',len(defined)==151)
check('all264_angle_slots_and113_missing',len(slots)==264 and int(slots.chi_defined.sum())==151)
expected_n={1:66,2:52,3:23,4:10}
for chi,n in expected_n.items():
    label='chi'+str(chi)
    part=defined[defined.chi==label]
    check(label+'_n',len(part)==n)
    for metric,field in [('fine_entropy','fine_entropy_four_origins_mean_nats'),('low_Fstar_width','width_four_origins_mean_degrees')]:
        xx=part.predicted_Hchi_nats_AA_centered.to_numpy()
        yy=part[field+'_AA_centered'].to_numpy()
        # Independent average-rank correlation of all saved residuals, including ties.
        rho=float(np.corrcoef(pd.Series(xx).rank(method='average'),pd.Series(yy).rank(method='average'))[0,1])
        r=chosen[(chosen.chi==label)&(chosen.metric==metric)].iloc[0]
        p=points[(points.panel=='a_chi_association')&(points.chi==chi)&(points.kind==metric)].iloc[0]
        check(label+'_'+metric+'_source_rho',abs(rho-float(r.rho_AA_centered))<1e-14)
        check(label+'_'+metric+'_plot_rho',float(p.x)==float(r.rho_AA_centered) and int(p.n_defined)==n)
        calculated=part[field]-part.groupby('amino_acid')[field].transform('mean')
        check(label+'_'+metric+'_AA_residuals',np.allclose(calculated,yy,rtol=0,atol=1e-14))
    calculated=part.predicted_Hchi_nats-part.groupby('amino_acid').predicted_Hchi_nats.transform('mean')
    check(label+'_joint_prediction_AA_residuals',np.allclose(calculated,part.predicted_Hchi_nats_AA_centered,rtol=0,atol=1e-14))
check('chi4_negative_width_retained',float(points[(points.panel=='a_chi_association')&(points.chi==4)&(points.kind=='low_Fstar_width')].x.iloc[0])<0)
audit=readj(OUT/(STEM+'.audit.json'))
check('association_full_axis',audit['association_axis_range']==[-1,1])
check('only_Q2_Q41_labelled',[(r['residue_id'],r['label_x_data']) for r in audit['a_left_Gln_labels']]==[(2,2.),(41,41.)])
check('unchanged_66_main_metrics',audit['frozen_overall_metrics']==readj(OLD/'source/fig7_all_residues_Gln_and_MD_free_energy.statistics.json'))
check('font_floor',audit['minimum_font_pt']>=7.5)
check('no_text_outside',audit['text_outside_canvas']==[])
im=Image.open(OUT/(STEM+'.png')).convert('RGB')
oldim=Image.open(OLD/'figures'/(OLDSTEM+'.png')).convert('RGB')
regions={'a_left_joint_curve_and_numeric_axis':(0,270,1380,640),'b_to_d_all_curves_labels_and_axes':(0,1000,min(im.width,oldim.width),im.height)}
for name,box in regions.items():
    check('pixel_identical_'+name,ImageChops.difference(im.crop(box),oldim.crop(box)).getbbox() is None)
source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'source').iterdir()) if p.is_file()}
for filename in source_hashes:
    original=OLD/filename
    if not original.exists():original=ROOT.parent/'rich_angle_figure_v30_20261014/experiment'/Path(filename).name
    if original.exists():check('source_unchanged_'+filename,sha(ROOT/filename)==sha(original))
info=Image.open(OUT/(STEM+'.png')).info
actualmm=[im.width/info['dpi'][0]*25.4,im.height/info['dpi'][1]*25.4]
audit.update(status='complete',visual_review='passed: full PNG directly viewed; all four chi layers, both metric symbols and n labels clear; chi4 negative width retained; no overlaps or clipping.',independent_check_count=len(checks),independent_checks=checks,pixel_identical_regions=regions,actual_export_png_pixels=list(im.size),actual_export_png_dpi=list(info['dpi']),actual_export_mm=actualmm,export_crop_policy='183x135mm configured canvas; inherited tight crop from v34. Changed right axis labels cause a small change of cropped outer width; original axes positions stay fixed.',source_local_copy_sha256=source_hashes,experiment_source_sha256={p.name:sha(p) for p in [EXP/'per_chi_correlations.csv',EXP/'figure_statistics.csv',EXP/'defined_angles_long.csv',EXP/'all_angle_slots_long.csv']})
(OUT/(STEM+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
statistics={'panel_a_left':{'n':66,'n_definition':'All valid ubiquitin residues, including three valid zero-reference Pro positions.','MAE':'mean absolute error of frozen joint-Hchi prediction vs original five-window mean reference, nats.','rho':'Spearman rank correlation.','same_as_v34':audit['frozen_overall_metrics']},'panel_a_right':{'units':'dimensionless Spearman rho','n_by_chi':expected_n,'source':'experiment/per_chi_correlations.csv; predictor=predicted_joint_Hchi; metric=fine_entropy or low_Fstar_width; field=rho_AA_centered','model_output':'One joint-Hchi value per residue from frozen full-model test fold 0. No single-angle model entropy or width predictions exist.','MD_entropy':'Natural-log Shannon entropy in 5-degree angle bins, mean of four origins 0,1.25,2.5,3.75 degrees, per chi.','MD_width':'Sum of all 5-degree bins with Fstar<=1, mean of the same four origins; degrees.','frame_scope':'Existing 0-100ns, 501 frames at 200ps intervals.','AA_control':'Within each chi layer, separately subtract same-amino-acid mean from model joint Hchi and MD metric; compute Spearman of all paired residuals. Singleton amino-acid types and constant groups remain zero residuals/ties; average ranks retained.','missing':'113 undefined slots retained in source table; not zero-imputed into correlations.','all_raw_rows':'151 defined residue-angle observations, 66/52/23/10 per layer.','intervals_tests':'No intervals, p-values, or new statistical tests; no independent MD repeats implied.','model_selection':'No new training, seed choice or site selection.'},'panels_b_d':{'same_as_v34':True,'rows':432,'n':'two residues x three angles x 72 five-degree bins','selection':'Q41/Q2 frozen-model minimum/maximum among all six GLN.','definition':'MD Fstar=-ln(p/pmax), origin 0 degrees; zero-count bins remain NaN and disconnected.','time_scope':'same501 frames as a-right, but b-d displays origin0 rather than four-origin average.'},'plot_rows':572,'independent_checks':len(checks),'source_traceability':audit['experiment_source_sha256'],'export':{'configured_mm':[183,135],'actual_mm':actualmm,'dpi':list(info['dpi'])}}
(OUT/(STEM+'.statistics.json')).write_text(json.dumps(statistics,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
report=f'''# 图7 v35 视觉与数据核验

已直接查看最终PNG。a左保留66个位点的全部预测/参考主曲线，只标Q2和Q41；Q41标签回实际41位置。a右全部显示χ1–χ4及n=66/52/23/10，两种既有MD指标均保留，χ4宽度负相关原样显示；采用完整−1到1相关系数横轴，蓝色实圆/灰色空菱形与上下轻微错位区分两指标，未加误差条。各文字、图例与点不重叠，无裁切，最小字体7.5 pt。

最终图中572行=左侧132行+右侧8行相关摘要+b–d432行。其余564条原数据逐值等同v34；a左主曲线区域与b–d区域逐像素相同。原151个有效逐角观察与113个未定义角槽均在experiment保存，独立重算逐角平均秩相关与图上八个系数一致。本图共有{len(checks)}项独立核验通过；MD逐角原点指标的独立原始轨迹复算另见experiment/audit.json。

右侧图是联合Hχ与单角MD指标的关联，模型没有单角预测值。每个χ层内分别减去相同氨基酸类型的均值，然后对全部残差计算Spearmanρ；单例和常数组以零残差/并列秩保留。细箱熵和宽度均由0–100ns的501帧、5°箱、四原点平均产生；b–d仍显示原点0°曲线。未进行训练、生产MD或新显著性检验。完整统计定义和溯源见statistics.json。

声明画布183×135mm，沿用原版tight外沿裁边，最终PNG为{im.width}×{im.height}、约300dpi，实际外沿约{actualmm[0]:.2f}×{actualmm[1]:.2f}mm。右侧轴标签使裁边宽度微变，原有轴位置、曲线及字体未变。同步输出PDF与可编辑文字SVG。
'''
(OUT/(STEM+'.QA_statistics_zh.md')).write_text(report,encoding='utf8')
plan=readj(ROOT/'figure_layout_plan.json');plan['status']='complete';plan['right_summary']['point_source']='experiment/per_chi_correlations.csv';(ROOT/'figure_layout_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':'complete','checks':len(checks),'points':len(points),'actual_png':list(im.size),'unchanged_564_original_rows':True,'all8_coefficients_verified':True},ensure_ascii=False))
