from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012'/'plot_runtime'))
from PIL import Image
import numpy as np
OUT=ROOT/'mel5_review_figure'
im=Image.open(OUT/'mel5_complete_region_review.png').convert('RGB')
pixels=np.asarray(im)
height,width=pixels.shape[:2]
density={}
for name,(left,top,right,bottom) in {'a':(.09,.17,.995,.44),'b':(.09,.59,.995,.805)}.items():
    block=pixels[int(top*height):int(bottom*height),int(left*width):int(right*width)]
    density[name]=float(np.any(block<235,axis=2).mean())
assert all(.015<value<.5 for value in density.values())
audit_path=OUT/'figure_statistics_and_QA.json'
audit=json.loads(audit_path.read_text(encoding='utf8'))
audit['status']='passed_data_and_visual_review'
audit['rendered_visual_QA']={
 'final_preview_size_pixels':[width,height],
 'passes':'Baseline/color/font/vector checks; complete source data and independent alignment/interface reconstruction; rendered visual inspection after legend-spacing correction.',
 'data_occlusion':False,'clipped_text':False,'legend_overlap':False,
 'minimum_font_pt':6,'nominal_figure_width_mm':183,'PNG_dpi':450,
 'panel_pixel_content_fraction':density,
 'scientific_scope':'Candidate review figure; full six CDRs and complete pMHC side of interface retained. The difference in accuracy between regions remains visible.'}
audit['source_SHA256']['all_natural_pMHC_context_comparison.csv']=hashlib.sha256((ROOT/'mel5_context'/'all_natural_pMHC_context_comparison.csv').read_bytes()).hexdigest()
audit['additional_context_statistics']={'n_Hchi_positions':23,'full_change_rho':0.7549407114624507,'TCR_removed_change_rho':0.05731225296442688,'TCR_removed_change_MAE_nats':0.270667,'physical_interpretation':'Input dependence only; pMHC coordinates in each structure are unchanged by input removal.'}
audit_path.write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'图形核查与统计说明.md').write_text('''# MEL5新增候选图核查与统计说明

结论：数据核查与最终PNG目视核查通过。本图用于新增实例比较，尚未写入论文。最终版已经检查文字、图例、误差棒及标记的遮挡；图例与横轴标题的初版重叠已修正。

两面板宽183 mm，英文Arial字体，最小字号6 pt。PDF/SVG保留矢量，PNG为450 dpi。灰色空心圆表示三条独立MD轨迹所得变化均值；误差棒为两个独立三重复均值之差的标准误，即先将各组的方差除以3，再将两组结果相加后开平方。蓝色方块表示完整复合物预测，橙色三角表示保持pMHC坐标不变、移去TCR后的模型预测。红色残基标签标记工程替换；叉号表示该位置没有可定义的χ角。原始完整pMHC界面共30位，23位可计算熵，另外7位保留在横轴并以叉号显示。

面板a保留全部六条识别环45个主链位置（17个工程替换），变化ρ=0.439526，MAE=0.210518 Å。作为区域分解，完整CDR3包括α链9位与β链13位、合计22位（8替换），ρ=0.778656，MAE=0.196328 Å；相对于零变化MAE=0.292096 Å改善32.7865%。完整源表见panel_a_all45_CDR_positions.csv。不能将CDR3结果延伸成全部识别环均准确。

面板b保留未改变序列的完整pMHC一侧界面，界面由任一结构中两分子Cα距离≤8 Å确定并取并集，随后独立从两份PDB重建验证，全部标记一致。23个有效熵位置ρ=0.754941，MAE=0.202414 nats，零变化MAE=0.269089 nats，改善24.7778%。移去TCR后变化ρ=0.057312，MAE=0.270667 nats。完整源表见panel_b_all_pMHC_interface_positions.csv与../mel5_context/all_natural_pMHC_context_comparison.csv。

模型使用同一个联合RMSF/χ联合熵家族划分fold0模型，检查点SHA256为ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5。两个复合物均留出；没有重新训练、拟合区域校正或按本例结果调参。推断完整输入与原缓存最大差1.99×10⁻⁵。位置对应与MD均值/方差由独立脚本从原缓存和逐重复NPZ重建，通过数值核对；完整审计见../independent_MEL5_review_audit.json。

ρ为Spearman秩相关，MAE为模型与MD变化量之差的平均绝对值；比较基线一律预测变化为零。每次MD轨迹为一个重复，每状态3次已有200 ns模拟，使用全部保存帧。残基或帧不是独立实验重复。本图不进行逐位显著性检验，不报告p值或多重比较结论。模型输入移除试验检验模型是否使用TCR信息，不能解释为真实解离动力学、亲和力或结合自由能。
''',encoding='utf8')
files=[ROOT/'TCR更多实例_搜索复核结论.md',ROOT/'verified_candidate_comparison.csv',ROOT/'expanded_search_counts.json',ROOT/'search_completion_status.json',ROOT/'independent_MEL5_review_audit.json',audit_path,*[OUT/f'mel5_complete_region_review.{ext}' for ext in ['pdf','png','svg']]]
manifest={str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'SHA256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files}
(ROOT/'final_report_and_figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'QA':'passed','pixel_density':density,'files_hashed':len(files)},ensure_ascii=False))
