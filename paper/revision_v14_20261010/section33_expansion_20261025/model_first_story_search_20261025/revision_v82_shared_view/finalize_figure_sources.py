from pathlib import Path
import json,hashlib
r=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');v=r/'revision_v82_shared_view';f=v/'figure3'
a=f/'structure_assets/render_manifest.json';d=json.loads(a.read_text(encoding='utf-8'));d.pop('global_N_lobe_alignment_RMSD_A',None);new=json.loads((v/'pymol/render_manifest.json').read_text(encoding='utf-8'));d['global_display_alignment']={k:val for k,val in new.items() if k not in ['assets','palette']};a.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
qa=json.loads((f/'figure_QA.json').read_text(encoding='utf-8'));qa['visual_review']='passed: root and independent numerical/visual review';qa['source_data_identical_to_v81']={}
for name in ['cam_md_state_contrast_all_common_residues.csv','cam_md_state_contrast_all_residues_five_blocks.csv']:
    values=(f/'source_data'/name).read_bytes();assert values==(r/'revision_v81_absolute_profiles/figure3/source_data'/name).read_bytes();qa['source_data_identical_to_v81'][name]=hashlib.sha256(values).hexdigest()
qa['display_geometry_audit']=new.get('geometry_audit',new.get('coordinate_audit',{}));qa['display_fit_RMSD_A']=new['fit_RMSD_A'];qa['numeric_correspondence_to_v81']='All 1364 plotted quantities identical; only panel letters and layout changed.'
(f/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
readme='''# 图3 v82 复现与统计说明

本版按用户要求，整体结构a使用共同的N叶螺旋I/IV参考系、相同正交相机、比例和共同裁边；原e局部放大紧邻a右侧并改为b。图内没有数值表及灰色说明脚注。方法、限制和关键绝对数值写入图注与3.3正文。a中的CaM序列相同，但构象确实不同，不能解释为只有伙伴肽链发生变化。

结构来自原模型输入的真实原子坐标。显示对齐只对整个IQCG复合物执行一个刚体变换，不分别扭转两叶瓣。24个Cα（成熟编号7–19、65–75）拟合RMSD为1.313123 Å。显示视角选择和刚体拟合不改变MD读数或模型预测。b沿用局部位点前后五残基N/CA/C拟合后的成对视角，完整连续肽链保留在场景中，镜头只局部取景。新的全结构PSE、输入PDB、自包含渲染脚本与原子几何审计位于../pymol。

c：145个主链位置与122个Hχ有效侧链位置，分别在MD参考及固定第0折Full预测中展示两态原值，暖赭=NHE1、蓝=IQCG；相同指标MD与Full共用刻度。无有效Hχ的Ala/Gly/Pro等残基留空断线。浅带为EF手环区，接触按跨链Cα≤8 Å定义。

d：完整8/21/3/90侧链接触分组的平均Hχ。e：全部九个Met。MD细线为同一条500 ns轨迹的五个连续100 ns窗口最小值到最大值，不是独立重复或置信区间。f：每态全部122侧链，共244个绝对Hχ对应点；不能用该绝对对应取代状态变化的预测检验，变化MAE/排序相关继续保留正文。

两条CaM记录均在该折训练中，本图为结构解释，不替代留出泛化评价。伙伴、构型及原始钙条件同时不同，不能作单因素归因。RMSF按完整CaM拟合，包含叶间运动；叶瓣独立拟合对照见源表与正文。

运行：python compose_fig3.py。依赖numpy、pandas、matplotlib、Pillow；source_data和structure_assets包含绘图所需数值与结构渲染。全1364个已绘数值与v81逐字段完全相同，CSV位于source_data/all_plotted_absolute_quantities.csv。输出为183×235 mm，450 dpi PNG、嵌入字体PDF及可编辑文字SVG。
'''
(f/'README_统计与复现.md').write_text(readme,encoding='utf-8')
print(json.dumps({'figure_QA':'passed','fit_RMSD_A':new['fit_RMSD_A'],'old_misleading_field_removed':True},ensure_ascii=False))