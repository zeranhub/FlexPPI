from pathlib import Path
import json, shutil
root=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025')
v=root/'revision_v81_absolute_profiles'; old=root/'revision_v80_redesign'; f=v/'figure3'; m=v/'manuscript'
assets=f/'structure_assets'; assets.mkdir(exist_ok=True)
manifest=json.loads((old/'pymol/render_manifest.json').read_text(encoding='utf-8'))
keys=['whole_NHE1','whole_IQCG','detail_focus_M72_NHE1','detail_focus_M72_IQCG','detail_focus_M109_NHE1','detail_focus_M109_IQCG']
for key in keys:
    item=manifest['assets'][key]; name=item.get('paired_trim',item)['file']; shutil.copy2(old/'pymol'/name,assets/name)
manifest['assets']={key:manifest['assets'][key] for key in keys}
(assets/'render_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
body=(f/'absolute_body.py').read_text(encoding='utf-8').replace("R=OUT.parent.parent/'revision_v80_redesign/pymol'", "R=OUT/'structure_assets'")
(f/'absolute_body.py').write_text(body,encoding='utf-8')
compose=(f/'compose_fig3.py').read_text(encoding='utf-8'); header=compose[:compose.index('# --- panel code ---')]
(f/'compose_fig3.py').write_text(header+'# --- panel code ---\n'+body,encoding='utf-8')
caption='图3：CaM两种识别构型的主链与侧链柔性。（a）NHE1与IQCG态的输入结构及各态全部122个侧链Hχ的MD参考与预测对应；结构深蓝、浅蓝分别表示CaM的N、C叶瓣，暖赭表示伙伴，青绿标出M72与M109。（b）145个共同主链位置的RMSF及122个有效侧链位置的Hχ，以暖赭、蓝色曲线分别表示NHE1、IQCG态；MD与Full预测上下对应，同一指标共用纵轴范围。浅色带表示EF手环区，下方标记输入结构中跨链Cα距离≤8 Å的接触位置；无有效Hχ的位置留空。（c）按两态输入接触划分的四个完整侧链组，两种颜色的点表示两态平均Hχ。（d）全部九个Met的两态Hχ。（e）M72与M109的结构环境及两态绝对读数。RMSF单位为Å，Hχ为nats；c、d中MD细线为同一轨迹五个连续100 ns时间段的最小至最大值，并非独立重复的置信区间。两条记录均用于固定第0折模型训练，本图用于结构解释。'
(f/'caption_zh.md').write_text(caption+'\n',encoding='utf-8')
plan=json.loads((old/'manuscript/v80_text_plan.json').read_text(encoding='utf-8'))
plan['output_docx']='E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v81_20261025_中文修订_图3双状态绝对柔性曲线.docx'
plan['replacements']['69']=plan['replacements']['69'].replace('以下变化统一定义为IQCG态减NHE1态（图3a）。','图3直接展示两态的绝对读数；以下变化统计统一定义为IQCG态减NHE1态。')
plan['replacements']['70']=plan['replacements']['70'].replace('（图3）。','（图3b）。')
plan['replacements']['73']='接触位置的重新布置与侧链变化的空间分配相对应。仅依输入结构，以跨链Cα距离不超过8 Å界定伙伴接触，将122个侧链完整划分为仅NHE1态接触、仅IQCG态接触、两态均接触和两态均不接触四组，分别含8、21、3和90个位点（图3c）。新进入IQCG接触范围的21个侧链，MD平均Hχ由1.776降至1.573 nats，完整模型由1.857降至1.634 nats；两态均在接触范围外的90个侧链，MD平均Hχ由1.625升至1.771 nats，完整模型由1.686升至1.765 nats。两态均接触的三个位置，MD与完整模型的平均Hχ分别由1.927降至1.750、由1.811降至1.589 nats。模型因而同时恢复了新接触区域的状态集中，以及其他区域保留或增加构象占据多样性的相反趋势；完整分组避免了仅选取少数相符位点的解释。'
plan['replacements']['74']=plan['replacements']['74'].replace('全部九个Met为这种分配提供了相同化学类型的具体落点。','全部九个Met为这种分配提供了相同化学类型的具体落点（图3d、e）。')
plan['notes']=['在v80全套修订基础上，将图3曲线、分组、Met及局部位点改为两态绝对读数，更新对应图注及3.3表述。','原始MD和冻结模型预测不变，变化量准确度统计继续用于检验状态变化预测。','图5、图S2及3.4沿用已核验的v80结果。']
(m/'v81_text_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
build=(old/'manuscript/build_v80_manuscript.py').read_text(encoding='utf-8').replace('v80_text_plan.json','v81_text_plan.json').replace('v80_source_audit.json','v81_source_audit.json').replace('v80_build_audit.json','v81_build_audit.json').replace('v80_final_text.txt','v81_final_text.txt')
(m/'build_v81_manuscript.py').write_text(build,encoding='utf-8')
readme='''# 图3：两态绝对柔性读数（v81）

NHE1用暖赭，IQCG用蓝色。b为RMSF/Hχ绝对值曲线，MD及Full上下排列，同一指标统一刻度；c为完整接触分组的两态均值，d为全部九个Met，e为两个关键位置的绝对读数。a为两态分别着色的MD与预测绝对Hχ对应，不能用该绝对对应代替对变化预测的检验；ΔHχ的MAE及排序相关仍保留正文。

使用v80完全相同的原始MD复算值与冻结Full第0折预测；没有重新训练或推理。两条记录均在该折训练中，故为结构解释案例，不是留出泛化证明。五个连续100 ns窗口来自各态同一条500 ns轨迹，细线表示最小至最大值，不是独立重复的置信区间。输入伙伴、构型与原始钙条件同时不同，不能建立单因素因果关系。Hχ缺失的Ala/Gly/Pro等位置保持断线；RMSF包含整链拟合下叶间运动贡献。

运行 `python compose_fig3.py`，需要numpy、pandas、matplotlib、Pillow。该目录内source_data与structure_assets已包含绘图所需数值及PyMOL渲染。完整可编辑PyMOL场景仍见v80源文件包。结构坐标不因这次更换纵轴而变化。

输出：450 dpi PNG、嵌入TrueType字体的PDF、可编辑文本SVG；原始物理画布183 × 245 mm。all_plotted_absolute_quantities.csv保留a–d绘制数值。figure_QA.json记录全位置数值和坐标范围检查。
'''
(f/'README_统计与复现.md').write_text(readme,encoding='utf-8')
print(json.dumps({'output':plan['output_docx'],'caption_characters':len(caption),'structural_assets':len(keys)},ensure_ascii=False))