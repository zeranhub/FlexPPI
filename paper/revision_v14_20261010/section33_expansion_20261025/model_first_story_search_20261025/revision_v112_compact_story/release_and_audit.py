from pathlib import Path
from hashlib import sha256
import json,shutil,zipfile
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

root=Path(__file__).resolve().parent;base=root.parents[2]
old=root.parent/'revision_v110_compact_five_panels'
pb_root=root.parent/'paperbanana_layout_v111_20261027'
pb=pb_root/'generation_retry_20261027_155023'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
assert digest(base/'paper_v98_20261027_中文修订_结构区域放大.docx')=='eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
for name in ['case_profiles_all97_with_masks.csv','case_SASA_all97_starting_complex.csv','standalone_all368_plotted_values.csv','same_residue_MD_and_Full_all824988.npz']:
    assert digest(root/name)==digest(old/name)
for p in (old/'pymol').rglob('*'):
    if p.is_file():assert digest(p)==digest(root/p.relative_to(old))
plotted=pd.read_csv(root/'figure3/all767_values_logical_and_display_panels.csv',float_precision='round_trip',dtype={'identity':str})
source=pd.read_csv(old/'figure3/main_original_panels_a_to_e_values.csv',float_precision='round_trip',dtype={'identity':str})
keys=['panel','kind','identity','series','chi']
for z in [plotted,source]:z['chi']=z.chi.fillna(-1).astype(int)
plotted=plotted.sort_values(keys).reset_index(drop=True);source=source.sort_values(keys).reset_index(drop=True)
pd.testing.assert_frame_equal(plotted[keys],source[keys],check_dtype=False)
assert np.array_equal(plotted.value,source.value) and len(plotted)==767
assert plotted.groupby('display_panel').size().to_dict()=={'c':117,'d':348,'e':302}
cells=pd.read_csv(root/'overlay_cells_x0_10_grid_50_35.csv')
for label,num in [('MD',820754),('Full',824082)]:
    assert cells[label+'_positions'].sum()==num
    assert (cells.loc[cells[label+'_positions']>0,label+'_opacity']>=.1).all()
    assert (cells.loc[cells[label+'_positions']==0,label+'_opacity']==0).all()
fonts=[];page=PdfReader(root/'figure3/fig3_PaperBanana_layout_native.pdf').pages[0]
for ref in page['/Resources']['/Font'].values():
    font=ref.get_object();fonts.append(str(font['/BaseFont']))
    for child in font.get('/DescendantFonts',[]):
        d=child.get_object()['/FontDescriptor'].get_object();assert any(k in d for k in ['/FontFile','/FontFile2','/FontFile3'])
assert all('Arial' in f for f in fonts)
with Image.open(root/'figure3/fig3_PaperBanana_layout_native.png') as im:im.verify()
qa_path=root/'figure3/figure_QA.json';qa=json.loads(qa_path.read_text(encoding='utf-8'))
qa.update(visual_review='passed: root and two independent reviewers opened final native PNG and actual PaperBanana A; no clipped labels or collisions',embedded_Arial_fonts=fonts,
          all767_values_independently_exact=True,original_structure_pixels_preserved=True,Word_unchanged=True)
qa_path.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')

# Completion-order filenames differ from the stable independent candidate IDs.
candidates=[]
for p in sorted(pb.glob('candidate_*_pipeline.json')):
    record=json.loads(p.read_text(encoding='utf-8'));number=p.name.split('_')[1]
    image=pb/f'candidate_{number}_final.png';assert image.exists()
    candidates.append({'candidate_id':int(record['candidate_id']),'saved_name':image.name,
       'image_sha256':digest(image),'record_path':str(p),'rounds':record['max_critic_rounds'],'eval_field':record['eval_image_field']})
candidates=sorted(candidates,key=lambda x:x['candidate_id'])
assert [c['candidate_id'] for c in candidates]==[0,1]
assert all(c['rounds']==2 for c in candidates)
status={'actual_PaperBanana_called':True,'completed':True,'exp_mode':'demo_full','retrieval':'auto',
  'source_run':str(pb),'first_attempt_402_preserved_in':'generation_01/failure.json','candidates':candidates,
  'selected_candidate_id':0,'native_implementation':str(root/'figure3/fig3_PaperBanana_layout_native.png'),
  'scientific_generated_pixels_used_in_native':False,'old_Figure3_supplied':False,'Figure2_style_reference_supplied':True,'Word_unchanged':True}
(pb_root/'run_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
(root/'PaperBanana_success_and_selection.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')

native_stem='figure3_v111_20261027_PaperBanana重新布局_真实数据'
concept_stems={0:'figure3_PaperBanana_v111_20261027_布局方案A',1:'figure3_PaperBanana_v111_20261027_布局方案B'}
for ext in ['png','pdf','svg']:shutil.copy2(root/'figure3'/f'fig3_PaperBanana_layout_native.{ext}',base/(native_stem+'.'+ext))
for c in candidates:shutil.copy2(pb/c['saved_name'],base/(concept_stems[c['candidate_id']]+'.png'))

caption='''图3，结构背景下的双尺度柔性预测及侧链占据关联。（a）相邻展示CaM–CaMKI（1MXE）和泛素–UBA（2QHO）的真实输入结构：分析链紫色、伙伴灰色、关键位点青色。（b）MD与Full预测的RMSF–联合Hχ分布透明叠加及边缘频率分布，使用同824,988共同有效位点、3,063记录；全量Spearman分别0.209、0.262。MD灰、Full蓝，所有非零格保留，无等密度轮廓。主视野0–10Å，MD/Full分别4,234/906个位点位于右侧范围外，未从统计删除；完整范围参考此前全范围图。重叠分布不等同逐残基准确率。（c、d）CaMKI294–318与泛素1–72的RMSF（蓝色左轴Å）及Hχ（赭色右轴nats）沿序列分布：Full实线、MD虚线；下方为对齐的起始复合物侧链SASA（青灰柱，Å²）。两例主轴均0–2.4Å/0–4.5nats，SASA均0–160Å²。CaMKI的RMSF/Hχ有效25/21位，泛素72/66位；原无定义位置保留断线，有效孤立点和Pro合法零值保留。SASA使用相同PyMOL3.0/1.4Å探针/dot_density4，全部97位均保留，是静态几何而非模型预测或MD均值；Gly无侧链导致SASA0不代表深埋。两例Full的RMSF/Hχ MAE分别0.115Å/0.254nats、0.040Å/0.080nats。（e）预测联合Hχ与MD单角细分箱占据熵在同氨基酸类型中心化后的关联，χ1–χ4全部观察数66/52/23/10、ρ0.50/0.45/0.52/0.13，四图横排；Q41菱形、Q2方形，与结构的Q41/Q62位置对照含义不同。所有数据与上一版保持一致，只更改版面与字母编号。
'''
(root/'figure3_caption_v111_zh.md').write_text(caption,encoding='utf-8')
readme='''# 图3 v111：实际PaperBanana方案A的真实数据实现

本次服务恢复后，真实调用本机PaperBanana完整demo_full：自动参考检索→Planner→Stylist→Visualizer→两轮Critic/Visualizer，成功生成两个4096×4096方案。没有把旧图3当输入，向Planner提供的是当前科学内容与论文图2风格参考。首轮曾HTTP402，本次成功结果与首次失败记录分开保存。当前候选完成顺序与独立方案号不同：A=candidate_id0，保存在candidate_2_final.png；B=candidate_id1，保存在candidate_1_final.png。不能把framework_final.png自动当A，它是先完成的B。

选择A的三排骨架：顶排相邻两结构＋紧凑总体密度图，中排等宽等高的两个案例双轴曲线＋SASA下轨，底排四个χ散点共一行。正式图为183×190mm、Arial9/10/11pt，沿用图2细轴线和颜色。a结构，b全量密度，cCaMKI，d泛素，e单角占据关联；原六Gln面板保持删除。结构、所有曲线、SASA和散点均由真实PyMOL图片及冻结数据原生绘制，未使用生成模型给出的分子外形、点坐标或数据曲线。与v110保留全部767标量严格相同，824988密度NPZ和有效mask完全相同，合法Pro零与A/G缺口保留；图注和数据来源见manifest及caption。

原始PaperBananaA/B只作为设计草案，不宜直接作为科学图：它们含错绘结构、错曲线、错误SASA双组柱、错散点和错刻度。正式图全部纠正，SASA只有单组静态值，ρ和n按真实数据保留。生成原图及每个阶段描述保存在相邻paperbanana_layout_v111_20261027/generation_retry_20261027_155023，最终图、输入哈希、原文描述及候选选择纳入本次结果包。

可复现：compose_fig3_v111.py调用native_components.py和既有case_profiles/plot_overlay；固定输入来自prepare_frozen_inputs.py。运行validate_components.py核验组件767数值，release_and_audit.py核验最终767标量、字体嵌入、全部非零密度格、Word哈希并导出打包。PDF/SVG保留原生矢量文本/曲线，分子为原始高分辨PyMOL光栅，PNG450dpi。没有训练、推理、新MD或Word修改。
'''
(root/'README_native_components.md').write_text(readme,encoding='utf-8')
(pb_root/'README_实际调用状态.md').write_text('''# 已成功恢复并完成PaperBanana新版面生成

首次generation_01在Retriever返回402，失败诊断保留。用户恢复服务后，本次generation_retry_20261027_155023成功完成两方案4K、最多两轮Critic的完整demo_full流程。run_status.json已更新为completed=true。A稳定candidate_id0对应candidate_2_final.png，B稳定candidate_id1对应candidate_1_final.png，文件号由完成顺序决定。

已按A的三排骨架制作真实数据版，路径见run_status.json。PB原图仅为布局草案，不能将生成的分子/曲线当真实科学内容。正式图使用原PyMOL图像与767标量、824988共同位点的原生绘图；Word未修改。旧输入/失败包保留作为历史记录，当前成功结果与原生实现见figure3_v111_20261027_PaperBanana布局对照.html及成功结果ZIP。
''',encoding='utf-8')

gallery=base/'figure3_v111_20261027_PaperBanana布局对照.html'
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图3 PaperBanana新版布局</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;background:white;color:#222;margin:28px}h1{font-size:25px}h2{font-size:20px}section{max-width:1150px;margin-bottom:30px}img{display:block;width:100%;height:auto}.native{max-width:950px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:24px}.grid figure{margin:0}a{color:#2166ac}@media(max-width:750px){.grid{grid-template-columns:1fr}}</style></head><body><h1>PaperBanana重新规划的图3</h1><p>已真实完成两套4K布局及两轮检查。采用A的三排结构：顶排结构与小型总体分布，中排两案例等宽曲线与SASA，底排四χ横排。下图已全部换回真实结构和原始数据，沿用图2字体与配色，Word未替换。</p><section class="native"><h2>采用A布局的真实数据版</h2><a href="NATIVE_IMAGE"><img src="NATIVE_IMAGE" alt="PaperBanana布局的真实数据实现"></a><p><a href="NATIVE_PDF">PDF</a> · <a href="NATIVE_SVG">SVG</a> · <a href="OLD_IMAGE">上一版v110</a></p></section><section><h2>PaperBanana实际生成的布局草案</h2><p>仅供比较版面，草案里的生成结构与数值有误，不作为论文数据图。</p><div class="grid"><figure><a href="A_IMAGE"><img src="A_IMAGE" alt="PaperBanana方案A"></a><figcaption>A：两案例等宽，已采用</figcaption></figure><figure><a href="B_IMAGE"><img src="B_IMAGE" alt="PaperBanana方案B"></a><figcaption>B：另一种宽度分配</figcaption></figure></div></section></body></html>'''
for key,value in {'NATIVE_IMAGE':native_stem+'.png','NATIVE_PDF':native_stem+'.pdf','NATIVE_SVG':native_stem+'.svg',
    'OLD_IMAGE':'figure3_v110_20261027_删d并收紧五面板.png','A_IMAGE':concept_stems[0]+'.png','B_IMAGE':concept_stems[1]+'.png'}.items():html=html.replace(key,value)
gallery.write_text(html,encoding='utf-8')

package=base/'figure3_v111_20261027_PaperBanana布局与真实数据完整包.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:z.write(p,'native/'+p.relative_to(root).as_posix())
    for name in ['scientific_content.txt','visual_intent.txt','run_from_scratch.py','input_manifest.json','run_status.json','manuscript_figure2_style_reference.png','prepare.py','README_实际调用状态.md']:
        z.write(pb_root/name,'PaperBanana_inputs/'+name)
    z.write(pb_root/'style_guides/neurips2025_diagram_style_guide.md','PaperBanana_inputs/style_guide.md')
    z.write(pb/'provenance.json','PaperBanana_output/provenance.json')
    for c in candidates:
        label='A' if c['candidate_id']==0 else 'B'
        z.write(pb/c['saved_name'],f'PaperBanana_output/{label}_generated_layout_concept.png')
        z.write(Path(c['record_path']),f'PaperBanana_output/{label}_pipeline.json')
with zipfile.ZipFile(package) as z:assert z.testzip() is None
print(json.dumps({'main':str(base/(native_stem+'.png')),'pdf':str(base/(native_stem+'.pdf')),'comparison':str(gallery),
    'actual_PaperBanana_completed':True,'adopted':'A / candidate_id0','source_values_unchanged':767,'Word_unchanged':True,
    'package':str(package),'package_bytes':package.stat().st_size},ensure_ascii=False,indent=2))
