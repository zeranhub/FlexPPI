from pathlib import Path
from hashlib import sha256
import json,shutil,zipfile
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

root=Path(__file__).resolve().parent;base=root.parents[2];old=root.parent/'revision_v109_case_SASA_tracks'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
for name in ['case_profiles.py','plot_overlay.py','case_profiles_all97_with_masks.csv','case_SASA_all97_starting_complex.csv',
             'standalone_all368_plotted_values.csv','same_residue_MD_and_Full_all824988.npz','direct_relationship_statistics.json']:
    assert digest(root/name)==digest(old/name)
for folder in ['pymol','source_inputs','sasa_audit','figure3/source_data']:
    for p in (old/folder).rglob('*'):
        if p.is_file():assert digest(p)==digest(root/p.relative_to(old))
previous=pd.read_csv(old/'figure3/main_original_panels_a_to_e_values.csv',float_precision='round_trip')
current=pd.read_csv(root/'figure3/main_original_panels_a_to_e_values.csv',float_precision='round_trip')
assert len(previous)==785 and len(current)==767
retained=previous[previous.panel.ne('d')].reset_index(drop=True).dropna(axis=1,how='all')
retained['identity']=retained.identity.astype(str)
current['identity']=current.identity.astype(str)
pd.testing.assert_frame_equal(retained,current.dropna(axis=1,how='all'),check_exact=True)
previous[previous.panel.eq('d')].to_csv(root/'removed_old_d_18_values.csv',index=False)
assert current.groupby('panel').size().to_dict()=={'b':117,'c':348,'e':302}
cells=pd.read_csv(root/'overlay_cells_x0_10_grid_50_35.csv')
for label,n in [('MD',820754),('Full',824082)]:
    assert cells[label+'_positions'].sum()==n
    assert (cells.loc[cells[label+'_positions']>0,label+'_opacity']>=.1).all()
    assert (cells.loc[cells[label+'_positions']==0,label+'_opacity']==0).all()
assert digest(base/'paper_v98_20261027_中文修订_结构区域放大.docx')=='eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
page=PdfReader(root/'figure3/fig3_dual_case_entropy.pdf').pages[0];fonts=[]
for ref in page['/Resources']['/Font'].values():
    font=ref.get_object();fonts.append(str(font['/BaseFont']))
    for child in font.get('/DescendantFonts',[]):
        desc=child.get_object()['/FontDescriptor'].get_object();assert any(k in desc for k in ['/FontFile','/FontFile2','/FontFile3'])
assert all('Arial' in f for f in fonts)
with Image.open(root/'figure3/fig3_dual_case_entropy.png') as im:im.verify()
qa_path=root/'figure3/figure_QA.json';qa=json.loads(qa_path.read_text(encoding='utf-8'))
assert qa['old_d_removed'] and qa['density']['contours'] is False
qa.update(visual_review='passed: root opened complete five-panel PNG; density legend and neighboring axes clear; no clipped labels or collisions',
          embedded_fonts=fonts,all767_retained_values_exact=True,source_arrays_unchanged=True,word_unchanged=True,
          old_density_width_mm=157,new_density_width_mm=60,old_height_mm=238,new_height_mm=181)
qa_path.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
caption='''布局更新：原d（六个Gln的Hχ/SASA环境对比）从主图删除，原f的RMSF–Hχ关系图移入右侧空位，现标d。a为双案例结构；b、c为CaMKI/泛素双指标曲线及起始复合物侧链SASA；d为全量MD/Full的RMSF–Hχ分布；e仍为χ1–χ4占据关联四图横排。a/b/c/e所有原读数和原定义保留。

d使用相同824,988个共同有效位点，覆盖3,063记录，全部位置的Spearman为MD0.209、Full0.262。灰色MD与浅蓝Full透明叠加，无轮廓线；视野内所有非零格保留，不设密度裁切阈值。上方和右侧保留边缘频率分布。主视野RMSF0–10Å，MD/Full仍分别4,234/906位置范围外，完整范围参考既有全范围版；相关系数按全量计算。为适配60mm宽轴，双方一起改成50×35共同网格；频数及两组共同最大值随网格改变，但原始坐标、mask、样本和ρ未改。透明度映射仍双方共用0.10+0.25*log1p(count)/log1p(shared_peak)，零格透明，最高0.35。密度重叠不表示逐残基预测准确率。

e的Q41/Q2沿用原单角占据关联定义，和结构中Q41/Q62的环境位置标记含义不同；删除原d不改变e的全部151个单角观察及中心化口径。b/c的SASA是同PyMOL参数的起始结构几何，非MD均值、非模型预测，详见sasa_audit和v109图注。Word未修改。
'''
(root/'caption_layout_revision_zh.md').write_text(caption,encoding='utf-8')
readme='''# v110：删旧d、缩窄原f，压成五面板

原d的六Gln环境对比从主图删除，原f移入右侧空位并顺序标为d；密度主轴从157mm宽缩为60mm宽、高31mm。a/b/c/e位置关系与数据保留，e仍四张横排。整图从183×238mm缩为183×181mm。没有用截图裁切替代重排，所有曲线/坐标/文本为原生矢量。

数据：主图原785项非密度读数只移除旧d的18项（12个Hχ、6个SASA），其余767逐项严格相同。移除值存removed_old_d_18_values.csv，原source_data完整保留以追溯。368曲线值、97SASA条值、151个χ观察和全量824988密度原坐标/mask完全未变，源缓存与v109哈希相同。新d双方共用50×35网格，密度统计仍视野内820754/824082、全量rho.208615155/.261904889；全部非零格画出，零格透明，无等值轮廓。网格较旧底行80×35改变，因此每格峰值不同属于显示分箱，不是预测变化。a分子图未重绘，b/c静态SASA未重算，Word未改。

复现：assemble_main.py从v109生成figure3/compose_fig3_v110.py；后者直接用包内冻结输入、case_profiles.py和plot_overlay.py重画。源inputs、PyMOL静态几何表/脚本和完整NPZ保留于包内。release_and_audit.py核对767值完全保留、全部非零密度格、字体嵌入及Word哈希，导出PNG450dpi/PDF/SVG并打包。图注见caption_layout_revision_zh.md。
'''
(root/'README_五面板紧凑版.md').write_text(readme,encoding='utf-8')
stem='figure3_v110_20261027_删d并收紧五面板'
for ext in ['png','pdf','svg']:shutil.copy2(root/'figure3'/('fig3_dual_case_entropy.'+ext),base/(stem+'.'+ext))
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图3 v110五面板</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;margin:28px;color:#222;background:white}h1{font-size:24px}section{max-width:1050px}img{display:block;width:100%;height:auto}a{color:#2166ac}</style></head><body><h1>图3 v110：删掉原d，原f缩到右侧</h1><p>原d环境对比已删除；原f移入右侧并重标为d，主轴宽157→60mm。两个案例保留双指标曲线和下方SASA，底部四χ仍横排。整图183×181mm，Word未改。</p><section><a href="IMAGE"><img src="IMAGE" alt="五面板紧凑图3"></a><p><a href="PDF">PDF</a> · <a href="SVG">SVG</a></p></section></body></html>'''
html=html.replace('IMAGE',stem+'.png').replace('href="PDF"','href="'+stem+'.pdf"').replace('href="SVG"','href="'+stem+'.svg"')
gallery=base/'figure3_v110_20261027_五面板紧凑预览.html';gallery.write_text(html,encoding='utf-8')
package=base/'figure3_v110_20261027_五面板紧凑_整图与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(root).as_posix())
with zipfile.ZipFile(package) as z:assert z.testzip() is None
print(json.dumps({'main':str(base/(stem+'.png')),'pdf':str(base/(stem+'.pdf')),'preview':str(gallery),'package':str(package),
                  'retained_values':767,'word_unchanged':True},ensure_ascii=False,indent=2))
