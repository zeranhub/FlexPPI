from pathlib import Path
from hashlib import sha256
import json,shutil,zipfile
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

root=Path(__file__).resolve().parent;base=root.parents[2]
old=root.parent/'revision_v103_direct_RMSF_Hchi'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
for name in ['same_residue_MD_and_Full_all824988.npz','direct_relationship_statistics.json','cohort_all3063_records.csv']:
    assert digest(root/name)==digest(old/name)
stats=json.loads((root/'direct_relationship_statistics.json').read_text(encoding='utf-8'))
assert digest(Path(stats['source']))==stats['source_sha256']
assert stats['common_valid_positions']==824988 and stats['records']==3063
for folder in ['figure3/source_data','pymol']:
    for p in (old/folder).rglob('*'):
        if p.is_file():assert digest(p)==digest(root/p.relative_to(old))
pd.testing.assert_frame_equal(pd.read_csv(old/'figure3/main_original_panels_a_to_e_values.csv',float_precision='round_trip'),
                             pd.read_csv(root/'figure3/main_original_panels_a_to_e_values.csv',float_precision='round_trip'),check_exact=True)
assert digest(base/'paper_v98_20261027_中文修订_结构区域放大.docx')=='eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
for p in root.glob('overlay_cells_*.csv'):
    z=pd.read_csv(p)
    for label in ['MD','Full']:
        expected=(824988 if '_x0_35_' in p.name else {'MD':820754,'Full':824082}[label])
        assert z[label+'_positions'].sum()==expected
    for label in ['MD','Full']:
        color_columns=[label+'_'+c for c in ['R','G','B','opacity']]
        assert ((z[color_columns]>=0)&(z[color_columns]<=1)).all().all()
        assert (z.loc[z[label+'_positions']>0,label+'_opacity']>=.10).all()
        assert (z[label+'_opacity']<=.35+1e-14).all()
        assert (z.loc[z[label+'_positions']==0,label+'_opacity']==0).all()
assert not list(root.glob('MD_contour_vertices_*.csv'))

for stem in ['RMSF_Hchi_overlay','RMSF_Hchi_overlay_full_range','figure3/fig3_dual_case_entropy']:
    page=PdfReader(root/(stem+'.pdf')).pages[0];fonts=[]
    for ref in page['/Resources']['/Font'].values():
        font=ref.get_object();fonts.append(str(font['/BaseFont']))
        for child in font.get('/DescendantFonts',[]):
            d=child.get_object()['/FontDescriptor'].get_object()
            assert any(k in d for k in ['/FontFile','/FontFile2','/FontFile3'])
    assert all('Arial' in f for f in fonts)
    with Image.open(root/(stem+'.png')) as im:im.verify()
    p=root/('figure3/figure_QA.json' if stem.startswith('figure3/') else stem+'_QA.json')
    qa=json.loads(p.read_text(encoding='utf-8'))
    density=qa.get('density',qa)
    assert density['contours'] is False and all(density['all_nonzero_cells_visible'].values())
    qa.update(visual_review='passed: root opened standalone central, full-range and complete main PNGs; no clipped labels or collisions',
              embedded_fonts=fonts,source_arrays_unchanged=True,word_unchanged=True)
    p.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')

caption='''f，MD参考与Full预测的RMSF–联合侧链Hχ分布。MD灰色密度层与Full浅蓝色密度层透明叠加，重叠处自然呈稍深的蓝灰色，不加任何等密度轮廓线。双方在当前视野内所有非零计数格均绘出，不设密度显示阈值，零格留白。两组均使用同824,988个共同有效残基位置，覆盖3,063条记录，采用相同网格、共同最大计数和相同对数计数到不透明度映射；非零映射下限0.10、最高0.35。上方、右侧为共同分箱和分母的边缘频率分布（MD灰虚线，Full蓝实线）。Spearman相关按全部824,988位置计算，MD为0.209、Full为0.262。主视野RMSF为0–10Å，MD/Full范围外分别4,234/906位置；完整0–35Å版包含双方全部位置。未作平滑、插值、拟合或抽样。有效Pro及合法Hχ零值保留，Ala/Gly按原标签定义排除。分布区域重叠不等同于同一残基预测逐点吻合；横带与残基类型、χ数量和离散状态定义有关，弱总体相关不代表严格独立。

正文可用：MD参考与Full预测均显示，主链位置波动与侧链取向多样性存在有限的整体关联，同一RMSF范围内仍可出现宽范围Hχ值。模型保留了这种分布特征，但总体相关程度略强于MD。因此，此图用于呈现两类动态摘要的对应关系。
'''
(root/'caption_and_claim_zh.md').write_text(caption,encoding='utf-8')
readme='''# v107：无轮廓线、柔和透明叠加

按用户要求删除f中所有等密度轮廓线，仅使用完整MD灰色hexbin层与完整Full浅蓝层透明叠加，重叠处稍深、略偏蓝灰。双方所有非零格均显示，MD右侧稀疏格保留。MD灰#62676B、Full蓝#2166AC，图例也使用与实际显示一致的浅色块。边缘分布仍为共享轴两条频率线。主图a–e原值及位置保留，e四张横排，f单图；Word未修改。

数据与划分：原五折OOF表variant=full，固定epoch40，对应留出折预测；两目标mask均有效、四坐标有限的824,988个位置，来自3,063条记录。原full总961,321行，共同mask排除136,333行；A/G无定义χ标签按原mask排除；44,822个有效Pro（MD Hχ为零44,673个）全部保留。没有新增训练、模型推理、MD或统计检验。统计为残基池化的描述Spearman，长记录权重较高，无P值、无置信区间。

四坐标/位置身份缓存NPZ、统计JSON和cohort表与v103哈希完全一致；其来源逐项原值审查已在v103完成，本次再次核对原CSV哈希。rho：MD0.20861515505308，Full0.2619048885568794，始终按全量算。各视图两组使用相同网格、边界，格计数总和与可见位置数完全相等。主视野0–10Å：MD820754可见、4234范围外；Full824082可见、906范围外。完整0–35Å双方各824988可见。

密度：双方非零格不透明度alpha=0.10+0.25*log1p(count)/log1p(两组共同peak)，最高0.35，零格为0；相较v106的最高0.50，整体更淡。共享正值下限保留稀疏MD尾格，不分别归一双方峰值。每个非零MD/Full格验证alpha>0、零格alpha0。MD先铺灰底，Full叠浅蓝，自然混色。没有等密度轮廓、核平滑、插值、拟合或抽样。两组原始格计数与每层RGBA另存CSV；没有90%概率边界或任何密度裁切。所有合法MD零Hχ值保留。每幅图内部双方同网格；不同尺寸/完整范围有不同网格，不跨网格作面积定量比较。分布重叠不等于逐残基准确率；显示深浅用于密度直观表达，重叠色不是第三类数据。上下边缘用相同hist箱、824988分母，MD灰虚线/Full蓝实线。

复现：plot_overlay.py生成独立图和全范围版，assemble_main.py生成figure3/compose_fig3_v107.py以调用原a–e原生绘图逻辑。已有缓存可直接运行plot_overlay.py与figure3/compose_fig3_v107.py重画，无需访问大CSV。PDF/SVG为原生矢量绘图，Arial嵌入，PNG450dpi。release_and_audit.py检验源值、网格计数、全部非零格alpha、无等密度轮廓、字体、导出及Word哈希并打包。
'''
(root/'README_叠加图与复现.md').write_text(readme,encoding='utf-8')

main='figure3_v107_20261027_无轮廓柔和叠加';single='figure3_f_v107_20261027_无轮廓柔和叠加';full='figure3_f_v107_20261027_无轮廓叠加完整范围'
for ext in ['png','pdf','svg']:
    shutil.copy2(root/'figure3'/('fig3_dual_case_entropy.'+ext),base/(main+'.'+ext))
    shutil.copy2(root/('RMSF_Hchi_overlay.'+ext),base/(single+'.'+ext))
    shutil.copy2(root/('RMSF_Hchi_overlay_full_range.'+ext),base/(full+'.'+ext))
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图3 v107：无轮廓柔和叠加</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;margin:28px;color:#222;background:white}h1{font-size:24px}section{max-width:1150px;margin-bottom:30px}img{display:block;width:100%;height:auto}.single{max-width:800px}a{color:#2166ac}</style></head><body><h1>图3 v107：去掉轮廓，降低不透明度</h1><p>浅灰MD与浅蓝Full透明叠加，重叠处自然稍深、略偏蓝灰。两组均使用0.10–0.35的共同不透明度映射；当前视野内所有非零格保留，包括右侧稀疏MD数据。数据与相关系数不变，Word未替换。</p><section class="single"><a href="SINGLE"><img src="SINGLE" alt="无轮廓的浅灰浅蓝密度叠加"></a><p><a href="FULL">双方全部位置的0–35 Å范围</a> · <a href="SPDF">独立图PDF</a></p></section><section><a href="MAIN"><img src="MAIN" alt="组图3"></a><p><a href="MPDF">组图3 PDF</a></p></section></body></html>'''
html=html.replace('SINGLE',single+'.png').replace('FULL',full+'.png').replace('SPDF',single+'.pdf').replace('MAIN',main+'.png').replace('MPDF',main+'.pdf')
gallery=base/'figure3_v107_20261027_无轮廓透明叠加预览.html';gallery.write_text(html,encoding='utf-8')
package=base/'figure3_v107_20261027_无轮廓透明叠加_整图与全量源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(root).as_posix())
with zipfile.ZipFile(package) as z:assert z.testzip() is None
print(json.dumps({'single':str(base/(single+'.png')),'main':str(base/(main+'.png')),'preview':str(gallery),
                  'package':str(package),'package_bytes':package.stat().st_size,'source_values_unchanged':True,'word_unchanged':True},ensure_ascii=False,indent=2))
