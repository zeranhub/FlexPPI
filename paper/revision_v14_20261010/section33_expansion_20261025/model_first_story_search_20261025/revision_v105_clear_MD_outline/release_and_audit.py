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
    assert ((z[['R','G','B','opacity']]>=0)&(z[['R','G','B','opacity']]<=1)).all().all()

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
    qa.update(visual_review='passed: root opened standalone central, full-range and complete main PNGs; no clipped labels or collisions',
              embedded_fonts=fonts,source_arrays_unchanged=True,word_unchanged=True)
    p.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')

caption='''f，MD参考与Full预测的RMSF–联合侧链Hχ分布。浅蓝色底图为Full预测逐格频数，深灰色轮廓为同网格MD参考频数；MD不填色，轮廓始终置于蓝色底图上方，避免重叠后混色。两组均使用同824,988个共同有效残基位置，覆盖3,063条记录。蓝色深浅按共同网格中两组最大计数确定的对数强度显示；MD虚线/实线分别是每格频率0.01%/0.1%的等值线，基于原始六边形格中心计数作线性插值，未作核平滑、拟合或抽样。这两个频率阈值不是包含某比例样本的累计概率区域。上方、右侧为共同分箱和分母的边缘频率轮廓（MD灰色虚线，Full蓝色实线）。Spearman相关按全部824,988位置计算，MD为0.209、Full为0.262。主视野RMSF为0–10Å，MD/Full范围外分别4,234/906位置；完整0–35Å版保留全部长尾。有效Pro及合法Hχ零值均保留，Ala/Gly按原标签定义排除。分布区域重叠不等同于同一残基预测逐点吻合；横带与残基类型、χ数量和离散状态定义有关，弱总体相关不代表严格独立。

正文可用：MD参考与Full预测均显示，主链位置波动与侧链取向多样性存在有限的整体关联，同一RMSF范围内仍可出现宽范围Hχ值。模型保留了这种分布特征，但总体相关程度略强于MD。因此，此图用于呈现两类动态摘要的对应关系。
'''
(root/'caption_and_claim_zh.md').write_text(caption,encoding='utf-8')
readme='''# v105：MD深灰轮廓与Full浅蓝密度

用户指出v104的暖灰与蓝色混色难看，且MD被淹没。本版取消双层填色与混色，使用Full浅蓝hexbin底图、上层MD深灰密度轮廓。颜色恢复为论文MD灰#62676B与Full蓝#2166AC。边缘分布仍为共享轴两条轮廓。主图a–e原值及位置保留，e四张横排，f单图；Word未修改。

数据与划分：原五折OOF表variant=full，固定epoch40，对应留出折预测；两目标mask均有效、四坐标有限的824,988个位置，来自3,063条记录。原full总961,321行，共同mask排除136,333行；A/G无定义χ标签按原mask排除；44,822个有效Pro（MD Hχ为零44,673个）全部保留。没有新增训练、模型推理、MD或统计检验。统计为残基池化的描述Spearman，长记录权重较高，无P值、无置信区间。

四坐标/位置身份缓存NPZ、统计JSON和cohort表与v103哈希完全一致；其来源逐项原值审查已在v103完成，本次再次核对原CSV哈希。rho：MD0.20861515505308，Full0.2619048885568794，始终按全量算。各视图两组使用相同网格、边界，格计数总和与可见位置数完全相等。主视野0–10Å：MD820754可见、4234范围外；Full824082可见、906范围外。完整0–35Å双方各824988可见。

密度：Full的透明强度=0.52*log1p(Full_count)/log1p(两组共同peak)，零计数透明。MD不填色，只以上层深灰轮廓显示；外层虚线每格频率0.01%，内层实线每格频率0.1%，分母为全部824988，对应计数阈值82.4988/824.988。轮廓由原始hexbin中心计数的三角网格作线性插值生成；没有核平滑、拟合或抽样，但存在常规等值线线性插值，不能声称完全无插值。两阈值是预设逐格频率，不是累计概率区间。原始两组格计数、Full颜色和MD轮廓全部顶点另存CSV。MD所有有效零值均保留，零附近的少量锯齿来自原始离散网格。每幅图内部双方同网格；主视野、完整范围和整图尺寸不同有相应网格数，不将跨网格的轮廓面积当定量比较。分布覆盖相同区域不等同于逐残基预测准确率。上下边缘用相同hist箱和824988分母，MD灰虚线/Full蓝实线。

复现：plot_overlay.py生成独立图和全范围版，assemble_main.py生成figure3/compose_fig3_v105.py以调用原a–e原生绘图逻辑。已有缓存可直接运行plot_overlay.py与figure3/compose_fig3_v105.py重画，无需访问大CSV。PDF/SVG为原生矢量绘图，Arial嵌入，PNG450dpi。release_and_audit.py检验源值、网格计数、字体、导出、Word哈希并打包。
'''
(root/'README_叠加图与复现.md').write_text(readme,encoding='utf-8')

main='figure3_v105_20261027_MD轮廓与预测密度';single='figure3_f_v105_20261027_MD轮廓与预测密度';full='figure3_f_v105_20261027_MD轮廓完整范围'
for ext in ['png','pdf','svg']:
    shutil.copy2(root/'figure3'/('fig3_dual_case_entropy.'+ext),base/(main+'.'+ext))
    shutil.copy2(root/('RMSF_Hchi_overlay.'+ext),base/(single+'.'+ext))
    shutil.copy2(root/('RMSF_Hchi_overlay_full_range.'+ext),base/(full+'.'+ext))
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图3 v105：MD轮廓与预测密度</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;margin:28px;color:#222;background:white}h1{font-size:24px}section{max-width:1150px;margin-bottom:30px}img{display:block;width:100%;height:auto}.single{max-width:800px}a{color:#2166ac}</style></head><body><h1>图3 v105：浅蓝预测密度与深灰MD轮廓</h1><p>MD轮廓置于预测密度之上，重叠处也保留MD颜色。外层虚线/内层实线分别对应MD每格频率0.01%/0.1%，并非累计概率区间。两组均为相同824,988个有效位置，相关系数与原版不变。Word未替换。</p><section class="single"><a href="SINGLE"><img src="SINGLE" alt="浅蓝预测密度与深灰MD轮廓"></a><p><a href="FULL">完整0–35 Å范围</a> · <a href="SPDF">独立图PDF</a></p></section><section><a href="MAIN"><img src="MAIN" alt="组图3"></a><p><a href="MPDF">组图3 PDF</a></p></section></body></html>'''
html=html.replace('SINGLE',single+'.png').replace('FULL',full+'.png').replace('SPDF',single+'.pdf').replace('MAIN',main+'.png').replace('MPDF',main+'.pdf')
gallery=base/'figure3_v105_20261027_MD轮廓叠加预览.html';gallery.write_text(html,encoding='utf-8')
package=base/'figure3_v105_20261027_MD轮廓叠加_整图与全量源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(root).as_posix())
with zipfile.ZipFile(package) as z:assert z.testzip() is None
print(json.dumps({'single':str(base/(single+'.png')),'main':str(base/(main+'.png')),'preview':str(gallery),
                  'package':str(package),'package_bytes':package.stat().st_size,'source_values_unchanged':True,'word_unchanged':True},ensure_ascii=False,indent=2))
