from pathlib import Path
from hashlib import sha256
from PIL import Image
from pypdf import PdfReader
from collections import Counter
import csv
import json
import shutil
import zipfile

root=Path(__file__).resolve().parent
base=root.parents[2]
old=root.parent/'revision_v98_structure_fill'
out=root/'figure3'
digest=lambda p:sha256(p.read_bytes()).hexdigest()

for folder in ['figure3/source_data','pymol']:
    for path in (old/folder).rglob('*'):
        if path.is_file():
            assert digest(path)==digest(root/path.relative_to(old)),str(path)
with (out/'source_data'/'all_plotted_values.csv').open(encoding='utf-8-sig',newline='') as f:
    rows=list(csv.DictReader(f))
assert len(rows)==853
counts=dict(Counter(row['panel'] for row in rows))
assert counts=={'b':75,'c':132,'d':12,'e':302,'f':332}
assert digest(out/'source_data'/'all_plotted_values.csv')=='233dcbf5e9222b248d46db9f5b377bde605c0905877209bfb8b375543b1a6b21'
assert digest(base/'paper_v98_20261027_中文修订_结构区域放大.docx')=='eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
fig2=base/'figure2_dense_20261024'/'fig2_dense.png'
assert digest(fig2)=='3695952ad99dca46f341ab4dc0f80871a4f3f1eaf24b5df9be9b1f182063f9e6'
pdf=PdfReader(out/'fig3_dual_case_entropy.pdf')
assert len(pdf.pages)==1
page=pdf.pages[0]
dimensions_mm=[float(page.mediabox.width)*25.4/72,float(page.mediabox.height)*25.4/72]
assert abs(dimensions_mm[0]-183)<.01 and abs(dimensions_mm[1]-202)<.01
fonts=[]
for ref in page['/Resources']['/Font'].values():
    font=ref.get_object()
    descendants=font.get('/DescendantFonts',[])
    descriptors=[]
    if '/FontDescriptor' in font:descriptors.append(font['/FontDescriptor'].get_object())
    for item in descendants:
        desc=item.get_object().get('/FontDescriptor')
        if desc:descriptors.append(desc.get_object())
    assert descriptors and all(any(k in d for k in ['/FontFile','/FontFile2','/FontFile3']) for d in descriptors)
    fonts.append(str(font.get('/BaseFont')))
assert all('Arial' in name for name in fonts),fonts
svg=(out/'fig3_dual_case_entropy.svg').read_text(encoding='utf-8')
assert '<text' in svg and svg.count('<image')==2
with Image.open(out/'fig3_dual_case_entropy.png') as im:
    raster_size=list(im.size)
    im.verify()
qa_path=out/'figure_QA.json'
qa=json.loads(qa_path.read_text(encoding='utf-8'))
qa['visual_review']='passed: root and independent reviewer opened full PNG and actual Figure2; no clipping/overlap or data truncation; coordinated fonts/strokes/palette'
qa['independent_data_review']='all source CSVs and every PyMOL asset byte-identical to v98; exact numerical equivalence confirmed independently'
qa['vector_audit']={'dimensions_mm':dimensions_mm,'embedded_fonts':fonts,'svg_editable_text':True,'svg_native_structure_raster_count':2,'raster_size':raster_size}
qa['existing_word_unchanged']=True
qa_path.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')

caption='''图3｜复合物上下文中的主链柔性、侧链局部环境与角度占据。
a，CaM–CaMKI（1MXE）与泛素–UBA（2QHO）的真实结构；被分析链为紫色，伙伴为浅灰，关键位点为青绿色。b，CaMKI残基294–318的绝对主链RMSF。Full与Intra使用同一结合态结构，分别采用完整复合物上下文与链内上下文；MD为参考。c，泛素全序列的联合侧链Hχ，保留原始无定义位点的断线与合法零值。b–d共享模型配色：Full蓝色、Intra棕色虚线、MD灰色；d中空心圆为MD、实心蓝圆为Full。d，全部6个Gln的侧链SASA与Hχ配对，竖线只连接同一残基的MD和Full结果；Q41/Q62突出同类型侧链在不同局部环境下的差异，六个位点并不支持SASA与Hχ的单调关系。e，按氨基酸类型均值校正后，同一联合Hχ预测与MD中χ1–χ4单角细箱占据熵的关联，分别保留66、52、23、10个位点；横轴复用同一联合预测，青绿色菱形和方形分别为Q41和Q2。ρ为Spearman相关系数。f，83条记录各自按氨基酸类型校正后的关联分布，分别以Full Hχ和Full RMSF对应MD细箱占据熵与低F*角宽；每组保留全部83个点，醒目标记与误差线分别为均值和既有bootstrap 95%区间。F*=−ln[p/pmax]为角度占据的无量纲相对势。图中所有数值和结构均沿用冻结源数据，未使用生成式概念稿中的分子或曲线。
'''
(root/'figure3_caption_zh.md').write_text(caption,encoding='utf-8')
statistics='''# 数值与复现说明

本次为B版空间布局与图2样式统一，不新增模型训练、推理、MD、筛样或显著性检验。已冻结模型结果和参考的协议、fold、seed、checkpoint等字段仍见各原始源表及原文；本次未改动这些字段。统计结果与上一版逐字段完全相同。

b是同一CaMKI肽段的25个残基，不是25个独立实验；MD、Full、Intra各25值，共75值。两模型MAE分别为0.1148592012Å、2.383203344Å，按25位点绝对误差均值计算。源列为cam_all25.csv的CaMKI_native_position、target_rmsf、pred_rmsf、intra_rmsf。

c是泛素72个序列位置中的66个有效位点，MD和Full各66值，共132值。原始6个A/G无定义位点仍断线；三个Pro的合法零参考保留。Full MAE为0.07953050649754886nats。源表ubiquitin_72_position_profile_with_original_mask.csv及其valid_hchi掩码完整未动。

d全部6个Gln各有1个MD和1个Full Hχ，共12纵值，另使用原始6个sidechain_SASA_A2横坐标。无回归、无因果推断，不新增显著性检验。源表ubiquitin_all6_Gln_Hchi_SASA_plot.csv。空心点为MD、蓝色实心点为Full；只有同一残基的竖向配对线。

e按chi1–chi4分组的位点数为66/52/23/10，总151个点、302个坐标。双方先减去各自的氨基酸类型均值。横轴是同一联合Hχ预测，纵轴为相应MD单角fine-bin entropy；不是四个独立预测头。Spearman rho沿用per_chi_correlations.csv：0.4975052713、0.4471864059、0.5247035573、0.1272727273（显示两位小数）。原始全部点来自defined_angles_long.csv。未画回归线或新增P值。

f的样本单位为record，每组83值，共4组332值。record内部先按氨基酸类型校正，再计算Full预测与MD占据描述量的Spearman rho。记录级点的横向抖动只用于显示、不参与统计。均值及bootstrap95%区间沿用retained_panel_statistics.csv；未重新采样或改变区间。细箱熵Hχ均值0.5098929841，CI0.4907300073–0.5274475608；RMSF均值0.3205688592，CI0.3052365755–0.3359174296。低F*角宽Hχ均值0.3995555843，CI0.3818306871–0.4160917128；RMSF均值0.2385388632，CI0.2206801593–0.2554976387。原始记录值来自retained_absolute_correlations_all498.csv，其中Full两预测量的四组用于主图，Intra源值仍随原表保留。

主图冻结长表all_plotted_values.csv有853行，哈希233dcbf5e9222b248d46db9f5b377bde605c0905877209bfb8b375543b1a6b21，与v98字节相同。SASA横坐标和均值区间也完全保留。两PyMOL原图、PSE、PDB及相机清单字节不变。

样式取自当前论文图2figure2_dense_20261024：Arial，轴标/刻度/图例9pt，普通标题10pt，面板字母11pt，轴线0.6pt，Full#2166AC实线1.25pt，Intra#B87936虚线1.1pt。按图2相同Word宽度167.6mm显示时，普通文字约8.24pt。图3原生画布183×202mm，PNG450dpi，PDF嵌入Arial字体，SVG保留可编辑文本；分子为原始真实栅格渲染。Word未替换。
'''
(root/'statistics_and_reproducibility.md').write_text(statistics,encoding='utf-8')

release_stem='figure3_v99_20261027_B布局_图2同款样式'
for ext in ['png','pdf','svg']:
    shutil.copy2(out/f'fig3_dual_case_entropy.{ext}',base/f'{release_stem}.{ext}')
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图2与图3统一样式对照</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;color:#151515;background:#fff;margin:28px}h1{font-size:24px;margin:0 0 12px}p{max-width:1100px}main{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:28px;align-items:start}figure{margin:0}figcaption{font-size:20px;font-weight:700;margin-bottom:12px}img{width:100%;height:auto;display:block}a{color:#2166ac}@media(max-width:760px){main{grid-template-columns:1fr}}</style></head><body><h1>图2与新图3：同宽样式对照</h1><p>图3采用B的空间安排，按论文当前图2统一Arial、字号、细轴线、蓝橙方法色及点线样式。所有结构与853个主图数值均来自真实冻结源文件。两图等宽展示，点击打开原始PNG。现有Word未替换。</p><main><figure><figcaption>论文当前图2</figcaption><a href="figure2_dense_20261024/fig2_dense.png"><img src="figure2_dense_20261024/fig2_dense.png" alt="当前图2"></a></figure><figure><figcaption>图3 v99：B布局、图2样式</figcaption><a href="NEWIMAGE"><img src="NEWIMAGE" alt="图3统一样式版"></a></figure></main></body></html>'''.replace('NEWIMAGE',release_stem+'.png')
gallery=base/'figure3_v99_20261027_与图2同宽对照.html'
gallery.write_text(html,encoding='utf-8')
readme='''# 图3 v99：B布局与图2统一样式

figure3/compose_fig3_v99.py是可直接运行的原生绘图脚本；同目录source_data包含完整冻结源表，publication_style.py为论文图2的真实共享配色模块。pymol为不变的真实结构资产。以安装numpy、pandas、matplotlib、Pillow的Python运行该脚本，即在figure3目录重建PNG/PDF/SVG。

本版采用PaperBanana B概念的a/b、c/d、e/f空间组织，未使用其合成分子、曲线、散点或标签。样式基准是当前v98 Word内图2的实际图像，已通过哈希匹配。figure3/figure_QA.json保存样式、字号、数值计数、裁边、嵌入字体和独立视觉审查结果；statistics_and_reproducibility.md给出统计单位与源列；figure3_caption_zh.md是包含共享图例及d空心/实心定义的图注草稿。本次不修改正文或Word。
'''
(root/'README.md').write_text(readme,encoding='utf-8')
zip_path=base/'figure3_v99_20261027_B布局_图2同款样式_源文件与核验.zip'
with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
    for path in sorted(root.rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts:
            archive.write(path,path.relative_to(root).as_posix())
with zipfile.ZipFile(zip_path) as archive:assert archive.testzip() is None
print(json.dumps({'png':str(base/(release_stem+'.png')),'pdf':str(base/(release_stem+'.pdf')),
 'svg':str(base/(release_stem+'.svg')),'comparison':str(gallery),'source_zip':str(zip_path),
 'font_names':fonts,'all_853_values_byte_identical':True,'all_structure_files_byte_identical':True,
 'word_unchanged':True,'source_zip_bytes':zip_path.stat().st_size},ensure_ascii=False,indent=2))
