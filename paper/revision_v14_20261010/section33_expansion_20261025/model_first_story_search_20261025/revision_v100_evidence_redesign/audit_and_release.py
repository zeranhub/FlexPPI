from pathlib import Path
from hashlib import sha256
from collections import Counter
from PIL import Image
from pypdf import PdfReader
import csv,json,shutil,zipfile,math

root=Path(__file__).resolve().parent;old=root.parent/'revision_v99_fig2_style_B';base=root.parents[2];out=root/'figure3'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
for folder in ['figure3/source_data','pymol']:
    for path in (old/folder).rglob('*'):
        if path.is_file():assert digest(path)==digest(root/path.relative_to(old)),str(path)
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
raw=read(out/'source_data'/'all_plotted_values.csv');main=read(out/'main_plotted_values_v100.csv')
assert len(raw)==853 and len(main)==859
def key(r):
    chi=r.get('chi','')
    return (r['panel'],r['kind'],r['identity'],r.get('series',''),str(int(float(chi))) if chi else '',r.get('metric',''))
reference={key(r):float(r['value']) for r in raw}
extra=[r for r in main if r['kind'] in ['static_sidechain_shielding_percent','sidechain_SASA_A2']]
original=[r for r in main if r not in extra]
assert len(original)==828 and len(extra)==31
for r in original:assert float(r['value'])==reference[key(r)]
missing=set(reference)-{key(r) for r in original}
assert len(missing)==25 and all(k[0]=='b' and k[3]=='Intra' for k in missing)
cam=read(out/'source_data'/'cam_all25.csv');gln=read(out/'source_data'/'ubiquitin_all6_Gln_Hchi_SASA_plot.csv')
shield={str(r['CaMKI_native_position']):100*float(r['sidechain_shielding_fraction']) for r in cam}
for r in cam:
    free=float(r['sidechain_SASA_without_CaM_same_coordinates_A2']);bound=float(r['sidechain_SASA_complex_A2'])
    assert abs((free-bound)/free-float(r['sidechain_shielding_fraction']))<1e-10
sasa={r['residue']:float(r['sidechain_SASA_A2']) for r in gln}
for r in extra:assert float(r['value'])==(shield if r['panel']=='b' else sasa)[r['identity']]
paired=read(out/'paired_record_comparison_v100.csv');wins={}
for metric in ['fine_entropy','lowF_width']:
    z=[r for r in paired if r['metric']==metric];assert len(z)==83
    wins[metric]=sum(float(r['Hchi_rho'])>float(r['RMSF_rho']) for r in z)
assert wins=={'fine_entropy':79,'lowF_width':78}
assert digest(base/'paper_v98_20261027_中文修订_结构区域放大.docx')=='eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
pdf=PdfReader(out/'fig3_dual_case_entropy.pdf');page=pdf.pages[0]
fonts=[]
for ref in page['/Resources']['/Font'].values():
    font=ref.get_object();fonts.append(str(font.get('/BaseFont')))
    for descendant in font.get('/DescendantFonts',[]):
        descriptor=descendant.get_object()['/FontDescriptor'].get_object()
        assert any(k in descriptor for k in ['/FontFile','/FontFile2','/FontFile3'])
assert all('Arial' in name for name in fonts)
for path in [out/'fig3_dual_case_entropy.png',out/'figS_context_and_record_intervals.png']:
    with Image.open(path) as im:im.verify()
qa_path=out/'figure_QA.json';qa=json.loads(qa_path.read_text(encoding='utf-8'))
qa.update(visual_review='passed: root and independent reviewer opened main and supplementary PNG; complete horizontal e row, intelligible b/d/f, no clipping/collisions',
 original_values_exactly_preserved_in_main=828,original_Intra_values_in_supplement=25,additional_preexisting_geometry_values=31,
 frozen_source_and_structure_bytes_unchanged=True,existing_word_unchanged=True,embedded_Arial_fonts=fonts,
 science_review='passed: static SASA shielding explicitly distinguished from dynamic binding effect; no monotonic SASA claim; f compares matched records without joint confidence intervals')
qa_path.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
caption='''图3｜主链运动、同类型侧链的环境差异及其角度占据关联。

a，CaM–CaMKI（1MXE）和泛素–UBA（2QHO）的真实结构。被分析链为紫色，伙伴为浅灰，关键位点为青绿色。b，上方为CaMKI残基294–318的MD参考与Full模型绝对RMSF，下方为起始结构的CaM侧链几何遮蔽比例。该比例以固定坐标删除CaM后的侧链SASA为基准，计算其与复合物中SASA的差，再除以固定坐标删除CaM后的SASA；不是游离态实验或MD，也不代表结合引起的动态变化。W303–M316区间以淡色标示。Full/Intra消融仅比较同一结合结构的不同模型上下文，移至补充稿，不能用来测量结合导致的柔性下降。c，泛素全部72序列位置的联合侧链Hχ；保留66个有效位点、原无定义位点的断线及合法Pro零值。b–d共享图例：灰色空心点为MD，蓝色实心点/线为Full。d，上方为全部六个Gln的MD/Full Hχ配对，按MD Hχ从低到高排列；下方在对应位置显示起始结构侧链SASA。Q41/Q62用于说明同类型侧链的不同局部环境和取向多样性，六位点不支持SASA与Hχ的单调关系。e，按氨基酸类型均值校正后，同一联合Hχ预测对MD χ1–χ4单角细箱占据熵的关联，四图横排，样本量66/52/23/10；Q41菱形、Q2方形。横轴复用同一联合预测，ρ为Spearman相关系数。f，同一条记录的Full RMSF占据关联作横轴，Full Hχ占据关联作纵轴，两个子图分别对应MD细箱占据熵和低F*角宽。每点是一条记录，每个子图83点；虚线为两者关联相等的参考线，79/83和78/83表示Hχ关联更高的记录数。两者关联均在记录内部控制氨基酸类型。F*=−ln[p/pmax]为无量纲相对占据势，其低值角宽不是结合自由能。f的均值及bootstrap95%区间保留于补充稿。

补充稿｜同一结合结构的上下文消融及记录级关联分布。a，完整25位MD、Full、Intra RMSF；Full和Intra的MAE分别为0.115Å和2.383Å。b，四组各83条记录的原始关联、均值及原有bootstrap95%区间；圆为Full Hχ，三角为Full RMSF。
'''
(root/'figure3_caption_zh.md').write_text(caption,encoding='utf-8')
report='''# v100 证据重排与统计说明

本版延续图2的Arial、字号、0.6pt轴线、蓝Full/灰MD等样式。修改来自用户对b/d/f意义和e排位的意见，不新增模型推理、MD或实验。

b的Intra回答模型是否利用上下文，不能证明结合诱导柔性改变。因此完整消融移到补图，主图改为50个MD/Full RMSF值和25个已有几何遮蔽值。RMSF绝对纵轴0起，Full MAE0.1148592012Å；Intra原MAE2.383203344Å未改。静态遮蔽严格由起始结构同坐标删除CaM前后的SASA计算，不应解释为动态因果。W303/M316遮蔽98.72%/97.50%；W303–M316的MD/Full平均RMSF0.553/0.513Å，余11位1.089/0.919Å，支持恢复中央低运动与两端较高运动的分布。

d的单位是同一泛素内6个Gln，不是6次独立实验。排序依MD Hχ，保留全部12纵值；下面6个SASA来自原表。MD排序Q41<Q40<Q31<Q62<Q49<Q2，Full仅交换Q40/Q31，因此不声称全排序正确。Q62暴露最高却非熵最高，不能写越暴露越灵活。没有跨残基的连线、拟合或新增P值。

e保留原151点和302个坐标、原掩码和相同轴范围；四图改为一行。样本量66/52/23/10，rho0.50/0.45/0.52/0.13，chi4弱关联如实保留。仍是同一个联合预测与不同单角参考比较。

f对相同sample_id进行配对，横轴Full RMSF的类型校正rho，纵轴Full Hχ的类型校正rho。每个目标83对，4组原始值共332个均未改。细箱熵79/83条、低F*角宽78/83条在相等参考线上方。这是记录数，非独立蛋白、复合物或湿实验数。均值0.510vs0.321、0.400vs0.239及原95%区间完整保留于补图，未构造二维联合CI，未新增显著性检验。

数据审计：全部source_data文件与v99字节相同，原853值清单保留；主图重用其中828个原始值，将25个Intra值移到补图，并从同一源表加入已有25个遮蔽值和6个SASA横坐标，因此新main_plotted_values_v100.csv为859行，不应误称原853值均仍在主图。所有PDB、PSE、PyMOL PNG和相机清单字节不变。新源表paired_record_comparison_v100.csv只把同一记录的原值并列并求差用于计数。

原生画布183×222mm，PNG450dpi，PDF嵌入Arial，SVG文本可编辑。Word未改。本版高度增加是因为e四格横排并让f保留独立完整比较；后续如需缩短，优先压缩面板间空隙，不删弱数据、不缩小至不可读字体。
'''
(root/'statistics_and_revision_notes.md').write_text(report,encoding='utf-8')
stem='figure3_v100_20261027_证据重排_e横排'
suppstem='figureS_v100_20261027_Intra与关联区间'
for ext in ['png','pdf','svg']:
    shutil.copy2(out/f'fig3_dual_case_entropy.{ext}',base/f'{stem}.{ext}')
    shutil.copy2(out/f'figS_context_and_record_intervals.{ext}',base/f'{suppstem}.{ext}')
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图3 v100 证据重排</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;color:#151515;background:#fff;margin:28px}h1{font-size:24px;margin:0 0 12px}p{max-width:1100px}main{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:28px;align-items:start}figure{margin:0}figcaption{font-size:20px;font-weight:700;margin-bottom:12px}img{width:100%;height:auto;display:block}a{color:#2166ac}@media(max-width:760px){main{grid-template-columns:1fr}}</style></head><body><h1>图3 v100：b/d/f证据重排，e四图横排</h1><p>b以MD/Full和已有静态遮蔽呈现中央低运动的结构环境；d将六个Gln的多样性与SASA对齐；f对同一记录的两预测量作直接比较。Intra和原均值区间移到补图，所有源数据保留。Word未改。点击图片打开原图。</p><main><figure><figcaption>v99</figcaption><a href="figure3_v99_20261027_B布局_图2同款样式.png"><img src="figure3_v99_20261027_B布局_图2同款样式.png" alt="v99"></a></figure><figure><figcaption>v100</figcaption><a href="MAIN"><img src="MAIN" alt="v100"></a></figure></main><p><a href="SUPP">查看补充稿：Intra与关联区间</a></p></body></html>'''.replace('MAIN',stem+'.png').replace('SUPP',suppstem+'.png')
gallery=base/'figure3_v100_20261027_证据重排对照.html';gallery.write_text(html,encoding='utf-8')
zip_path=base/'figure3_v100_20261027_证据重排_主图补图与源文件.zip'
with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
    for path in sorted(root.rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts:archive.write(path,path.relative_to(root).as_posix())
with zipfile.ZipFile(zip_path) as archive:assert archive.testzip() is None
print(json.dumps({'main':str(base/(stem+'.png')),'pdf':str(base/(stem+'.pdf')),'supplement':str(base/(suppstem+'.png')),'comparison':str(gallery),'source_package':str(zip_path),'source_files_unchanged':True,'word_unchanged':True,'main_original_values':828,'added_preexisting_geometry_values':31,'Intra_values_moved_to_supplement':25,'paired_counts':wins},ensure_ascii=False,indent=2))
