from pathlib import Path
from hashlib import sha256
import json,shutil,zipfile
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

root=Path(__file__).resolve().parent;base=root.parents[2];old=root.parent/'revision_v108_dual_axis_case_profiles'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
manifest=json.loads((root/'case_profile_manifest.json').read_text(encoding='utf-8'))
for source in manifest['sources']:
    assert digest(Path(source['original']))==digest(Path(source['frozen']))==source['sha256']
assert digest(root/'same_residue_MD_and_Full_all824988.npz')==digest(old/'same_residue_MD_and_Full_all824988.npz')
assert digest(root/'plot_overlay.py')==digest(old/'plot_overlay.py')
assert digest(root/'overlay_cells_x0_10_grid_80_35.csv')==digest(old/'overlay_cells_x0_10_grid_80_35.csv')
for p in (old/'pymol').rglob('*'):
    if p.is_file():assert digest(p)==digest(root/p.relative_to(old))
assert digest(base/'paper_v98_20261027_中文修订_结构区域放大.docx')=='eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
data=pd.read_csv(root/'case_profiles_all97_with_masks.csv',float_precision='round_trip')
drawn=pd.read_csv(root/'standalone_all368_plotted_values.csv',float_precision='round_trip')
assert len(data)==97 and len(drawn)==368
for (case,metric,series),z in drawn.groupby(['case','metric','series']):
    mask='valid_rmsf' if metric=='RMSF' else 'valid_hchi';unit='A' if metric=='RMSF' else 'nats'
    reference=data[data.case.eq(case)&data[mask]].sort_values('residue_id');z=z.sort_values('residue_id')
    assert np.array_equal(z.residue_id,reference.residue_id)
    assert np.array_equal(z.value,reference[series+'_'+metric+'_'+unit])
previous=pd.read_csv(old/'figure3/main_original_panels_a_to_e_values.csv',float_precision='round_trip')
current=pd.read_csv(root/'figure3/main_original_panels_a_to_e_values.csv',float_precision='round_trip')
assert len(previous)==713 and len(current)==785
previous=previous[previous.kind.ne('static_sidechain_shielding_percent')]
keys=[c for c in previous.columns if c!='value']
matched=previous.merge(current,on=keys,suffixes=('_old','_new'),validate='one_to_one')
assert len(matched)==688
precision_delta=float(abs(matched.value_old-matched.value_new).max());assert precision_delta<=1e-14
sasa=pd.read_csv(root/'case_SASA_all97_starting_complex.csv',float_precision='round_trip')
assert len(sasa)==97 and not sasa.duplicated(['case','residue_id']).any()
for case,panel in [('CaMKI','b'),('Ubiquitin','c')]:
    expected=sasa[sasa.case.eq(case)].sort_values('residue_id')
    plotted=current[current.panel.eq(panel)&current.kind.eq('static_sidechain_SASA_A2')].sort_values('identity',key=lambda s:s.astype(int))
    assert np.array_equal(expected.residue_id,plotted.identity.astype(int))
    assert np.array_equal(expected.sidechain_SASA_A2,plotted.value)
    assert expected.sidechain_SASA_A2.between(0,160).all()
assert not current.kind.eq('static_sidechain_shielding_percent').any()
for stem in ['case_RMSF_Hchi_dual_axis','figure3/fig3_dual_case_entropy']:
    page=PdfReader(root/(stem+'.pdf')).pages[0];fonts=[]
    for ref in page['/Resources']['/Font'].values():
        font=ref.get_object();fonts.append(str(font['/BaseFont']))
        for child in font.get('/DescendantFonts',[]):
            desc=child.get_object()['/FontDescriptor'].get_object();assert any(k in desc for k in ['/FontFile','/FontFile2','/FontFile3'])
    assert all('Arial' in f for f in fonts)
    with Image.open(root/(stem+'.png')) as im:im.verify()
    qp=root/('figure3/figure_QA.json' if stem.startswith('figure3/') else 'case_profiles_QA.json')
    qa=json.loads(qp.read_text(encoding='utf-8'));qa.update(visual_review='passed: root and two independent reviewers opened standalone and complete main PNGs; axes, labels, masks and singleton points visible',
       embedded_fonts=fonts,profile_values_exact_from_frozen_sources=True,f_density_unchanged=True,word_unchanged=True,
       historical_CSV_precision_max_difference=precision_delta)
    qp.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')

caption='''b、c，CaMKI调控片段（1MXE，原生294–318）与泛素（2QHO，链A1–72）的主链RMSF、联合侧链Hχ及侧链溶剂暴露。上层蓝色左轴为RMSF（Å），赭色右轴为Hχ（nats）；Full实线、MD虚线。下层青灰色柱图为相同残基的起始复合物侧链溶剂可及面积SASA（Å²），两层共用残基横轴，两例SASA均为0–160Å²。SASA在精确起始完整复合物中计算，统一PyMOL3.0.0、探针1.4Å、dot_density4，排除氢原子及侧链定义中的N/CA/C/O/OXT/OT1/OT2；它是静态结构环境，不是模型输出，也不是MD时间平均值。CaMKI25位SASA按统一方法重算，泛素72位沿用与d面板同源的既有计算。SASA不受Hχ有效性mask限制，故全部97位保留；Gly无侧链重原子导致侧链SASA为0，不代表深埋。两例曲线均用RMSF0–2.4Å/Hχ0–4.5nats固定零起点，无Hχ定义位置两来源均留空，小点保留有效孤立位置，合法Pro零Hχ保留。CaMKI RMSF/Hχ有效25/21位，泛素72/66位；Full MAE分别0.115Å/0.254nats和0.040Å/0.080nats。CaMKI为家族留出fold1、泛素为chemical_identity_corrected_v1原记录留出fold0，均seed7、epoch40对应留出折单模型，未新增训练、模型预测或MD。左右轴单位不同，曲线交叉不表示数值相等；曝光与柔性并排呈现不建立因果关系。
'''
(root/'caption_dual_axis_profiles_zh.md').write_text(caption,encoding='utf-8')
readme='''# v109：双指标曲线下增加对齐的溶剂暴露柱图

主图b、c保留同残基横轴的双纵轴曲线，并在正下方增加青灰色侧链SASA柱图。原b的CaM遮蔽百分比已换成侧链溶剂可及面积，c新增同口径柱图。上下层严格对齐，仅底层显示残基编号。SASA两例统一0–160Å²；RMSF0–2.4Å/Hχ0–4.5nats仍不变。主曲线368有效读数、mask、合法零和孤立点均不变。a/d/e/f数据保持；e/f向下移8mm为SASA轨道留空间。Word未修改。

SASA来源及统一：泛素all119_input_environment.csv的A链72位与d的六Gln同源；原CaM SASA由MDTraj计算，因此本版用PyMOL3.0.0在精确原始CaM复合物重算25位，采用与泛素一致的1.4Å探针/dot_density4，复合物全部重原子参与遮挡，侧链排除N CA C O OXT OT1 OT2。新几何表与输入哈希、PyMOL版本、参数及脚本位于sasa_audit/。只做静态几何，不是模型推理/MD；SASA不受Hχ mask限制，97位置全部绘出。泛素G10/G35/G47/G53无侧链重原子，侧链SASA0表示几何定义，不是深埋；A28/A46仍有SASA柱。case_SASA_all97_starting_complex.csv保存全部条值。不要将本版b下方继续写成CaM遮蔽比例，不是同坐标删除CaM的ΔSASA。

来源：source_inputs/cam_all25.csv四列与原家族fold1源一致，25RMSF/21Hχ有效，原生294–318，A299/306/309/311 Hχ无定义。泛素源source_inputs/original_frozen_Full_A003R.npz与label NPZ，chain A72位，target[:,0/1]和mean[:,0/1]，valid[:,0/1]原mask；72RMSF/66Hχ有效，G10/A28/G35/A46/G47/G53 Hχ留空，P19/P37/P38的MD合法零及Full近零原值保留。case_profiles_all97_with_masks.csv保存97位置四项原读数及mask，standalone_all368_plotted_values.csv保存全部368有效曲线读数。模型来源为现成对应留出折单模型，CaMKI familyfold1与Ub原fold0不是同一套划分；详细哈希和MAE在case_profile_manifest.json。无新模型或MD计算。

准确度：CaMKI RMSF MAE0.1148592012Å、Hχ0.254417929857nats；泛素RMSF0.039814019369Å、Hχ0.079530506498nats。原v108主图713值中25个旧遮蔽值退出本版主图，余688值逐项保留；97个SASA值加入，总785值。f脚本、NPZ和计数CSV均与v108逐字节一致，保持无轮廓浅色透明叠加。

复现：运行case_profiles.py制作独立双案例曲线及SASA轨道，assemble_main.py生成figure3/compose_fig3_v109.py，后者绘制整图。所有冻结源和新几何输入/脚本包含于包内；已有缓存可直接画，不需重跑几何、模型或MD。PDF/SVG原生矢量，Arial嵌入，PNG450dpi。release_and_audit.py检查368曲线值、688旧值保留、97SASA条值、源哈希、f不变及字体并打包。
'''
(root/'README_双指标曲线与复现.md').write_text(readme,encoding='utf-8')
main='figure3_v109_20261027_双指标曲线与溶剂暴露';single='figure3_bc_v109_20261027_双指标曲线下加SASA'
for ext in ['png','pdf','svg']:
    shutil.copy2(root/'figure3'/('fig3_dual_case_entropy.'+ext),base/(main+'.'+ext))
    shutil.copy2(root/('case_RMSF_Hchi_dual_axis.'+ext),base/(single+'.'+ext))
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图3 v109：曲线下方增加SASA</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;margin:28px;color:#222;background:white}h1{font-size:24px}section{max-width:1150px;margin-bottom:30px}img{display:block;width:100%;height:auto}a{color:#2166ac}</style></head><body><h1>图3 v109：每个案例下方增加溶剂暴露柱图</h1><p>上方蓝色RMSF、赭色Hχ曲线；下方青灰色为起始复合物侧链SASA，残基逐位对齐。两例SASA统一PyMOL3.0、1.4Å探针、dot_density4；柱图不是模型预测或MD平均值。CaMKI原遮蔽比例轨道已换成SASA。Word未改。</p><section><a href="SINGLE"><img src="SINGLE" alt="曲线下方的侧链SASA柱图"></a><p><a href="SPDF">双案例曲线与SASA PDF</a></p></section><section><a href="MAIN"><img src="MAIN" alt="组图3"></a><p><a href="MPDF">组图3 PDF</a></p></section></body></html>'''
html=html.replace('SINGLE',single+'.png').replace('SPDF',single+'.pdf').replace('MAIN',main+'.png').replace('MPDF',main+'.pdf')
gallery=base/'figure3_v109_20261027_曲线与SASA预览.html';gallery.write_text(html,encoding='utf-8')
package=base/'figure3_v109_20261027_曲线与SASA_整图与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(root).as_posix())
with zipfile.ZipFile(package) as z:assert z.testzip() is None
print(json.dumps({'single':str(base/(single+'.png')),'main':str(base/(main+'.png')),'preview':str(gallery),'package':str(package),
                  'original_valid_case_values':368,'SASA_values':97,'main_values':785,'word_unchanged':True},ensure_ascii=False,indent=2))
