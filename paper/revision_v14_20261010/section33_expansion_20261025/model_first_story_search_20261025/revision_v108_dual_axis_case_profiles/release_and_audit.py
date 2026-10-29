from pathlib import Path
from hashlib import sha256
import json,shutil,zipfile
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

root=Path(__file__).resolve().parent;base=root.parents[2];old=root.parent/'revision_v107_soft_alpha_overlay'
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
assert len(previous)==527 and len(current)==713
keys=[c for c in previous.columns if c!='value']
matched=previous.merge(current,on=keys,suffixes=('_old','_new'),validate='one_to_one')
assert len(matched)==527
precision_delta=float(abs(matched.value_old-matched.value_new).max());assert precision_delta<=1e-14
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

caption='''b、c，CaMKI调控片段（1MXE，原生294–318）与泛素（2QHO，链A1–72）的主链RMSF和联合侧链Hχ沿序列分布。共用残基横轴，蓝色左轴表示RMSF（Å），赭色右轴表示Hχ（nats）；实线为Full预测，虚线为MD参考。两案例均使用RMSF0–2.4Å、Hχ0–4.5nats固定范围，两纵轴都从零开始。无Hχ定义的残基两来源均留空，不跨缺口插值；小点确保有效孤立位置也可见，合法Pro零Hχ保留。CaMKI的RMSF/Hχ分别有25/21个有效位置，泛素分别72/66个。Full的RMSF/Hχ MAE分别为CaMKI0.115Å/0.254nats、泛素0.040Å/0.080nats。CaMKI使用现成家族留出fold1，泛素使用原chemical_identity_corrected_v1记录留出fold0，均为seed7、epoch40对应留出折单模型，未新增训练、预测或MD。b下方CaM遮蔽比例仍为起始结构同坐标删除CaM的静态SASA几何指标。颜色和线型分别编码指标与来源；左右轴单位不同，曲线交叉不代表数值相等，峰谷对应也不用于作因果解释。
'''
(root/'caption_dual_axis_profiles_zh.md').write_text(caption,encoding='utf-8')
readme='''# v108：每个案例同时展示主链RMSF与侧链Hχ

主图b、c改为同残基横轴的双纵轴曲线：蓝色左轴RMSF（Å）、赭色右轴Hχ（nats），与图2a指标色一致；Full实线、MD虚线。两例都用RMSF0–2.4Å/Hχ0–4.5nats，从0起，不通过单独缩放制造对应趋势。Hχ小点保留无定义位置两侧的孤立有效点，不跨缺口连线。两指标单位不同，曲线交叉或峰谷接近不表示相等、因果耦合或数值可互换。

来源：source_inputs/cam_all25.csv四列与原家族fold1源一致，25RMSF/21Hχ有效，原生294–318，A299/306/309/311 Hχ无定义。泛素源source_inputs/original_frozen_Full_A003R.npz与label NPZ，chain A72位，target[:,0/1]和mean[:,0/1]，valid[:,0/1]原mask；72RMSF/66Hχ有效，G10/A28/G35/A46/G47/G53 Hχ留空，P19/P37/P38的MD合法零及Full近零原值保留。case_profiles_all97_with_masks.csv保存97位置四项原读数及mask，standalone_all368_plotted_values.csv保存全部368有效曲线读数。模型来源为现成对应留出折单模型，CaMKI familyfold1与Ub原fold0不是同一套划分；详细哈希和MAE在case_profile_manifest.json。无新模型或MD计算。

准确度：CaMKI RMSF MAE0.1148592012Å、Hχ0.254417929857nats；泛素RMSF0.039814019369Å、Hχ0.079530506498nats。CaMKI Hχ的准确度应使用这里21有效位置的值，不能沿用不同范围旧值。主图原527项科学读数保留，新增42个CaMKI Hχ值及144个泛素RMSF值，共713项；旧导出CSV在少数泛素Hχ末位存在最大约4.44e-16读写舍入，本次严格恢复原NPZ精度。a分子结构、d同类侧链环境、e四χ横排、f无轮廓浅色叠加均保留。f脚本、NPZ和计数CSV均与v107逐字节一致。

复现：运行case_profiles.py制作独立双案例曲线，assemble_main.py生成figure3/compose_fig3_v108.py，后者绘制整图。冻结源均包含于包内；已有case_profiles_all97缓存可直接重画，f已有全量NPZ可复现。PDF/SVG原生矢量，Arial嵌入，PNG450dpi。Word未改。release_and_audit.py检查全部368曲线值、527旧值保留、源哈希、Pro/mask、f不变和字体，并打包。
'''
(root/'README_双指标曲线与复现.md').write_text(readme,encoding='utf-8')
main='figure3_v108_20261027_双指标案例曲线';single='figure3_bc_v108_20261027_每例RMSF与Hchi同图'
for ext in ['png','pdf','svg']:
    shutil.copy2(root/'figure3'/('fig3_dual_case_entropy.'+ext),base/(main+'.'+ext))
    shutil.copy2(root/('case_RMSF_Hchi_dual_axis.'+ext),base/(single+'.'+ext))
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图3 v108：每例双指标曲线</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;margin:28px;color:#222;background:white}h1{font-size:24px}section{max-width:1150px;margin-bottom:30px}img{display:block;width:100%;height:auto}a{color:#2166ac}</style></head><body><h1>图3 v108：每个案例同时看RMSF与Hχ</h1><p>蓝色左轴RMSF，赭色右轴Hχ；实线Full，虚线MD。两个案例使用相同纵轴范围，保留原始无定义位置的断线和有效孤立点。组图b、c已替换，f保留无轮廓透明叠加，Word未修改。</p><section><a href="SINGLE"><img src="SINGLE" alt="每个案例的双指标曲线"></a><p><a href="SPDF">双案例曲线PDF</a></p></section><section><a href="MAIN"><img src="MAIN" alt="组图3"></a><p><a href="MPDF">组图3 PDF</a></p></section></body></html>'''
html=html.replace('SINGLE',single+'.png').replace('SPDF',single+'.pdf').replace('MAIN',main+'.png').replace('MPDF',main+'.pdf')
gallery=base/'figure3_v108_20261027_双指标曲线预览.html';gallery.write_text(html,encoding='utf-8')
package=base/'figure3_v108_20261027_双指标曲线_整图与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(root).as_posix())
with zipfile.ZipFile(package) as z:assert z.testzip() is None
print(json.dumps({'single':str(base/(single+'.png')),'main':str(base/(main+'.png')),'preview':str(gallery),'package':str(package),
                  'original_valid_case_values':368,'main_values':713,'word_unchanged':True},ensure_ascii=False,indent=2))
