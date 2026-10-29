from pathlib import Path
from hashlib import sha256
from PIL import Image
from pypdf import PdfReader
import json,shutil,zipfile
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parent;base=root.parents[2];old=root.parent/'revision_v100_evidence_redesign'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
stats=json.loads((root/'direct_relationship_statistics.json').read_text(encoding='utf-8'))
assert digest(Path(stats['source']))==stats['source_sha256']=='293cef1437db032d2dde6cb176869e43f8ee4af27a7358cab20242e2da94b96d'
assert stats['common_valid_positions']==824988 and stats['records']==3063
for folder in ['figure3/source_data','pymol']:
    for p in (old/folder).rglob('*'):
        if p.is_file():assert digest(p)==digest(root/p.relative_to(old))
old_values=pd.read_csv(old/'figure3'/'main_plotted_values_v100.csv',float_precision='round_trip')
new_values=pd.read_csv(root/'figure3'/'main_original_panels_a_to_e_values.csv',float_precision='round_trip')
pd.testing.assert_frame_equal(old_values[old_values.panel.ne('f')].reset_index(drop=True).dropna(axis=1,how='all'),new_values.dropna(axis=1,how='all'),check_exact=True)
assert digest(base/'paper_v98_20261027_中文修订_结构区域放大.docx')=='eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
checks=[]
for stem in ['RMSF_Hchi_direct_density','RMSF_Hchi_direct_density_full_range','figure3/fig3_dual_case_entropy']:
    p=root/(stem+'.pdf');page=PdfReader(p).pages[0];fonts=[]
    for ref in page['/Resources']['/Font'].values():
        font=ref.get_object();fonts.append(str(font.get('/BaseFont')))
        for child in font.get('/DescendantFonts',[]):
            d=child.get_object()['/FontDescriptor'].get_object();assert any(k in d for k in ['/FontFile','/FontFile2','/FontFile3'])
    assert all('Arial' in f for f in fonts)
    with Image.open(root/(stem+'.png')) as im:im.verify()
    q=root/('figure3/figure_QA.json' if stem.startswith('figure3') else stem+'_QA.json')
    qa=json.loads(q.read_text(encoding='utf-8'));qa.update(visual_review='passed: root and independent reviewer opened density views; main root reviewed a–f; no clipped text or collision',embedded_Arial_fonts=fonts,all_values_exact_from_frozen_OOF=True,word_unchanged=True)
    q.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');checks.append(stem)

caption='''f，主链RMSF与联合侧链Hχ的直接关系。左为MD参考，右为对应Full折外预测；使用相同824,988个共同有效残基位置，覆盖3,063条记录。深色表示六边形格中的位置数较多，两图共用对数计数色标、相同横纵轴和网格；上方与右侧为相应边缘频率分布。图示主视野为RMSF0–10Å；MD/预测分别有4,234（0.513%）/906（0.110%）个位点位于右侧范围外，完整0–35Å图随补充稿提供。相关系数始终由全部824,988个位点计算，分别为0.209与0.262，未按显示范围筛选。有效Pro及其合法零Hχ保留，Ala/Gly按原χ标签无定义掩码排除。横带与氨基酸类型、侧链χ数量及固定区间状态定义有关，因此整体弱相关不等同于两指标相互独立。

建议正文：MD中主链位置波动与侧链联合状态多样性具有有限的整体关联，且同一RMSF范围内可出现宽范围Hχ值；二者并非一一对应。Full预测保留了这种总体分布特征，但整体耦合程度略强于MD。该图用于说明两类动态摘要的关系，不据此声称严格信息独立、完全不可替代或恢复了精确耦合机制。
'''
(root/'figure3_f_caption_and_claim_zh.md').write_text(caption,encoding='utf-8')
readme='''# v103：回到直接RMSF–Hχ关系

按用户旧图重画左右六边形密度+上/右边缘分布：MD RMSF–MD联合Hχ，对应Full预测RMSF–预测联合Hχ。主图a–e保持v100的真实数据、结构、样式和e四图横排；f只替换为上述直接关系。v101/v102的复杂配对分析未并入本版主图。

数据来自现有五折OOF逐残基表，仅variant=full，两目标mask均有效、四个坐标均有限。MD/Full两图使用完全相同824,988位置、3,063记录；原full表961,321行，其中136,333因共同mask排除。A/G无定义标签排除，44,822个有效Pro保留，其中MD Hχ零值44,673个。原表及共同mask对应信息由same_residue_MD_and_Full_all824988.npz和cohort_all3063_records.csv保存，四个原始浮点坐标与原CSV逐项相同。

统计：全部位置的Spearman MD0.2086151551、Full0.2619048886；是按残基池化的描述值，较长记录权重较高，不是独立残基实验检验。没有新增P值或回归。按氨基酸类型中心化后rho分别0.3451963753/0.5051875285，同记录同类型中心化0.3613679095/0.4886433120，仅用于解释边界，未塞进主图。横带包含残基类型与χ数量先验，不能将弱总体rho说成信息独立。模型的类型内耦合更强，不宜称精确恢复耦合强度。

显示：主视野沿用户参考0–10Å。MD超界4,234/824,988，Full906/824,988，统计仍按全量算。完整0–35Å图双方都显示全部位置。密度计数逐格验算，双方色标一致，对数显示格内计数；边缘直方图用全部824,988共同分母，边缘坐标范围一致。没有平滑或抽样。完整区间图作为长尾核验，不能将主视野图误写成展示全部RMSF范围。

文件：plot_direct_density.py生成独立双图并提供draw_pair；assemble_main.py生成调用原a–e生产逻辑的新组图脚本；figure3/compose_fig3_v103.py可重建整图。NPZ与统计JSON已有缓存可复现绘图，无需新模型推理、训练或MD。PDF/SVG是原生矢量，Arial嵌入/文本可编辑。Word未修改。
'''
(root/'README_直接关系与数值说明.md').write_text(readme,encoding='utf-8')
mainstem='figure3_v103_20261027_直接RMSF_Hchi关系'
pairstem='figure3_f_v103_20261027_RMSF_Hchi直接关系'
fullstem='figure3_f_v103_20261027_RMSF_Hchi完整范围'
for ext in ['png','pdf','svg']:
    shutil.copy2(root/'figure3'/f'fig3_dual_case_entropy.{ext}',base/(mainstem+'.'+ext))
    shutil.copy2(root/f'RMSF_Hchi_direct_density.{ext}',base/(pairstem+'.'+ext))
    shutil.copy2(root/f'RMSF_Hchi_direct_density_full_range.{ext}',base/(fullstem+'.'+ext))
html='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>图3直接关系版</title><style>body{font:17px/1.7 Arial,"Microsoft YaHei",sans-serif;color:#151515;background:white;margin:28px}h1{font-size:24px}img{display:block;width:100%;height:auto}section{max-width:1250px;margin-bottom:32px}a{color:#2166ac}</style></head><body><h1>图3 v103：RMSF–Hχ直接关系</h1><p>两图采用相同824,988个共同有效OOF位置；相关系数按全部位置计算。主视野沿旧参考为0–10Å，完整范围图单独提供。复杂匹配方案未放进主图，e仍为四图横排，Word未替换。</p><section><a href="PAIR"><img src="PAIR" alt="RMSF与Hchi直接关系"></a><p><a href="FULL">查看完整0–35Å范围</a> · <a href="PDF">下载整图PDF</a></p></section><section><a href="MAIN"><img src="MAIN" alt="整张图3"></a></section></body></html>'''.replace('PAIR',pairstem+'.png').replace('FULL',fullstem+'.png').replace('PDF',mainstem+'.pdf').replace('MAIN',mainstem+'.png')
gallery=base/'figure3_v103_20261027_直接关系预览.html';gallery.write_text(html,encoding='utf-8')
package=base/'figure3_v103_20261027_直接关系_整图与全量源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:archive.write(p,p.relative_to(root).as_posix())
with zipfile.ZipFile(package) as archive:assert archive.testzip() is None
print(json.dumps({'main':str(base/(mainstem+'.png')),'pair':str(base/(pairstem+'.png')),'full_range':str(base/(fullstem+'.png')),'pdf':str(base/(mainstem+'.pdf')),'gallery':str(gallery),'package':str(package),'package_bytes':package.stat().st_size,'word_unchanged':True},ensure_ascii=False,indent=2))
