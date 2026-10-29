from pathlib import Path
from hashlib import sha256
from PIL import Image
from pypdf import PdfReader
import json,shutil,zipfile

root=Path(__file__).resolve().parent;base=root.parents[2]
digest=lambda p:sha256(p.read_bytes()).hexdigest()
plan=json.loads((root/'analysis_plan_before_results.json').read_text(encoding='utf-8'))
for filename,expected in plan['inputs'].items():assert digest(Path(filename))==expected
assert digest(base/'paper_v98_20261027_中文修订_结构区域放大.docx')=='eb1bc676be591345f29ab45c2b4f8f8ef9b94ea2cba6c64d74d0781d41ad4e92'
fonts=[]
for ref in PdfReader(root/'f_direct_matched_Hchi_concept.pdf').pages[0]['/Resources']['/Font'].values():
    fonts.append(str(ref.get_object().get('/BaseFont')))
assert all('Arial' in f for f in fonts)
with Image.open(root/'f_direct_matched_Hchi_concept.png') as im:im.verify()
qa=json.loads((root/'plot_QA.json').read_text(encoding='utf-8'));qa.update(visual_review='passed: root opened PNG; all record lines retained, axes/readouts explicit, titles separated from mean annotations',independent_scientific_audit='passed: deterministic matching, no outcome selection, no repeated residue, equal AA then record weighting; one group hand-verified',embedded_Arial_fonts=fonts)
(root/'plot_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
report='''# 新f方案：相近主链运动，仍可区分侧链占据

建议将原“RMSF与Hχ两个相关系数的散点比较”改成相近主链运动下的直接MD对照。图的两个横轴都是预测Hχ较低/较高；左图纵轴是MD参考主链RMSF，右图是MD单角细箱占据熵。每条灰线是一条记录内先对氨基酸类型等权后的结果，蓝线为83条记录等权均值及来源分层bootstrap95%区间。左图几乎水平、右图整体上升，直接回答模型识别了什么。两轴单位不同，不比较两条线的视觉坡度或高度大小。

本次是新增的事后分析设计，使用冻结的原模型/MD数值；未训练模型、未推理、未运行新MD。匹配主阈值0.10Å在计算匹配结果之前写入analysis_plan_before_results.json；同时预先固定0.05/0.20Å敏感性和预测RMSF匹配检查，没有按结果选择最佳阈值。

匹配规则：同一记录、同一种氨基酸内，按MD参考RMSF和residue_index排序，从小到大以相邻残基作不重叠匹配。若两者RMSF差超过阈值，跳过较小者再试；满足则两残基各使用一次。配对完成之后，仅用冻结Full Hχ预测决定低/高方向，不设Hχ差异阈值，不按MD侧链结果挑成功样本。主分析有13,968对，共27,936个位点，输入33,493个位点中5,557个没有进入配对；没有整条记录被排除，83条均有配对，最少5对、最多336对。来自ATLAS30、DynaRepo23、iTrace30，记录并非83个独立家族。旧相关分析要求同类型至少4个位点，新配对只需至少2个，因此不是原33,142位点集合。

每个氨基酸组先对所有配对取均值，再对记录内氨基酸组等权，最后对83条记录等权。没有把13,968个相关配对当成独立实验做置信区间。区间为来源分层的4,000次整记录bootstrap，seed20261027，同一批索引用于低、高和差值；不是按家族分层的泛化区间。

主结果：低/高预测Hχ组的MD参考RMSF均值1.115768/1.118007Å，差0.002239Å，95%区间0.000450–0.004059Å；两组近似相同，但不能声称完全无差异。MD侧链细箱占据熵2.530946/2.688912nats/χ，差0.157966，95%区间0.143817–0.172513，83条中82条为正。低F*角宽36.4413/40.4155°/χ，差3.97423°，95%区间3.40266–4.52747，83条中78条为正；此图单独保存为角范围验证。

敏感性：RMSF阈值0.05/0.20Å时，侧链熵均差0.158950/0.152660，所有83条记录仍保留。改为匹配模型预测RMSF且阈值0.10Å时，MD侧链熵差0.161620，95%区间0.148001–0.174955；MD参考RMSF差−0.005900Å，95%区间−0.027829–0.012498。所有预先规定的6组检查完整在summary_all_prespecified_analyses.csv，不仅报告正面主阈值。

需要保留的边界：original_RMSF取自冻结模型target，侧链占据来自replica1缓存；不能写成完全同一副本的严格匹配实验。参考RMSF只是近似匹配，不能声称排除了所有主链影响；这里支持的是条件接近时的区分能力，不证明严格信息独立、因果机制、完全不可替代或跨家族泛化。若要证明在预测任务上RMSF无法替代Hχ，还需预先规定的RMSF-only与RMSF+Hχ交叉验证增益。

建议正文句：在同一记录内、主链参考RMSF相近的同类残基之间，较高的模型Hχ仍对应较高的MD侧链角度占据熵，提示模型能够识别仅由主链运动幅度难以直接看出的侧链多样性差异。

本次只交付新f的实际数据原型与分析包，未替换v100整图或Word。字体、线重、颜色沿用论文图2。
'''
(root/'README_新f呈现与统计边界.md').write_text(report,encoding='utf-8')
stem='figure3_f_v101_20261027_相近主链_侧链差异'
for ext in ['png','pdf','svg']:shutil.copy2(root/f'f_direct_matched_Hchi_concept.{ext}',base/f'{stem}.{ext}')
suppstem='figure3_f_v101_20261027_角度范围验证'
for ext in ['png','pdf','svg']:shutil.copy2(root/f'f_direct_matched_angular_range.{ext}',base/f'{suppstem}.{ext}')
manifest={'main':str(base/(stem+'.png')),'reference_RMSF_tolerance_A':.1,'matched_pairs':13968,'records':83,'existing_main_and_word_unchanged':True,
          'outputs':[{'file':p.name,'sha256':digest(p)} for p in sorted(root.iterdir()) if p.is_file() and p.name!='release_direct_concept.py']}
(root/'release_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
archivepath=base/'figure3_f_v101_20261027_相近主链匹配分析与源图.zip'
with zipfile.ZipFile(archivepath,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
    for p in sorted(root.iterdir()):
        if p.is_file():archive.write(p,p.name)
with zipfile.ZipFile(archivepath) as archive:assert archive.testzip() is None
print(json.dumps({'image':str(base/(stem+'.png')),'pdf':str(base/(stem+'.pdf')),'width_validation':str(base/(suppstem+'.png')),'archive':str(archivepath),'bytes':archivepath.stat().st_size,'word_unchanged':True},ensure_ascii=False,indent=2))
