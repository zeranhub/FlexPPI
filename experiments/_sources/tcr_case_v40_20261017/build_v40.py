"""Integrate a validated same-pMHC TCR comparison into the Chinese manuscript."""
from pathlib import Path
from copy import deepcopy
import json,hashlib,re,zipfile
from collections import Counter
from docx import Document
from docx.shared import Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.text.paragraph import Paragraph
from docx.oxml import OxmlElement
from lxml import etree
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v39_20261017_中文修订_读者表述与英文图表_三线表.docx'
OUT=BASE/'paper_v40_20261017_中文修订_TCR识别环境变化实例_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='61c88902753e1b0a826f2f3750a3dbdb40e1e57a8dbf08c5d60fa8d559721f48'
d=Document(SOURCE);paras=list(d.paragraphs);changes=[]
def settext(p,t):
    assert not p._p.xpath('.//m:oMath')
    old=p.text;nodes=p._p.xpath('.//w:t');assert nodes
    nodes[0].text=t
    for n in nodes[1:]:n.text=''
    changes.append({'before':old,'after':t})
def insert_before(anchor,text,template):
    el=OxmlElement('w:p');anchor._p.addprevious(el);p=Paragraph(el,anchor._parent)
    if template._p.pPr is not None:el.append(deepcopy(template._p.pPr))
    r=p.add_run(text)
    if template.runs and template.runs[0]._r.rPr is not None:r._r.insert(0,deepcopy(template.runs[0]._r.rPr))
    p.paragraph_format.page_break_before=False
    return p
body=paras[94];head=paras[91];caption=paras[102]
# Shift only the existing final Results section and its figure references.
for p in paras:
    for n in p._p.xpath('.//w:t'):
        if n.text:n.text=n.text.replace('图6','图7')
    if p.text.startswith(('3.5 ','3.5.')):settext(p,re.sub(r'^3\.5','3.6',p.text))
settext(paras[2],paras[2].text.replace('KIX与泛素–UBA实例分别检验活动区域排序和同类侧链状态差异。','KIX与泛素–UBA实例检验活动区域排序和同类侧链状态差异；同一肽–MHC被SB27和SB47识别时，肽的主链波动变化与MD的排序相关为0.945，中央Q7的侧链联合熵增幅也与MD一致。'))
settext(paras[9],paras[9].text.replace('KIX和泛素–UBA实例分别对照活动区域定位与同类型侧链差异。','KIX和泛素–UBA实例对照活动区域定位与同类型侧链差异；ITrace中的SB27/SB47比较进一步检验相同肽–MHC在不同TCR识别环境下的动态变化。'))
settext(paras[49],paras[49].text.replace('结构实例分别考察主链活动区域的定位与同类型残基侧链朝向多样性的区分。','结构实例考察主链活动区域定位、TCR识别环境改变下的动态响应，以及同类型残基侧链朝向多样性的区分。'))
method1='在244条ITrace记录中，按肽与MHC序列匹配得到124对具有不同TCR的结构组合，比较同一肽位置的预测变化与MD变化。SB27/SB47兼具已报道的识别方式差异和清晰的局部预测响应，用于展开分析[28]。SB27（2AK4）与SB47（4JRY）的模型输入具有相同的13肽LPEPLPQGQLTAY和MHC序列；预测采用2.5节全组成链家族划分的第0折模型，该折训练来自ATLAS与DynaRepo。每种复合物使用ITrace的三条200 ns重复轨迹，保留0–200 ns的全部1,001–1,003个保存帧，按2.2节重新计算肽链RMSF与Hχ并核对原标签。'
method2='TCR实例中的变化定义为SB47减SB27：在13个肽位置比较RMSF变化，在全部11个具有χ角的位置比较Hχ变化。MD参考为各组三条重复轨迹的均值；均值差标准误由两组重复间的样本方差分别除以三、相加后开方得到。为量化TCR邻域对预测的贡献，另从每个输入中移除TCR的D、E链，保留肽、MHC及β2微球蛋白的全部原坐标，使用同一模型重新预测。'
insert_before(paras[51],method1,body);insert_before(paras[51],method2,body)
result=[
'3.5 相同肽–MHC在不同TCR识别环境下的动态变化',
'SB27和SB47识别同一条HLA-B35呈递的病毒13肽，但结合方式明显不同：SB27跨越肽中央的隆起区，SB47的识别位置则向肽N端偏移[28]。两种输入中肽与MHC的序列相同，因此，这一比较检验模型能否把不同的TCR结合环境对应到同一条肽的局部动态变化（图6a）。',
'在全部13个肽位置上，模型与MD的RMSF变化排序相关为0.945，变化量的平均绝对误差为0.093 Å，13个位点的变化方向均与三重复MD的均值同向（图6c）。模型识别了中央P6–L10区域波动增强与N端L1波动降低的不同响应。其中Q7的MD RMSF由0.392升至0.742 Å，模型由0.850升至1.449 Å。主链输出既定位了主要变化区域，也重现了沿肽序列的增减分布。',
'侧链变化在Q7上提供了第二项相符的动态读数。该位点的MD联合Hχ由1.683升至2.757 nats，模型由1.726升至2.825 nats，预测增量1.099 nats与MD增量1.074 nats接近（图6b、d）。SB47的三条重复轨迹均具有更高的Q7联合熵，其范围为2.621–2.928 nats，高于SB27的1.571–1.883 nats。全肽11个侧链有效位置的对照同时显示，Q9的熵增幅未被模型重现（图6d）；该实例的侧链变化证据主要来自Q7。',
'保持各自结合构象中的肽–MHC坐标并移除TCR后，模型对Q7的Hχ预测分别为2.680与2.663 nats，两种结构的预测差由1.099变为−0.018 nats（图6b）。在相同肽与MHC化学组成下，完整TCR邻域使模型区分了Q7的侧链状态多样性。结合全肽主链变化及三条MD重复，这一实例支持模型将识别方式差异转换为具体肽位点的动态响应。'
]
for i,t in enumerate(result):insert_before(paras[91],t,head if i==0 else body)
figp=insert_before(paras[91],'',body)
figp.alignment=WD_ALIGN_PARAGRAPH.CENTER
figp.paragraph_format.keep_with_next=True
figp.add_run().add_picture(str(ROOT/'figure6/fig6_tcr_case.png'),width=Inches(6.4))
cap='图6：不同TCR识别环境下同一肽–MHC的动态变化预测。（a）SB27（2AK4）与SB47（4JRY）的输入结构，按MHC肽结合区Cα对齐并采用相同视角与比例；展示MHC结合区、TCR可变区与完整13肽，Q7为红色。（b）Q7联合Hχ：空心圆为三条MD重复，灰色点线连接重复均值；蓝色为完整复合物预测，橙色为保留各自肽–MHC坐标、移除TCR后的同一模型预测。（c）全部13个肽位置的RMSF变化。（d）全部11个侧链有效位置的Hχ变化，G8与A12无χ角，灰色叉号仅标记无值位置。c、d均取SB47减SB27，MD点为三重复均值差，误差线为均值差标准误；蓝色方点为模型预测差。红色底纹标出Q7。'
newcap=insert_before(paras[91],cap,caption);newcap.paragraph_format.keep_with_next=False
(ROOT/'figure6/caption_zh.md').write_text(cap,encoding='utf8')
settext(paras[106],paras[106].text.replace('KIX中的共同高波动环区和泛素中的同类型位点分别给出具体验证；','KIX中的共同高波动环区、SB27/SB47识别环境改变后的全肽主链响应及Q7侧链熵增幅，以及泛素中的同类型位点给出具体验证；'))
settext(paras[108],paras[108].text.replace('KIX与泛素分别展示活动区域定位和同类型侧链多样性区分。','KIX与泛素展示活动区域定位和同类型侧链多样性区分；SB27/SB47实例进一步重现相同肽–MHC在不同TCR环境下的主链变化与Q7侧链熵升高。'))
lastref=next(p for p in paras if p.text.startswith('[27]'))
ref='[28] Yu Chih Liu, John J. Miles, Michelle A. Neller, Emma Gostick, David A. Price, Anthony W. Purcell, James McCluskey, Scott R. Burrows, Jamie Rossjohn, and Stephanie Gras. Highly divergent T-cell receptor binding modes underlie specific recognition of a bulged viral peptide bound to a human leukocyte antigen class I molecule. Journal of Biological Chemistry, 288(22):15442–15454, 2013. doi:10.1074/jbc.M112.447185.'
insert_before(lastref._p.getnext() is not None and Paragraph(lastref._p.getnext(),lastref._parent) or paras[-1],ref,lastref)
# Paragraphs inserted before existing section: keep body together only as required by headings.
d.save(OUT)
alltext='\n'.join(p.text for p in d.paragraphs)
(QA/'v40_content.txt').write_text('\n'.join(f'P{i}: {p.text}' for i,p in enumerate(d.paragraphs)),encoding='utf8')
assert len(d.tables)==10 and len([p for p in d.paragraphs if re.match(r'^图[1-7]：',p.text)])==7
assert len([p for p in d.paragraphs if re.match(r'^\[\d+\]',p.text)])==28
assert not any(x in alltext for x in ['冻结','发布权重','发布模型'])
audit={'status':'built_pending_render_QA','source':str(SOURCE),'source_sha256':sha(SOURCE),'document':str(OUT),'document_sha256':sha(OUT),'new_experiment':'SB27 (2AK4) vs SB47 (4JRY), same peptide/MHC sequence, three200ns ITrace replicas per structure','new_figure':str(ROOT/'figure6/fig6_tcr_case.png'),'changes':changes,'inserted_methods':[method1,method2],'inserted_results':result,'caption':cap,'reference':ref,'existing_tables':10,'new_training':0,'new_MD':0,'new_same_model_inferences':4}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(str(OUT));print(sha(OUT))
