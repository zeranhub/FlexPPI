from pathlib import Path
from zipfile import ZipFile
from lxml import etree
import json
s=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');v=s/'revision_v87_general_story_compact';old=s/'revision_v86_dense_evidence'
plan=json.loads((v/'narrative/v87_paragraph_plan.json').read_text(encoding='utf-8'));repl=plan['paragraph_replacements']
background=plan['insert_before_original_paragraph_69'].replace('连续片段的运动范围与同一片段内部的位置差异','连续片段的主链波动与不同位置的侧链状态差异')
repl['69']='以钙调蛋白（CaM）与CaM依赖性蛋白激酶I（CaMKI）的识别复合物为例。CaM具有N端和C端两个叶瓣；CaMKI在未被激活时，可由自身调控片段限制激酶活性，该片段的CaM结合区与自抑制区部分重叠[29,30]。复合物1MXE显示，CaM的C端和N端叶瓣分别识别CaMKI294–318调控肽的W303与M316，形成跨越14个残基的双锚点包裹构型[29]（图3a、b）。已有功能实验发现，C端叶瓣能够紧密结合调控肽，激酶激活仍需要N端叶瓣参与[29]。因此，本例进一步考察：两个叶瓣对主链及侧链状态的局部约束，能否为维持自抑制解除后的调控构型提供动态解释？'
repl['72']='图3：主链与侧链柔性的互补解析及结构对应。a–e采用1MXE来源的结合态输入；Full/Intra为全组成链家族留出第1折的完整/链内模型，分别用蓝圆实线与暖赭方点虚线表示。（a）CaM–CaMKI整体结构，CaM叶瓣与连接段为灰色，调控肽为紫色，关键侧链为青色。（b）与a同旋转方向的W303–C端叶瓣、M316–N端叶瓣局部放大。（c）全部25个主链与21个有效侧链的RMSF/Hχ绝对预测值，序列与曲线共享坐标；浅底色与细条标出W303–M316核心，四个无χ角的Ala留空。（d）由全部25×144残基对距离汇总各肽位置接触的CaM残基数，按N端叶瓣、C端叶瓣及连接段堆叠，阈值为4.5 Å，含零接触位置；完整3600对矩阵随源数据提供。（e）全部21个有效侧链按类型排列的Hχ配对值及CaM静态遮蔽比例；遮蔽以相同肽坐标保留/去除CaM时的侧链可及面积计算，不能代替真实游离集合。（f）自抑制CaMKI结构1A06及M316–L45/S99同视角局部放大，虚线表示重原子距离而非氢键，与a、b交叉对照不同接触环境。1A06缺失317–318不补建；该静态对照不表示观测到状态转变，1MXE也不是完整激酶激活态。MD核查见图S3，方法与解释边界见补充材料D。'
oldplan=json.loads((old/'manuscript/v86_text_plan.json').read_text(encoding='utf-8'))['paragraph_replacements']
repl['161']=oldplan['161'].replace('图3d保留25个CaMKI位置与144个CaM位置的全部3600个组合，以各残基对间最小重原子距离区分≤4.0 Å、4.0–4.5 Å与>4.5 Å；','图3d以25个CaMKI位置与144个CaM位置的全部3600个最小重原子距离，统计每个肽位置在4.5 Å内接触的CaM残基数，按叶瓣与连接段汇总并保留零接触位置；完整距离矩阵及≤4.0 Å、4.0–4.5 Å分级源图随源数据提供。')
repl['162']=oldplan['162'].replace('图3f的调控片段回接触和激活效应依据','正文的调控片段回接触与激活假说依据')
repl['175']=json.loads((s/'revision_v85_1MXE/manuscript/v85_text_plan.json').read_text(encoding='utf-8'))['paragraph_replacements']['175'].replace('图3f的激活及回接触假说','3.3节的激活及回接触假说')
(v/'narrative/paragraph_replacements.json').write_text(json.dumps({'paragraph_replacements':repl,'insert_before_original_paragraph_69':background},ensure_ascii=False,indent=2),encoding='utf-8')
code=(old/'build_v86_manuscript.py').read_text(encoding='utf-8')
code=code.replace("revision_v86_dense_evidence'","revision_v87_general_story_compact'")
code=code.replace("SOURCE=R/'paper_v85_20261026_中文修订_1MXE机制假说与图3重绘.docx';OUTPUT=R/'paper_v86_20261026_中文修订_图3完整接触与侧链环境.docx'","SOURCE=R/'paper_v86_20261026_中文修订_图3完整接触与侧链环境.docx';OUTPUT=R/'paper_v87_20261026_中文修订_互补柔性背景与图3凝练.docx'")
code=code.replace("fig3_1MXE_dense_evidence.png","fig3_joint_flexibility_compact.png")
code=code.replace("Figure 3 1MXE complete contacts and sidechain environments","Figure 3 Complementary backbone and sidechain flexibility")
code=code.replace("True 1MXE and 1A06 structures, complete contact map, all valid sidechains and absolute model profiles","True 1MXE and 1A06 structural contacts, all25 lobe counts, all21 sidechains, absolute profiles")
needle="parts['word/document.xml']=etree.tostring(root,encoding='UTF-8',xml_declaration=True,standalone=True)"
ins="""background=json.loads((V/'narrative/paragraph_replacements.json').read_text(encoding='utf-8'))['insert_before_original_paragraph_69']
intro=deepcopy(pars[69])
for key in list(intro.attrib):
 if key.endswith('paraId') or key.endswith('textId'):del intro.attrib[key]
set_text(intro,background);body.insert(body.index(pars[69]),intro)
assert len(body.findall('w:p',NS))==186
"""
code=code.replace(needle,ins+needle)
code=code.replace("finaltext={str(i):''.join(p.xpath('.//w:t/text()',namespaces=NS)) for i,p in enumerate(pars)}","finaltext={str(i):''.join(p.xpath('.//w:t/text()',namespaces=NS)) for i,p in enumerate(pars)}")
code=code.replace("v86_text_plan.json","v87_text_plan.json").replace("v86_final_text.txt","v87_final_text.txt").replace("v86_build_audit.json","v87_build_audit.json")
code=code.replace("'paragraph_replacements':repl}","'paragraph_replacements':repl,'insert_before_original_paragraph_69':background}")
code=code.replace("'\\n'.join(finaltext.values())","'\\n'.join(finaltext[str(i)]+('\\n'+background if i==68 else '') for i in range(len(pars)))")
code=code.replace("'\\n\\n'.join(finaltext[str(i)] for i in [68,69,70,73,74,75])","'\\n\\n'.join([finaltext['68'],background]+[finaltext[str(i)] for i in [69,70,73,74,75]])")
code=code.replace("all_other_paragraphs_identical=True","all_other_original_paragraphs_identical=True,inserted_background_before_original_paragraph=69")
code=code.replace("revision_v85_1MXE/render","revision_v86_dense_evidence/render")
(v/'build_v87_manuscript.py').write_text(code,encoding='utf-8')
print('Prepared final text and OOXML builder; no Word file written yet.')