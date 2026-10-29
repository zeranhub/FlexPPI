from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree
from copy import deepcopy
from PIL import Image
import json,hashlib,shutil
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');V=R/'section33_expansion_20261025/model_first_story_search_20261025/revision_v85_1MXE';M=V/'manuscript';M.mkdir(exist_ok=True)
SOURCE=R/'paper_v83_20261025_中文修订_图3模型柔性再分配与关键区域.docx';OUTPUT=R/'paper_v85_20261026_中文修订_1MXE机制假说与图3重绘.docx'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','pic':'http://schemas.openxmlformats.org/drawingml/2006/picture','pr':'http://schemas.openxmlformats.org/package/2006/relationships'}
def q(prefix,name):return '{'+NS[prefix]+'}'+name
def serialize(x):return etree.tostring(x,encoding='UTF-8',xml_declaration=True,standalone=True)
def signature(x):return etree.tostring(x,method='c14n')
def set_text(p,text):
 old=p.findall('w:r',NS);style=deepcopy(old[0].find('w:rPr',NS)) if old and old[0].find('w:rPr',NS) is not None else None
 for node in list(p):
  if node.tag!=q('w','pPr'):p.remove(node)
 r=etree.SubElement(p,q('w','r'))
 if style is not None:r.append(style)
 t=etree.SubElement(r,q('w','t'));t.text=text
 if text.startswith(' ') or text.endswith(' '):t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
def set_setting(p,name,value):
 pr=p.find('w:pPr',NS)
 if pr is None:pr=etree.Element(q('w','pPr'));p.insert(0,pr)
 n=pr.find('w:'+name,NS)
 if n is None:n=etree.SubElement(pr,q('w',name))
 n.set(q('w','val'),str(value))
with ZipFile(SOURCE) as z:parts={i.filename:z.read(i.filename) for i in z.infolist()}
root=etree.fromstring(parts['word/document.xml']);body=root.find('w:body',NS);pars=body.findall('w:p',NS);assert len(pars)==185
before=[signature(p) for p in pars];tables=[signature(t) for t in body.findall('w:tbl',NS)]
repl=json.loads((V/'narrative/paragraph_replacements_draft.json').read_text(encoding='utf-8'))['paragraph_replacements']
repl['69']=repl['69'].replace('与激酶自身的接触[30]','与激酶自身的接触[29,30]')
repl['161']=repl['161'].replace('作为完整互补集合','作为完整对照集合').replace('。 起始结构','。起始结构')
repl['162']+=' 完整体系145个有效侧链的Hχ MAE为Full 0.193119、Intra 0.228446 nats；其原值Spearman相关分别为0.928208和0.932082，原值排序并非所有统计口径均由Full占优。输入结构及对应提供的蛋白轨迹不含Ca离子，本例不用于量化钙效应。'
repl['72']=repl['72'].replace('蓝色与暖赭色在全部定量面板中保持相同含义','蓝色圆点/实线与暖赭方点/虚线在全部定量面板中保持相同含义').replace('与a相邻的局部结构放大','与a相邻、保持同一旋转方向的局部结构放大').replace('CaM为灰色、CaMKI调控肽为紫色','CaM两个叶瓣为深浅灰色、CaMKI调控肽为紫色，关键侧链为青色').replace('假说中的回接触与激活效应尚未由本例模型或MD直接测量','虚线箭头为提出的后果，回接触与激活效应尚未由本例模型或MD直接测量')
repl['175']=repl['175'].replace('MAE与Spearman相关ρ按相同有效位置计算','全体系RMSF的Full/Intra MAE为0.189/0.794 Å、Spearman相关为0.776/0.299，Hχ的MAE为0.193/0.228 nats、相关为0.928/0.932').replace('（a）CaMKI','黑色点线为MD，蓝色圆点实线为Full，暖赭方点虚线为Intra。（a）CaMKI')
repl['128']='[29] J. A. Clapperton, S. R. Martin, S. J. Smerdon, S. J. Gamblin, and P. M. Bayley. Structure of the complex of calmodulin with the target sequence of calmodulin-dependent protein kinase I: studies of the kinase activation mechanism. Biochemistry, 41(50):14669–14679, 2002. doi:10.1021/bi026660t.'
repl['129']='[30] J. Goldberg, A. C. Nairn, and J. Kuriyan. Structural basis for the autoinhibition of calcium/calmodulin-dependent protein kinase I. Cell, 84(6):875–887, 1996. doi:10.1016/S0092-8674(00)81066-1.'
# Tighten the newly revised discussion/conclusion to prevent a stranded CPU-timing paragraph.
repl['92']='完整复合物环境使联合柔性预测能够定位具体的动态约束。CaM–CaMKI家族留出实例恢复了双锚点间的低主链波动核心，并保留同类型侧链的位置差异。结合C端叶瓣能够紧密结合、激活仍需N端叶瓣的实验，可提出N端叶瓣约束调控片段、帮助维持自抑制解除构型的假说；当前短肽复合物尚未直接检验回接触或活性。MEL5亲和成熟实例将HLA Q72的侧链状态集中与未替换TCR βGly54的主链波动降低联系到同一接触片区。起始结构的原子接近及重复轨迹中接触与氢键几何占据的增加，共同支持局部稳定化。联合读数可为已有功能实验提供有残基定位的动态解释和检验方向。'
repl['95']='FlexPPI以共享主链几何联合预测逐残基RMSF与侧链联合状态熵，提供互补的骨架波动和侧链状态读数。体系留出、外部方法与同类型评价支持预测优势，家族隔离对照支持完整结构环境的贡献。CaM–CaMKI实例为调控片段构型维持提出可检验假说，MEL5实例将联合变化联系到跨界面接触网络；83条记录与泛素留出预测进一步支持Hχ对实际角度覆盖及占据多样性的表征。'
for i,t in repl.items():set_text(pars[int(i)],t)
rels=etree.fromstring(parts['word/_rels/document.xml.rels']);rm={n.get('Id'):n.get('Target') for n in rels}
figures=[]
for index,imgfile,name in [(71,V/'figure3/fig3_1MXE_mechanistic_hypothesis.png','Figure 3 1MXE CaMKI flexibility and mechanistic hypothesis'),(174,V/'supplement_1MXE/figS3_1MXE_validation.png','Figure S3 1MXE family-held-out predictions and MD validation')]:
 p=pars[index];rid=p.xpath('.//a:blip/@r:embed',namespaces=NS);assert len(rid)==1
 assert len(root.xpath('.//a:blip[@r:embed="'+rid[0]+'"]',namespaces=NS))==1
 target='word/'+rm[rid[0]];parts[target]=imgfile.read_bytes();im=Image.open(imgfile);cx=round(164*36000);cy=round(cx*im.height/im.width)
 for node in p.xpath('.//wp:extent|.//a:xfrm/a:ext',namespaces=NS):node.set('cx',str(cx));node.set('cy',str(cy))
 for node in p.xpath('.//wp:docPr|.//pic:cNvPr',namespaces=NS):node.set('name',name);node.set('descr',name)
 set_setting(p,'keepNext',1);set_setting(pars[index+1],'keepNext',0);set_setting(pars[index+1],'keepLines',1)
 figures.append(dict(paragraph=index,target=target,file=str(imgfile),width_mm=164,height_mm=cy/36000,sha256=hashlib.sha256(parts[target]).hexdigest()))
# Preserve paragraphs and tables outside the authorized narrative/image replacements.
changed=set(map(int,repl))|{71,72,174,175}
for i,p in enumerate(pars):
 if i not in changed:assert signature(p)==before[i],i
assert tables==[signature(t) for t in body.findall('w:tbl',NS)]
all_text='\n'.join(''.join(p.xpath('.//w:t/text()',namespaces=NS)) for p in pars)
for old in ['NHE1','IQCG','2YGG','4LZX','CaM训练内案例','CaM两种识别构型','两条CaM轨迹']:
 assert old not in all_text,old
assert '机制假说' in all_text and '尚未直接' in all_text and '家族' in all_text
parts['word/document.xml']=serialize(root)
with ZipFile(OUTPUT,'w',ZIP_DEFLATED) as z:
 for name,data in parts.items():z.writestr(name,data)
with ZipFile(OUTPUT) as z:assert z.testzip() is None
(M/'v85_text_plan.json').write_text(json.dumps({'source_docx':str(SOURCE),'output_docx':str(OUTPUT),'paragraph_replacements':repl},ensure_ascii=False,indent=2),encoding='utf-8')
(M/'v85_final_text.txt').write_text(all_text,encoding='utf-8')
for index,path in [(72,V/'figure3/caption_zh.md'),(175,V/'supplement_1MXE/caption_zh.md')]:path.write_text(repl[str(index)]+'\n',encoding='utf-8')
(V/'narrative/3.3_定稿.md').write_text('\n\n'.join(repl[str(i)] for i in [68,69,70,73,74,75]),encoding='utf-8')
audit=dict(source=str(SOURCE),output=str(OUTPUT),modified_paragraphs=sorted(changed),protected_other_paragraphs=True,all_six_tables_identical=True,other_figures_unchanged=True,figures=figures,old_case_text_removed=True,docx_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest())
(M/'v85_build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
render=V/'render';render.mkdir(exist_ok=True)
for f in ['render_docx_word_wrapper.py','export_word_readonly.ps1']:shutil.copy2(V.parent/'revision_v83_model_story/render'/f,render/f)
print(json.dumps(audit,ensure_ascii=False))