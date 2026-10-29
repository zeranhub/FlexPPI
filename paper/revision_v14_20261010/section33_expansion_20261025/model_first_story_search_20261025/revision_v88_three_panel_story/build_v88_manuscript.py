from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree
from copy import deepcopy
from PIL import Image
import json,hashlib,shutil
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');V=R/'section33_expansion_20261025/model_first_story_search_20261025/revision_v88_three_panel_story';M=V/'manuscript';M.mkdir(exist_ok=True)
SOURCE=R/'paper_v87_20261026_中文修订_互补柔性背景与图3凝练.docx';OUTPUT=R/'paper_v88_20261026_中文修订_图3三面板主线.docx'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','pic':'http://schemas.openxmlformats.org/drawingml/2006/picture','pr':'http://schemas.openxmlformats.org/package/2006/relationships'}
def q(prefix,name):return '{'+NS[prefix]+'}'+name
def signature(x):return etree.tostring(x,method='c14n')
def set_text(p,text):
 old=p.findall('w:r',NS);style=deepcopy(old[0].find('w:rPr',NS)) if old and old[0].find('w:rPr',NS) is not None else None
 for n in list(p):
  if n.tag!=q('w','pPr'):p.remove(n)
 r=etree.SubElement(p,q('w','r'))
 if style is not None:r.append(style)
 etree.SubElement(r,q('w','t')).text=text
with ZipFile(SOURCE) as z:parts={i.filename:z.read(i.filename) for i in z.infolist()}
original_parts=parts.copy();root=etree.fromstring(parts['word/document.xml']);body=root.find('w:body',NS);pars=body.findall('w:p',NS);assert len(pars)==186
before=[signature(p) for p in pars];tables=[signature(t) for t in body.findall('w:tbl',NS)]
repl=json.loads((V/'narrative/paragraph_replacements.json').read_text(encoding='utf-8'))['paragraph_replacements']
for i,t in repl.items():set_text(pars[int(i)],t)
p=pars[72];rels=etree.fromstring(parts['word/_rels/document.xml.rels']);rm={n.get('Id'):n.get('Target') for n in rels}
rid=p.xpath('.//a:blip/@r:embed',namespaces=NS);assert len(rid)==1
img=V/'figure3/fig3_three_panel_story.png';target='word/'+rm[rid[0]];assert target=='word/media/image5.png';parts[target]=img.read_bytes()
im=Image.open(img);cx=round(164*36000);cy=round(cx*im.height/im.width)
for n in p.xpath('.//wp:extent|.//a:xfrm/a:ext',namespaces=NS):n.set('cx',str(cx));n.set('cy',str(cy))
for n in p.xpath('.//wp:docPr|.//pic:cNvPr',namespaces=NS):n.set('name','Figure 3 Complementary backbone and sidechain flexibility');n.set('descr','Three panels: full structural location, M316 static contact contexts, complete absolute flexibility profiles')
changed=set(map(int,repl))|{72}
for i,p in enumerate(pars):
 if i not in changed:assert signature(p)==before[i],i
assert tables==[signature(t) for t in body.findall('w:tbl',NS)]
parts['word/document.xml']=etree.tostring(root,encoding='UTF-8',xml_declaration=True,standalone=True)
changed_parts=[k for k in parts if parts[k]!=original_parts[k]];assert set(changed_parts)=={'word/document.xml',target}
with ZipFile(OUTPUT,'w',ZIP_DEFLATED) as z:
 for name,data in parts.items():z.writestr(name,data)
with ZipFile(OUTPUT) as z:assert z.testzip() is None
finaltext={str(i):''.join(p.xpath('.//w:t/text()',namespaces=NS)) for i,p in enumerate(pars)}
(M/'v88_text_plan.json').write_text(json.dumps({'source_docx':str(SOURCE),'output_docx':str(OUTPUT),'paragraph_replacements':repl,'index_basis':'Actual v87 source paragraph numbers'},ensure_ascii=False,indent=2),encoding='utf-8')
(M/'v88_final_text.txt').write_text('\n'.join(finaltext.values()),encoding='utf-8')
(V/'narrative/3.3_定稿.md').write_text('\n\n'.join(finaltext[str(i)] for i in [68,69,70,71,74,75,76]),encoding='utf-8')
(V/'figure3/caption_zh.md').write_text(repl['73']+'\n',encoding='utf-8')
audit=dict(source=str(SOURCE),output=str(OUTPUT),modified_paragraphs=sorted(changed),changed_zip_parts=changed_parts,all_other_original_paragraphs_identical=True,all_six_tables_identical=True,all_other_images_identical=True,figure_mm=[164,cy/36000],figure_sha256=hashlib.sha256(parts[target]).hexdigest(),docx_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest())
(M/'v88_build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
render=V/'render';render.mkdir(exist_ok=True)
for f in ['render_docx_word_wrapper.py','export_word_readonly.ps1']:shutil.copy2(V.parent/'revision_v87_general_story_compact/render'/f,render/f)
print(json.dumps(audit,ensure_ascii=False))