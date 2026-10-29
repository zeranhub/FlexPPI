"""Replace Figure 2 and the approved targeted prose while preserving native parts."""
from pathlib import Path
from copy import deepcopy
import hashlib,json,re,zipfile
from lxml import etree
from PIL import Image

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v44_20261018_中文修订_预测一致性与高动态识别_三线表.docx'
OUT=BASE/'paper_v45_20261018_中文修订_紧凑多层模型评价图_三线表.docx'
FIG=ROOT/'figure2/fig2_dense_capability.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='4156791ba646d8137754386ea387b920aaaf1bdac041614414a17028915ffcab'
assert FIG.is_file()
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
 'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
 'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
 'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
 'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
with zipfile.ZipFile(SOURCE) as z:
    infos=z.infolist();parts={n.filename:z.read(n.filename) for n in infos}
root=etree.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
old_nodes=list(body);txt=lambda p:''.join(p.xpath('.//w:t/text()',namespaces=NS))
def find(prefix):
    matches=[p for p in body if p.tag=='{'+NS['w']+'}p' and txt(p).startswith(prefix)]
    assert len(matches)==1,(prefix,len(matches))
    return matches[0]
changes=[];insertions=[]
def set_text(p,value):
    assert not p.xpath('.//m:oMath',namespaces=NS)
    ts=p.xpath('.//w:t',namespaces=NS);assert ts
    previous=txt(p)
    if previous==value:return
    changes.append({'source_body_index':old_nodes.index(p),'before':previous,'after':value})
    ts[0].text=value
    for n in ts[1:]:n.text=''
def set_property(p,name,value):
    pr=p.find('w:pPr',NS)
    if pr is None:pr=etree.Element('{'+NS['w']+'}pPr');p.insert(0,pr)
    el=pr.find('w:'+name,NS)
    if el is None:el=etree.SubElement(pr,'{'+NS['w']+'}'+name)
    el.set('{'+NS['w']+'}val',str(value))

proposal=json.loads((ROOT/'text_review/manuscript_proposal.json').read_text(encoding='utf8'))
assert proposal['source_sha256']==sha(SOURCE)
for item in proposal['replacements']:
    node=find(item['prefix']);assert txt(node)==item['old']
    set_text(node,item['new'])
template=find('评价先在每条记录的有效残基上计算')
assert not template.xpath('.//m:oMath|.//w:commentRangeStart|.//w:commentReference|.//w:bookmarkStart',namespaces=NS)
for item in proposal['insertions']:
    anchor=find(item['anchor_prefix']);node=etree.Element('{'+NS['w']+'}p')
    pr=template.find('w:pPr',NS)
    if pr is not None:node.append(deepcopy(pr))
    run=etree.SubElement(node,'{'+NS['w']+'}r')
    rpr=template.find('w:r/w:rPr',NS)
    if rpr is not None:run.append(deepcopy(rpr))
    etree.SubElement(run,'{'+NS['w']+'}t').text=item['text']
    assert item['position'] in ('before','after')
    if item['position']=='before':anchor.addprevious(node)
    else:anchor.addnext(node)
    insertions.append(item)
cap=find('图2：');fig=cap.getprevious();assert fig.xpath('.//a:blip',namespaces=NS)
set_text(cap,(ROOT/'figure2/caption_zh.md').read_text(encoding='utf8').strip())
blip=fig.xpath('.//a:blip',namespaces=NS);assert len(blip)==1
rid=blip[0].get('{'+NS['r']+'}embed')
rels=etree.fromstring(parts['word/_rels/document.xml.rels']);rm={n.get('Id'):n.get('Target') for n in rels}
media='word/'+rm[rid];assert media=='word/media/image2.png'
inline=blip[0].xpath('ancestor::wp:inline',namespaces=NS)[0]
width=round(6.6*914400)
with Image.open(FIG) as im:height=round(width*im.height/im.width)
for extent in inline.xpath('./wp:extent|.//a:xfrm/a:ext',namespaces=NS):
    extent.set('cx',str(width));extent.set('cy',str(height))
for prop in inline.xpath('./wp:docPr',namespaces=NS):
    prop.set('name','Figure 2 compact comparative prediction evidence')
    prop.set('descr','Record-weighted prediction agreement and four-method trends, whole-record and within-type precision-recall, and source-stratified average precision')
set_property(fig,'keepNext',1);set_property(cap,'keepLines',1);set_property(cap,'keepNext',0)
alltext='\n'.join(txt(n) for n in body)
assert len(root.xpath('.//w:tbl',namespaces=NS))==10
assert len(root.xpath('.//m:oMath',namespaces=NS))==36
assert len(root.xpath('.//a:blip',namespaces=NS))==6
assert len(re.findall(r'^\[\d+\]',alltext,re.M))==28
assert len(re.findall(r'^图[1-5]：',alltext,re.M))==5
for token in ['发布权重','发布模型','冻结','SB27','SB47']:assert token not in alltext
updates={media:FIG.read_bytes(),'word/document.xml':etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)}
with zipfile.ZipFile(OUT,'w') as z:
    for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
with zipfile.ZipFile(OUT) as z:assert z.testzip() is None
paragraphs=[p for p in body if p.tag=='{'+NS['w']+'}p']
(QA/'v45_content.txt').write_text('\n'.join(f'P{i}: {txt(p)}' for i,p in enumerate(paragraphs)),encoding='utf8')
audit={'status':'built_pending_render_QA','source':str(SOURCE),'source_sha256':sha(SOURCE),
 'document':str(OUT),'document_sha256':sha(OUT),'changes':changes,'insertions':insertions,
 'changed_zip_parts':list(updates),'figure2':str(FIG),'figure2_sha256':sha(FIG),
 'width_EMU':width,'height_EMU':height,'native_tables':10,'native_equations':36,'references':28,
 'main_figures':5,'supplementary_figures':1,'new_MD':0,'new_training':0,'new_inference':0}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
