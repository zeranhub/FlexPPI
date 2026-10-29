"""v114: replace finalized Figure3 and synchronize its four-panel references.
All scientific narrative and results remain inherited from v113.
"""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree
from copy import deepcopy
from PIL import Image
import hashlib,json,re
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
STORY=BASE/'section33_expansion_20261025/model_first_story_search_20261025'
HERE=Path(__file__).resolve().parent
SOURCE=BASE/'paper_v113_20261027_中文修订_图3对齐与图例精修.docx'
TARGET=BASE/'paper_v114_20261027_中文修订_双案例合并与泛素视角.docx'
FIGURE=BASE/'figure3_v114_20261027_双案例合并与结构视角.png'
PRIOR=STORY/'revision_v113_axis_legend/doc_edit'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
W='{'+NS['w']+'}'
def sha(data):return hashlib.sha256(data).hexdigest()
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def set_text(p,new):
    assert not p.xpath('.//w:fldChar|.//w:fldSimple|.//w:drawing',namespaces=NS)
    runs=p.findall(W+'r')
    rpr=next((deepcopy(r.find(W+'rPr')) for r in runs if r.find(W+'rPr') is not None),None)
    for child in list(p):
        if child.tag!=W+'pPr':p.remove(child)
    r=etree.SubElement(p,W+'r')
    if rpr is not None:r.append(rpr)
    etree.SubElement(r,W+'t').text=new
assert FIGURE.is_file(),FIGURE
source_sha=sha(SOURCE.read_bytes())
qa=json.loads((PRIOR/'render_QA.json').read_text(encoding='utf-8'))
assert qa['visual_status']=='all30pagespassed'
assert qa['sha256']==source_sha
old_hashes={p.name:sha(p.read_bytes()) for p in (PRIOR/'render').glob('page-*.png')}
assert len(old_hashes)==30
with ZipFile(SOURCE) as z:
    items=z.infolist();payload={i.filename:z.read(i.filename) for i in items}
root=etree.fromstring(payload['word/document.xml']);original=deepcopy(root)
ps=root.xpath('//w:body/w:p',namespaces=NS)
changes=[]
for index,p in enumerate(ps):
    old=text(p);new=old
    if old.startswith('图3：'):
        new=new.replace('（c、d）CaMKI294–318与泛素1–72的RMSF','（c）左为CaMKI294–318，右为泛素1–72，两例均显示RMSF')
        new=new.replace('（浅绿色柱，Å²）','（青灰柱，Å²）')
        new=new.replace('（e）同一联合Hχ预测','（d）同一联合Hχ预测')
        assert '（c、d）' not in new and '（e）' not in new
    elif index in (70,71):
        new=new.replace('图3c','图3c左')
    elif index==74:
        new=new.replace('图3d','图3c右')
    elif index==75:
        new=new.replace('图3a、d','图3a及c右')
    elif index in (76,164,182):
        new=new.replace('图3e','图3d')
    elif index==155:
        new=new.replace('图3a、c中的1MXE预测','图3a及c左中的1MXE预测')
    elif index==165:
        new=new.replace('图3c、d两例','图3c两例')
        new=new.replace('图3c的CaMKI','图3c左侧CaMKI').replace('图3d的泛素','图3c右侧泛素').replace('图3c的SASA','图3c左侧的SASA')
    elif index==168:
        new=new.replace('预测见图3d','预测见图3c右').replace('主文图3a、d、e','主文图3a、c右及d')
    if new!=old:
        set_text(p,new);changes.append({'paragraph_index':index,'before':old,'after':new})
expected_changed=[70,71,73,74,75,76,155,164,165,168,182]
assert [x['paragraph_index'] for x in changes]==expected_changed
alltext='\n'.join(text(p) for p in ps)
assert '图3e' not in alltext
assert '（c、d）CaMKI' not in alltext
assert '浅绿色柱' not in alltext
assert '3.4 TCR–pMHC识别中的界面柔性与亲和力差异' in alltext
# Exactly the same image relationship and physical width.
figures=[p for p in ps if p.xpath('.//a:blip[@r:embed="rId13"]',namespaces=NS)]
assert len(figures)==1
figp=figures[0]
rels=etree.fromstring(payload['word/_rels/document.xml.rels'])
rel=[r for r in rels if r.get('Id')=='rId13'][0]
image_part='word/'+rel.get('Target');assert image_part=='word/media/image5.png'
new_image=FIGURE.read_bytes()
with Image.open(FIGURE) as im:iw,ih=im.size;im.verify()
extent=figp.xpath('.//wp:extent',namespaces=NS)[0]
cx=int(extent.get('cx'));old_cy=int(extent.get('cy'));cy=round(cx*ih/iw)
extent.set('cy',str(cy))
for x in figp.xpath('.//a:xfrm/a:ext',namespaces=NS):x.set('cx',str(cx));x.set('cy',str(cy))
if cy!=old_cy:changes.append({'paragraph_index':ps.index(figp),'image_height_changed_to_emu':cy})
payload[image_part]=new_image
payload['word/document.xml']=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
    for item in items:z.writestr(item,payload[item.filename])
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
    changed_parts=[n for n in a.namelist() if a.read(n)!=b.read(n)]
    assert set(changed_parts)=={'word/document.xml',image_part}
    assert b.testzip() is None
    assert sha(b.read(image_part))==sha(new_image)
oldps=original.xpath('//w:body/w:p',namespaces=NS)
assert len(oldps)==len(ps)
actual_changed=[i for i,(a,b) in enumerate(zip(oldps,ps)) if etree.tostring(a)!=etree.tostring(b)]
assert actual_changed==sorted([x['paragraph_index'] for x in changes])
fields=lambda r:[etree.tostring(p) for p in r.xpath('//w:fldChar|//w:fldSimple',namespaces=NS)]
assert fields(original)==fields(root)
refs=lambda r:[text(p) for p in r.xpath('//w:body/w:p',namespaces=NS) if re.match(r'^\[\d+\]',text(p))]
assert refs(original)==refs(root)
assert len(refs(root))==30
assert sha(SOURCE.read_bytes())==source_sha
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':source_sha,'target_sha256':sha(TARGET.read_bytes()),'figure':str(FIGURE),'figure_sha256':sha(new_image),'figure_pixels':[iw,ih],'on_page_size_mm':[cx/36000,cy/36000],'changed_zip_parts':changed_parts,'changed_paragraph_indices':actual_changed,'preserved_fields':len(fields(root)),'preserved_references':len(refs(root)),'paragraph_changes':changes,'validation':'ZIP/XML checks passed; pending render and changed-page visual QA'}
(HERE/'document_change_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(HERE/'v113_prior_render_hashes.json').write_text(json.dumps(old_hashes,indent=2),encoding='utf-8')
(HERE/'v113_prior_render_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='paragraph_changes'},ensure_ascii=False,indent=2))


