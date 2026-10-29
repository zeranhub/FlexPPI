"""Update only Figure3 image and its caption from v112 to v113.
Run using bundled Codex Python after the new scientific figure is finalized.
"""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree
from copy import deepcopy
from PIL import Image
import hashlib,json
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
STORY=BASE/'section33_expansion_20261025/model_first_story_search_20261025'
HERE=Path(__file__).resolve().parent
SOURCE=BASE/'paper_v112_20261027_中文修订_双尺度柔性与角度占据整合.docx'
TARGET=BASE/'paper_v113_20261027_中文修订_图3对齐与图例精修.docx'
FIGURE=BASE/'figure3_v113_20261027_轴标题对齐与双色图例.png'
PRIOR_QA=STORY/'revision_v112_compact_story/doc_edit/render_QA.json'
PRIOR_RENDER=STORY/'revision_v112_compact_story/doc_edit/render'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
W='{'+NS['w']+'}'
def sha(data):return hashlib.sha256(data).hexdigest()
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
assert FIGURE.is_file(),FIGURE
source_sha=sha(SOURCE.read_bytes())
prior_qa=json.loads(PRIOR_QA.read_text(encoding='utf-8'))
assert prior_qa['visual_status']=='all30pagespassed'
assert prior_qa['sha256']==source_sha
prior_hashes={p.name:sha(p.read_bytes()) for p in PRIOR_RENDER.glob('page-*.png')}
assert len(prior_hashes)==30
with ZipFile(SOURCE) as z:
    items=z.infolist();payload={i.filename:z.read(i.filename) for i in items}
root=etree.fromstring(payload['word/document.xml'])
original=deepcopy(root)
ps=root.xpath('//w:body/w:p',namespaces=NS)
captions=[p for p in ps if text(p).startswith('图3：')]
assert len(captions)==1
caption=captions[0];old_caption=text(caption)
assert '残基横向间距相同' in old_caption
assert '青灰柱' in old_caption
new_caption=old_caption.replace('两例纵轴范围一致，残基横向间距相同。下方为共用残基横轴的起始复合物侧链SASA（青灰柱，Å²）','两例纵轴范围一致；每例曲线与SASA共用残基横轴。下方为起始复合物侧链SASA（浅绿色柱，Å²）')
assert new_caption!=old_caption
assert '残基横向间距相同' not in new_caption
assert '青灰柱' not in new_caption
runs=caption.findall(W+'r')
rpr=next((deepcopy(r.find(W+'rPr')) for r in runs if r.find(W+'rPr') is not None),None)
for child in list(caption):
    if child.tag!=W+'pPr':caption.remove(child)
r=etree.SubElement(caption,W+'r')
if rpr is not None:r.append(rpr)
etree.SubElement(r,W+'t').text=new_caption
# The figure's relation and media filename stay unchanged.
figures=[p for p in ps if p.xpath('.//a:blip[@r:embed="rId13"]',namespaces=NS)]
assert len(figures)==1
figp=figures[0]
rels=etree.fromstring(payload['word/_rels/document.xml.rels'])
rel=[r for r in rels if r.get('Id')=='rId13'][0]
image_part='word/'+rel.get('Target')
assert image_part=='word/media/image5.png'
new_image=FIGURE.read_bytes()
with Image.open(FIGURE) as im:iw,ih=im.size;im.verify()
extent=figp.xpath('.//wp:extent',namespaces=NS)[0]
cx=int(extent.get('cx'));old_cy=int(extent.get('cy'));cy=round(cx*ih/iw)
extent.set('cy',str(cy))
for x in figp.xpath('.//a:xfrm/a:ext',namespaces=NS):x.set('cx',str(cx));x.set('cy',str(cy))
payload[image_part]=new_image
payload['word/document.xml']=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
    for item in items:z.writestr(item,payload[item.filename])
# Every paragraph other than Figure3 and its caption must be byte-identical XML.
oldps=original.xpath('//w:body/w:p',namespaces=NS)
assert len(oldps)==len(ps)
changed_indices=[i for i,(a,b) in enumerate(zip(oldps,ps)) if etree.tostring(a)!=etree.tostring(b)]
expected=sorted([ps.index(caption),ps.index(figp)])
assert set(changed_indices).issubset(set(expected)),(changed_indices,expected)
assert ps.index(caption) in changed_indices
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
    changed_parts=[n for n in a.namelist() if a.read(n)!=b.read(n)]
    assert set(changed_parts)=={'word/document.xml',image_part}
    assert b.testzip() is None
    assert sha(b.read(image_part))==sha(new_image)
oldfields=[etree.tostring(p) for p in original.xpath('//w:fldChar|//w:fldSimple',namespaces=NS)]
newfields=[etree.tostring(p) for p in root.xpath('//w:fldChar|//w:fldSimple',namespaces=NS)]
assert oldfields==newfields
assert sha(SOURCE.read_bytes())==source_sha
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':source_sha,'target_sha256':sha(TARGET.read_bytes()),'figure':str(FIGURE),'figure_sha256':sha(new_image),'figure_pixels':[iw,ih],'on_page_size_mm':[cx/36000,cy/36000],'prior_on_page_size_mm':[cx/36000,old_cy/36000],'changed_zip_parts':changed_parts,'changed_paragraph_indices':changed_indices,'preserved_fields':len(oldfields),'caption_before':old_caption,'caption_after':new_caption,'validation':'ZIP/XML invariants passed; pending render and changed-page visual review'}
(HERE/'document_change_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(HERE/'v112_prior_render_hashes.json').write_text(json.dumps(prior_hashes,indent=2),encoding='utf-8')
(HERE/'v112_prior_render_QA.json').write_text(json.dumps(prior_qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ('caption_before','caption_after')},ensure_ascii=False,indent=2))

