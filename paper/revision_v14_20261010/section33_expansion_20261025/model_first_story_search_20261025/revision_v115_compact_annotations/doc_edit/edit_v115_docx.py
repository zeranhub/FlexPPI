"""v115: finalized Figure3 image plus truthful annotation wording only."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree
from copy import deepcopy
from PIL import Image
import hashlib,json,re
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
STORY=BASE/'section33_expansion_20261025/model_first_story_search_20261025'
HERE=Path(__file__).resolve().parent
SOURCE=BASE/'paper_v114_20261027_中文修订_双案例合并与泛素视角.docx'
TARGET=BASE/'paper_v115_20261027_中文修订_图3紧凑标注精修.docx'
FIGURE=BASE/'figure3_v115_20261027_紧凑标注与统一细线.png'
PRIOR=STORY/'revision_v114_combined_profiles/doc_edit'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
W='{'+NS['w']+'}'
def sha(data):return hashlib.sha256(data).hexdigest()
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
assert FIGURE.is_file(),FIGURE
source_sha=sha(SOURCE.read_bytes())
prior_qa=json.loads((PRIOR/'render_QA.json').read_text(encoding='utf-8'))
assert prior_qa['visual_status']=='all30pagespassed' and prior_qa['sha256']==source_sha
prior_hashes={p.name:sha(p.read_bytes()) for p in (PRIOR/'render').glob('page-*.png')}
assert len(prior_hashes)==prior_qa['page_count']
with ZipFile(SOURCE) as z:
    items=z.infolist();payload={i.filename:z.read(i.filename) for i in items}
root=etree.fromstring(payload['word/document.xml']);original=deepcopy(root)
ps=root.xpath('//w:body/w:p',namespaces=NS)
caption=[p for p in ps if text(p).startswith('图3：')]
assert len(caption)==1
caption=caption[0];old_caption=text(caption)
new_caption=old_caption.replace('青绿色标记W303/M316与Q41/Q62','青绿色位点名称及标记对应W303/M316与Q41/Q62').replace('青绿色位置标记对应W303/M316与Q41/Q62','细竖虚线指示W303/M316与Q41/Q62')
assert new_caption!=old_caption
assert '细竖虚线指示W303/M316与Q41/Q62' in new_caption
assert not caption.xpath('.//w:fldChar|.//w:fldSimple|.//w:drawing',namespaces=NS)
runs=caption.findall(W+'r')
rpr=next((deepcopy(r.find(W+'rPr')) for r in runs if r.find(W+'rPr') is not None),None)
for child in list(caption):
    if child.tag!=W+'pPr':caption.remove(child)
r=etree.SubElement(caption,W+'r')
if rpr is not None:r.append(rpr)
etree.SubElement(r,W+'t').text=new_caption
figp=[p for p in ps if p.xpath('.//a:blip[@r:embed="rId13"]',namespaces=NS)]
assert len(figp)==1
figp=figp[0]
rels=etree.fromstring(payload['word/_rels/document.xml.rels'])
rel=[r for r in rels if r.get('Id')=='rId13'][0]
part='word/'+rel.get('Target');assert part=='word/media/image5.png'
new_image=FIGURE.read_bytes()
with Image.open(FIGURE) as im:iw,ih=im.size;im.verify()
extent=figp.xpath('.//wp:extent',namespaces=NS)[0]
cx=int(extent.get('cx'));old_cy=int(extent.get('cy'));cy=round(cx*ih/iw)
extent.set('cy',str(cy))
for x in figp.xpath('.//a:xfrm/a:ext',namespaces=NS):x.set('cx',str(cx));x.set('cy',str(cy))
payload[part]=new_image
payload['word/document.xml']=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
    for item in items:z.writestr(item,payload[item.filename])
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
    changed_parts=[n for n in a.namelist() if a.read(n)!=b.read(n)]
    assert set(changed_parts)=={'word/document.xml',part}
    assert b.testzip() is None
    assert sha(b.read(part))==sha(new_image)
oldps=original.xpath('//w:body/w:p',namespaces=NS)
assert len(oldps)==len(ps)
changed_indices=[i for i,(a,b) in enumerate(zip(oldps,ps)) if etree.tostring(a)!=etree.tostring(b)]
expected=sorted([ps.index(caption)]+([ps.index(figp)] if cy!=old_cy else []))
assert changed_indices==expected
fields=lambda r:[etree.tostring(p) for p in r.xpath('//w:fldChar|//w:fldSimple',namespaces=NS)]
refs=lambda r:[text(p) for p in r.xpath('//w:body/w:p',namespaces=NS) if re.match(r'^\[\d+\]',text(p))]
assert fields(original)==fields(root)
assert refs(original)==refs(root)
assert sha(SOURCE.read_bytes())==source_sha
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':source_sha,'target_sha256':sha(TARGET.read_bytes()),'figure':str(FIGURE),'figure_sha256':sha(new_image),'figure_pixels':[iw,ih],'on_page_size_mm':[cx/36000,cy/36000],'prior_on_page_size_mm':[cx/36000,old_cy/36000],'changed_zip_parts':changed_parts,'changed_paragraph_indices':changed_indices,'preserved_fields':len(fields(root)),'preserved_references':len(refs(root)),'caption_before':old_caption,'caption_after':new_caption,'validation':'ZIP/XML invariants passed; pending render and changed-page visual QA'}
(HERE/'document_change_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(HERE/'v114_prior_render_hashes.json').write_text(json.dumps(prior_hashes,indent=2),encoding='utf-8')
(HERE/'v114_prior_render_QA.json').write_text(json.dumps(prior_qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ('caption_before','caption_after')},ensure_ascii=False,indent=2))

