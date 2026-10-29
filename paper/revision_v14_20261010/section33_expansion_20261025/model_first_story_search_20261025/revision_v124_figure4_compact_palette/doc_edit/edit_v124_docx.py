# Replace Figure4 and set its new aspect ratio; only its caption adds colour and SPR provenance.
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree
from io import BytesIO
from PIL import Image
from xml.sax.saxutils import escape
import hashlib,json,re
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
STORY=BASE/'section33_expansion_20261025/model_first_story_search_20261025'
HERE=Path(__file__).resolve().parent
SOURCE=BASE/'paper_v123_20261027_中文修订_图3橙色位点与黑色SASA.docx'
TARGET=BASE/'paper_v124_20261027_中文修订_图4紧凑蓝橙灰.docx'
FIGURE=BASE/'figure4_v124_20261027_紧凑蓝橙灰.png'
PRIOR=STORY/'revision_v123_orange_sites_black_sasa/doc_edit'
PART='word/media/image41.png'
sha=lambda x:hashlib.sha256(x).hexdigest()
source_sha=sha(SOURCE.read_bytes())
qa=json.loads((PRIOR/'render_QA.json').read_text(encoding='utf-8'))
assert qa['visual_status']=='all30pagespassed' and qa.get('sha256',qa.get('document_sha256'))==source_sha
prior_hashes={p.name:sha(p.read_bytes()) for p in (PRIOR/'render').glob('page-*.png')}
assert len(prior_hashes)==30
with ZipFile(SOURCE) as z:
 items=z.infolist();payload={i.filename:z.read(i.filename) for i in items}
original_xml=payload['word/document.xml'];xml=original_xml
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
root=etree.fromstring(xml)
captions=[p for p in root.findall('.//w:body/w:p',ns) if ''.join(p.xpath('.//w:t/text()',namespaces=ns)).startswith('图4：')]
assert len(captions)==1
cap=captions[0];old_caption=''.join(cap.xpath('.//w:t/text()',namespaces=ns))
new_caption=old_caption.replace('成熟复合物概览及Q72–β53/β54局部结构；','成熟复合物概览及Q72–β53/β54局部结构；TCRβ为蓝色，TCRα与HLA为灰色，肽与Q72为橙色，局部原子随所属残基着色并标注名称；解离常数来自SPR实验[28]；')
assert new_caption!=old_caption and len(cap.xpath('.//w:t',namespaces=ns))==1
old_cap_bytes=escape(old_caption).encode();new_cap_bytes=escape(new_caption).encode()
assert xml.count(old_cap_bytes)==1
xml=xml.replace(old_cap_bytes,new_cap_bytes)
rels=etree.fromstring(payload['word/_rels/document.xml.rels']);rid=[r.get('Id') for r in rels if r.get('Target')=='media/image41.png'];assert len(rid)==1
blocks=list(re.finditer(rb'<wp:inline\b[^>]*>.*?</wp:inline>',xml,re.S))
blocks=[m for m in blocks if ('r:embed="'+rid[0]+'"').encode() in m.group()];assert len(blocks)==1
m=blocks[0];block=m.group();inline=etree.fromstring(block,parser=etree.XMLParser(recover=True)) if False else None
old_extent=root.xpath('//wp:inline[.//a:blip/@r:embed="'+rid[0]+'"]',namespaces=ns)[0].find('wp:extent',ns)
cx=int(old_extent.get('cx'));old_cy=int(old_extent.get('cy'));new_cy=round(cx*152/183)
old_attr=('cy="'+str(old_cy)+'"').encode();new_attr=('cy="'+str(new_cy)+'"').encode();assert block.count(old_attr)==2
new_block=block.replace(old_attr,new_attr);xml=xml[:m.start()]+new_block+xml[m.end():]
# Exactly reverse the three intentional XML edits to recover original bytes.
assert xml.replace(new_cap_bytes,old_cap_bytes).replace(new_block,block)==original_xml
new_root=etree.fromstring(xml)
old_paragraphs=[''.join(p.xpath('.//w:t/text()',namespaces=ns)) for p in root.findall('.//w:body/w:p',ns)]
new_paragraphs=[''.join(p.xpath('.//w:t/text()',namespaces=ns)) for p in new_root.findall('.//w:body/w:p',ns)]
assert len(old_paragraphs)==len(new_paragraphs)
changed_text=[i for i,(a,b) in enumerate(zip(old_paragraphs,new_paragraphs)) if a!=b]
assert len(changed_text)==1 and old_paragraphs[changed_text[0]]==old_caption
assert len(root.xpath('//w:fldSimple|//w:instrText',namespaces=ns))==len(new_root.xpath('//w:fldSimple|//w:instrText',namespaces=ns))
new_image=FIGURE.read_bytes()
with Image.open(BytesIO(new_image)) as im:new_size=im.size;im.verify()
payload[PART]=new_image;payload['word/document.xml']=xml
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
 for item in items:z.writestr(item,payload[item.filename])
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
 assert a.namelist()==b.namelist() and b.testzip() is None
 changed_parts=[n for n in a.namelist() if a.read(n)!=b.read(n)]
 assert set(changed_parts)=={'word/document.xml',PART}
 assert sha(b.read('word/media/image5.png'))==sha(a.read('word/media/image5.png'))
assert sha(SOURCE.read_bytes())==source_sha
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':source_sha,'target_sha256':sha(TARGET.read_bytes()),'figure':str(FIGURE),'figure_sha256':sha(new_image),'figure_pixels':new_size,'changed_zip_parts':changed_parts,'unchanged_other_ZIP_parts':True,'unchanged_Figure3':True,'changed_caption_paragraph':changed_text[0],'caption_before':old_caption,'caption_after':new_caption,'drawing_extent_EMU':{'cx':cx,'old_cy':old_cy,'new_cy':new_cy},'drawing_extent_mm':{'width':cx/36000,'old_height':old_cy/36000,'new_height':new_cy/36000},'all_other_body_text_and_fields_identical':True,'validation':'Pending final render and page inspection.'}
(HERE/'document_change_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(HERE/'v123_prior_render_hashes.json').write_text(json.dumps(prior_hashes,indent=2),encoding='utf-8')
(HERE/'v123_prior_render_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
(HERE.parent/'caption_zh.md').write_text(new_caption+'\n',encoding='utf-8')
print(json.dumps({k:report[k] for k in ['target','target_sha256','changed_zip_parts','drawing_extent_mm','all_other_body_text_and_fields_identical']},ensure_ascii=False))
