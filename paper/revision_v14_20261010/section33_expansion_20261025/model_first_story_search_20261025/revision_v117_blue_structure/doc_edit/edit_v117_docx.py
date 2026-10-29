"""Update the genuine Figure3 image and one colour word in its caption."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from io import BytesIO
from copy import deepcopy
from lxml import etree
from PIL import Image
import json,hashlib
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
STORY=BASE/'section33_expansion_20261025/model_first_story_search_20261025'
HERE=Path(__file__).resolve().parent
SOURCE=BASE/'paper_v116_20261027_中文修订_图3图例外置.docx'
TARGET=BASE/'paper_v117_20261027_中文修订_图3蓝色结构与SASA同色.docx'
FIGURE=BASE/'figure3_v117_20261027_蓝色结构与SASA同色.png'
PRIOR=STORY/'revision_v116_external_density_legend/doc_edit'
PART='word/media/image5.png'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
def sha(x):return hashlib.sha256(x).hexdigest()
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
source_sha=sha(SOURCE.read_bytes())
qa=json.loads((PRIOR/'render_QA.json').read_text(encoding='utf-8'))
assert qa['visual_status']=='all30pagespassed' and qa['sha256']==source_sha
old_hashes={p.name:sha(p.read_bytes()) for p in (PRIOR/'render').glob('page-*.png')}
assert len(old_hashes)==30
with ZipFile(SOURCE) as z:
    items=z.infolist();payload={i.filename:z.read(i.filename) for i in items}
root=etree.fromstring(payload['word/document.xml']);original=deepcopy(root)
paragraphs=root.xpath('//w:body/w:p',namespaces=NS)
caption=next(p for p in paragraphs if text(p).startswith('图3：'))
before=text(caption);count=0
for t in caption.xpath('.//w:t',namespaces=NS):
    if t.text and '紫色为分析链' in t.text:
        t.text=t.text.replace('紫色为分析链','蓝色为分析链');count+=1
assert count==1
after=text(caption)
assert before.replace('紫色为分析链','蓝色为分析链')==after
new_image=FIGURE.read_bytes()
with Image.open(BytesIO(payload[PART])) as im:old_size=im.size
with Image.open(BytesIO(new_image)) as im:new_size=im.size;im.verify()
assert old_size==new_size
payload[PART]=new_image
payload['word/document.xml']=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
    for item in items:z.writestr(item,payload[item.filename])
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
    changed=[n for n in a.namelist() if a.read(n)!=b.read(n)]
    assert set(changed)=={'word/document.xml',PART},changed
    assert b.testzip() is None
oldps=original.xpath('//w:body/w:p',namespaces=NS)
changed_p=[i for i,(a,b) in enumerate(zip(oldps,paragraphs)) if etree.tostring(a)!=etree.tostring(b)]
assert changed_p==[paragraphs.index(caption)]
assert original.xpath('//w:fldChar/@w:fldCharType',namespaces=NS)==root.xpath('//w:fldChar/@w:fldCharType',namespaces=NS)
assert sha(SOURCE.read_bytes())==source_sha
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':source_sha,'target_sha256':sha(TARGET.read_bytes()),
    'figure':str(FIGURE),'figure_sha256':sha(new_image),'image_pixels':new_size,
    'changed_zip_parts':changed,'changed_paragraph_indices':changed_p,
    'caption_before':before,'caption_after':after,'other_paragraphs_and_ZIP_parts_identical':True,
    'validation':'ZIP/XML checks passed; pending rendered changed-page review'}
(HERE/'document_change_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(HERE/'v116_prior_render_hashes.json').write_text(json.dumps(old_hashes,indent=2),encoding='utf-8')
(HERE/'v116_prior_render_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if not k.startswith('caption')},ensure_ascii=False))
