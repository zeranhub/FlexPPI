"""v119 replaces Figure3's image part only, retaining all document XML."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from io import BytesIO
from PIL import Image
import hashlib,json
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
STORY=BASE/'section33_expansion_20261025/model_first_story_search_20261025'
HERE=Path(__file__).resolve().parent
SOURCE=BASE/'paper_v118_20261027_中文修订_图3蓝橙灰统一配色.docx'
TARGET=BASE/'paper_v119_20261027_中文修订_图3对齐图2标准色.docx'
FIGURE=BASE/'figure3_v119_20261027_图2标准蓝橙统一.png'
PRIOR=STORY/'revision_v118_blue_orange/doc_edit'
PART='word/media/image5.png'
def sha(data):return hashlib.sha256(data).hexdigest()
assert FIGURE.is_file(),FIGURE
source_sha=sha(SOURCE.read_bytes())
qa=json.loads((PRIOR/'render_QA.json').read_text(encoding='utf-8'))
assert qa['visual_status']=='all30pagespassed' and qa['sha256']==source_sha
prior_hashes={p.name:sha(p.read_bytes()) for p in (PRIOR/'render').glob('page-*.png')}
assert len(prior_hashes)==30
with ZipFile(SOURCE) as z:
    items=z.infolist();payload={i.filename:z.read(i.filename) for i in items}
old_image=payload[PART];new_image=FIGURE.read_bytes()
with Image.open(BytesIO(old_image)) as im:old_size=im.size
with Image.open(BytesIO(new_image)) as im:new_size=im.size;im.verify()
assert old_size==new_size,(old_size,new_size)
payload[PART]=new_image
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
    for item in items:z.writestr(item,payload[item.filename])
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
    assert a.namelist()==b.namelist()
    changed_parts=[n for n in a.namelist() if a.read(n)!=b.read(n)]
    assert changed_parts==[PART]
    assert a.read('word/document.xml')==b.read('word/document.xml')
    assert b.testzip() is None
    assert sha(b.read(PART))==sha(new_image)
assert sha(SOURCE.read_bytes())==source_sha
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':source_sha,'target_sha256':sha(TARGET.read_bytes()),'figure':str(FIGURE),'figure_sha256':sha(new_image),'figure_pixels':new_size,'changed_zip_parts':changed_parts,'all_XML_identical':True,'all_other_ZIP_parts_identical':True,'validation':'OnlyFigure3imagepartchanged;pendingrenderandchangedpagevisualQA'}
(HERE/'document_change_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(HERE/'v118_prior_render_hashes.json').write_text(json.dumps(prior_hashes,indent=2),encoding='utf-8')
(HERE/'v118_prior_render_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
