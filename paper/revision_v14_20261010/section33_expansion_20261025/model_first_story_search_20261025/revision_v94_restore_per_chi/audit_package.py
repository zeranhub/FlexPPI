from pathlib import Path
import json,hashlib,zipfile,csv
from pypdf import PdfReader
from lxml import etree
R=Path(__file__).resolve().parent;OLD=R.parent/'revision_v93_dual_case_entropy';BASE=R.parents[2];O=R/'audit';O.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
DOC=BASE/'paper_v94_20261026_中文修订_恢复逐角占据关联.docx';SRC=BASE/'paper_v93_20261026_中文修订_双案例与构象占据整合.docx'
assert sha(SRC)=='6c603bfff03d033610d228854c05c0cb65bce4ceee3e9f9dae2bf8b3433aab44'
render=R/'render/release_01';pdf=render/(DOC.stem+'.pdf');reader=PdfReader(pdf);assert len(reader.pages)==29
changed=[i for i in range(1,30) if sha(render/f'page-{i}.png')!=sha(OLD/f'render/release_05/page-{i}.png')];assert changed==[11,12,13,29]
with zipfile.ZipFile(SRC) as z:before={n:z.read(n) for n in z.namelist()}
with zipfile.ZipFile(DOC) as z:
    assert z.testzip() is None
    after={n:z.read(n) for n in z.namelist()}
assert {n for n in before if before[n]!=after[n]}=={'word/document.xml','word/media/image5.png','word/media/image40.png'}
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
root=etree.fromstring(after['word/document.xml']);old=etree.fromstring(before['word/document.xml'])
assert [etree.tostring(x,method='c14n') for x in root.findall('.//w:tbl',ns)]==[etree.tostring(x,method='c14n') for x in old.findall('.//w:tbl',ns)]
texts=[''.join(x.xpath('.//w:t/text()',namespaces=ns)) for x in root.findall('.//w:body/w:p',ns)];full='\n'.join(texts)
assert 'χ1–χ3呈中等正相关' in full and '图3f' in full and '图S4c' in full
assert not any(x in full for x in ['图5：','图5d','图 S4c、d'])
assert json.loads((O/'figure_numerical_audit.json').read_text(encoding='utf-8'))['status']=='passed'
assert '图3：' in reader.pages[10].extract_text().replace(' ','') and '图S4：' in reader.pages[28].extract_text().replace(' ','')
assert not any('\ufffd' in p.extract_text() for p in reader.pages)
audit=dict(status='passed',docx_sha256=sha(DOC),source_unchanged=True,pages=29,changed_pages=changed,unchanged25_pages_inherit_previous_visual_QA=True,root_actually_viewed_all4_changed_pages=True,main_figure_and_caption_page=11,S4_and_caption_page=29,all_six_tables_unchanged=True,other44_media_unchanged=True,source_value_audit='figure_numerical_audit.json')
(O/'final_document_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
files=sorted(p for p in R.rglob('*') if p.is_file() and not {'render','__pycache__'}.intersection(p.relative_to(R).parts) and p.name!='archive_contents.csv' and p.suffix!='.zip')
manifest=O/'archive_contents.csv'
with manifest.open('w',encoding='utf-8-sig',newline='') as out:
    w=csv.writer(out);w.writerow(['relative_path','bytes','sha256'])
    for p in files:w.writerow([p.relative_to(R).as_posix(),p.stat().st_size,sha(p)])
files.append(manifest);archive=BASE/'paper_v94_20261026_恢复逐角关联_源文件与审计.zip'
if archive.exists():
    assert sha(archive)=='d660bcbb084db14f4d66feda7bd340d6cab667f960495141b3017fd3081c720f','Refuse to overwrite an unrelated archive'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,Path(R.name)/p.relative_to(R))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps(dict(**audit,archive_files=len(files),archive_bytes=archive.stat().st_size,archive_sha256=sha(archive)),ensure_ascii=False))
