from pathlib import Path
from datetime import datetime, timezone
import csv, hashlib, json, zipfile
from xml.etree import ElementTree as ET
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
assert sha(ROOT/'analysis_plan.json')=='fb676c2ee8fc1a0c4c02cfef100c920051368bb5ba30b7c3f57edfa7abf836f6'
review=read(ROOT/'final_case_group_review.json')
assert review['status']=='numeric_verification_passed_with_scientific_scope_restrictions'
for item in review['sources']:
    assert sha(item['path'])==item['sha256'],item['path']
manifest_path=ROOT/'case_figure/figure_manifest.json'
figure=read(manifest_path)
for item in figure['source_files']:
    assert sha(item['path'])==item['sha256']
rows=list(csv.DictReader((ROOT/'case_figure/figure_exact_source.csv').open(encoding='utf-8-sig',newline='')))
assert {k:sum(r['panel']==k for r in rows) for k in 'abcd'}=={'a':4,'b':40,'c':135,'d':90}
assert [sum(r['panel']=='d' and r['contrast']==contrast and r['same_AA_valid_Hchi']=='True' for r in rows) for contrast in ['TCR maturation','Peptide G4A','Peptide I7A']]==[23,23,22]
pdf_path=ROOT/'case_figure/MEL5_opposite_affinity_case_series.pdf'
pdf=PdfReader(pdf_path)
font_audit=[]
for _,reference in pdf.pages[0]['/Resources']['/Font'].get_object().items():
    f=reference.get_object()
    if '/DescendantFonts' in f:
        f=f['/DescendantFonts'][0].get_object()
    descriptor=f.get('/FontDescriptor')
    assert descriptor is not None
    descriptor=descriptor.get_object()
    embedded=any(key in descriptor for key in ['/FontFile','/FontFile2','/FontFile3'])
    assert embedded
    font_audit.append(dict(font=str(f.get('/BaseFont')),embedded=embedded))
figure.update(status='passed',visual_QA=dict(status='passed',image='MEL5_opposite_affinity_case_series.png',checks=['All40 peptide values and all3 fixed contrast rows retained','All45CDR and30pMHC positions retained with same-AA chemical masking','G4A P10 and Q72 adverse model output patterns visible','English labels readable with no clipping or overlap; editable SVG and embedded-font PDF verified','No MD line or fitted affinity direction; exact source table archived']),pdf_font_embedding_audit=font_audit,script_sha256=sha(ROOT/'plot_case_series.py'))
manifest_path.write_text(json.dumps(figure,ensure_ascii=False,indent=2),encoding='utf8')
document=ROOT/'局部柔性下降与相反亲和结果_案例分析_三线表_20261022.docx'
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with zipfile.ZipFile(document) as package:
    docxml=ET.fromstring(package.read('word/document.xml'))
    stylexml=ET.fromstring(package.read('word/styles.xml'))
    assert not docxml.findall('.//w:pBdr',ns)
    assert not stylexml.findall('.//w:pBdr',ns)
    assert len(docxml.findall('.//w:tbl',ns))==1
    assert len(docxml.findall('.//w:drawing',ns))==1
    text=''.join(e.text or '' for e in docxml.findall('.//w:t',ns))
    assert all(t in text for t in ['60,000','0.1305','0.6605','0.7281','0.821','0.016','58.7','未修改论文'])
render=ROOT/'document_QA/case_analysis_render.pdf'
rendered=PdfReader(render)
assert len(rendered.pages)==4
assert all((ROOT/f'document_QA/page-{i}.png').exists() for i in range(1,5))
docqa=dict(status='passed',document_sha256=sha(document),pages=4,all_pages_visually_inspected=True,three_line_native_tables=1,title_border_removed=True,english_figure_embedded=True,no_clipping_or_overlap=True,renderer='Packaged render_docx.py attempted; failed because soffice.exe unavailable. Native Word read-only PDF export + bundled PDFium PNG rendering completed.',main_manuscript_modified=False)
(ROOT/'document_QA/visual_and_structure_audit.json').write_text(json.dumps(docqa,ensure_ascii=False,indent=2),encoding='utf8')
paper=ROOT.parent.parent/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'
assert sha(paper)=='5f79ddf7d657dc9b782c9fb432e373a1e7c170c1c88d019e528fbf73c366c222'
audit=dict(status='complete',date_local='2026-10-22',finished_utc=datetime.now(timezone.utc).isoformat(),exploratory=True,primary_case_group='MEL5 four observed complex states; fixed3contrasts',scientific_verdict=review['scientific_verdict'],independent_numeric_review=review['status'],source_values_checked=review['source_values_checked'],main_manuscript_sha256=sha(paper),main_manuscript_modified=False,new_MD=0,new_training=0,new_inference=0,figure_QA='passed',document_QA='passed',files=[])
files=[p for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and 'document_QA' not in p.parts and p.name!='delivery_audit.json']
for p in files:audit['files'].append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
(ROOT/'delivery_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
files.append(ROOT/'delivery_audit.json')
files.append(ROOT/'document_QA/visual_and_structure_audit.json')
archive=ROOT.parent/'相反亲和结果_完整案例分析_20261022.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as package:
    for p in files:package.write(p,str(Path(ROOT.name)/p.relative_to(ROOT)))
with zipfile.ZipFile(archive) as package:
    assert package.testzip() is None
    for p in files:
        key=str(Path(ROOT.name)/p.relative_to(ROOT)).replace('\\','/')
        assert hashlib.sha256(package.read(key)).hexdigest()==sha(p)
print(json.dumps(dict(status='complete',document=str(document),pages=4,archive=str(archive),archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),main_manuscript_unchanged=True),ensure_ascii=False,indent=2))
