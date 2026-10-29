from pathlib import Path
from zipfile import ZipFile
from lxml import etree
from pypdf import PdfReader
from PIL import Image
import json,hashlib,re,sys

HERE=Path(__file__).resolve().parent
BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
DOC=BASE/'paper_v93_20261026_中文修订_双案例与构象占据整合.docx'
SOURCE=BASE/'paper_v90_20261026_中文修订_图3数据呈现与机制推断.docx'
RENDER=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'render/release_05'
PDF=RENDER/(DOC.stem+'.pdf')
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
sha=lambda data:hashlib.sha256(data).hexdigest()
sig=lambda x:etree.tostring(x,method='c14n')
with ZipFile(DOC) as z:
    root=etree.fromstring(z.read('word/document.xml'));body=root.find('w:body',ns)
    ps=body.findall('w:p',ns)
    text=[''.join(p.xpath('.//w:t/text()',namespaces=ns)) for p in ps]
    full='\n'.join(text)
    media={n:sha(z.read(n)) for n in z.namelist() if n.startswith('word/media/')}
    finalparts={n:z.read(n) for n in z.namelist()}
with ZipFile(SOURCE) as z:
    sourceparts={n:z.read(n) for n in z.namelist()}
    sr=etree.fromstring(z.read('word/document.xml'))
    source_tables=[sig(t) for t in sr.find('w:body',ns).findall('w:tbl',ns)]
assert source_tables==[sig(t) for t in body.findall('w:tbl',ns)]
changed=[n for n,d in finalparts.items() if sourceparts.get(n)!=d]
assert set(changed)=={'word/document.xml','word/media/image5.png','word/media/image40.png'}
unchanged_media=[n for n in media if sourceparts[n]==finalparts[n]]
assert len(unchanged_media)==44
maincaps=[t for t in text if re.match(r'^图\d+[：:]',t)]
suppcaps=[t for t in text if re.match(r'^图S\d+[：:]',t)]
assert [re.match(r'^图(\d+)',t).group(1) for t in maincaps]==['1','2','3','4']
assert [re.match(r'^图S(\d+)',t).group(1) for t in suppcaps]==['1','2','3','4']
assert not re.search(r'图\s*5',full)
assert len(body.xpath('.//a:blip',namespaces=ns))==8
headings=[t for t in text if re.match(r'^3\.[345] ',t)]
assert headings==[
 '3.3 典型识别体系中的主链约束与侧链环境差异',
 '3.4 侧链预测与构象占据及相对自由能角区的对应',
 '3.5 TCR–pMHC识别中的界面柔性与亲和力差异']
reader=PdfReader(PDF);pagetext=[p.extract_text() for p in reader.pages]
pngs=[RENDER/f'page-{i+1}.png' for i in range(len(reader.pages))]
assert all(p.exists() for p in pngs)
norm=lambda s:re.sub(r'\s+','',s)
pdftext=''.join(norm(t) for t in pagetext)
assert '图5' not in pdftext
assert '\ufffd' not in pdftext
anchors={s:[i+1 for i,t in enumerate(pagetext) if norm(s) in norm(t)] for s in headings+[maincaps[2][:35],maincaps[3][:30],suppcaps[-1][:35]]}
textfile=HERE/'final_render_text.txt'
textfile.write_text('\n\n'.join(f'PAGE {i+1}\n{t}' for i,t in enumerate(pagetext)),encoding='utf-8')
audit={
 'docx':str(DOC),'docx_sha256':sha(DOC.read_bytes()),'source_docx_sha256':sha(SOURCE.read_bytes()),
 'render_folder':str(RENDER),'render_pdf':str(PDF),'render_pdf_sha256':sha(PDF.read_bytes()),
 'page_count':len(reader.pages),'png_count':len(pngs),'page_image_dimensions':sorted(set(Image.open(p).size for p in pngs)),
 'source_table_xml_identical':True,'source_table_count':len(source_tables),'unchanged_media_count':len(unchanged_media),
 'changed_zip_parts':changed,'active_inline_pictures':8,'main_figure_numbers':[1,2,3,4],'supplement_figure_numbers':[1,2,3,4],
 'main_figure5_references':0,'section_order':headings,'heading_and_caption_start_pages':anchors,
 'figure3_sha256':media['word/media/image5.png'],'supplement4_sha256':media['word/media/image40.png'],
 'replacement_glyph_found':False,
 'visual_review':{'status':'passed' if sha(DOC.read_bytes())=='6c603bfff03d033610d228854c05c0cb65bce4ceee3e9f9dae2bf8b3433aab44' else 'pending fresh human review',
 'root_pages':'Root actually viewed release04 pages1-12 and release05 changed pages1,2,12,13,29; same-hash pages inherit review.',
 'independent_document_agent_pages':'Document agent viewed release04 pages11-29. Release05 changed only1,2,12,13,29; root reviewed all changed pages.',
 'scope':'No new clipping, overlap, broken tables or detached captions. Source figure-page whitespace outside revised sections retained.'},
 'page_png_sha256':{p.name:sha(p.read_bytes()) for p in pngs}}
(HERE/'v93_final_document_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['docx_sha256','page_count','changed_zip_parts','unchanged_media_count','heading_and_caption_start_pages']},ensure_ascii=False))
