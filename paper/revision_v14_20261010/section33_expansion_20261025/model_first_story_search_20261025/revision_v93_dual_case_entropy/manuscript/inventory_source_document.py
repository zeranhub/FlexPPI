from pathlib import Path
from zipfile import ZipFile
import json, re, hashlib
from lxml import etree
from docx import Document
from pypdf import PdfReader

BASE = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
SOURCE = BASE / 'paper_v90_20261026_中文修订_图3数据呈现与机制推断.docx'
OUT = BASE / r'section33_expansion_20261025\model_first_story_search_20261025\revision_v93_dual_case_entropy\manuscript'
OUT.mkdir(parents=True, exist_ok=True)
NS = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main', 'a':'http://schemas.openxmlformats.org/drawingml/2006/main', 'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships', 'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}

def txt(el):
    return ''.join(el.xpath('.//w:t/text()', namespaces=NS))

doc=Document(SOURCE)
with ZipFile(SOURCE) as z:
    root=etree.fromstring(z.read('word/document.xml'))
    relroot=etree.fromstring(z.read('word/_rels/document.xml.rels'))
    rels={r.get('Id'):{'target':r.get('Target'),'type':r.get('Type'),'mode':r.get('TargetMode')} for r in relroot}
    body=root.find('w:body',NS)
    paragraphs=[]; media=[]; tables=[]; pidx=0
    for bidx,el in enumerate(body):
        tag=etree.QName(el).localname
        if tag=='p':
            p=doc.paragraphs[pidx]
            images=[]
            for b in el.xpath('.//a:blip', namespaces=NS):
                rid=b.get('{'+NS['r']+'}embed') or b.get('{'+NS['r']+'}link')
                r=rels.get(rid,{})
                images.append({'relationship_id':rid,**r})
            wext=el.xpath('.//wp:extent',namespaces=NS)
            paragraph={'paragraph_index':pidx,'paragraph_number':pidx+1,'body_child_index':bidx,
                'xpath':f'/w:document/w:body/w:p[{pidx+1}]','style_name':p.style.name,'text':txt(el),
                'images':images,'image_extents_emu':[dict(x.attrib) for x in wext],
                'page_break_before':p.paragraph_format.page_break_before,
                'keep_with_next':p.paragraph_format.keep_with_next,
                'explicit_page_breaks':len(el.xpath('.//w:br[@w:type="page"]',namespaces=NS)),
                'field_instructions':el.xpath('.//w:instrText/text()',namespaces=NS),
                'bookmark_names':el.xpath('.//w:bookmarkStart/@w:name',namespaces=NS),
                'raw_xml_sha256':hashlib.sha256(etree.tostring(el)).hexdigest()}
            paragraphs.append(paragraph);pidx+=1
        elif tag=='tbl':
            tables.append({'table_index':len(tables),'body_child_index':bidx,'text':txt(el),'n_rows':len(el.findall('w:tr',NS))})
    for name in z.namelist():
        if name.startswith('word/media/'):
            data=z.read(name)
            target=OUT/'source_media'/Path(name).name
            target.parent.mkdir(exist_ok=True)
            target.write_bytes(data)
            media.append({'zip_member':name,'size':len(data),'sha256':hashlib.sha256(data).hexdigest(),'extracted_to':str(target)})
    all_fields=root.xpath('//w:instrText/text()',namespaces=NS)
    comments=[n for n in z.namelist() if re.search(r'(comment|footnote|endnote)',n)]

pattern=re.compile(r'(?:图\s*S?\d+(?:[a-z](?:[—–\-、,，][a-z])?)?|(?:Figure|Fig\.?)\s*S?\d+|3\.[345]|补充材料|补充图|Supplement)',re.I)
refs=[]
for p in paragraphs:
    matches=[m.group() for m in pattern.finditer(p['text'])]
    if matches: refs.append({'paragraph_index':p['paragraph_index'],'body_child_index':p['body_child_index'],'text':p['text'],'matches':matches})
inv={'source':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'index_convention':'zero-based paragraph_index=python-docx Document.paragraphs index; body_child_index includes tables and sectPr; xpath p ordinal is one-based',
    'paragraph_count':len(paragraphs),'body_child_count':len(body),'table_count':len(tables),'image_count':len(media),
    'paragraphs':paragraphs,'tables':tables,'media':media,'cross_reference_occurrences':refs,
    'field_instructions':all_fields,'comment_or_note_parts':comments,
    'sections':[{'page_width':s.page_width,'page_height':s.page_height,'left_margin':s.left_margin,'right_margin':s.right_margin,'top_margin':s.top_margin,'bottom_margin':s.bottom_margin} for s in doc.sections]}
baseline_pdf = BASE / r'section33_expansion_20261025\model_first_story_search_20261025\revision_v90_data_led_story\render\release_01\paper_v90_20261026_中文修订_图3数据呈现与机制推断.pdf'
if baseline_pdf.exists():
    pdf=PdfReader(baseline_pdf)
    pagetext=[p.extract_text() for p in pdf.pages]
    clean=lambda s:re.sub(r'\s+','',s)
    cleanpages=[clean(t) for t in pagetext]
    for p in paragraphs:
        if len(p['text'])>=12:
            needle=clean(p['text'])[:36]
            p['baseline_pdf_start_pages']=[i+1 for i,s in enumerate(cleanpages) if needle in s]
    inv['baseline_render']={'pdf':str(baseline_pdf),'pages':len(pdf.pages),'sha256':hashlib.sha256(baseline_pdf.read_bytes()).hexdigest(),'note':'Existing v90 Word read-only render; page mapping from normalized first 36 text characters, pages are baseline only.'}
    (OUT/'baseline_pdf_text.txt').write_text('\n\n'.join(f'PAGE {i+1}\n{t}' for i,t in enumerate(pagetext)),encoding='utf-8')
(OUT/'document_inventory.json').write_text(json.dumps(inv,ensure_ascii=False,indent=2),encoding='utf-8')
lines=[]
for p in paragraphs:
    lines.append(f"[{p['paragraph_index']:03d} | body {p['body_child_index']:03d} | {p['style_name']}] {p['text']}")
    for im in p['images']:lines.append('  IMAGE '+json.dumps(im,ensure_ascii=False))
(OUT/'source_paragraphs_with_indices.txt').write_text('\n\n'.join(lines),encoding='utf-8')
print(json.dumps({k:inv[k] for k in ('source_sha256','paragraph_count','body_child_count','table_count','image_count')},ensure_ascii=False))
print('Inventory and extracted media written.')
