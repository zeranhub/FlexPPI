from pathlib import Path
import hashlib,json,zipfile
from lxml import etree
from docx import Document

OUT=Path('E:/Project/FlexPPI/experiments/_sources/completed_results_20261012')
DOC=Path('E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v19_20261012_中文修订_补充实验完成_三线表.docx')
document=Document(DOC)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
full=lambda element:''.join(element.xpath('.//w:t/text() | .//m:t/text()',namespaces=ns))
paragraphs=[]
for index,paragraph in enumerate(document.paragraphs):
    text=''.join(paragraph._p.xpath('.//w:t/text()'))
    math=''.join(paragraph._p.xpath('.//m:t/text()'))
    paragraphs.append(dict(index=index,text=text,math=math,style=paragraph.style.name))
tables=[]
for index,table in enumerate(document.tables):
    rows=[[cell.text for cell in row.cells] for row in table.rows]
    tables.append(dict(index=index,rows=rows))
with zipfile.ZipFile(DOC) as archive:
    xml=etree.fromstring(archive.read('word/document.xml'))
    body=xml.find('w:body',ns)
    blocks=[]
    p_index=0;t_index=0
    for node in body:
        if node.tag.endswith('}p'):
            blocks.append(dict(type='paragraph',index=p_index,text=full(node)))
            p_index+=1
        elif node.tag.endswith('}tbl'):
            blocks.append(dict(type='table',index=t_index,rows=tables[t_index]['rows']))
            t_index+=1
    media={name:hashlib.sha256(archive.read(name)).hexdigest() for name in archive.namelist() if name.startswith('word/media/')}
report=dict(docx=str(DOC),sha256=hashlib.sha256(DOC.read_bytes()).hexdigest(),paragraphs=paragraphs,tables=tables,blocks=blocks,media_sha256=media)
(OUT/'v19_content_extracted.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
lines=[]
for block in blocks:
    if block['type']=='paragraph':
        if block['text']:lines.append('P'+str(block['index'])+': '+block['text'])
    else:
        for index,row in enumerate(block['rows']):lines.append('T'+str(block['index'])+'R'+str(index)+': '+' | '.join(row))
(OUT/'v19_content_extracted.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(dict(sha256=report['sha256'],paragraphs=len(paragraphs),tables=len(tables),media=len(media)),indent=2))
