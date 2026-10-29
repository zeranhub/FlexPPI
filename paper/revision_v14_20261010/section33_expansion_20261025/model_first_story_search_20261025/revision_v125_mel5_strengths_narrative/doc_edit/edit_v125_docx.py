from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree
from hashlib import sha256
from xml.sax.saxutils import escape
import json,re
ROOT=Path(__file__).resolve().parent.parent
PLAN=json.loads((ROOT/'narrative_changes.json').read_text(encoding='utf-8'))
SOURCE=Path(PLAN['source']);TARGET=Path(PLAN['target']);sha=lambda b:sha256(b).hexdigest()
with ZipFile(SOURCE) as z:items=z.infolist();parts={i.filename:z.read(i.filename) for i in items}
original=parts['word/document.xml'];ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};root=etree.fromstring(original);ps=root.findall('.//w:body/w:p',ns)
pattern=re.compile(rb'<w:p(?:\s[^>]*)?>.*?</w:p>',re.S);blocks=list(pattern.finditer(original))
replacements=[];draft=[]
for index,entry in PLAN['paragraph_replacements'].items():
 old=entry['before'];new=entry['after'];assert ''.join(ps[int(index)].xpath('.//w:t/text()',namespaces=ns))==old
 ob=escape(old).encode('utf-8');nb=escape(new).encode('utf-8');found=[m for m in blocks if ob in m.group()];assert len(found)==1,(index,len(found));m=found[0];assert m.group().count(ob)==1
 replacements.append((m.start(),m.end(),m.group().replace(ob,nb)));draft.append((int(index),new))
for index,entries in PLAN['paragraphs_inserted_after'].items():
 anchor=''.join(ps[int(index)].xpath('.//w:t/text()',namespaces=ns));ab=escape(anchor).encode('utf-8');found=[m for m in blocks if ab in m.group()];assert len(found)==1;m=found[0]
 chunks=[]
 for j,entry in enumerate(entries):
  ti=int(entry['template']);template=''.join(ps[ti].xpath('.//w:t/text()',namespaces=ns));tb=escape(template).encode('utf-8');f=[x for x in blocks if tb in x.group()];assert len(f)==1
  chunk=f[0].group();assert chunk.count(tb)==1
  chunk=chunk.replace(tb,escape(entry['text']).encode('utf-8'))
  for attr in [b'w14:paraId',b'w14:textId']:
   new_id=sha((index+':'+str(j)+':'+attr.decode()).encode())[:8].upper().encode()
   chunk=re.sub(attr+rb'="[^"]*"',attr+b'="'+new_id+b'"',chunk)
  chunks.append(chunk);draft.append((int(index)+.1+.01*j,entry['text']))
 existing=next((k for k,r in enumerate(replacements) if r[0]==m.start()),None)
 if existing is not None:
  a,b,payload=replacements[existing];replacements[existing]=(a,b,payload+b''.join(chunks))
 else:replacements.append((m.start(),m.end(),m.group()+b''.join(chunks)))
xml=original
for start,end,new in sorted(replacements,reverse=True):xml=xml[:start]+new+xml[end:]
newroot=etree.fromstring(xml);nps=newroot.findall('.//w:body/w:p',ns)
assert len(nps)==len(ps)+sum(len(e) for e in PLAN['paragraphs_inserted_after'].values())
# Exact expected body paragraph sequence and untouched non-paragraph blocks.
expected=[]
for i,p in enumerate(ps):
 t=''.join(p.xpath('.//w:t/text()',namespaces=ns));expected.append(PLAN['paragraph_replacements'].get(str(i),{}).get('after',t))
 expected.extend(e['text'] for e in PLAN['paragraphs_inserted_after'].get(str(i),[]))
assert [''.join(p.xpath('.//w:t/text()',namespaces=ns)) for p in nps]==expected
for query in ['//w:fldSimple','//w:instrText','//w:drawing','//w:tbl','//w:sectPr','//w:footnoteReference','//w:commentRangeStart','//w:commentRangeEnd']:
 a=root.xpath(query,namespaces=ns);b=newroot.xpath(query,namespaces=ns);assert len(a)==len(b),(query,len(a),len(b))
 assert all(etree.tostring(x)==etree.tostring(y) for x,y in zip(a,b)),query
for term in ['遗憾的是','仅','仍明显落后','效果有限','试过','失败了','模型揭示氢键','模型预测氢键']:
 assert not any(term in text for _,text in draft),term
parts['word/document.xml']=xml
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
 for i in items:z.writestr(i,parts[i.filename])
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
 assert a.namelist()==b.namelist() and b.testzip() is None
 changed=[n for n in a.namelist() if a.read(n)!=b.read(n)];assert changed==['word/document.xml'],changed
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':sha(SOURCE.read_bytes()),'target_sha256':sha(TARGET.read_bytes()),'changed_zip_parts':changed,'body_paragraph_count_before':len(ps),'body_paragraph_count_after':len(nps),'replaced_paragraphs':len(PLAN['paragraph_replacements']),'inserted_paragraphs':sum(len(e) for e in PLAN['paragraphs_inserted_after'].values()),'all_figures_media_tables_fields_comments_and_other_parts_identical':True,'all_edited_text_passes_user_expression_scan':True,'style_and_drawing_XML_preserved':True,'validation':'Pending rendering and visual verification.'}
(ROOT/'doc_edit/document_change_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
text='\n\n'.join(t for i,t in sorted(draft) if 78<=i<85)
(ROOT/'MEL5_正文改写.md').write_text(text+'\n',encoding='utf-8')
(ROOT/'modified_text_all.md').write_text('\n\n'.join(str(i)+'\n'+t for i,t in sorted(draft))+'\n',encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
