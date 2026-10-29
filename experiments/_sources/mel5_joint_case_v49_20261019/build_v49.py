"""Replace the authorized TCR case, preserving every unrelated native Word part."""
from pathlib import Path
from copy import deepcopy
import hashlib,json,re,zipfile
from lxml import etree
from PIL import Image

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v48_20261019_中文修订_KRAS亲和力机制与绝对RMSF_三线表.docx'
OUT=BASE/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'
FIG=ROOT/'figure4/fig4_mel5_joint.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='952c789ef6271f88f7ac619ac78291f71536c53f3a7bdeddb00eaebdcf003726'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
 'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
 'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
 'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
 'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
with zipfile.ZipFile(SOURCE) as z:
    infos=z.infolist();parts={p.filename:z.read(p.filename) for p in infos}
root=etree.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
original_body=list(body)
txt=lambda n:''.join(n.xpath('.//w:t/text()',namespaces=NS))
pars=[p for p in body if p.tag=='{'+NS['w']+'}p']
old={i:txt(p) for i,p in enumerate(pars)}
contract=json.loads((ROOT/'prose_v49.json').read_text(encoding='utf8'))
new={int(i):t for i,t in contract['replacements'].items()}
for i,(a,b) in contract['substring_replacements'].items():
    i=int(i);assert old[i].count(a)==1,(i,a)
    new[i]=old[i].replace(a,b)
new[144]=(ROOT/'literature_review/reference28.txt').read_text(encoding='utf8').strip()
assert new[144].startswith('[28]') and 'Madura' in new[144] and 'M113.464560' in new[144]
changes=[]
for i,text in sorted(new.items()):
    p=pars[i];assert old[i]!=text,i
    assert not p.xpath('.//m:oMath|.//w:fldChar|.//w:hyperlink',namespaces=NS),i
    ts=p.xpath('.//w:t',namespaces=NS);assert ts
    ts[0].text=text
    for t in ts[1:]:t.text=''
    changes.append({'source_paragraph_index':i,'source_body_index':original_body.index(p),
                    'before':old[i],'after':text,'element':p})

# One additional connected paragraph develops the residue-level mechanism.
p=etree.Element('{'+NS['w']+'}p')
pr=pars[94].find('w:pPr',NS)
if pr is not None:p.append(deepcopy(pr))
r=etree.SubElement(p,'{'+NS['w']+'}r')
rpr=pars[94].find('w:r/w:rPr',NS)
if rpr is not None:r.append(deepcopy(rpr))
etree.SubElement(r,'{'+NS['w']+'}t').text=contract['insert_after_94']
pars[94].addnext(p)
inserted=[{'body_index':list(body).index(p),'after_source_paragraph_index':94,'text':txt(p)}]

cap=pars[92];fig=cap.getprevious()
blips=fig.xpath('.//a:blip',namespaces=NS);assert len(blips)==1
rid=blips[0].get('{'+NS['r']+'}embed')
rels=etree.fromstring(parts['word/_rels/document.xml.rels'])
media='word/'+{n.get('Id'):n.get('Target') for n in rels}[rid]
assert media=='word/media/image41.png',media
inline=blips[0].xpath('ancestor::wp:inline',namespaces=NS)[0]
width=round(6.6*914400)
with Image.open(FIG) as im:height=round(width*im.height/im.width)
for e in inline.xpath('./wp:extent|.//a:xfrm/a:ext',namespaces=NS):
    e.set('cx',str(width));e.set('cy',str(height))
for pr in inline.xpath('./wp:docPr',namespaces=NS):
    pr.set('name','Figure 4 MEL5 affinity maturation and paired dynamic readouts')
    pr.set('descr','MEL5 and alpha24beta17 structures and SPR affinity; absolute model RMSF across all 45 CDR positions; absolute Hchi at all 23 valid pMHC interface positions; existing MD contact occupancy at Q72-beta53/beta54')
def prop(node,name,value):
    pr=node.find('w:pPr',NS)
    if pr is None:pr=etree.Element('{'+NS['w']+'}pPr');node.insert(0,pr)
    e=pr.find('w:'+name,NS)
    if e is None:e=etree.SubElement(pr,'{'+NS['w']+'}'+name)
    e.set('{'+NS['w']+'}val',str(value))
prop(fig,'keepNext',1);prop(cap,'keepLines',1);prop(cap,'keepNext',0)

# Keep the complete argument together before its full-page figure.  Placing a
# 170-mm figure immediately after the introductory paragraph stranded most of
# that page as white space in the first Word render.
p.addnext(fig)
fig.addnext(cap)
inserted[0]['body_index']=list(body).index(p)

alltext='\n'.join(txt(n) for n in body)
for token in ['KRAS','7OW5','7OW6','JDIa41b1','Poole','M113.464560']:
    assert (token in alltext)==(token=='M113.464560'),token
assert len(root.xpath('.//w:tbl',namespaces=NS))==10
assert len(root.xpath('.//m:oMath',namespaces=NS))==36
assert len(root.xpath('.//a:blip',namespaces=NS))==6
assert len(re.findall(r'^\[\d+\]',alltext,re.M))==28
updates={media:FIG.read_bytes(),'word/document.xml':etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)}
with zipfile.ZipFile(OUT,'w') as z:
    for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
with zipfile.ZipFile(OUT) as z:assert z.testzip() is None
for c in changes:c['body_index']=list(body).index(c.pop('element'))
(ROOT/'figure4/caption_zh.md').write_text(new[92]+'\n',encoding='utf8')
finalpars=[n for n in body if n.tag=='{'+NS['w']+'}p']
(QA/'v49_content.txt').write_text('\n'.join(f'P{i}: {txt(n)}' for i,n in enumerate(finalpars)),encoding='utf8')
audit={'status':'built_pending_render_QA','source':str(SOURCE),'source_sha256':sha(SOURCE),
 'document':str(OUT),'document_sha256':sha(OUT),'figure4':str(FIG),'figure4_sha256':sha(FIG),
 'prose_contract_sha256':sha(ROOT/'prose_v49.json'),
 'caption_body_index':list(body).index(cap),'prose_changes':changes,'inserted_body_nodes':inserted,
 'layout_note':'Figure 4 and its caption moved together after all four Section 3.4 prose paragraphs to avoid a mostly empty preceding page.',
 'source_to_final_body_index':{str(i):list(body).index(n) for i,n in enumerate(original_body)},
 'editable_source_body_indices':[original_body.index(fig)]+[c['source_body_index'] for c in changes],
 'changed_zip_parts':list(updates),'width_EMU':width,'height_EMU':height,
 'native_tables':10,'native_equations':36,'references':28,'references_intentionally_replaced':[28],
 'main_figures':5,'supplementary_figures':1,'new_MD':0,'new_training':0,'new_inference':0,
 'existing_MD_reanalysis':{'trajectories':6,'frames':6012,'source':'tcr_joint_case_search_20261019/residue_story/trajectory_link'}}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
