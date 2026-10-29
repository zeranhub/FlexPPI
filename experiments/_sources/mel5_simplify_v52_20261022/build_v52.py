"""Remove the rejected local molecular close-ups and synchronize the v52 manuscript."""
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import re
import zipfile
from lxml import etree
from PIL import Image

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent.parent
QA = ROOT / 'manuscript_QA'; QA.mkdir(exist_ok=True)
SOURCE = BASE / 'paper_v51_20261022_中文修订_MEL5跨界面残基机制与图4重绘_三线表.docx'
OUT = BASE / 'paper_v52_20261022_中文修订_MEL5图4去除局部结构图_三线表.docx'
FIG = ROOT / 'figure4/fig4_mel5_compact.png'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE) == '2f99757f7063cc6506b25720dbe20a97197ac9d1e55eef0926d8837187f75e8d'
NS = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
      'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
with zipfile.ZipFile(SOURCE) as z:
    infos = z.infolist(); parts = {p.filename:z.read(p.filename) for p in infos}
root = etree.fromstring(parts['word/document.xml'])
body = root.find('w:body', NS)
original_nodes = [etree.tostring(p) for p in body]
pars = [p for p in body if p.tag == '{'+NS['w']+'}p']
txt = lambda p: ''.join(p.xpath('.//w:t/text()',namespaces=NS))
old = {i:txt(p) for i,p in enumerate(pars)}
assert old[89].startswith('3.4 MEL5') and old[95].startswith('图4：')
new = {}
opening='将两个复合物以HLA对齐并从同一视角观察，可直接比较受体改造前后的局部配对关系（图4a）。'
assert old[92].count(opening)==1
new[92]=old[92].replace(opening,'为解释这两个位于界面两侧的模型响应，进一步检查了βCDR2附近的局部接触。')
method='跨界面局部结构沿用来源于3HG1与4JFF的MD准备起始坐标，以HLA第1–180位共同Cα原子对齐两个复合物，并用相同视角和尺度显示。图中标注Q72 OE1–β53主链N和Q72 NE2–βGly54主链O的原子间距。'
assert old[52].count(method)==1
new[52]=old[52].replace(method,'局部原子间距由来源于3HG1与4JFF的MD准备起始坐标计算，分别测量Q72 OE1–β53主链N及Q72 NE2–βGly54主链O的距离。')
new[95]='图4：FlexPPI定位MEL5亲和成熟后的界面双动态响应。（a）复合物概览、文献[28]的SPR亲和力，以及界面两侧HLA Q72的Hχ与TCR βGly54的RMSF两态模型读数。两种TCR结合相同的HLA-A*02:01与改良MART-1十肽ELAGIGILTV。（b）六个完整CDR中全部45个对齐位点的两态模型RMSF，识别环之间分隔，突出βCDR2及Gly54。橙色为MEL5，紫色为α24β17；编号统一到MEL5，配对与共同编号规则见方法。'
changes=[]
for i,text in sorted(new.items()):
    p=pars[i]
    assert text!=old[i]
    assert not p.xpath('.//m:oMath|.//w:fldChar|.//w:hyperlink',namespaces=NS)
    ts=p.xpath('.//w:t',namespaces=NS); assert ts
    ts[0].text=text
    for t in ts[1:]: t.text=''
    changes.append(dict(paragraph_index=i,body_index=list(body).index(p),before=old[i],after=text))
fig=pars[94]
rid=fig.xpath('.//a:blip',namespaces=NS)[0].get('{'+NS['r']+'}embed')
rels=etree.fromstring(parts['word/_rels/document.xml.rels'])
media='word/'+{e.get('Id'):e.get('Target') for e in rels}[rid]
assert media=='word/media/image41.png'
inline=fig.xpath('.//wp:inline',namespaces=NS)[0]
width=round(6.6*914400)
with Image.open(FIG) as im:height=round(width*im.height/im.width)
for e in inline.xpath('./wp:extent|.//a:xfrm/a:ext',namespaces=NS):
    e.set('cx',str(width)); e.set('cy',str(height))
for prop in inline.xpath('./wp:docPr',namespaces=NS):
    prop.set('name','Figure 4 MEL5 compact paired model readouts')
    prop.set('descr','Two panels: complex locator, published SPR affinity and compact paired Q72 Hchi / betaGly54 RMSF values; all45CDR model RMSF curves. Both local molecular close-ups removed.')
alltext='\n'.join(txt(p) for p in body)
case='\n'.join(txt(pars[i]) for i in range(89,94))
assert not any(t in alltext for t in ['图4c','图4d','（图4c）'])
assert not any(t in case for t in ['预测失败','10个位点','G4A','I7A','移除TCR','2.141','2.245'])
assert all(t in case for t in ['5.85','2.94','10.58','3.12','1.792','1.375','1.118','0.747','Q72','βGly54','23个'])
counts=dict(native_tables=len(root.xpath('.//w:tbl',namespaces=NS)),
            native_equations=len(root.xpath('.//m:oMath',namespaces=NS)),
            images=len(root.xpath('.//a:blip',namespaces=NS)),
            references=len(re.findall(r'^\[\d+\]',alltext,re.M)))
assert counts==dict(native_tables=10,native_equations=36,images=6,references=28)
updates={media:FIG.read_bytes(),'word/document.xml':etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)}
with zipfile.ZipFile(OUT,'w') as z:
    for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
editable={c['body_index'] for c in changes}|{list(body).index(fig)}
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None
    changed=[p for p in parts if z.read(p)!=parts[p]]
    assert set(changed)==set(updates)
    reread=etree.fromstring(z.read('word/document.xml')).find('w:body',NS)
    assert len(reread)==len(original_nodes)
    for i,p in enumerate(reread):
        if i not in editable: assert etree.tostring(p)==original_nodes[i],i
(ROOT/'figure4/caption_zh.md').write_text(new[95]+'\n',encoding='utf-8')
(ROOT/'prose_v52.json').write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8')
(QA/'v52_content.txt').write_text('\n'.join(f'P{i}: {txt(p)}' for i,p in enumerate(pars)),encoding='utf-8')
audit=dict(status='built_pending_render_QA',source=str(SOURCE),source_sha256=sha(SOURCE),
           document=str(OUT),document_sha256=sha(OUT),figure=str(FIG),figure_sha256=sha(FIG),
           changed_zip_parts=changed,edited_body_indices=sorted(editable),prose_changes=changes,
           body_nodes_unchanged_outside_authorized_indices=True,source_preserved=True,
           figure_dimensions_EMU=dict(cx=width,cy=height),new_training=0,new_inference=0,new_MD=0,**counts)
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
