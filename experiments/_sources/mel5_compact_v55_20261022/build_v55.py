"""Replace Figure4 with the user-selected four-panel compact composite; synchronize its prose."""
from pathlib import Path
import hashlib,json,re,zipfile
from lxml import etree
from PIL import Image

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx'
OUT=BASE/'paper_v55_20261022_中文修订_MEL5精简联合图4_三线表.docx'
FIG=ROOT/'figure4/fig4_mel5_compact_joint.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='90b66a4029e22b0f67e52d3176839e547a8265690e5ff21d385ab7b5fecd1027'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
with zipfile.ZipFile(SOURCE) as z:
    infos=z.infolist();parts={i.filename:z.read(i.filename) for i in infos}
root=etree.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
original_nodes=[etree.tostring(p) for p in body]
pars=[p for p in body if p.tag=='{'+NS['w']+'}p']
txt=lambda p:''.join(p.xpath('.//w:t/text()',namespaces=NS))
old={i:txt(p) for i,p in enumerate(pars)}
assert old[89].startswith('3.4 MEL5') and old[95].startswith('图4：')
new={}
new[91]=old[91].replace('（图4d）','（图4d下部）')
new[93]=old[93].replace('（图4e）','（图4d上部）').replace('（图4f）','（图4d下部）')
new[95]='图4：MEL5体系中受体亲和成熟与肽替换的残基动态响应。（a）成熟复合物概览及Q72–β53/β54局部结构；虚线为Q72 OE1–β53 N和NE2–β54 O的原子间距。（b）四态全部10个肽位置的预测RMSF；P4为识别区域，P10为C端HLA锚定位点。（c）亲和成熟前后全部45个CDR位点的RMSF。（d）三种对照的联合动态响应：上部为CDR ΔRMSF，下部为pMHC界面ΔHχ；两项分别标示单位和色标。三行依次为4JFF−3HG1、4JFD−4JFF和4JFE−4JFF，蓝、红表示降低、升高，灰格表示Hχ无定义或氨基酸类型不同。界面为四态最近Cα距离≤8 Å的并集，共30位，三行Hχ有效数为23、23、22。TCR采用MEL5共同编号；各态均由同一模型对已有复合物结构预测，状态定义见方法。'
changes=[]
for i,s in sorted(new.items()):
    p=pars[i];assert s!=old[i]
    assert not p.xpath('.//m:oMath|.//w:fldChar|.//w:hyperlink',namespaces=NS)
    ts=p.xpath('.//w:t',namespaces=NS);assert ts
    ts[0].text=s
    for t in ts[1:]:t.text=''
    changes.append(dict(paragraph_index=i,body_index=list(body).index(p),before=old[i],after=s))
# The full composite and its complete caption share one manuscript page.
cap=pars[95]
for r in cap.findall('w:r',NS):
    rp=r.find('w:rPr',NS)
    if rp is None:rp=etree.SubElement(r,'{'+NS['w']+'}rPr')
    for tag in ['sz','szCs']:
        e=rp.find('w:'+tag,NS)
        if e is None:e=etree.SubElement(rp,'{'+NS['w']+'}'+tag)
        e.set('{'+NS['w']+'}val','21')
pp=cap.find('w:pPr',NS)
if pp is None:pp=etree.Element('{'+NS['w']+'}pPr');cap.insert(0,pp)
sp=pp.find('w:spacing',NS)
if sp is None:sp=etree.SubElement(pp,'{'+NS['w']+'}spacing')
sp.set('{'+NS['w']+'}line','260');sp.set('{'+NS['w']+'}lineRule','auto')
fig=pars[94]
rid=fig.xpath('.//a:blip',namespaces=NS)[0].get('{'+NS['r']+'}embed')
rels=etree.fromstring(parts['word/_rels/document.xml.rels'])
media='word/'+{e.get('Id'):e.get('Target') for e in rels}[rid]
assert media=='word/media/image41.png'
inline=fig.xpath('.//wp:inline',namespaces=NS)[0]
width=round(6.6*914400)
with Image.open(FIG) as im:height=round(width*im.height/im.width)
for e in inline.xpath('./wp:extent|.//a:xfrm/a:ext',namespaces=NS):e.set('cx',str(width));e.set('cy',str(height))
for e in inline.xpath('./wp:docPr',namespaces=NS):
    e.set('name','Figure 4 MEL5 compact four-panel joint response')
    e.set('descr','Four panels: complex and Q72-beta53/54 atom-pair distances; four-state peptide RMSF; all45 CDR absolute RMSF; unified upper CDR RMSF and lower pMHC Hchi heatmaps for three contrasts. No ranked bar panel. Orange MEL5, purple alpha24beta17, blue G4A, sage I7A. Blue-white-rose signed changes; grey ineligible Hchi.')
alltext='\n'.join(txt(p) for p in body)
case='\n'.join(txt(pars[i]) for i in range(89,94))
assert not any(s in case for s in ['预测失败','10个位点','移除TCR','5.46','8.27','54对'])
assert all(s in case for s in ['图4b','图4c','图4d上部','图4d下部','G4A','I7A','5.85','10.58','2.94','3.12'])
assert not any(s in alltext for s in ['图4e','图4f','（e、f）','p表示肽'])
counts=dict(native_tables=len(root.xpath('.//w:tbl',namespaces=NS)),native_equations=len(root.xpath('.//m:oMath',namespaces=NS)),images=len(root.xpath('.//a:blip',namespaces=NS)),references=len(re.findall(r'^\[\d+\]',alltext,re.M)))
assert counts==dict(native_tables=10,native_equations=36,images=6,references=28)
updates={media:FIG.read_bytes(),'word/document.xml':etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)}
with zipfile.ZipFile(OUT,'w') as z:
    for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
editable={c['body_index'] for c in changes}|{list(body).index(fig)}
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None
    changed=[p for p in parts if z.read(p)!=parts[p]];assert set(changed)==set(updates)
    reread=etree.fromstring(z.read('word/document.xml')).find('w:body',NS)
    assert len(reread)==len(original_nodes)
    for i,p in enumerate(reread):
        if i not in editable:assert etree.tostring(p)==original_nodes[i],i
(ROOT/'figure4/caption_zh.md').write_text(new[95]+'\n',encoding='utf-8')
(ROOT/'prose_v55.json').write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8')
(QA/'v55_content.txt').write_text('\n'.join(f'P{i}: {txt(p)}' for i,p in enumerate(pars)),encoding='utf-8')
audit=dict(status='built_pending_render_QA',source=str(SOURCE),source_sha256=sha(SOURCE),document=str(OUT),document_sha256=sha(OUT),figure=str(FIG),figure_sha256=sha(FIG),changed_zip_parts=changed,edited_body_indices=sorted(editable),prose_changes=changes,body_nodes_unchanged_outside_authorized_indices=True,source_preserved=True,figure_dimensions_EMU=dict(cx=width,cy=height),caption_font_pt=10.5,new_training=0,new_inference=0,new_MD=0,**counts)
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
