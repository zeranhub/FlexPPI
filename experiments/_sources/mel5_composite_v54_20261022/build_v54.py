"""Replace Figure4 with the user-selected six-panel composite; synchronize its prose."""
from pathlib import Path
import hashlib,json,re,zipfile
from lxml import etree
from PIL import Image

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v53_20261022_中文修订_MEL5残基接触与双动态联合图_三线表.docx'
OUT=BASE/'paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx'
FIG=ROOT/'figure4/fig4_mel5_composite.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='d35736b5165c9e12525f414ffa2a2c966bfde17d01c7247310bff1834ebe9852'
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
new[52]='TCR实例采用ITrace中的MEL5、亲和成熟变体α24β17及两种肽替换复合物。3HG1与4JFF分别为MEL5和α24β17结合HLA-A*02:01呈递的改良MART-1十肽ELAGIGILTV（ELA）；4JFD与4JFE分别为α24β17识别ELA第4位Gly→Ala（G4A）及第7位Ile→Ala（I7A）的复合物[28]。预测采用2.5节全组成链家族划分的第0折完整模型，训练集合为ATLAS和DynaRepo的2,430条记录，未使用ITrace。四态分别输入各自已有复合物的氨基酸类型及N、Cα、C坐标。通过全链序列对齐，将MEL5六个互补决定区（CDR）的边界投影到其他状态，保留全部45个对应位点及四态全部10个肽位置。三种差值依次为4JFF−3HG1、4JFD−4JFF和4JFE−4JFF。pMHC界面取四态中与TCR最近Cα距离≤8 Å的位置并集，共30位；Hχ仅在前后氨基酸相同且侧链χ角有定义时作差，三种对照有效位点数分别为23、23和22。图文TCR编号统一到MEL5，故成熟TCR的βPhe53和βGly54分别对应其输入E54与E55。局部结构原子间距由3HG1与4JFF的MD准备起始坐标计算，测量Q72 OE1–β53主链N及Q72 NE2–βGly54主链O。'
new[89]='3.4 MEL5体系中跨界面接触与双动态响应的残基解析'
new[90]=old[90]+'进一步纳入成熟TCR识别G4A和I7A肽的复合物，以比较受体改造与抗原替换在同一体系中的动态响应，四态肽曲线见图4b。'
new[91]=old[91].replace('（图4b）','（图4c）').replace('在对侧pMHC界面的23个有效且氨基酸类型相同的位点中，','在对侧pMHC界面的23个有效且氨基酸类型相同的位点中（图4d），')
new[92]='亲和成熟后的βCDR2包含Gly52→Pro52、Ile53→Phe53等替换，原研究指出Ile53→Phe53带来了额外MHC接触[28]。在来源于3HG1与4JFF的起始坐标中，Q72侧链OE1与β53主链N的距离由5.85缩短至2.94 Å，Q72侧链NE2与βGly54主链O的距离由10.58缩短至3.12 Å（图4a）。Q72酰胺两端由此靠近相邻两处TCR骨架，将β53的接触重排与未替换Gly54联系到同一局部网络。模型在该片区同时预测到Q72侧链状态集中和Gly54主链波动降低，为界面两侧共同参与的局部稳定化提供动态解释。'
new[93]=old[93].replace('该实例展示了模型将主链与侧链输出共同用于定位受体改造后局部动态响应的能力。','在α24β17背景引入G4A或I7A后，模型将新增主链响应主要定位至CDR2α和CDR3α（图4e），并解析出HLA K68、R75和Q72等位点幅度不同的侧链响应（图4f）。这一组比较把受体改造与抗原替换产生的动态变化区分到肽、识别环和HLA侧链的具体位置。')
new[95]='图4：MEL5体系中受体亲和成熟与肽替换的残基动态响应。（a）成熟复合物概览及Q72–β53/β54局部结构；虚线为Q72 OE1–β53 N和NE2–β54 O的原子间距。（b）四态全部10个肽位置的预测RMSF；P4为识别区域，P10为C端HLA锚定位点。（c）亲和成熟前后全部45个CDR位点的RMSF。（d）23个有效界面位点的Hχ变化排序，p表示肽。（e、f）三种对照的CDR ΔRMSF与pMHC界面ΔHχ，依次为4JFF−3HG1、4JFD−4JFF和4JFE−4JFF；蓝、红表示降低、升高，灰格表示Hχ无定义或氨基酸类型不同。界面为四态最近Cα距离≤8 Å的并集，共30位，三行Hχ有效数为23、23、22。TCR采用MEL5共同编号；各态均由同一模型对已有复合物结构预测，状态定义见方法。'
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
    e.set('name','Figure 4 MEL5 six-panel structural and dynamic response composite')
    e.set('descr','Six panels: complex and Q72-beta53/54 atom-pair distances; four-state peptide RMSF; all45 CDR absolute RMSF; all23 valid interface Hchi changes; three-contrast complete CDR and pMHC interface heatmaps. Orange MEL5, purple alpha24beta17, blue G4A, sage I7A. Blue-white-rose signed changes; grey ineligible Hchi.')
alltext='\n'.join(txt(p) for p in body)
case='\n'.join(txt(pars[i]) for i in range(89,94))
assert not any(s in case for s in ['预测失败','10个位点','移除TCR','5.46','8.27','54对'])
assert all(s in case for s in ['图4b','图4c','图4d','图4e','图4f','G4A','I7A','5.85','10.58','2.94','3.12'])
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
(ROOT/'prose_v54.json').write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8')
(QA/'v54_content.txt').write_text('\n'.join(f'P{i}: {txt(p)}' for i,p in enumerate(pars)),encoding='utf-8')
audit=dict(status='built_pending_render_QA',source=str(SOURCE),source_sha256=sha(SOURCE),document=str(OUT),document_sha256=sha(OUT),figure=str(FIG),figure_sha256=sha(FIG),changed_zip_parts=changed,edited_body_indices=sorted(editable),prose_changes=changes,body_nodes_unchanged_outside_authorized_indices=True,source_preserved=True,figure_dimensions_EMU=dict(cx=width,cy=height),caption_font_pt=10.5,new_training=0,new_inference=0,new_MD=0,**counts)
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
