"""Integrate the residue-pair contact map and aligned model tracks into v53."""
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
SOURCE = BASE / 'paper_v52_20261022_中文修订_MEL5图4去除局部结构图_三线表.docx'
OUT = BASE / 'paper_v53_20261022_中文修订_MEL5残基接触与双动态联合图_三线表.docx'
FIG = ROOT / 'figure4/fig4_mel5_contact_map.png'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE) == 'bcf8a06760532e63b330341f315fc075acc8640f4b69e46c968e2538a74b61eb'
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
method='局部原子间距由来源于3HG1与4JFF的MD准备起始坐标计算，分别测量Q72 OE1–β53主链N及Q72 NE2–βGly54主链O的距离。'
assert old[52].count(method)==1
new[52]=old[52].replace(method,'局部接触分析沿用来源于3HG1与4JFF的MD准备起始坐标，计算HLA第65–76位与TCR β48–56位的最短重原子间距。在任一态与该TCR片段距离≤6 Å且Hχ有效的六个HLA位置上，保留与九个TCR位置的全部54对距离；残基对距≤4 Å定义为近接触。图4a按这些位置对齐两态模型Hχ与RMSF。')
new[92]='为解释这两个位于界面两侧的模型响应，进一步将局部接触关系与两端的动态预测对齐（图4a）。亲和成熟后的βCDR2包含Gly52→Pro52、Ile53→Phe53等替换，原研究指出Ile53→Phe53带来了额外MHC接触[28]。在包含六个HLA邻接位点与TCR β48–56位的局部窗口中，Q72保留了与β51位的近接触，并新增与β52、β53、β54及β56位的近接触；Q72与β53、βGly54的最短重原子间距分别由5.46和8.27 Å缩短至2.94和3.12 Å，而与β55位的距离由3.94变为4.76 Å。这一接触重排将Q72侧链与βCDR2及相邻骨架联系到同一局部网络。模型同时预测到对侧Q72侧链状态集中和该受体片段的主链活动收紧，为界面两侧共同参与的局部稳定化提供动态解释。'
new[95]='图4：FlexPPI将MEL5亲和成熟后的界面接触与两侧动态响应对应到具体残基。（a）中部每格为一对HLA–TCR残基，橙色圆与紫色方块分别表示MEL5和α24β17的最短重原子距离：实心为≤4 Å近接触，空心为4–6 Å，小点为>6 Å。完整显示六个局部Hχ有效HLA位置与TCR β48–56位的全部54对距离；左侧Hχ预测与HLA行对齐，上方RMSF预测与TCR列对齐。突出Q72行及β53/β54列，Q72行数字为两态距离（Å），TCR位点下方标注氨基酸替换。局部位点按方法中的几何规则确定，复合物图仅作定位。（b）六个完整CDR中全部45个对齐位点的两态模型RMSF，识别环之间分隔，突出βCDR2及Gly54。颜色与形状在两部分保持一致，编号统一到MEL5；结构距离取自已有MD准备起始坐标，模型输入为氨基酸类型及N、Cα、C坐标。'
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
    prop.set('name','Figure 4 MEL5 residue contact map with aligned dynamic predictions')
    prop.set('descr','Two panels: complete54pair HLA-TCR contact map with row-aligned Hchi and column-aligned RMSF predictions; all45CDR absolute model RMSF. Colored paired markers encode closest heavy-atom distance; no local molecular close-ups or text cards.')
alltext='\n'.join(txt(p) for p in body)
case='\n'.join(txt(pars[i]) for i in range(89,94))
assert not any(t in alltext for t in ['图4c','图4d','（图4c）'])
assert not any(t in case for t in ['预测失败','10个位点','G4A','I7A','移除TCR','2.141','2.245'])
assert all(t in case for t in ['5.46','2.94','8.27','3.12','1.792','1.375','1.118','0.747','Q72','βGly54','23个'])
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
(ROOT/'prose_v53.json').write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8')
(QA/'v53_content.txt').write_text('\n'.join(f'P{i}: {txt(p)}' for i,p in enumerate(pars)),encoding='utf-8')
audit=dict(status='built_pending_render_QA',source=str(SOURCE),source_sha256=sha(SOURCE),
           document=str(OUT),document_sha256=sha(OUT),figure=str(FIG),figure_sha256=sha(FIG),
           changed_zip_parts=changed,edited_body_indices=sorted(editable),prose_changes=changes,
           body_nodes_unchanged_outside_authorized_indices=True,source_preserved=True,
           figure_dimensions_EMU=dict(cx=width,cy=height),new_training=0,new_inference=0,new_MD=0,**counts)
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
