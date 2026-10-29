"""Integrate the authorized MEL5 two-panel interface revision into v50."""
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
SOURCE = BASE / 'paper_v50_20261022_中文修订_MEL5残基接触网络与图4精简_三线表.docx'
OUT = BASE / 'paper_v51_20261022_中文修订_MEL5跨界面残基机制与图4重绘_三线表.docx'
FIG = ROOT / 'figure4/fig4_mel5_interface.png'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE) == '9c05035e9af5ba5a7c2291f29ef939be7bc5e2d676d3e30c8d6c3e59cc36aaf8'
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
new = {
89:'3.4 MEL5亲和成熟中跨界面接触与双动态响应的残基解析',
90:'MART-1（又称Melan-A）是黑色素瘤相关抗原，其肽段由HLA-A*02:01呈递并被T细胞识别。MEL5是识别这一复合物的天然TCR，α24β17是经体外亲和成熟获得的变体，两者结合相同的MART-1来源改良十肽ELAGIGILTV。Madura等通过表面等离子体共振（SPR）测得，亲和成熟后解离常数由18 μM降至0.6 nM，亲和力增强约30,000倍，主要伴随解离减慢；晶体结构中总体结合姿态和肽构象相近，TCR与MHC的接触数则由56增至97，βCDR2参与的pMHC接触由6增至38[28]。这一同抗原、不同受体的对照，用于检验FlexPPI能否将接触重排定位到界面两侧具体残基的主链与侧链动态响应（图4）。',
91:'在六个识别环的全部45个对齐位点中，模型将主链波动降低集中定位到β链识别环（图4b）。βCDR2的五个位点预测RMSF均降低，均值由0.896降至0.636 Å；βCDR3的13个位点也均降低，均值由0.801降至0.648 Å。在对侧pMHC界面的23个有效且氨基酸类型相同的位点中，HLA Q72的预测Hχ降幅最大，由1.792降至1.375 nats，表示其三个χ角的联合状态占据更集中；与之相对、未发生替换的TCR βGly54预测RMSF由1.118降至0.747 Å。两个输出由此共同指向Q72–β53/β54跨界面片区（图4a）。',
92:'将两个复合物以HLA对齐并从同一视角观察，可直接比较受体改造前后的局部配对关系（图4a）。亲和成熟后的βCDR2包含Gly52→Pro52、Ile53→Phe53等替换，原研究指出Ile53→Phe53带来了额外MHC接触[28]。在本研究用于预测的结构中，Q72侧链OE1与β53主链N的距离由5.85缩短至2.94 Å，Q72侧链NE2与βGly54主链O的距离由10.58缩短至3.12 Å。Q72侧链的酰胺两端由此靠近相邻两处TCR骨架，将β53的接触重排与未替换Gly54连接到同一局部网络。模型在这一片区同时预测到Q72侧链状态集中和Gly54主链波动降低，为界面两侧共同参与的局部稳定化提供动态解释。',
93:'已有ITrace三条独立重复轨迹进一步支持这一接触网络的稳定化：Q72–βGly54近接触的平均帧占比由24.4%升至97.4%，Q72两端同时与β53/β54主链满足氢键几何的平均帧占比由0%升至37.1%；Q72 Hχ和βGly54 RMSF的三重复均值也均降低。结合原研究的SPR和晶体结构结果，FlexPPI把βCDR2接触增加所伴随的动态响应落实到对侧Q72与相邻未替换Gly54，支持该片区更稳定地贴合HLA、参与维持跨界面接触的机制解释。该实例展示了模型将主链与侧链输出共同用于定位受体改造后局部动态响应的能力。',
95:'图4：FlexPPI定位MEL5亲和成熟后的跨界面接触与双动态响应。（a）左侧为复合物定位与文献[28]的SPR亲和力；中央和右侧为MEL5与α24β17经HLA对齐后的同视角、同尺度局部结构，突出HLA Q72、TCR β52、βIle53/βPhe53及未替换βGly54。虚线连接Q72 OE1与β53主链N、Q72 NE2与βGly54主链O，并标注原子间距；Q72 Hχ与βGly54 RMSF的模型读数直接标在对应残基附近。（b）六个完整CDR中全部45个对齐位点的两态模型RMSF，识别环之间分隔，突出βCDR2及Gly54。橙色为MEL5，紫色为α24β17；HLA及Q72以青色表示，编号统一到MEL5。局部结构沿用来源于3HG1和4JFF的MD准备起始坐标，预测输入为氨基酸类型及N、Cα、C坐标；配对与共同编号规则见方法。',
}
sentence='为检验Q72预测对受体环境的响应，在保留各自pMHC坐标和序列的条件下从输入移除TCR两条链，并比较该位点的模型Hχ。'
assert old[52].count(sentence) == 1
new[52] = old[52].replace(sentence, '跨界面局部结构沿用来源于3HG1与4JFF的MD准备起始坐标，以HLA第1–180位共同Cα原子对齐两个复合物，并用相同视角和尺度显示。图中标注Q72 OE1–β53主链N和Q72 NE2–βGly54主链O的原子间距。')
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
    prop.set('name','Figure 4 MEL5 paired cross-interface residue contacts')
    prop.set('descr','Two panels: shared-view and shared-scale MEL5/alpha24beta17 interface structures, Q72-beta53/beta54 distances and direct FlexPPI readouts; all45CDR model RMSF curves. No panel c or d.')
alltext='\n'.join(txt(p) for p in body)
case='\n'.join(new[i] for i in range(89,94))
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
(ROOT/'prose_v51.json').write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8')
(QA/'v51_content.txt').write_text('\n'.join(f'P{i}: {txt(p)}' for i,p in enumerate(pars)),encoding='utf-8')
audit=dict(status='built_pending_render_QA',source=str(SOURCE),source_sha256=sha(SOURCE),
           document=str(OUT),document_sha256=sha(OUT),figure=str(FIG),figure_sha256=sha(FIG),
           changed_zip_parts=changed,edited_body_indices=sorted(editable),prose_changes=changes,
           body_nodes_unchanged_outside_authorized_indices=True,source_preserved=True,
           figure_dimensions_EMU=dict(cx=width,cy=height),new_training=0,new_inference=0,new_MD=0,**counts)
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
