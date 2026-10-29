"""Focus the authorized MEL5 argument and replace Figure4, preserving other DOCX parts."""
from pathlib import Path
from copy import deepcopy
import hashlib,json,re,zipfile
from lxml import etree
from PIL import Image

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'
OUT=BASE/'paper_v50_20261022_中文修订_MEL5残基接触网络与图4精简_三线表.docx'
FIG=ROOT/'figure4/fig4_mel5_focus.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='5f79ddf7d657dc9b782c9fb432e373a1e7c170c1c88d019e528fbf73c366c222'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
with zipfile.ZipFile(SOURCE) as z:
    infos=z.infolist();parts={p.filename:z.read(p.filename) for p in infos}
root=etree.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
original_body=list(body);pars=[p for p in body if p.tag=='{'+NS['w']+'}p']
txt=lambda p:''.join(p.xpath('.//w:t/text()',namespaces=NS))
old={i:txt(p) for i,p in enumerate(pars)}
assert old[89].startswith('3.4 MEL5') and old[95].startswith('图4：')
new={
 89:'3.4 MEL5亲和成熟中Q72–βGly54接触片区的双动态解析',
 90:'MART-1（又称Melan-A）是黑色素瘤相关抗原，其肽段由HLA-A*02:01呈递并被T细胞识别。MEL5是识别这一复合物的天然TCR，α24β17是经体外亲和成熟获得的变体，两者结合相同的MART-1来源改良十肽ELAGIGILTV。Madura等通过表面等离子体共振（SPR）测得，其解离常数由18 μM降至0.6 nM，亲和力增强约30,000倍，主要伴随解离减慢；晶体结构中整体结合姿态和肽构象相近，TCR与MHC的接触数则由56增至97，βCDR2参与的pMHC接触由6增至38[28]。这一同抗原、不同受体的对照，用于考察FlexPPI能否将界面重排联系到具体残基的主链与侧链动态响应（图4a）。',
 91:'对全部六个识别环的45个对齐位点进行比较，模型将主链波动降低集中定位到β链识别环（图4b）。βCDR2的五个位点预测RMSF均降低，均值由0.896降至0.636 Å；βCDR3的13个位点也均降低，均值由0.801降至0.648 Å。βCDR2的这一响应与原研究中该区域接触显著增加相对应[28]。在对侧pMHC界面，23个有效且氨基酸类型相同的位点中有17个预测Hχ降低，其中HLA Q72的降幅最大（图4c），由此将侧链响应进一步聚焦到βCDR2邻近的接触片区。',
 92:'在该片区，模型预测Q72的Hχ由1.792降至1.375 nats，表示其三个χ角的联合状态占据更集中；与之相邻且未发生替换的TCR βGly54，其预测RMSF由1.118降至0.747 Å。保留各自pMHC坐标并从输入中移除TCR后，Q72的预测Hχ由2.141变为2.245 nats，完整复合物中的降低不再出现，说明该预测响应使用了对侧受体环境（图4c）。局部结构给出了这两项读数的空间联系：亲和成熟后的βCDR2包含Gly52→Pro52、Ile53→Phe53等替换，Q72侧链酰胺的两端分别朝向βPhe53的主链N和βGly54的主链O，使侧链状态集中与相邻受体骨架的波动降低对应到同一跨界面接触网络（图4a）。',
 93:'已有ITrace三条独立重复轨迹进一步支持该接触网络的稳定化：Q72–βGly54近接触的平均帧占比由24.4%升至97.4%，Q72两端同时与β53/β54主链满足氢键几何的平均帧占比由0%升至37.1%；Q72 Hχ和βGly54 RMSF的三重复均值也均降低。原研究通过晶体结构指出，βIle53→Phe53带来了额外MHC接触[28]；FlexPPI则将这一接触重排对应到对侧Q72的侧链状态与未替换Gly54的主链活动，为βCDR2更稳定地贴合HLA、以及SPR所见的解离减慢提供了局部动态解释。该实例展示了联合输出将受体改造的响应定位到相邻未替换残基及其跨界面接触网络的能力。',
 95:'图4：FlexPPI联合定位MEL5亲和成熟后的跨界面动态响应。橙色为MEL5，紫色为α24β17。（a）复合物概览、文献[28]报道的SPR亲和力与Q72–βPhe53/βGly54局部接触网络，旁列Q72 Hχ与βGly54 RMSF的模型变化量。（b）六个完整CDR中45个对齐位点的两态模型RMSF，识别环之间分隔。（c）左侧为全部23个有效且氨基酸类型相同的pMHC界面位点的模型Hχ变化，按变化量排序，突出Q72，p表示肽位点；右侧为完整复合物与保留各自pMHC坐标、仅从输入移除TCR时Q72的两态模型Hχ。变化量均为α24β17减去MEL5；原子结构用于展示接触位置，预测输入为残基类型及N、Cα、C坐标；界面、共同编号及输入对照定义见方法。',
}
replacements={
 2:('MEL5亲和成熟实例中，模型将HLA Q72侧链取向的集中与未替换TCR βGly54主链波动的降低对应到同一接触片区，为识别界面稳定化提供残基层面的动态解释。','MEL5亲和成熟实例将HLA Q72的侧链联合状态集中与未替换TCR βGly54的主链波动降低联系到Q72–βPhe53/βGly54接触网络，为局部界面稳定化提供残基层面的动态解释。'),
 52:('另在保留各自pMHC全部坐标的条件下移除TCR链，比较相同位点的模型输出。','为检验Q72预测对受体环境的响应，在保留各自pMHC坐标和序列的条件下从输入移除TCR两条链，并比较该位点的模型Hχ。'),
 111:('MEL5亲和成熟实例进一步将HLA Q72的取向集中与未替换βGly54的主链波动降低对应到同一接触网络，展示两种读数联合提出残基层面机制解释的能力。','MEL5实例进一步将HLA Q72的侧链联合状态集中与未替换βGly54的主链波动降低联系到Q72–βPhe53/βGly54接触网络，展示两类输出共同定位受体改造后局部动态响应的能力。'),
 113:('MEL5亲和成熟实例将未改变类型的HLA Q72与TCR βGly54联系到同一接触片区，通过侧链取向与骨架波动的联合响应，为识别界面稳定化提供残基层面的动态解释。','MEL5亲和成熟实例将HLA Q72的侧链联合状态集中与TCR βGly54的主链波动降低联系到Q72–βPhe53/βGly54接触网络，为识别界面稳定化提供残基层面的动态解释。'),
}
for i,(first,second) in replacements.items():
    assert old[i].count(first)==1,(i,first)
    new[i]=old[i].replace(first,second)
changes=[]
for i,text in sorted(new.items()):
    p=pars[i];assert text!=old[i]
    assert not p.xpath('.//m:oMath|.//w:fldChar|.//w:hyperlink',namespaces=NS),i
    ts=p.xpath('.//w:t',namespaces=NS);assert ts
    ts[0].text=text
    for t in ts[1:]:t.text=''
    changes.append(dict(paragraph_index=i,body_index=original_body.index(p),before=old[i],after=text))
caption=pars[95];fig=pars[94]
assert len(fig.xpath('.//a:blip',namespaces=NS))==1
rid=fig.xpath('.//a:blip',namespaces=NS)[0].get('{'+NS['r']+'}embed')
rels=etree.fromstring(parts['word/_rels/document.xml.rels'])
media='word/'+{e.get('Id'):e.get('Target') for e in rels}[rid]
assert media=='word/media/image41.png'
inline=fig.xpath('.//wp:inline',namespaces=NS)[0]
width=round(6.6*914400)
with Image.open(FIG) as im:height=round(width*im.height/im.width)
for e in inline.xpath('./wp:extent|.//a:xfrm/a:ext',namespaces=NS):e.set('cx',str(width));e.set('cy',str(height))
for prop in inline.xpath('./wp:docPr',namespaces=NS):
    prop.set('name','Figure 4 MEL5 Q72 beta53 beta54 contact network')
    prop.set('descr','Three panels: real molecular network and paired model changes; complete45CDR RMSF; all23interface Hchi changes and fixed-coordinate TCR-removal Q72 input contrast. No panel d.')
alltext='\n'.join(txt(p) for p in body)
case_text='\n'.join(new[i] for i in range(89,94))
assert '图4d' not in alltext
assert not any(token in case_text for token in ['预测失败','10个位点','G4A','I7A','4JFD','4JFE'])
assert all(token in case_text for token in ['Q72','βPhe53','βGly54','23个','2.141','2.245','37.1%'])
assert len(root.xpath('.//w:tbl',namespaces=NS))==10
assert len(root.xpath('.//m:oMath',namespaces=NS))==36
assert len(root.xpath('.//a:blip',namespaces=NS))==6
assert len(re.findall(r'^\[\d+\]',alltext,re.M))==28
updates={media:FIG.read_bytes(),'word/document.xml':etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)}
with zipfile.ZipFile(OUT,'w') as z:
    for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None
    changed=[p for p in parts if z.read(p)!=parts[p]]
    assert set(changed)==set(updates)
    rr=etree.fromstring(z.read('word/document.xml'));bb=rr.find('w:body',NS)
    editable={c['body_index'] for c in changes}|{original_body.index(fig)}
    assert len(bb)==len(original_body)
    for index,(before,after) in enumerate(zip(original_body,bb)):
        if index not in editable:
            # original_body contains mutated editable nodes only; all other nodes are untouched.
            assert etree.tostring(before)==etree.tostring(after),index
contract=dict(replacements={str(i):text for i,text in new.items()},main_claim='FlexPPI links HLAQ72 joint-sidechain-state concentration and unchanged TCRbetaGly54 backbone fluctuation reduction to the Q72-betaPhe53/betaGly54 contact network.',figure_panels=['a spatial network and model response','b all45CDR absolute RMSF','c all23interface Hchi changes and Q72 TCR-input contrast'],removed_panels=['old d separate MD occupancy'],no_new_training_or_MD_or_inference=True)
(ROOT/'prose_v50.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'figure4/caption_zh.md').write_text(new[95]+'\n',encoding='utf8')
(QA/'v50_content.txt').write_text('\n'.join(f'P{i}: {txt(p)}' for i,p in enumerate(pars)),encoding='utf8')
audit=dict(status='built_pending_render_QA',source=str(SOURCE),source_sha256=sha(SOURCE),document=str(OUT),document_sha256=sha(OUT),figure=str(FIG),figure_sha256=sha(FIG),changed_zip_parts=changed,body_nodes_unchanged_outside_authorized_indices=True,edited_body_indices=sorted(editable),prose_changes=changes,native_tables=10,native_equations=36,images=6,references=28,source_preserved=True,figure_dimensions_EMU=dict(cx=width,cy=height),new_inference=0,new_MD=0,new_training=0)
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False,indent=2))
