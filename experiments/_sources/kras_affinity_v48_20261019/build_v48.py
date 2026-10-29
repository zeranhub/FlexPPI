"""Replace Figure 4 and its connected case argument while preserving native parts."""
from pathlib import Path
import hashlib, json, re, zipfile
from lxml import etree
from PIL import Image

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v47_20261019_中文修订_模型框架图重绘_三线表.docx'
OUT=BASE/'paper_v48_20261019_中文修订_KRAS亲和力机制与绝对RMSF_三线表.docx'
FIG=ROOT/'figure4/fig4_kras_affinity.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='7b37df1a8726609844681969befeb15e8f77a8ac342ef2a59d7345d05320dd7f'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
 'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
 'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
 'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
 'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
with zipfile.ZipFile(SOURCE) as z:
    infos=z.infolist();parts={p.filename:z.read(p.filename) for p in infos}
root=etree.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
txt=lambda n:''.join(n.xpath('.//w:t/text()',namespaces=NS))
paragraphs=[p for p in body if p.tag=='{'+NS['w']+'}p']
old={i:txt(p) for i,p in enumerate(paragraphs)}
new={}
new[2]=old[2].replace('KRAS G12D实例中，模型与MD的全肽主链波动变化排序相关为0.806，未突变位点的变化量误差比零变化基线低40.2%。','KRAS G12D实例中，模型将单点抗原替换伴随的主链约束增强定位到第6位及邻近未突变位置，与文献报道的HLA介导识别机制相呼应。')
new[9]=old[9].replace('ITrace中的KRAS G12D对照进一步检验模型能否预测单点抗原替换伴随的周围肽段动态变化。','ITrace中的KRAS G12D对照进一步分析模型预测的局部动态变化与文献结合实验及相互作用网络的联系。')
new[52]='TCR实例采用ITrace中的KRAS野生型和G12D复合物7OW5与7OW6。两者包含相同的工程TCR JDIa41b1与HLA-A*11:01，十肽序列分别为VVVGAGGVGK与VVVGADGVGK，仅肽第6位由Gly替换为Asp，对应KRAS蛋白的G12D变化[28]。预测采用2.5节全组成链家族划分的第0折完整模型，该折训练集合由ATLAS和DynaRepo的2,430条记录组成，未使用ITrace。将两个复合物的原始预测按全部10个肽位置配对，以Å为单位分别绘制野生型与G12D的绝对RMSF曲线。'
new[53]='补充的MD分析使用ITrace为每种复合物提供的三条200 ns独立重复轨迹，保留0–200 ns的全部1,001–1,003个保存帧，按2.2节逐链对齐后重新计算肽链RMSF，并核对原标签。按相同肽位置比较两种复合物的重复轨迹结果，考察模型定位区域的动态趋势。'
new[89]='3.4 KRAS G12D新抗原识别中的局部主链稳定化'
new[90]='KRAS是一种调控细胞生长信号的小GTP酶，在GTP结合的活化态与GDP结合的失活态之间切换。G12D将KRAS第12位的甘氨酸替换为天冬氨酸，使其更难通过GTP水解回到失活状态，可推动持续的致癌信号。细胞内KRAS降解产生的肽段可由I类人类白细胞抗原（HLA）呈递到细胞表面，其中携带G12D的肽可成为T细胞受体（TCR）识别的肿瘤新抗原。Poole等通过表面等离子体共振（SPR）测得，经亲和力优化的同一JDIa41b1 TCR对HLA-A*11:01呈递的KRAS野生型与G12D十肽，解离常数分别为3.0 μM与743±18 pM，相差约4,000倍；X射线晶体结构却显示两种复合物的整体结合姿态相近[28]。我们以这组对照分析，模型能否从相近的静态结构中定位与抗原识别差异相关的局部动态响应（图4a）。'
new[92]='图4：FlexPPI预测KRAS G12D替换后的肽段主链动态响应。（a）JDIa41b1 TCR识别HLA-A*11:01呈递的KRAS野生型（7OW5）与G12D（7OW6）肽段的输入结构，以及按MHC肽结合区Cα对齐的肽段叠合；两种复合物使用相同视角与比例，标出肽段第6位G→D替换，并列出文献[28]通过SPR测得的解离常数。（b）模型分别根据两个复合物结构预测的全部10个肽位置RMSF绝对值，橙线为野生型，紫线为G12D，底纹标出替换位点。'
new[93]='两条模型预测曲线显示，动态差异集中在肽中央及C端一侧，而N端V1–V3基本重合（图4b）。替换位G6/D6的预测RMSF从1.080降至0.718 Å，未替换的邻位G7从1.246降至0.904 Å；A5与V8也分别由0.885和1.052 Å降至0.680和0.886 Å。响应因此延伸到氨基酸未改变的相邻骨架，其中G7位于TCR的CDR3β接触区域[28]。模型定位的约束增强覆盖了突变位点与已有识别接触位置，显示其能够从复合物几何中识别局部替换伴随的周围动态变化。ITrace独立重复轨迹中，G12D肽中央区域也呈较低波动，为这一局部响应提供补充佐证。'
new[94]='原研究的结构分析为这一预测模式提供了具体的相互作用解释：D6侧链埋入HLA结合槽，新增与HLA R114的盐桥及与Q70的氢键，而两种复合物的直接TCR接触大体相似[28]。变温SPR显示，G12D结合具有更有利的焓贡献；作者进一步对每种复合物进行25条4 ns MD模拟，并通过MM/PBSA能量分解将亲和力差异主要联系到增强的静电相互作用[28]。D6在HLA内部获得额外约束，模型则将较低主链波动定位到D6及邻近的G7等识别骨架位置：这一空间对应把新增的肽–HLA相互作用与TCR接触区域的动态变化连在一起，为文献提出的HLA介导抗原选择性机制提供了逐位置的动态支持。'
new[110]=old[110].replace('KRAS G12D单点替换后的肽段主链波动响应','KRAS G12D中与HLA内相互作用相呼应的局部主链响应')
new[112]=old[112].replace('KRAS G12D实例进一步表明，模型能够从两种静态复合物中预测单点替换伴随的肽段主链波动差异，并在未突变位点上重现主要响应。','KRAS G12D实例进一步表明，模型能够从两种相近的静态复合物中定位替换位点及邻近未突变位置的主链约束增强，为HLA介导的抗原识别机制提供局部动态支持。')
changes=[]
for i,text in new.items():
    p=paragraphs[i]
    assert text!=old[i],i
    assert not p.xpath('.//m:oMath|.//w:fldChar|.//w:hyperlink',namespaces=NS),i
    ts=p.xpath('.//w:t',namespaces=NS);assert ts
    ts[0].text=text
    for t in ts[1:]:t.text=''
    changes.append({'paragraph_index':i,'body_index':list(body).index(p),'before':old[i],'after':text})
cap=paragraphs[92];fig=cap.getprevious()
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
    pr.set('name','Figure 4 KRAS antigen recognition and predicted RMSF')
    pr.set('descr','WT and G12D TCR-HLA structures with published SPR affinities, and model-predicted absolute peptide RMSF at all ten positions')
def prop(p,name,value):
    pr=p.find('w:pPr',NS)
    if pr is None:pr=etree.Element('{'+NS['w']+'}pPr');p.insert(0,pr)
    e=pr.find('w:'+name,NS)
    if e is None:e=etree.SubElement(pr,'{'+NS['w']+'}'+name)
    e.set('{'+NS['w']+'}val',str(value))
prop(fig,'keepNext',1);prop(cap,'keepLines',1);prop(cap,'keepNext',0)
alltext='\n'.join(txt(p) for p in body)
assert len(root.xpath('.//w:tbl',namespaces=NS))==10
assert len(root.xpath('.//m:oMath',namespaces=NS))==36
assert len(root.xpath('.//a:blip',namespaces=NS))==6
assert len(re.findall(r'^\[\d+\]',alltext,re.M))==28
assert '低估了这些位置的变化幅度' not in alltext
assert '均值差标准误' not in txt(cap)
updates={media:FIG.read_bytes(),'word/document.xml':etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)}
with zipfile.ZipFile(OUT,'w') as z:
    for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
with zipfile.ZipFile(OUT) as z:assert z.testzip() is None
(ROOT/'figure4/caption_zh.md').write_text(new[92]+'\n',encoding='utf8')
(QA/'v48_content.txt').write_text('\n'.join(f'P{i}: {txt(p)}' for i,p in enumerate(paragraphs)),encoding='utf8')
audit={'status':'built_pending_render_QA','source':str(SOURCE),'source_sha256':sha(SOURCE),
 'document':str(OUT),'document_sha256':sha(OUT),'figure4':str(FIG),'figure4_sha256':sha(FIG),
 'caption_body_index':list(body).index(cap),'prose_changes':changes,
 'editable_body_indices':[list(body).index(fig)]+[item['body_index'] for item in changes],
 'changed_zip_parts':list(updates),'width_EMU':width,'height_EMU':height,
 'native_tables':10,'native_equations':36,'references':28,'main_figures':5,'supplementary_figures':1,
 'new_MD':0,'new_training':0,'new_inference':0}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
