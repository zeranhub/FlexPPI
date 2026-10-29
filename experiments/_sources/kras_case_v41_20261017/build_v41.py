"""Replace the manuscript case locally, retaining every other DOCX package part."""
from pathlib import Path
import hashlib, json, re, zipfile
from lxml import etree
from PIL import Image

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent.parent
QA = ROOT / 'manuscript_QA'
QA.mkdir(exist_ok=True)
SOURCE = BASE / 'paper_v40_20261017_中文修订_TCR识别环境变化实例_三线表.docx'
OUT = BASE / 'paper_v41_20261017_中文修订_KRAS新抗原动态响应_三线表.docx'
FIGURE = ROOT / 'figure6/fig6_kras_case.png'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE) == 'bd9851a0f6dabdeb6f4a7c2fbbe9fe0dd20ec0820a156578db7f624cd990f5ca'
assert FIGURE.exists()
ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
      'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
with zipfile.ZipFile(SOURCE) as z:
    original = {i.filename:z.read(i.filename) for i in z.infolist()}
    infos = z.infolist()
root = etree.fromstring(original['word/document.xml'])
paragraphs = root.xpath('./w:body/w:p', namespaces=ns)
text = lambda p: ''.join(p.xpath('.//w:t/text()',namespaces=ns))
changes = []
def find_start(prefix):
    found = [p for p in paragraphs if text(p).startswith(prefix)]
    assert len(found)==1, (prefix, len(found))
    return found[0]
def set_text(p,new):
    assert not p.xpath('.//m:oMath',namespaces=ns)
    nodes = p.xpath('.//w:t',namespaces=ns)
    assert nodes
    old = text(p)
    nodes[0].text = new
    for n in nodes[1:]: n.text = ''
    changes.append({'paragraph':paragraphs.index(p),'before':old,'after':new})
def replace_phrase(prefix, old, new):
    p=find_start(prefix)
    assert old in text(p),old
    set_text(p,text(p).replace(old,new))

replace_phrase('主链位置波动与侧链取向分布',
 '同一肽–MHC被SB27和SB47识别时，肽的主链波动变化与MD的排序相关为0.945，中央Q7的侧链联合熵增幅也与MD一致。',
 'KRAS G12D实例中，模型与MD的全肽主链波动变化排序相关为0.806，未突变位点的变化量误差比零变化基线低40.2%。')
replace_phrase('我们整合ATLAS、DynaRepo与ITrace',
 'ITrace中的SB27/SB47比较进一步检验相同肽–MHC在不同TCR识别环境下的动态变化。',
 'ITrace中的KRAS G12D对照进一步检验模型能否预测单点抗原替换伴随的周围肽段动态变化。')
replace_phrase('结构实例考察主链活动区域定位',
 'TCR识别环境改变下的动态响应', '单点抗原替换后的肽段动态响应')

prose=json.loads((ROOT/'text_review/replacements.json').read_text(encoding='utf8'))
# Adapt the independently reviewed prose to the final two-panel figure and connected argument.
prose['result_paragraphs'][0]=prose['result_paragraphs'][0].replace(
 '这一对照因此提供了具体的动态预测问题：仅一个肽位点改变时，模型能否从相近的静态结构中辨认整条肽的运动差异（图6a）？',
 '我们比较两种复合物中整条肽的主链波动，以检验模型能否从相近的静态结构中辨认单点替换伴随的动态差异（图6a）。')
prose['result_paragraphs'][1]=prose['result_paragraphs'][1].replace('（图6b、c）','（图6b）')
prose['result_paragraphs'][2]=prose['result_paragraphs'][2].replace('（图6c、d）','（图6b）').replace(
 '两种复合物各三条MD轨迹中，7个位点的RMSF取值范围完全分开。正面证据由未替换的整段肽共同提供，表明模型能够将抗原局部变化对应到邻近及其它未改变位置的主链运动响应。',
 '其中7个位点在G12D的三条MD轨迹中均比野生型三条轨迹具有更低的RMSF。这些未替换位置的响应表明，模型能够将抗原局部变化对应到周围肽段的主链运动差异。')
prose['method2']=prose['method2'].replace('为检验完整复合物环境的作用，采用相同家族划分和训练集合的链内模型，在相同位置比较预测变化。','')
# This case is supported by the Full predictions and replicated MD. The cached chain-only
# control has poor absolute calibration for these short peptides and is not interpreted here.
prose['result_paragraphs']=prose['result_paragraphs'][:3]
set_text(find_start('在244条ITrace记录中，按肽与MHC序列匹配'),prose['method1'])
set_text(find_start('TCR实例中的变化定义为SB47减SB27'),prose['method2'])
set_text(find_start('3.5 相同肽–MHC'),prose['result_heading'])
old_starts=['SB27和SB47识别同一条','在全部13个肽位置上，模型与MD',
            '侧链变化在Q7上提供了','保持各自结合构象中的肽–MHC坐标并移除TCR后']
assert len(prose['result_paragraphs'])==3
for start,new in zip(old_starts,prose['result_paragraphs']): set_text(find_start(start),new)
removed=find_start(old_starts[-1])
changes.append({'paragraph':paragraphs.index(removed),'before':text(removed),'after':'','removed':True})
removed.getparent().remove(removed)
paragraphs.remove(removed)

caption=('图6：KRAS G12D单点替换后的肽段主链动态响应预测。（a）JDIa41b1 TCR识别HLA-A*11:01呈递的KRAS野生型（7OW5）与G12D（7OW6）肽段的输入结构，以及按MHC肽结合区Cα对齐的肽段叠合；两种复合物使用相同视角与比例，标出肽段第6位G→D替换。（b）全部10个肽位置的RMSF变化，统一取G12D减野生型；MD为每种状态三条200 ns重复轨迹的均值差，误差线为均值差标准误。模型预测保留全部位置，底纹标出突变位点；负值表示G12D肽段的波动较低。图内分别报告全部10个位点及其余9个氨基酸未变位点的变化排序相关与MAE，零变化基线将各位点变化均设为0。')
caption_p=find_start('图6：不同TCR识别环境下')
set_text(caption_p,caption)
# Keep the complete caption on the same page as its figure.
ppr=caption_p.find('w:pPr',namespaces=ns)
assert ppr is not None
keep=ppr.find('w:keepLines',namespaces=ns)
if keep is None:keep=etree.SubElement(ppr,'{'+ns['w']+'}keepLines')
keep.set('{'+ns['w']+'}val','1')
(ROOT/'figure6/caption_zh.md').write_text(caption,encoding='utf8')
replace_phrase('两类输出对应不同的筛选问题',
 'SB27/SB47识别环境改变后的全肽主链响应及Q7侧链熵增幅',
 'KRAS G12D单点替换后的肽段主链波动响应')
replace_phrase('FlexPPI从静态主链几何联合估计逐残基RMSF',
 'SB27/SB47实例进一步重现相同肽–MHC在不同TCR环境下的主链变化与Q7侧链熵升高。',
 'KRAS G12D实例进一步表明，模型能够从两种静态复合物中预测单点替换伴随的肽段主链波动差异，并在未突变位点上重现主要响应。')
set_text(find_start('[28]'),prose['reference28'])

# Reuse the original image part and relationship; only its binary contents and aspect ratio change.
blips=root.xpath('.//a:blip[@r:embed="rId52"]',namespaces=ns)
assert len(blips)==1
inline=blips[0].xpath('ancestor::wp:inline',namespaces=ns)[0]
width=round(6.4*914400)
with Image.open(FIGURE) as im: height=round(width*im.height/im.width)
for ext in inline.xpath('./wp:extent|.//a:xfrm/a:ext',namespaces=ns):
    ext.set('cx',str(width));ext.set('cy',str(height))
for prop in inline.xpath('./wp:docPr',namespaces=ns):
    prop.set('name','Figure 6 KRAS G12D peptide dynamic response')
    prop.set('descr','Wild-type and KRAS G12D TCR-pMHC structures and all ten peptide RMSF changes compared with three-replica molecular dynamics.')

alltext='\n'.join(text(p) for p in root.xpath('.//w:p',namespaces=ns))
for forbidden in ['SB27','SB47','Q7','R7','2AK4','4JRY','LPEPLPQGQLTAY','M112.447185','冻结','发布权重','发布模型']:
    assert forbidden not in alltext,forbidden
assert len(root.xpath('.//m:oMath',namespaces=ns))==36
assert len(root.xpath('.//w:tbl',namespaces=ns))==10
assert len(root.xpath('.//a:blip',namespaces=ns))==8
assert len(re.findall(r'^\[\d+\]',alltext,re.M))==28
updates={'word/document.xml':etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),
         'word/media/image41.png':FIGURE.read_bytes()}
with zipfile.ZipFile(OUT,'w') as z:
    for info in infos:z.writestr(info,updates.get(info.filename,original[info.filename]))
with zipfile.ZipFile(OUT) as z: assert z.testzip() is None
(QA/'v41_content.txt').write_text('\n'.join(f'P{i}: {text(p)}' for i,p in enumerate(paragraphs)),encoding='utf8')
audit={'status':'built_pending_render_QA','source':str(SOURCE),'source_sha256':sha(SOURCE),
       'document':str(OUT),'document_sha256':sha(OUT),'changed_zip_parts':list(updates),
       'changes':changes,'figure':str(FIGURE),'figure_sha256':sha(FIGURE),
       'existing_native_tables':10,'existing_native_equations':36,'new_training':0,'new_MD':0,
       'primary_model_checkpoint_sha256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',
       'primary_model_training_ITrace_records':0,'MD_reference':'three existing 200ns ITrace replicas per state; full saved frames; recomputed labels',
       'delta':'7OW6 G12D minus 7OW5 wild-type','figure_dimensions_EMU':[width,height]}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
