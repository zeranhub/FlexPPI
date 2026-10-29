"""Add concise, source-checked biological introductions to the three protein cases."""
from pathlib import Path
import hashlib,json,zipfile
from lxml import etree

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v41_20261017_中文修订_KRAS新抗原动态响应_三线表.docx'
OUT=BASE/'paper_v42_20261017_中文修订_蛋白背景与实例动机_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='58a15ba884180814dc562d276a3fb1c0e321b95a73cc302e0b0e68a4dafc0a68'
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as z:
    infos=z.infolist();original={i.filename:z.read(i.filename) for i in infos}
root=etree.fromstring(original['word/document.xml'])
paragraphs=root.xpath('./w:body/w:p',namespaces=ns)
text=lambda p:''.join(p.xpath('.//w:t/text()',namespaces=ns))
changes=[]
def find(prefix):
    matches=[p for p in paragraphs if text(p).startswith(prefix)]
    assert len(matches)==1,(prefix,len(matches))
    return matches[0]
def replace(p,new):
    assert not p.xpath('.//m:oMath',namespaces=ns)
    nodes=p.xpath('.//w:t',namespaces=ns);assert nodes
    changes.append({'paragraph':paragraphs.index(p),'before':text(p),'after':new})
    nodes[0].text=new
    for n in nodes[1:]:n.text=''

kix=('转录共激活因子CBP（CREB结合蛋白）通过结合多种转录调控蛋白参与基因表达调控，KIX是其介导这些相互作用的小型结构域。'
     '转录因子c-Myb与转录调控蛋白MLL的激活片段可同时结合KIX相对两侧的不同表面[22]。'
     '我们在无伙伴、仅c-Myb、仅MLL及双伙伴四种组成中，将留出模型的逐残基RMSF与同一KIX链87个位点的MD结果对应，'
     '检验模型能否定位连接相邻螺旋的L12活动环区（图5a–d）。')
replace(find('KIX转录共激活结构域在不同表面结合'),kix)

kras_background=('KRAS是一种调控细胞生长信号的小GTP酶，在GTP结合的活化态与GDP结合的失活态之间切换。'
                 'G12D将KRAS第12位的甘氨酸替换为天冬氨酸，使其更难通过GTP水解回到失活状态，可推动持续的致癌信号。'
                 '细胞内KRAS降解产生的肽段可由I类人类白细胞抗原（HLA）呈递到细胞表面，其中携带G12D的肽可成为T细胞受体（TCR）识别的肿瘤新抗原。')
p=find('JDIa41b1 TCR对HLA-A*11:01呈递')
old=text(p)
new=old.replace('JDIa41b1 TCR对HLA-A*11:01呈递的KRAS野生型与G12D肽',
                '经亲和力优化的JDIa41b1 TCR对HLA-A*11:01呈递的KRAS野生型与G12D十肽')
replace(p,kras_background+new)

ubiquitin=('泛素是由76个氨基酸组成的小蛋白，通过共价连接到其它蛋白参与蛋白质降解等过程。'
           'E3泛素连接酶EDD的N端含有泛素相关结构域（UBA），可通过疏水接触和氢键识别泛素表面[23]。'
           '我们在该复合物中比较泛素不同位置的同种侧链，检验模型能否辨认局部环境对应的取向多样性。')
p=find('泛素–EDD UBA复合物提供同类型侧链的留出验证')
replace(p,text(p).replace('泛素–EDD UBA复合物提供同类型侧链的留出验证[23]。',ubiquitin))

# Identify the same objects unambiguously in the technical methods without repeating the biology.
p=find('结构实例考察主链活动区域定位')
replace(p,text(p).replace('我们对KIX转录共激活结构域开展','我们对CBP的KIX结构域开展'))
p=find('侧链实例使用泛素–EDD UBA复合物')
replace(p,text(p).replace('侧链实例使用泛素–EDD UBA复合物2QHO',
                         '侧链实例使用泛素与EDD的泛素相关结构域（UBA）形成的复合物2QHO'))

# Place the unchanged KRAS figure after the biological setup so the page does not
# end with a large gap. Quantitative model results then continue below the figure.
caption=find('图6：KRAS G12D单点替换后的')
figure=caption.getprevious()
assert figure.xpath('.//a:blip',namespaces=ns)
anchor=find('在全部10个肽位置上，完整模型')
anchor.addprevious(figure)
anchor.addprevious(caption)

newxml=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with zipfile.ZipFile(OUT,'w') as z:
    for i in infos:z.writestr(i,newxml if i.filename=='word/document.xml' else original[i.filename])
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None
    changed=[i for i in original if z.read(i)!=original[i]]
assert changed==['word/document.xml']
assert len(root.xpath('.//m:oMath',namespaces=ns))==36
assert len(root.xpath('.//w:tbl',namespaces=ns))==10
assert len(root.xpath('.//a:blip',namespaces=ns))==8
alltext='\n'.join(text(p) for p in root.xpath('.//w:p',namespaces=ns))
for oldcase in ['SB27','SB47','2AK4','4JRY','Q7','发布权重','发布模型','冻结']:assert oldcase not in alltext
assert len(changes)==5
audit={'status':'built_pending_render_QA','source':str(SOURCE),'source_sha256':sha(SOURCE),
       'document':str(OUT),'document_sha256':sha(OUT),'changed_zip_parts':changed,
       'changes':changes,'references_reused':[22,23,28],
       'layout_adjustment':'Unchanged Figure6 and caption moved directly after the KRAS biological setup, before quantitative results.',
       'unchanged':'all8figures,10native3line_tables,36nativeequations,originalcommentsandanchors,allnumericalresults',
       'new_MD':0,'new_training':0,'new_inference':0}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
(QA/'v42_content.txt').write_text('\n'.join(f'P{i}: {text(p)}' for i,p in enumerate(paragraphs)),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
