"""Integrate the two approved source-data figure merges, retaining native Word parts."""
from pathlib import Path
from copy import deepcopy
import hashlib,json,re,zipfile
from lxml import etree
from PIL import Image

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v42_20261017_中文修订_蛋白背景与实例动机_三线表.docx'
OUT=BASE/'paper_v43_20261018_中文修订_能力与构象合并图_三线表.docx'
F2=ROOT/'figure2/fig2_model_capability.png'
F5=ROOT/'figure5/fig5_joint_readouts_conformations.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='0d56dd11e3abab0b7e12bbd75ef33c352d949bc305212c01f9eb5509541df887'
assert F2.is_file() and F5.is_file()
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
with zipfile.ZipFile(SOURCE) as z:
    infos=z.infolist();original={i.filename:z.read(i.filename) for i in infos}
root=etree.fromstring(original['word/document.xml']);body=root.find('w:body',ns)
original_body=list(body)
text=lambda p:''.join(p.xpath('.//w:t/text()',namespaces=ns))
changes=[];deleted=[]
def find(prefix):
    matches=[p for p in body if p.tag=='{'+ns['w']+'}p' and text(p).startswith(prefix)]
    assert len(matches)==1,(prefix,len(matches))
    return matches[0]
def set_text(p,new):
    assert not p.xpath('.//m:oMath',namespaces=ns)
    nodes=p.xpath('.//w:t',namespaces=ns);assert nodes
    old=text(p)
    if old==new:return
    changes.append({'source_body_index':original_body.index(p),'before':old,'after':new})
    nodes[0].text=new
    for n in nodes[1:]:n.text=''
def drop(p):
    assert not p.xpath('.//m:oMath|.//w:commentRangeStart|.//w:commentReference|.//w:bookmarkStart',namespaces=ns)
    deleted.append({'source_body_index':original_body.index(p),'text':text(p)})
    body.remove(p)
def paragraph_property(p,name,value):
    pr=p.find('w:pPr',ns)
    if pr is None:
        pr=etree.Element('{'+ns['w']+'}pPr');p.insert(0,pr)
    el=pr.find('w:'+name,ns)
    if el is None:el=etree.SubElement(pr,'{'+ns['w']+'}'+name)
    el.set('{'+ns['w']+'}val',value)

proposal=json.loads((ROOT/'text_review/manuscript_integration_proposal.json').read_text(encoding='utf8'))
assert proposal['source_sha256']==sha(SOURCE)
# Capture the original illustration/caption nodes before any numbering changes.
caps={n:find('图'+str(n)+'：') for n in range(1,8)}
figs={n:caps[n].getprevious() for n in range(1,8)}
assert all(p.xpath('.//a:blip',namespaces=ns) for p in figs.values())
for item in proposal['paragraph_replacements_outside_new35']:
    p=find(item['find_paragraph_start']);assert text(p)==item['old']
    set_text(p,item['new'])

# Remove the obsolete diagnostic paragraph and the two obsolete active figure blocks.
drop(find('图3同时展示模型在不同动态幅度下的误差'))
drop(figs[3]);drop(caps[3])
for prefix in ['3.3 主链与侧链读数的互补组合','全体折外预测的配对分析检验两个输出是否提供不同的局部信息',
               '同记录同类型中，低RMSF／高Hχ与高RMSF／低Hχ分别包含']:
    drop(find(prefix))
drop(figs[4]);drop(caps[4])
set_text(find('3.4 KIX结构域中主链活动区域的定位'),'3.3 KIX结构域中主链活动区域的定位')
set_text(find('3.5 KRAS G12D单点抗原变化的主链动态响应'),'3.4 KRAS G12D单点抗原变化的主链动态响应')

# Reorganize the final Results section around the population-to-case figure narrative.
prefixes=['3.6 侧链动态预测的构象对应与群体验证',
          '侧链联合熵的预测能力需要对应到实际取向占据',
          '3.6.2 侧链预测对低自由能取向范围与角度占据的验证',
          '在83条记录的33,142个同类型有效位置上',
          '按侧链预测在各记录同类型内形成高、低半组',
          '对同一组位置，完整模型Hχ、链内Hχ及完整模型RMSF',
          '3.6.1 泛素同类型侧链多样性的预测与轨迹对照',
          '泛素是由76个氨基酸组成的小蛋白',
          'Q41与Q62均为具有三个侧链角的谷氨酰胺',
          '在泛素各χ角层内控制氨基酸类型后']
nodes=[find(p) for p in prefixes]
newparas=proposal['new_section_35_ordered_paragraphs'];assert len(newparas)==len(nodes)
for node,item in zip(nodes,newparas):set_text(node,item['text'])
anchor=find('4 讨论')
for node in nodes+[figs[7],caps[7]]:anchor.addprevious(node)

# The retained molecular figures keep all data and visual content; only their numbering changes.
set_text(caps[5],text(caps[5]).replace('图5：','图3：',1))
set_text(caps[6],text(caps[6]).replace('图6：','图4：',1))
set_text(caps[2],(ROOT/'figure2/caption_zh.md').read_text(encoding='utf8').strip())
set_text(caps[7],(ROOT/'figure5/caption_zh.md').read_text(encoding='utf8').strip())
assert text(caps[2]).startswith('图2：') and text(caps[7]).startswith('图5：')

relationships=etree.fromstring(original['word/_rels/document.xml.rels'])
rels={n.get('Id'):n.get('Target') for n in relationships}
updates={};image_log=[]
def replace_image(fig,cap,path,number):
    blips=fig.xpath('.//a:blip',namespaces=ns);assert len(blips)==1
    rid=blips[0].get('{'+ns['r']+'}embed');part='word/'+rels[rid]
    assert part in original
    inline=blips[0].xpath('ancestor::wp:inline',namespaces=ns)[0]
    width=round(6.6*914400)
    with Image.open(path) as im:height=round(width*im.height/im.width)
    for extent in inline.xpath('./wp:extent|.//a:xfrm/a:ext',namespaces=ns):
        extent.set('cx',str(width));extent.set('cy',str(height))
    for prop in inline.xpath('./wp:docPr',namespaces=ns):
        prop.set('name',f'Figure {number} merged quantitative evidence')
        prop.set('descr','Model capability' if number==2 else 'Joint dynamic readouts and MD conformational validation')
    paragraph_property(fig,'keepNext','1');paragraph_property(cap,'keepLines','1')
    paragraph_property(cap,'keepNext','0')
    pr=fig.find('w:pPr',ns);ind=pr.find('w:ind',ns)
    if ind is None:ind=etree.SubElement(pr,'{'+ns['w']+'}ind')
    for key in ['left','right','firstLine']:ind.set('{'+ns['w']+'}'+key,'0')
    updates[part]=path.read_bytes()
    image_log.append({'figure':number,'part':part,'source':str(path),'sha256':sha(path),'width_EMU':width,'height_EMU':height})
replace_image(figs[2],caps[2],F2,2)
replace_image(figs[7],caps[7],F5,5)
for prop in figs[6].xpath('.//wp:docPr',namespaces=ns):
    prop.set('name','Figure 4 KRAS G12D peptide dynamic response')

alltext='\n'.join(text(p) for p in body if p.tag=='{'+ns['w']+'}p')
assert not re.search(r'图[67]',alltext)
assert not re.search(r'^3\.6',alltext,re.M)
assert len(re.findall(r'^图[1-5]：',alltext,re.M))==5
assert len(re.findall(r'^\[\d+\]',alltext,re.M))==28
assert len(root.xpath('.//a:blip',namespaces=ns))==6
assert len(root.xpath('.//w:tbl',namespaces=ns))==10
assert len(root.xpath('.//m:oMath',namespaces=ns))==36
for token in ['SB27','SB47','Q7','发布权重','发布模型','冻结']:assert token not in alltext
# Byte-preserve every other OOXML part, including all original comments and unused image parts.
updates['word/document.xml']=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with zipfile.ZipFile(OUT,'w') as z:
    for info in infos:z.writestr(info,updates.get(info.filename,original[info.filename]))
with zipfile.ZipFile(OUT) as z:assert z.testzip() is None
document_paragraphs=[p for p in body if p.tag=='{'+ns['w']+'}p']
(QA/'v43_content.txt').write_text('\n'.join(f'P{i}: {text(p)}' for i,p in enumerate(document_paragraphs)),encoding='utf8')
audit={'status':'built_pending_render_QA','source':str(SOURCE),'source_sha256':sha(SOURCE),
       'document':str(OUT),'document_sha256':sha(OUT),'changed_zip_parts':list(updates),
       'changes':changes,'deleted_nodes':deleted,'image_replacements':image_log,
       'main_figures':5,'supplementary_figures':1,'native_tables':10,'native_equations':36,'references':28,
       'figure_map':proposal['main_figure_map'],'new_MD':0,'new_training':0,'new_inference':0,
       'retained_original_source_images':'Unreferenced original Figure3/4 parts retained inside the ZIP; they are not active document figures.',
       'results_order':['3.1','3.2','3.3 KIX','3.4 KRAS','3.5','3.5.1 population','3.5.2 ubiquitin']}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','changed_zip_parts']},ensure_ascii=False))
