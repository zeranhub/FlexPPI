"""Apply the agreed single-mutation-to-residue-dynamics narrative without changing results."""
from pathlib import Path
from copy import deepcopy
from zipfile import ZipFile,ZIP_DEFLATED
from hashlib import sha256
from lxml import etree
import json
ROOT=Path(__file__).resolve().parent
PAPER=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
SOURCE=PAPER/'paper_v127_20261028_中文修订_TCR结合模型对照.docx'
TARGET=PAPER/'paper_v128_20261028_中文修订_单点突变与残基动态叙事.docx'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
q=lambda n:'{'+NS['w']+'}'+n
with ZipFile(SOURCE) as z:
 items=z.infolist();parts={x.filename:z.read(x.filename) for x in items}
assert sha256(SOURCE.read_bytes()).hexdigest()=='7ed98727b33f4658611846651f68017639cd6eab509594ed4384b22d92a54dcf'
doc=etree.fromstring(parts['word/document.xml']);body=doc.find('w:body',NS)
ps=body.findall('w:p',NS)
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def find(prefix):
 nodes=[p for p in ps if text(p).startswith(prefix)]
 assert len(nodes)==1,(prefix,len(nodes))
 return nodes[0]
def paragraph(value,template):
 p=etree.Element(q('p'))
 pp=template.find('w:pPr',NS)
 if pp is not None:p.append(deepcopy(pp))
 r=etree.SubElement(p,q('r'));rp=template.find('w:r/w:rPr',NS)
 if rp is not None:r.append(deepcopy(rp))
 etree.SubElement(r,q('t')).text=value
 return p
changes=[]
def replace(prefix,new):
 old=find(prefix);p=paragraph(new,old);body.replace(old,p)
 changes.append({'before':text(old),'after':new})
 return p

motivation=(
 '单个残基替换即可对TCR–pMHC亲和力产生显著影响。在α24β17体系中，G4A与I7A使SPR解离常数由0.6 nM分别升至36与31 μM，亲和力下降约60,000与51,700倍，而两种替换态仍有可测结合[28]。'
 'ERGO-II、NetTCR-2.2与ImRex对这些替换的结合评分响应存在差异：NetTCR-2.2与ImRex默认版给出评分下降的方向，ERGO-II两套权重则未给出该方向[33–35]。'
 '这些评分不直接表征KD，模型还存在不同程度的同家族训练覆盖，因此本例的评分对照不能作为严格家族外泛化的证据。完整评分、训练覆盖与不结合对照见表S8、图S8及补充材料G，通用绝对亲和力对照见补充材料F。'
)
transition=(
 '即使整体评分能够提示结合倾向的变化，单独使用该评分仍难以判断哪些残基的动态性质发生改变，以及这些改变如何联系到界面相互作用的重新组织。'
 '为进一步分析这一问题，FlexPPI提供逐残基的主链波动与侧链构象多样性读数，并结合重复分子动力学模拟，检验肽替换相关的识别环动态响应及接触重分配。'
 '这些分析将整体结合差异联系到具体残基和局部约束，为理解突变影响亲和力的可能机制提供可检验线索。'
)
a=replace('肽单点替换提供了亲和力',motivation)
a.addnext(paragraph(transition,a));changes.append({'inserted_after':'3.4 single-residue substitution motivation','after':transition})
old_disc=text(find('MEL5家族留出实例'))
prefix=old_disc.split(' 肽替换的完整核查显示')[0]
assert prefix.endswith('后续功能研究。')
discussion=(prefix+
 '肽替换的核查进一步显示，I7A的识别环主链变化与MD具有区域性对应，G4A的恢复则较弱。'
 '单点替换引起的亲和力巨变提示，整体结合差异仍需联系到具体位置的相互作用与动态组织；FlexPPI与重复MD的结合为这一分析提供残基分辨的证据，并支持提出关于局部约束和接触重分配的机制假设。'
 '结合态柔性不能直接换算为总结合熵，原研究强调的溶剂热力学贡献仍需独立检验。'
)
replace('MEL5家族留出实例',discussion)
training=(
 'NetTCR统一使用pan模式，不将ELA的肽特异模型与突变肽的pan模型混用，也不加TCRbase缩放。'
 '公开训练表含野生MEL5–ELA的完整六CDR输入，也含成熟受体双CDR3–ELA记录，后者的CDR1/2仍为野生序列。'
 '预测前发现这两组双CDR3相关的11行均位于分区0，主结果固定为t.0.v.1–4四个官方模型均值，其训练及验证不使用分区0；全20模型为敏感性对照。'
 '该处理仅排除已确认重叠，不保证所有近缘家族隔离。ImRex按官方README默认选择downsample，并保留full版；其训练含野生MEL5–ELA阳性，成熟CDR3β未见，但最近训练序列为野生MEL5，两条15残基序列中11位相同。'
 'ERGO-II的VDJdb训练中也存在野生β–ELA记录。三方法已核查的公开训练表中均未发现G4A或I7A肽；精确未见不能证明家族隔离，同家族覆盖也不能确定评分变化的来源。'
 '因此，本例仅作带训练覆盖限制的回顾性对照，不作为统一家族隔离的公平性能基准。'
)
replace('NetTCR统一使用pan模式',training)
interpretation=(
 '下述方向均指各配置内部的结合评分变化，不表示预测了KD的变化幅度。ImRex默认版对四核心态均给出高于0.5的分数，且G4A与I7A分数均低于成熟ELA；'
 '但I5A（实验NB）和YLEPGPVTA负对照的评分分别为0.631与0.611，在固定0.5参考切分下判为结合。full版对I5A为0.232，对I7A却高于ELA。'
 'NetTCR主结果对野生MEL5与成熟ELA的评分为0.962与0.219，其排序与实验亲和成熟方向不同；I7A阳性与I5A阴性的评分为0.078与0.076。'
 'ERGO-II两权重的G4A/I7A分数均高于ELA。作者推理接口主要输出连续分数，0.5只作未调参且未经本SPR实验校准的参考，不能把阈值结果直接作为方法普遍失效的证据。'
)
replace('识别倾向、结合效应量和结合态动态分别解释',interpretation)
parts['word/document.xml']=etree.tostring(doc,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
 for item in items:z.writestr(item,parts[item.filename])
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
 assert b.testzip() is None
 changed=[n for n in a.namelist() if a.read(n)!=b.read(n)]
 assert changed==['word/document.xml']
 assert a.namelist()==b.namelist()
 olddoc=etree.fromstring(a.read('word/document.xml'))
 for xpath in ['//w:tbl','//w:fldSimple','//w:instrText','//w:commentRangeStart','//w:commentRangeEnd','//w:sectPr','//w:drawing']:
  oldnodes=olddoc.xpath(xpath,namespaces=NS);newnodes=doc.xpath(xpath,namespaces=NS)
  assert len(oldnodes)==len(newnodes),xpath
  assert all(etree.tostring(x)==etree.tostring(y) for x,y in zip(oldnodes,newnodes)),xpath
 media=[n for n in a.namelist() if n.startswith('word/media/')]
 assert len(media)==50 and all(a.read(n)==b.read(n) for n in media)
 # Existing complete dynamics and contact paragraphs must remain literally unchanged.
 for pre in ['肽替换的动态对应','重复MD将这些结合态差异']:
  orig=text(find(pre));matches=doc.xpath('//w:p[starts-with(string(.),"'+pre+'")]',namespaces=NS)
  assert any(text(p)==orig for p in body.findall('w:p',NS))
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),'target_sha256':sha256(TARGET.read_bytes()).hexdigest(),'changed_original_parts':changed,'original_media_preserved':len(media),'all_figures_tables_fields_comments_and_sections_unchanged':True,'existing_dynamics_and_contact_evidence_unchanged':True,'changes':changes,'render_validation':'pending'}
(ROOT/'revision_changes.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'3.4叙事与讨论修订.md').write_text('3.4 单点替换与残基动态分析\n\n'+motivation+'\n\n'+transition+'\n\n讨论中的对应段落\n\n'+discussion+'\n',encoding='utf-8')
print(json.dumps({'target':str(TARGET),'paragraph_changes':len(changes),'preserved_media':len(media),'changed_package_parts':changed},ensure_ascii=False))
