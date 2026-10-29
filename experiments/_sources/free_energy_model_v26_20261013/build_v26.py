from pathlib import Path
from collections import Counter
from copy import deepcopy
import hashlib,json,zipfile,re
import pandas as pd
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from docx.shared import Inches
from docx.text.paragraph import Paragraph
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v25_20261013_中文修订_图7预测对照明确化_三线表.docx'
OUT=BASE/'paper_v26_20261013_中文修订_自由能范围模型验证_三线表.docx'
FIG=ROOT/'figures/fig7_model_free_energy_and_occupancy_validation.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='523290c8297ba2034cf249bd07c57256cdc1920cb0f9530b4acd21d0e6a9335b'
d=Document(SOURCE);P={i:p for i,p in enumerate(d.paragraphs)}
changes=[];removed=[]
def setp(p,text):
    assert not p._p.xpath('.//m:oMath')
    changes.append({'before':p.text,'after':text})
    nodes=p._p.xpath('.//w:t');assert nodes
    nodes[0].text=text;nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for n in nodes[1:]:n.text=''
def change(i,text):setp(P[i],text)
def remove(p):
    assert not p._p.xpath('.//w:commentRangeStart | .//w:commentRangeEnd | .//w:commentReference | .//m:oMath')
    removed.append(p.text);p._p.getparent().remove(p._p)
def after(p,text):
    new=deepcopy(P[104]._p)
    for child in list(new):
        if child.tag!=qn('w:pPr'):new.remove(child)
    pp=Paragraph(new,d._body);pp.add_run(text)
    pp.paragraph_format.page_break_before=False;pp.paragraph_format.keep_with_next=False
    p._p.addnext(new);return pp
def summary(name,mode='replica1'):
    f=pd.read_csv(ROOT/'experiment'/f'{mode}_{name}.csv')
    return f[(f.threshold==1)&(f.scope=='all')&(f.dataset=='all')]
c=summary('correlation_summary').set_index('predictor');h=summary('high_low_summary').set_index('predictor')
diff=summary('paired_correlation_difference_summary').set_index('contrast')
width=c.loc['full_hchi'];rmsf=c.loc['full_rmsf'];intra=c.loc['intrachain_hchi'];high=h.loc['full_hchi']
assert int(width.n_records)==int(high.n_records)==83 and int(high.positive_records)==83
f3=lambda x:f'{x:.3f}'
f2=lambda x:f'{x:.2f}'

change(2,P[2].text.replace('泛素的多角分布与经验投影自由能展开这一预测差异的具体取向形态。',f'在同一队列中，侧链输出与MD低相对自由能角区总宽度的排序相关为{f3(width["mean"])}，按预测确定的高值组实际角区平均宽{f2(high["mean"])}°/角，将读数进一步对应到主占据取向范围。'))
change(9,P[9].text.replace('进一步固定83条记录的折外输出，检验侧链预测能否识别实际连续角度占据的差异，并用泛素的角度密度与经验投影自由能解释预测对应的构象形态。','进一步固定83条记录的折外输出，检验模型排序及其选出的高低位置组能否区分MD低相对自由能角区总宽度与细粒度角度占据。'))
change(56,'式(5)给出侧链监督目标与联合类别相对占据的联系。为检验模型估计值能否识别主要访问的取向范围，我们从MD逐角经验投影自由能提取低值角区总宽度，再与冻结折外预测对应。高占据角度具有较低的F*(q)[3,24,25]；角区宽度强调主占据范围，细箱熵同时响应弱占据区域及各箱占比，两类轨迹读数共同用于评价预测的局部动态含义。')
change(59,'低相对自由能角区总宽度按每个可定义χ角计算。以5°等宽箱统计MD概率密度p̂(q)，令F*(q)=−ln[p̂(q)/p̂max]，最高密度位置为零；将F*≤1的全部角度箱宽相加，等价于统计密度至少为最高值1/e的角区总长度，相互分离的合格箱均计入。对0°、1.25°、2.5°和3.75°四个箱格原点分别计算后平均，再对该残基全部可定义角平均，单位为°/角；角度按圆周处理。主评价阈值固定为1，0.5和2的同协议结果随源数据提供。')
change(60,'第二个轨迹指标为平均细箱角度占据熵：对每个可定义χ角计算5°细箱Shannon熵，先对相同四个箱格原点平均，再对该残基全部角平均，单位为nats/角。两指标沿用完全相同的位置和观测窗口。记录内同一种氨基酸至少具有四个共同有效位置时参与评价，保留33,142个位点；分别从各轨迹指标及Full Hχ、Intra Hχ、Full RMSF折外预测减去记录内同类型均值，在相同位置上计算Spearman相关，再对记录等权平均。')
change(61,P[61].text.replace('计算两组真实细箱熵之差','计算两组真实角区总宽度或细箱熵之差'))
change(96,'侧链联合熵的预测能力需要对应到实际取向占据。本节先通过泛素的完整留出预测与同类型位点对照，检验模型是否识别不同的侧链状态多样性，再以83条记录的冻结折外输出检验其排序及预测选出的高低位置组是否具有不同的MD低自由能角区范围和细粒度占据。')
change(99,P[99].text.replace('（图7a）','（图6g）'))
remove(P[100])
change(103,'3.6.2 冻结侧链预测对低自由能取向范围与角度占据的验证')
gap=diff.loc['full_hchi_minus_full_rmsf']
new1=after(P[103],f'在83条记录的33,142个同类型有效位置上，完整模型Hχ预测与MD低相对自由能角区总宽度的记录平均相关为{f3(width["mean"])}，链内Hχ为{f3(intra["mean"])}，完整模型RMSF为{f3(rmsf["mean"])}（图7a）。侧链输出相对主链输出的配对相关增量为{f3(gap["mean"]) }，95%重采样区间为{f3(gap.bootstrap_95_low)}–{f3(gap.bootstrap_95_high)}。同记录同类型的控制使这一比较检验相同化学种类中不同结构位置的取向范围，而各输出均在相同位置上计分。')
half=pd.read_csv(ROOT/'experiment/replica1_high_low_per_record.csv')
half=half[(half.threshold==1)&(half.scope=='all')&(half.predictor=='full_hchi')]
after(new1,f'仅按完整模型Hχ预测确定高、低半组后，低组与高组的真实低自由能角区总宽度分别为{f2(half.low_MD_mean.mean())}和{f2(half.high_MD_mean.mean())}°/角，两组之差为{f2(high["mean"])}°/角，95%区间为{f2(high.bootstrap_95_low)}–{f2(high.bootstrap_95_high)}，83条记录均为正值（图7c）。链内侧链预测与完整模型主链预测确定的组差分别为{f2(h.loc["intrachain_hchi","mean"])}和{f2(h.loc["full_rmsf","mean"])}°/角。该检验中位置由静态结构的预测值选择，真实MD范围用于评价所选位置，支持侧链输出筛选主占据取向覆盖更广的同类型位点。')
change(104,P[104].text.replace('（图8a）','（图7b）').replace('在83条记录的33,142个同类型有效位置上，','对同一组位置，'))
change(105,P[105].text.replace('（图8b）','（图7d）').replace('因此，模型从静态几何给出的侧链高低值能够用于筛选同类型中角度占据更分散的位置；泛素的自由能投影展示了这一群体联系在具体旋转轴上的形态。','低自由能角区总宽度检验主占据取向范围，细箱熵检验整体占据多样性；两项结果共同支持模型从静态几何筛选同类型中取向范围更广、占据更分散的位置。'))

caption=P[102];pic=caption._p.getprevious();assert pic.xpath('.//w:drawing')
for child in list(pic):
    if child.tag!=qn('w:pPr'):pic.remove(child)
pp=Paragraph(pic,d._body);pp.alignment=WD_ALIGN_PARAGRAPH.CENTER
pp.paragraph_format.page_break_before=True;pp.paragraph_format.keep_with_next=True
pp.add_run().add_picture(str(FIG),width=Inches(6.2))
setp(caption,(ROOT/'figures/fig7_model_free_energy_and_occupancy_validation.caption_zh.md').read_text(encoding='utf8').strip())
P[105]._p.addnext(pic);pic.addnext(caption._p)
fig8=P[107];pic8=fig8._p.getprevious();assert pic8.xpath('.//w:drawing')
assert not pic8.xpath('.//w:commentRangeStart | .//w:commentReference')
pic8.getparent().remove(pic8);remove(fig8)
change(110,'侧链输出的物理意义由预测对照和预测驱动的MD检验建立。泛素留出预测识别类型均值之外的位置差异，并在Q41与Q62上保留联合熵的高低顺序。83条记录中，侧链预测与低自由能角区总宽度及细箱占据熵的同类型相关均高于主链预测；仅按侧链预测确定的高值位置组在全部记录中具有更广的主占据取向范围和更高的细箱占据熵。这使读数的用途从对固定标签的回归，延伸到优先筛选具有较广实际取向范围的同类型位点。')
next(p for p in d.paragraphs if p.text=='4 讨论').paragraph_format.page_break_before=True
change(113,P[113].text.replace('冻结折外预测对连续角度占据的群体验证进一步支持侧链输出筛选实际取向差异的能力，泛素自由能投影解释这种差异的具体形态。','冻结折外预测对MD低自由能角区范围及连续角度占据的群体验证，进一步支持侧链输出筛选取向范围更广、占据更分散位点的能力。'))
# Supplement E retains all previous quantitative evidence and adds new repeat result.
ac=summary('correlation_summary','available_replicas').set_index('predictor')
ah=summary('high_low_summary','available_replicas').set_index('predictor')
e=next(p for p in d.paragraphs if p.text.startswith('共同窗口中的33,493个位点有351个'))
setp(e,e.text+f'低自由能角区总宽度采用相同位置，使用全部可用重复的指标均值时，完整侧链输出的记录平均相关为{f3(ac.loc["full_hchi","mean"])}，仅按该预测确定的实际高低组差为{f2(ah.loc["full_hchi","mean"])}°/角；逐角、阈值及角数分层结果随源数据提供。')

d.save(OUT)
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    entries=[(info,a.read(info.filename) if info.filename=='word/comments.xml' else b.read(info.filename)) for info in b.infolist()]
tmp=OUT.with_suffix('.preserved.tmp')
with zipfile.ZipFile(tmp,'w') as z:
    for info,data in entries:z.writestr(info,data)
tmp.replace(OUT)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    anchors=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    math=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    tables=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//w:tbl',namespaces=ns)]
    assert anchors(ax)==anchors(bx) and math(ax)==math(bx) and len(math(bx))==36
    assert tables(ax)==tables(bx) and a.read('word/comments.xml')==b.read('word/comments.xml')
    assert len(bx.xpath('.//a:blip',namespaces=ns))==8
src=Document(SOURCE);dst=Document(OUT)
media=lambda doc:[hashlib.sha256(doc.part.related_parts[p.xpath('.//a:blip')[0].get(qn('r:embed'))].blob).hexdigest() for p in doc._element.body.xpath('./w:p') if p.xpath('.//a:blip')]
ma,mb=media(src),media(dst)
assert len(mb)==8 and sum(s in ma for s in mb)==7 and mb.count(sha(FIG))==1
text='\n'.join(p.text for p in dst.paragraphs)
assert '图8' not in text and len(dst.tables)==11
assert [p.text for p in src.paragraphs if re.match(r'^\[\d+\]',p.text)]==[p.text for p in dst.paragraphs if re.match(r'^\[\d+\]',p.text)]
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_document':str(SOURCE),'source_sha256':sha(SOURCE),'main_figures':7,'supplementary_figures':1,'native_three_line_tables':11,'all405tablecellsandtableformatpreserved':True,'math_elements':36,'math_and_comments_preserved':True,'other7activeimagesexact':True,'new_figure7_sha256':sha(FIG),'changes':changes,'removed_paragraphs':removed,'review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
lines=[]
for i,p in enumerate(dst.paragraphs):
    t=''.join(p._p.xpath('.//w:t/text() | .//m:t/text()'))
    if t:lines.append(f'P{i}: {t}')
for i,t in enumerate(dst.tables):
    for j,row in enumerate(t.rows):lines.append(f'T{i}R{j}: '+' | '.join(c.text for c in row.cells))
(QA/'v26_content.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')
oldroot=ROOT.parent/'figure7_model_validation_v25_20261013'
for source,target in [('export_v25_word.ps1','export_v26_word.ps1'),('render_v25.py','render_v26.py')]:
    text=(oldroot/source).read_text(encoding='utf8').replace(SOURCE.name,OUT.name).replace('v25_word','v26_word').replace('render_v25','render_v26').replace('v25_page_map','v26_page_map')
    (ROOT/target).write_text(text,encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k not in {'changes','removed_paragraphs'}},ensure_ascii=False))
