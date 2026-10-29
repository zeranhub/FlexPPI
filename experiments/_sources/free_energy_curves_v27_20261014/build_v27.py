from pathlib import Path
from collections import Counter
from copy import deepcopy
import hashlib,json,zipfile,re
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from docx.shared import Inches
from docx.text.paragraph import Paragraph
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v26_20261013_中文修订_自由能范围模型验证_三线表.docx'
OUT=BASE/'paper_v27_20261014_中文修订_模型选择与自由能曲线验证_三线表.docx'
FIG=ROOT/'figures/fig7_prediction_selected_MD_free_energy_curves.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='6b497baa0688b587bb1cc5be9eaa665f883da2536a688a2730dbea50d949b6e4'
d=Document(SOURCE);changes=[];removed=[]
def find(prefix):return next(p for p in d.paragraphs if p.text.startswith(prefix))
def setp(p,text):
    assert not p._p.xpath('.//m:oMath')
    changes.append({'before':p.text,'after':text})
    nodes=p._p.xpath('.//w:t');assert nodes
    nodes[0].text=text;nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for n in nodes[1:]:n.text=''
def remove(p):
    assert not p._p.xpath('.//w:commentRangeStart | .//w:commentRangeEnd | .//w:commentReference | .//m:oMath')
    removed.append(p.text);p._p.getparent().remove(p._p)
def after(p,text):
    template=find('在83条记录的33,142个同类型有效位置上，完整模型Hχ预测与MD低')
    new=deepcopy(template._p)
    for child in list(new):
        if child.tag!=qn('w:pPr'):new.remove(child)
    pp=Paragraph(new,d._body);pp.add_run(text)
    pp.paragraph_format.page_break_before=False;pp.paragraph_format.keep_with_next=False
    p._p.addnext(new);return pp

method=find('侧链实例使用泛素–EDD')
setp(method,'侧链实例使用泛素–EDD UBA复合物2QHO对应的DynaRepo记录A003R[23]，输入完整两链结构，采用第0折留出预测及主评价标签。前5,000帧分为五个连续1,000帧窗口，联合状态熵按窗口计算后取均值；状态占据率与累积覆盖合并全部帧。谷氨酰胺三个χ角各分为三个等宽区间，共27类联合状态。图6显示各位点占据率最高三类的实际MD代表帧，以位点9 Å内局部主链对齐。图7仅按冻结Hχ预测在泛素六个有效谷氨酰胺位点中排序，选择最低Q41与最高Q2；逐角曲线统一取0–100 ns、200 ps间隔的501帧，原五窗口联合熵按主标签协议核对。')
wm=find('低相对自由能角区总宽度按每个可定义χ角计算。')
setp(wm,wm.text.replace('主评价阈值固定为1，0.5和2的同协议结果随源数据提供。','主统计阈值固定为1；范围曲线在E=0至3、间隔0.1处计算W(E)，展示随截止值放宽的角区总长度。'))
sampling=find('位置筛选仅使用各预测读数：')
setp(sampling,sampling.text+'范围曲线沿用主统计的同一组配对抽样索引，报告各阈值处的区间。')

case=find('Q41与Q62均为具有三个侧链角的谷氨酰胺。')
setp(case,'Q41与Q62均为具有三个侧链角的谷氨酰胺，留出模型的联合Hχ预测分别为2.000和2.684 nats，MD参考为1.890和2.631 nats（图6g）。MD中覆盖90%观测帧分别需要七类和十二类联合状态（图6f、h）。模型保持这两个同类型位置的高低差异，并对应到实际取向组合的多样性。')
after(case,'进一步以模型预测确定曲线验证位点：在泛素六个有效谷氨酰胺位置中，Hχ预测最低为Q41（2.000 nats），最高为Q2（2.990 nats），原MD联合熵参考分别为1.890和3.009 nats。图7a–c展示这两个预测端点的实际逐角自由能投影。两个位置的χ1主占据区均在−60°附近；Q2的χ2同时覆盖反式及−60°邻域，χ3的低值角区也更广。在共同501帧窗口和四个箱格原点平均下，Q41与Q2的逐角低自由能总宽度分别为27.50和64.17°/角。模型先确定高低位置，真实曲线检验所选位置的取向范围。')

corr=find('在83条记录的33,142个同类型有效位置上，完整模型Hχ预测与MD低')
setp(corr,'在83条记录的33,142个同类型有效位置上，完整模型Hχ预测与MD低相对自由能角区总宽度的记录平均相关为0.400，链内Hχ为0.391，完整模型RMSF为0.239。侧链输出相对主链输出的配对相关增量为0.161，95%重采样区间为0.144–0.179。各输出采用相同位置和记录内同类型控制，检验相同化学种类中不同结构位置的取向范围。')
halves=find('仅按完整模型Hχ预测确定高、低半组后，低组与高组的真实低自由能')
setp(halves,'按同一预测排序依据在各记录同类型内形成高、低半组，测量实际MD范围随自由能阈值放宽的覆盖曲线（图7d）。E=1时，低组与高组分别覆盖35.41和43.17°/角，差值为7.76°/角，95%区间为7.26–8.26，83条记录均为正值。链内侧链输出与完整模型主链输出确定的组差分别为7.54和5.09°/角。位置组由原静态结构预测确定，MD曲线及宽度用于评价这些位置，支持侧链输出筛选主占据取向覆盖更广的同类型位点。')
entropy=find('对同一组位置，完整模型Hχ预测与真实平均细箱')
setp(entropy,'对同一组位置，完整模型Hχ、链内Hχ及完整模型RMSF与真实平均细箱角度占据熵的记录平均相关分别为0.510、0.493和0.321；完整侧链输出相对主链输出的配对增量为0.189，95%区间为0.172–0.206。仅按完整侧链预测确定的高、低半组真实细箱熵差为0.293 nats/角，95%区间为0.277–0.307，83条记录均为正值；链内侧链与完整主链输出的组差分别为0.281和0.198 nats/角（表S7及配套数据）。低自由能角区总宽度检验主占据范围，细箱熵检验整体占据多样性；两项结果支持侧链输出比主链输出更好地区分同类型位置的实际占据差异。')
remove(find('仅按完整模型Hχ预测确定高、低半组时，真实平均细箱'))

caption=find('图7：');pic=caption._p.getprevious();assert pic.xpath('.//w:drawing')
for child in list(pic):
    if child.tag!=qn('w:pPr'):pic.remove(child)
pp=Paragraph(pic,d._body);pp.alignment=WD_ALIGN_PARAGRAPH.CENTER
pp.paragraph_format.page_break_before=True;pp.paragraph_format.keep_with_next=True
pp.add_run().add_picture(str(FIG),width=Inches(6.2))
setp(caption,(ROOT/'figures/fig7_prediction_selected_MD_free_energy_curves.caption_zh.md').read_text(encoding='utf8').strip())
discussion=find('侧链输出的物理意义由预测对照和预测驱动的MD检验建立。')
setp(discussion,'侧链输出的物理意义由预测驱动的MD检验建立。泛素留出预测识别类型均值之外的位置差异；按预测选出的谷氨酰胺高低端点，其逐角MD曲线呈现不同的主占据范围。83条记录中，同类型侧链预测与低自由能角区总宽度及细箱占据熵的联系均高于主链预测；预测高值位置组在全部记录中具有更广的主占据取向范围和更高的细箱熵。具体曲线说明模型选择了怎样的位置，群体范围曲线检验相同排序依据在更多蛋白上的表现，支持优先筛选具有较广实际取向范围的同类型位点。')
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
assert len(dst.tables)==11
assert [p.text for p in src.paragraphs if re.match(r'^\[\d+\]',p.text)]==[p.text for p in dst.paragraphs if re.match(r'^\[\d+\]',p.text)]
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_document':str(SOURCE),'source_sha256':sha(SOURCE),'last_delivered_document':str(BASE/'paper_v25_20261013_中文修订_图7预测对照明确化_三线表.docx'),'main_figures':7,'supplementary_figures':1,'three_line_tables':11,'all405tablecellsandnativeformatpreserved':True,'math_elements':36,'math_and_comments_preserved':True,'other7activeimagesexact':True,'new_figure7_sha256':sha(FIG),'changes':changes,'removed':removed,'review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
lines=[]
for i,p in enumerate(dst.paragraphs):
    t=''.join(p._p.xpath('.//w:t/text() | .//m:t/text()'))
    if t:lines.append(f'P{i}: {t}')
for i,t in enumerate(dst.tables):
    for j,row in enumerate(t.rows):lines.append(f'T{i}R{j}: '+' | '.join(c.text for c in row.cells))
(QA/'v27_content.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')
oldroot=ROOT.parent/'free_energy_model_v26_20261013'
for source,target in [('export_v26_word.ps1','export_v27_word.ps1'),('render_v26.py','render_v27.py')]:
    text=(oldroot/source).read_text(encoding='utf8').replace(SOURCE.name,OUT.name).replace('v26_word','v27_word').replace('render_v26','render_v27').replace('v26_page_map','v27_page_map')
    (ROOT/target).write_text(text,encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k not in {'changes','removed'}},ensure_ascii=False))
