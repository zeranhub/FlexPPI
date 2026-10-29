"""Insert only verified, completed analyses; preserve the original manuscript."""
from pathlib import Path
import hashlib,json,zipfile
from copy import deepcopy
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

root=Path(__file__).resolve().parent
source=root.parent/'paper_v15_20261010_中文修订_三线表.docx'
output=root.parent/'paper_v16_20261010_中文修订_已完成分析_三线表.docx'
doc=Document(source)

def replace_text(paragraph,text):
    # Keep native comment anchors and paragraph properties.
    for child in list(paragraph._p):
        if child.tag==qn('w:r') and not child.findall('.//'+qn('w:commentReference')):
            paragraph._p.remove(child)
    paragraph.add_run(text)

def after(paragraph,text,style='Normal'):
    new=doc.add_paragraph(text,style)
    paragraph._p.addnext(new._p)
    return new

paragraphs=list(doc.paragraphs)
replace_text(paragraphs[29],paragraphs[29].text.replace('其校准与应用边界在讨论中说明','其覆盖率诊断与事后分组校准见补充材料C，应用边界在讨论中说明'))
replace_text(paragraphs[93],
 '监督标签还受轨迹长度、重复或窗口设计、对齐方式和侧链分箱影响；温度、溶剂、小分子配体状态和模拟时长均未作为模型输入，因此对于相同主链结构而观测条件不同的体系，不能给出条件特异的动态预测。统一轨迹处理、增加独立重复，以及联合学习与独立单目标训练的比较仍需补充。冻结测试预测的方差诊断显示，模型原生95%高斯区间对RMSF和Hχ的覆盖率分别为51.76%和85.19%。测试折内按序列关联分量互换拟合的事后校准将覆盖率提高至96.30%和94.84%，但RMSF平均区间宽度由0.876 Å增至6.469 Å（补充材料C）。这表明原生方差尚未校准，改善覆盖率也伴随明显的精度代价，不能将方差头直接解释为可靠的实验置信区间。')
after(paragraphs[16],'逐记录的PDB来源注释及其适用范围见补充材料A。')
replace_text(paragraphs[132],'A 数据来源与标签观测')
after(paragraphs[133],
 '为补充数据的生物学背景，我们通过RCSB PDB官方Data API获取输入结构的条目分类、分子实体描述、EC编号与UniProt对应[29]。3,063条记录对应2,751个PDB条目，其中2,748个条目可解析，覆盖3,060条记录；其余3条记录保留未知状态。逐记录注释与原始API响应随复现材料提供。PDB条目分类是来源提供的结构描述，可能包含多种标签，不等同于经过人工验证的统一功能分类，也不作为训练输入或据此推断功能预测能力。')
ref=next(p for p in doc.paragraphs if p.text.startswith('[28]'))
after(ref,'[29] RCSB Protein Data Bank. Data API documentation. https://data.rcsb.org/index.html. Accessed 10 October 2026.','Body Text')

content=(root/'uncertainty_baseline/uncertainty_results_zh.md').read_text(encoding='utf-8').split('\n\n')
doc.add_page_break()
doc.add_paragraph('C 方差诊断与事后分组校准','Heading 2')
methods=next(x for x in content if x.startswith('补充方法可用文字：')).split('：',1)[1]
doc.add_paragraph(methods)
for paragraph in content[:3]:doc.add_paragraph(paragraph)
caption=doc.add_paragraph('表S4：冻结五折测试预测的总体95%高斯区间覆盖率与平均宽度。覆盖率及宽度均先在记录内计算，再按记录等权汇总；数值为五折均值 ± 样本标准差。','Caption')
caption.paragraph_format.keep_with_next=True
caption.paragraph_format.page_break_before=True
raw=json.loads((root/'uncertainty_baseline/figures/table_uncertainty_95_compact.json').read_text(encoding='utf-8'))
table=doc.add_table(rows=1,cols=4)
table.autofit=False
widths=[1.00,1.80,1.75,1.75]
for column,width in zip(table.columns,widths):column.width=Inches(width)
for cell,label,width in zip(table.rows[0].cells,['指标','区间方案','覆盖率（%）','平均宽度'],widths):
    cell.text=label;cell.width=Inches(width)
for row in raw['rows']:
    if row[1]!='总体':continue
    for label,cov,width in [('原生模型方差',row[2],row[5]),('模型加标签变异',row[3],row[6]),('事后分组校准',row[4],row[7])]:
        cells=table.add_row().cells
        for cell,text,w in zip(cells,[row[0].split('（')[0],label,cov,width],widths):cell.text=text;cell.width=Inches(w)
for number,row in enumerate(table.rows):
    row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
    if number==0:row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
    for cell in row.cells:
        for paragraph in cell.paragraphs:
            paragraph.alignment=WD_ALIGN_PARAGRAPH.CENTER
            paragraph.paragraph_format.space_after=Pt(4)
            paragraph.paragraph_format.space_before=Pt(4)
            for run in paragraph.runs:
                run.font.size=Pt(10);run.font.bold=number==0
                run.font.name='Times New Roman'
                run._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'宋体')
doc.add_paragraph('RMSF宽度单位为Å，Hχ宽度单位为nats。标签变异不是均值标准误。两侧校准记录均未参与对应折模型训练，每个评价记录只使用另一侧关联分量拟合的系数；该事后诊断不改变均值预测，也不提供分布无关的覆盖保证。')
doc.add_page_break()
picture=doc.add_paragraph()
picture.paragraph_format.space_after=Pt(4)
picture.paragraph_format.keep_with_next=True
picture.add_run().add_picture(str(root/'uncertainty_baseline/figures/figS_uncertainty_calibration_corrected_v1.png'),width=Inches(6.30))
figure_caption=next(x for x in content if x.startswith('补充图图注：')).split('：',1)[1]
figure_caption=figure_caption.replace('图 S（编号由主稿统一）。','图S14：')
doc.add_paragraph(figure_caption,'Caption')

# Explicit cell borders avoid Word's cell/table border precedence differences.
def add_border(parent,side,value,size=0):
    e=OxmlElement('w:'+side);e.set(qn('w:val'),value)
    if size:e.set(qn('w:sz'),str(size));e.set(qn('w:color'),'000000');e.set(qn('w:space'),'0')
    parent.append(e)
for table in doc.tables:
    prop=table._tbl.tblPr
    for element in list(prop):
        if element.tag in {qn('w:tblBorders'),qn('w:shd')}:prop.remove(element)
    borders=OxmlElement('w:tblBorders')
    for side in ['top','left','bottom','right','insideH','insideV']:
        add_border(borders,side,'single' if side in {'top','bottom'} else 'nil',12 if side in {'top','bottom'} else 0)
    prop.append(borders)
    for i,row in enumerate(table.rows):
        for cell in row.cells:
            prop=cell._tc.get_or_add_tcPr()
            for element in list(prop):
                if element.tag in {qn('w:tcBorders'),qn('w:shd')}:prop.remove(element)
            borders=OxmlElement('w:tcBorders')
            for side in ['top','left','bottom','right','insideH','insideV']:
                outer=(i==0 and side=='top') or (i==len(table.rows)-1 and side=='bottom')
                header=i==0 and side=='bottom'
                add_border(borders,side,'single' if outer or header else 'nil',12 if outer else 6 if header else 0)
            prop.append(borders)
doc.save(output)

def archive_parts(path):
    with zipfile.ZipFile(path) as z:return {name:z.read(name) for name in z.namelist()}
before,after_parts=archive_parts(source),archive_parts(output)
images=[name for name in before if name.startswith('word/media/')]
assert all(before[name]==after_parts[name] for name in images)
assert before['word/comments.xml']==after_parts['word/comments.xml']
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
audit=dict(input=str(source),output=str(output),input_sha256=sha(source),output_sha256=sha(output),
    tables=len(doc.tables),old_images_preserved=len(images),new_images=1,native_comments_preserved=True,
    inserted_completed_analyses=['frozen OOF uncertainty diagnostics and group crossfit calibration','official PDB metadata annotation'],
    not_yet_inserted=['new single target training','all-chain family training','new external purge','new MD triplicates','new baselines','harmonized labels'])
(root/'completed_results_integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
