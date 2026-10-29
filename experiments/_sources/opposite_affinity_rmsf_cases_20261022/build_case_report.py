"""Create a separate Chinese case-analysis report; the main paper is read-only."""
from pathlib import Path
import re
from docx import Document
from docx.shared import Cm, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_ALIGN_VERTICAL
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'局部柔性下降与相反亲和结果_案例分析_20261022.md'
OUTPUT=ROOT/'局部柔性下降与相反亲和结果_案例分析_三线表_20261022.docx'

def font(run,size=None,bold=None):
    run.font.name='Arial'
    run._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'宋体')
    if size: run.font.size=Pt(size)
    if bold is not None: run.bold=bold
    run.font.color.rgb=RGBColor(0,0,0)

def add_hyperlink(paragraph,text,url):
    from docx.opc.constants import RELATIONSHIP_TYPE as RT
    rid=paragraph.part.relate_to(url,RT.HYPERLINK,is_external=True)
    link=OxmlElement('w:hyperlink');link.set(qn('r:id'),rid)
    rr=OxmlElement('w:r');pr=OxmlElement('w:rPr')
    f=OxmlElement('w:rFonts');f.set(qn('w:ascii'),'Arial');f.set(qn('w:hAnsi'),'Arial');f.set(qn('w:eastAsia'),'宋体');pr.append(f)
    color=OxmlElement('w:color');color.set(qn('w:val'),'000000');pr.append(color)
    underline=OxmlElement('w:u');underline.set(qn('w:val'),'single');pr.append(underline)
    rr.append(pr);tt=OxmlElement('w:t');tt.text=text;rr.append(tt);link.append(rr);paragraph._p.append(link)

def add_text(paragraph,text,size=None):
    position=0
    for match in re.finditer(r'\[([^\]]+)\]\((https?://[^)]+)\)',text):
        if match.start()>position:font(paragraph.add_run(text[position:match.start()]),size)
        add_hyperlink(paragraph,match.group(1),match.group(2));position=match.end()
    if position<len(text):font(paragraph.add_run(text[position:]),size)

def table(doc,rows):
    t=doc.add_table(rows=1,cols=len(rows[0]));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=False
    widths=[3.7,4.4,3.0,2.4,2.5]
    for i,cell in enumerate(t.rows[0].cells):cell.width=Cm(widths[i])
    for i,row in enumerate(rows):
        cells=t.rows[0].cells if i==0 else t.add_row().cells
        for j,(cell,value) in enumerate(zip(cells,row)):
            cell.width=Cm(widths[j]);cell.vertical_alignment=WD_ALIGN_VERTICAL.CENTER
            p=cell.paragraphs[0];p.paragraph_format.space_after=Pt(5);p.paragraph_format.space_before=Pt(5);p.paragraph_format.line_spacing=1.15
            p.alignment=WD_ALIGN_PARAGRAPH.LEFT if j<3 else WD_ALIGN_PARAGRAPH.CENTER
            add_text(p,value,9)
            for run in p.runs:run.bold=i==0
            borders=OxmlElement('w:tcBorders')
            for side in ['top','left','bottom','right','insideH','insideV']:
                e=OxmlElement(f'w:{side}');e.set(qn('w:val'),'nil');borders.append(e)
            if i==0:
                for side,sz in [('top','12'),('bottom','6')]:
                    e=OxmlElement(f'w:{side}');e.set(qn('w:val'),'single');e.set(qn('w:sz'),sz);e.set(qn('w:color'),'000000');borders.append(e)
            if i==len(rows)-1:
                e=OxmlElement('w:bottom');e.set(qn('w:val'),'single');e.set(qn('w:sz'),'12');e.set(qn('w:color'),'000000');borders.append(e)
            cell._tc.get_or_add_tcPr().append(borders)
        cant=OxmlElement('w:cantSplit');t.rows[i]._tr.get_or_add_trPr().append(cant)
    header=OxmlElement('w:tblHeader');t.rows[0]._tr.get_or_add_trPr().append(header)
    doc.add_paragraph().paragraph_format.space_after=Pt(0)

def main():
    doc=Document();s=doc.sections[0]
    s.page_width=Cm(21);s.page_height=Cm(29.7)
    s.top_margin=Cm(1.8);s.bottom_margin=Cm(1.8);s.left_margin=Cm(1.8);s.right_margin=Cm(1.8)
    normal=doc.styles['Normal'];normal.font.name='Arial';normal.font.size=Pt(10.5)
    normal.element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'宋体')
    normal.paragraph_format.line_spacing=1.25;normal.paragraph_format.space_after=Pt(7)
    for name in ['Title','Heading 1','Heading 2']:
        style=doc.styles[name];style.font.name='Arial';style.font.color.rgb=RGBColor(0,0,0)
        style.element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'黑体')
        style.paragraph_format.keep_with_next=True
    doc.styles['Title'].font.size=Pt(19);doc.styles['Heading 1'].font.size=Pt(13)
    doc.styles['Heading 1'].paragraph_format.space_before=Pt(10)
    lines=SOURCE.read_text(encoding='utf8').splitlines();i=0;table_count=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        if line.startswith('# '):
            p=doc.add_paragraph(line[2:],style='Title')
        elif line.startswith('## '):
            if line=='## 第4位识别面变化被第10位最大值掩盖':
                doc.add_page_break()
                p=doc.add_paragraph('MEL5四态案例对照图',style='Heading 1')
                p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER
                p.add_run().add_picture(str(ROOT/'case_figure/MEL5_opposite_affinity_case_series.png'),width=Cm(17.1))
                p=doc.add_paragraph('图1　同一MEL5谱系的亲和力与完整空间动态预测。a为文献SPR；b为四态全10肽位点模型RMSF；c为全部45个CDR位置的三组变化；d为四态界面并集30个位点中同类型有效Hχ变化，灰色表示不适合比较。变化均为后态减参照态，蓝色降低、红色升高。图件用于比较空间预测模式，不表示已建立亲和力分类器；MD支持及不一致结果见正文。')
                p.paragraph_format.line_spacing=1.1
                for run in p.runs:font(run,9)
                doc.add_page_break()
            doc.add_paragraph(line[3:],style='Heading 1')
        elif line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                values=[v.strip() for v in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r':?-+:?',v) for v in values):rows.append(values)
                i+=1
            table(doc,rows);table_count+=1;continue
        else:
            p=doc.add_paragraph();add_text(p,line)
            p.paragraph_format.widow_control=True
        i+=1
    footer=s.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.CENTER
    r=footer.add_run();field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');r._r.addnext(field)
    doc.core_properties.title='局部柔性下降与相反亲和结果的案例分析'
    doc.core_properties.subject='MEL5四态与A6 G4A的结构动力学和结合实验分析'
    doc.core_properties.author=''
    # Remove the template's inherited blue Title rule and any paragraph borders.
    for style in doc.styles:
        for border in list(style.element.findall('.//' + qn('w:pBdr'))):
            border.getparent().remove(border)
    for paragraph in doc.paragraphs:
        for border in list(paragraph._p.findall('.//' + qn('w:pBdr'))):
            border.getparent().remove(border)
    doc.save(OUTPUT)
    print(OUTPUT);print(f'Native three-line tables: {table_count}; main manuscript unchanged.')

if __name__=='__main__':main()
