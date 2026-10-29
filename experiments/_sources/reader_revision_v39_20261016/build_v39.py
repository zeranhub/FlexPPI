from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib, json, re, zipfile, struct
from lxml import etree
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.table import _Cell
from docx.text.paragraph import Paragraph

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent.parent
QA = ROOT/'manuscript_QA'; QA.mkdir(exist_ok=True)
SOURCE = BASE/'paper_v38_20261016_中文修订_图7标注说明与图例调整_三线表.docx'
OUT = BASE/'paper_v39_20261017_中文修订_读者表述与英文图表_三线表.docx'
d = Document(SOURCE)
paras = list(d.paragraphs)
changes=[]
sha=lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

def setp(i, text):
    p=paras[i]; before=p.text
    assert not p._p.xpath('.//m:oMath'), i
    nodes=p._p.xpath('.//w:t'); assert nodes, i
    nodes[0].text=text
    for n in nodes[1:]: n.text=''
    if before!=text: changes.append({'original_paragraph':i,'before':before,'after':text})

def replace_nodes(i, old, new):
    p=paras[i]; before=p.text; count=0
    for n in p._p.xpath('.//w:t'):
        if old in (n.text or ''):
            count+=n.text.count(old); n.text=n.text.replace(old,new)
    assert count==1,(i,old,count)
    changes.append({'original_paragraph':i,'before':before,'after':p.text})

setp(2, paras[2].text.replace('在ThermoPCD外部集合中，完整模型的Hχ误差相对链内模型降低9.4%。','').replace('冻结折外预测','折外预测'))
setp(9, '我们整合ATLAS、DynaRepo与ITrace的轨迹标签，在五折体系留出评价中比较完整模型、链内模型、静态特征与类型均值对照[16,17,18]。整体及类型内成绩检验几何输入和局部环境的预测信息，并与已有RMSF预测方法比较；ThermoPCD数据用于补充评价另一轨迹来源中的预测表现。KIX和泛素–UBA实例分别对照活动区域定位与同类型侧链差异。在83条记录中，将折外预测与MD低相对自由能角区总宽度及细粒度角度占据对应，检验模型排序及其选出的高低位置组能否区分实际取向范围与占据多样性。各项评价共同检验两个结构读数的准确性及其在局部动力学分析中的含义。')
replace_nodes(23,'Pro按冻结标签协议使用两个角','Pro使用χ1和χ2两个角')
setp(40, '外部评价采用ThermoPCD的32条抗体–抗原复合物轨迹记录，涵盖生理与发热温度范围[20,21]。纳入记录在建立标签和计分前确定，均有可追溯的初始拓扑、原始轨迹、温度及随机种子信息。训练集合按上述全组成链序列与体系关联关系排除相关记录；MMseqs2搜索敏感度为7.5，其余关联条件与家族划分一致。')
setp(41, '排除与外部集合关联的633条训练记录后，完整模型、链内模型、静态特征回归与类型均值对照均使用其余2,430条记录。完整与链内采用种子7训练40轮的单模型。32条外部记录属于同一个关联家族，结果逐记录等权汇总。')
replace_nodes(44,'BackFlip2采用作者发布的2.1代码、提交887c6a7及发布权重，对相同静态结构逐链推断，提供RMSF对照。','BackFlip-2采用2.1版本（代码提交887c6a7），对相同静态结构逐链预测RMSF。')
replace_nodes(44,'发布模型保留原训练权重；静态特征及类型均值使用本文各对应训练集合拟合。','')
setp(45, 'PLANET-MD[11]比较同时使用序列和结构特征的full配置，以及仅使用序列的sequence-only配置。各条链独立预测后合并对应记录；采用默认300 K与float32精度，RMSF输出由nm转换为Å。RMSF-net[9]采用仅结构输入的PDB-only配置，将完整蛋白结构转换成4 Å分辨率、1.5 Å网格的原子序数加权高斯模拟密度，评价其相对柔性分数的排序及记录内标准化剖面误差。这两类方法与ANM、GNM均覆盖3,063条记录中的961,321个RMSF有效位置。已有方法直接使用作者提供的模型；静态特征和类型均值对照使用本文对应训练集合拟合。')
setp(49, paras[49].text.replace('模型固定使用第1折留出权重，相关记录A00AQ位于该折测试集合。','预测由第1折模型生成，该折训练不包含KIX相关记录A00AQ。'))
setp(51, paras[51].text.replace('冻结预测','预测').replace('每个χ层内，预测及MD指标分别减去相同氨基酸类型均值，再用于散点和Spearman相关计算。','在每个χ角层内，分别减去同一种氨基酸在预测值与MD指标中的平均值，再用于散点和Spearman相关计算。'))
setp(54, paras[54].text.replace('冻结预测','模型预测'))
setp(56, paras[56].text.replace('冻结折外预测','折外预测'))
setp(58, paras[58].text.replace('83条冻结记录','83条记录'))
setp(65, '表1：主链波动与侧链联合状态熵的预测表现。A为3,063条记录的五折体系留出评价，B为ThermoPCD的32条外部记录。各项先逐记录计算，再对记录等权平均。MAE越低越好，Spearman相关ρ越高越好。')
setp(67, '在已有RMSF方法中，BackFlip-2和PLANET-MD (full)的平均绝对误差分别为0.537和0.568 Å，低于FlexPPI的0.577 Å；两者的Spearman相关也较高。PLANET-MD (full)与FlexPPI的RMSE分别为1.186和0.792 Å，反映两者在较大误差上的差异。共同位置、弹性网络及来源分层比较见补充材料C。')
setp(85, '在ThermoPCD的32条外部记录中，完整模型Hχ MAE为0.208 nats，相对链内模型和静态特征回归分别降低9.4%和1.5%，对应排序相关由0.884和0.886提高至0.895（表1B）。RMSF MAE为0.474 Å，与两项对照的0.471 Å接近。该集合中，完整结构输入的收益主要体现在Hχ。')
paras[71]._p.addnext(paras[85]._p)
setp(91, paras[91].text.replace('固定留出模型','留出模型'))
setp(96, paras[96].text.replace('冻结折外输出','折外预测'))
setp(98, paras[98].text.replace('参考与预测的排序相关为0.979；分别减去训练折的氨基酸类型均值后，相关仍为0.868（图7a）。','参考与预测的排序相关为0.979（图7a）。将预测值与MD参考值均减去相应氨基酸在训练集中的平均Hχ后，相关仍为0.868。'))
setp(101, '3.6.2 侧链预测对低自由能取向范围与角度占据的验证')
setp(103, paras[103].text.replace('冻结侧链预测','侧链预测'))
setp(106, paras[106].text.replace('折线为原始Hχ；扣除类型均值后的ρ由预测与MD参考分别减去训练折中该氨基酸类型的同一个均值后计算。','折线与所标ρ均对应原始Hχ。').replace('第0折冻结完整模型','第0折完整模型'))
setp(112, paras[112].text.replace('冻结折外预测','折外预测'))
setp(157, '弹性网络对照采用ProDy 2.6.1的ANM与GNM，由完整输入的Cα构建，距离截断为15 Å、弹簧常数为1，保留最低20个非零模态；对平方波动之和开方，得到相对柔性分数[4,27]。相同961,321个有效RMSF位置上，各项先逐记录计算，再等权平均（表S4）。')
setp(158, paras[158].text.replace('发布模型与相对柔性基线','RMSF预测方法与相对柔性基线'))
setp(167, paras[167].text.replace('冻结预测','折外预测'))
setp(170, paras[170].text.replace('固定留出模型','留出模型'))

# Remove the independent external section after retaining its complete results in Table 1B.
removed=[]
for i in (84,86,87,88,89):
    p=paras[i]
    assert not p._p.xpath('.//w:commentRangeStart|.//w:commentRangeEnd|.//w:commentReference|.//m:oMath'),i
    removed.append({'original_paragraph':i,'text':p.text,'images':len(p._p.xpath('.//a:blip'))})
    p._p.getparent().remove(p._p)
original_tables=list(d.tables)
source_main=[[c.text for c in r.cells] for r in original_tables[0].rows]
source_ext=[[c.text for c in r.cells] for r in original_tables[1].rows]
assert not original_tables[1]._tbl.xpath('.//w:commentRangeStart|.//w:commentRangeEnd|.//w:commentReference|.//m:oMath')
original_tables[1]._tbl.getparent().remove(original_tables[1]._tbl)

# Native Word table: methods aligned at the left, related metrics grouped at the top.
table=d.add_table(rows=2, cols=5)
original_tables[0]._tbl.addprevious(table._tbl)
original_tables[0]._tbl.getparent().remove(original_tables[0]._tbl)
table.alignment=WD_TABLE_ALIGNMENT.CENTER; table.autofit=False
widths=[2.8,.85,.85,.85,.85]
for c,w in zip(table.columns,widths): c.width=Inches(w)
table.cell(0,0).merge(table.cell(1,0)).text='方法'
table.cell(0,1).merge(table.cell(0,2)).text='RMSF (Å)'
table.cell(0,3).merge(table.cell(0,4)).text='Hχ (nats)'
for j,t in enumerate(['MAE ↓','ρ ↑','MAE ↓','ρ ↑'],1): table.cell(1,j).text=t

name_map={
 '氨基酸类型均值':'氨基酸类型均值', '静态特征模型':'静态特征回归',
 'RMSF-net PDB-only':'RMSF-net (PDB-only)',
 'PLANET seq 发布权重':'PLANET-MD (sequence)', 'PLANET full 发布权重':'PLANET-MD (full)',
 'BackFlip2发布模型':'BackFlip-2','仅链内模型':'FlexPPI (Intra)','完整模型':'FlexPPI (Full)'}
by_a={r[0]:r[1:] for r in source_main[1:]}; by_b={r[0]:r[1:] for r in source_ext[1:]}
order_a=['氨基酸类型均值','静态特征模型','RMSF-net PDB-only','PLANET seq 发布权重','PLANET full 发布权重','BackFlip2发布模型','仅链内模型','完整模型']
order_b=['氨基酸类型均值','静态特征模型','BackFlip2发布模型','仅链内模型','完整模型']
panels=[]; numeric_audit=[]
for title,order,lookup in [('A  五折体系留出（3,063条记录）',order_a,by_a),('B  ThermoPCD（32条外部记录）',order_b,by_b)]:
    row=table.add_row(); row.cells[0].merge(row.cells[-1]).text=title
    panels.append(len(table.rows)-1)
    for name in order:
        row=table.add_row(); values=lookup[name]
        texts=[name_map[name]]+['—' if x in ('无','不适用') else x for x in values]
        for c,t in zip(row.cells,texts):c.text=t
        numeric_audit.append({'panel':title[0],'method':name_map[name],'source_values':values,'display_values':texts[1:]})

tblpr=table._tbl.tblPr
tblpr.find(qn('w:tblW')).set(qn('w:type'),'dxa')
tblpr.find(qn('w:tblW')).set(qn('w:w'),str(int(6.2*1440)))
for key in ('tblBorders','tblCellMar','tblCaption'):
    for e in list(tblpr):
        if e.tag==qn('w:'+key):tblpr.remove(e)
mar=OxmlElement('w:tblCellMar')
for side,val in [('top',40),('bottom',40),('left',90),('right',90)]:
    e=OxmlElement('w:'+side);e.set(qn('w:w'),str(val));e.set(qn('w:type'),'dxa');mar.append(e)
tblpr.append(mar)
def add_border(parent,side,val,size=0):
    b=OxmlElement('w:'+side); b.set(qn('w:val'),val)
    if size:b.set(qn('w:sz'),str(size));b.set(qn('w:color'),'000000')
    parent.append(b)
for ri,row in enumerate(table.rows):
    trpr=row._tr.get_or_add_trPr()
    no=OxmlElement('w:cantSplit');trpr.append(no)
    if ri<2:
        repeat=OxmlElement('w:tblHeader');trpr.append(repeat)
    physical_cells=[]; grid_column=0
    for tc in row._tr.tc_lst:
        physical_cells.append((grid_column,_Cell(tc,table)))
        grid_column+=tc.grid_span
    for ci,cell in physical_cells:
        cell.width=Inches(sum(widths[ci:ci+cell._tc.grid_span]))
        cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
        cp=cell._tc.get_or_add_tcPr()
        for el in list(cp):
            if el.tag in (qn('w:tcBorders'),qn('w:shd')):cp.remove(el)
        borders=OxmlElement('w:tcBorders')
        for side in ('top','bottom','left','right','insideH','insideV'):
            outer=(ri==0 and side=='top') or (ri==len(table.rows)-1 and side=='bottom')
            header=(ri==1 and side=='bottom')
            add_border(borders,side,'single' if outer or header else 'nil',12 if outer else 6 if header else 0)
        cp.append(borders)
        for p in cell.paragraphs:
            p.alignment=WD_ALIGN_PARAGRAPH.LEFT if ci==0 or ri in panels else WD_ALIGN_PARAGRAPH.CENTER
            pf=p.paragraph_format;pf.first_line_indent=Inches(0);pf.left_indent=Inches(0);pf.right_indent=Inches(0)
            pf.space_before=Pt(5 if ri in panels else 2);pf.space_after=Pt(3 if ri in panels else 2)
            pf.line_spacing=1.0;pf.keep_with_next=ri<len(table.rows)-1;pf.keep_together=True
            if ci==0 and ri not in panels and ri>1 and cell.text in ('RMSF-net (PDB-only)','FlexPPI (Intra)'):
                pf.space_before=Pt(6)
            for run in p.runs:
                run.font.name='Times New Roman';run.font.size=Pt(10.5)
                run._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'宋体')
                run.bold=ri<2 or ri in panels or (ci==0 and cell.text=='FlexPPI (Full)')

note_el=OxmlElement('w:p');table._tbl.addnext(note_el)
note=Paragraph(note_el,d._body)
note.text='注：Full为完整结构模型，Intra为链内模型。B中本文模型及基础对照使用排除关联记录后的2,430条训练记录；已有方法采用作者提供的模型。—表示不提供该目标或其绝对幅值。'
note.paragraph_format.first_line_indent=Pt(0);note.paragraph_format.space_before=Pt(4);note.paragraph_format.space_after=Pt(8)
note.paragraph_format.line_spacing=1.0
for r in note.runs:r.font.size=Pt(9.5)

# Keep method names consistent in the supplementary comparison as well.
replacement_names={'BackFlip2发布模型':'BackFlip-2','PLANET full 发布权重':'PLANET-MD (full)','PLANET seq 发布权重':'PLANET-MD (sequence)','RMSF-net PDB-only':'RMSF-net (PDB-only)'}
for t in d.tables[1:]:
    for row in t.rows:
        for cell in row.cells:
            for node in cell._tc.xpath('.//w:t'):
                for old,new in replacement_names.items():node.text=(node.text or '').replace(old,new)

# Close numbering gaps left by the removed external figure and section.
for p in d.paragraphs:
    before=p.text
    for n in p._p.xpath('.//w:t'):
        n.text=re.sub(r'图([67])',lambda m:'图'+str(int(m.group(1))-1),n.text or '')
        n.text=re.sub(r'3\.([56])(?=\.|\s)',lambda m:'3.'+str(int(m.group(1))-1),n.text)
        n.text=n.text.replace('BackFlip2','BackFlip-2')
    if p.text!=before:changes.append({'renumber_or_name':True,'before':before,'after':p.text})

# Figure replacement manifest keys retain their original numbers for unambiguous provenance.
manifest=ROOT/'figure_manifest.json'
image_replacements=[]
if manifest.exists():
    for key,item in json.loads(manifest.read_text(encoding='utf8')).items():
        old_no=int(key) if key!='S1' else 'S1'
        new_no=old_no-1 if isinstance(old_no,int) and old_no>=6 else old_no
        prefix='图'+str(new_no)+'：'
        cap=next(p for p in d.paragraphs if p.text.startswith(prefix))
        img=cap._p.getprevious()
        assert len(img.xpath('.//a:blip'))==1,(key,prefix)
        # Replace the media part directly to retain all drawing dimensions and anchors.
        blip=img.xpath('.//a:blip')[0]
        part=d.part.related_parts[blip.get(qn('r:embed'))]
        path=Path(item['png']);part._blob=path.read_bytes()
        width_px,height_px=struct.unpack('>II',part._blob[16:24])
        extent=img.xpath('.//wp:extent')[0]
        new_height=round(int(extent.get('cx'))*height_px/width_px)
        extent.set('cy',str(new_height))
        for shape_extent in img.xpath('.//a:xfrm/a:ext'):shape_extent.set('cy',str(new_height))
        image_replacements.append({'original_figure':key,'new_figure':new_no,'png':str(path),'sha256':sha(path)})

d.save(OUT)
# python-docx is conservative with unknown parts; retain the comments XML byte for byte.
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    files=[(inf,a.read(inf.filename) if inf.filename=='word/comments.xml' else b.read(inf.filename)) for inf in b.infolist()]
temp=OUT.with_suffix('.preserve.tmp')
with zipfile.ZipFile(temp,'w') as z:
    for inf,data in files:z.writestr(inf,data)
temp.replace(OUT)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    ax=etree.fromstring(a.read('word/document.xml'));bx=etree.fromstring(b.read('word/document.xml'))
    anchors=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    maths=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    assert anchors(ax)==anchors(bx)
    assert a.read('word/comments.xml')==b.read('word/comments.xml')
    assert maths(ax)==maths(bx) and len(maths(bx))==36
    assert len(bx.xpath('.//w:tbl',namespaces=ns))==10
    assert len(bx.xpath('.//a:blip',namespaces=ns))==7
    assert all(not t.xpath('.//w:commentRangeStart|.//m:oMath') for t in (original_tables[0]._tbl,original_tables[1]._tbl))
current=Document(OUT)
fulltext='\n'.join(p.text for p in current.paragraphs)+'\n'+'\n'.join(c.text for t in current.tables for r in t.rows for c in r.cells)
assert not any(term in fulltext for term in ('冻结','发布权重','发布模型','扣除氨基酸类型均值后','3.6 ','图7：','表2：'))
(QA/'v39_content.txt').write_text('\n'.join(f'P{i}: '+p.text for i,p in enumerate(current.paragraphs) if p.text),encoding='utf8')
audit={'source':str(SOURCE),'source_sha256':sha(SOURCE),'output':str(OUT),'output_sha256':sha(OUT),'changes':changes,'removed_external_paragraphs':removed,'external_values_fully_retained':True,'table1_numeric_audit':numeric_audit,'figures':image_replacements,'math_elements_preserved':36,'comments_and_anchors_preserved':True,'main_figures':6,'supplementary_figures':1,'native_tables':10,'render_status':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k not in ('changes','removed_external_paragraphs','table1_numeric_audit','figures')},ensure_ascii=False))
