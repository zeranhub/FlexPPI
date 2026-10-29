"""Create a separate, native Word case atlas with all 32 publication-sized figures."""
from pathlib import Path
from collections import Counter
import json,hashlib,re,zipfile,math,csv
from docx import Document
from docx.shared import Mm,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from PIL import Image

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
OUT=BASE/'FlexPPI_32个潜在案例图册_背景证据与局限_20261021.docx'
QA=ROOT/'document_QA';QA.mkdir(exist_ok=True)
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf8')

def font(r,size=12,bold=False,east='宋体'):
    r.font.name='Times New Roman';r.font.size=Pt(size);r.bold=bold;r.font.color.rgb=RGBColor(0,0,0)
    rf=r._element.get_or_add_rPr().get_or_add_rFonts()
    for k,v in [('ascii','Times New Roman'),('hAnsi','Times New Roman'),('eastAsia',east)]:rf.set(qn('w:'+k),v)
    return r

doc=Document()
sec=doc.sections[0]
sec.page_width=Mm(210);sec.page_height=Mm(297)
sec.top_margin=Mm(15);sec.bottom_margin=Mm(15)
sec.left_margin=Mm(13.5);sec.right_margin=Mm(13.5)
sec.header_distance=Mm(4);sec.footer_distance=Mm(6)
for name in ['Normal','Title','Heading 1','Heading 2','Caption']:
    s=doc.styles[name]
    s.font.name='Times New Roman';s.font.color.rgb=RGBColor(0,0,0)
    s._element.get_or_add_rPr().get_or_add_rFonts().set(qn('w:eastAsia'),'宋体' if name in ['Normal','Caption'] else '黑体')
    p=s.paragraph_format;p.line_spacing=1.2;p.space_after=Pt(7)
    p.widow_control=True
    pr=s._element.find(qn('w:pPr'))
    if pr is not None:
        for tag in ['pBdr','shd']:
            e=pr.find(qn('w:'+tag))
            if e is not None:pr.remove(e)
doc.styles['Normal'].font.size=Pt(12)
doc.styles['Title'].font.size=Pt(18);doc.styles['Title'].font.bold=True
doc.styles['Heading 1'].font.size=Pt(14);doc.styles['Heading 1'].font.bold=True
doc.styles['Heading 2'].font.size=Pt(12);doc.styles['Heading 2'].font.bold=True
doc.styles['Caption'].font.size=Pt(10)

def para(text='',size=12,lead=None,style=None,before=False):
    p=doc.add_paragraph(style=style)
    p.paragraph_format.page_break_before=before
    p.paragraph_format.first_line_indent=Pt(24 if style is None and size==12 else 0)
    if lead:font(p.add_run(lead+' '),size,True)
    font(p.add_run(text),size)
    return p

def heading(text,before=False,level=1):
    p=doc.add_paragraph(style=f'Heading {level}')
    p.paragraph_format.page_break_before=before
    p.paragraph_format.keep_with_next=True
    font(p.add_run(text.replace('：',' ').replace(':',' ')),14 if level==1 else 12,True,'黑体')
    return p

def bookmark(p,name,idx):
    a=OxmlElement('w:bookmarkStart');a.set(qn('w:id'),str(idx));a.set(qn('w:name'),name)
    b=OxmlElement('w:bookmarkEnd');b.set(qn('w:id'),str(idx))
    p._p.insert(1 if p._p.find(qn('w:pPr')) is not None else 0,a);p._p.append(b)

def internal_link(p,label,anchor):
    h=OxmlElement('w:hyperlink');h.set(qn('w:anchor'),anchor)
    r=OxmlElement('w:r');pr=OxmlElement('w:rPr')
    color=OxmlElement('w:color');color.set(qn('w:val'),'000000');pr.append(color)
    sz=OxmlElement('w:sz');sz.set(qn('w:val'),'20');pr.append(sz);r.append(pr)
    t=OxmlElement('w:t');t.text=label;r.append(t);h.append(r);p._p.append(h)

def table(headers,rows,widths,size=10):
    t=doc.add_table(rows=1,cols=len(headers));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=False
    for c,w in zip(t.columns,widths):c.width=Mm(w)
    pr=t._tbl.tblPr
    old=pr.find(qn('w:tblBorders'))
    if old is not None:pr.remove(old)
    borders=OxmlElement('w:tblBorders')
    for key,val,sz in [('top','single','12'),('bottom','single','12'),('left','nil','0'),('right','nil','0'),('insideH','nil','0'),('insideV','nil','0')]:
        e=OxmlElement('w:'+key);e.set(qn('w:val'),val);e.set(qn('w:sz'),sz);e.set(qn('w:color'),'000000');borders.append(e)
    pr.append(borders)
    margin=OxmlElement('w:tblCellMar')
    for key,val in [('top',70),('bottom',70),('left',80),('right',80)]:
        e=OxmlElement('w:'+key);e.set(qn('w:w'),str(val));e.set(qn('w:type'),'dxa');margin.append(e)
    pr.append(margin)
    for j,h in enumerate(headers):t.rows[0].cells[j].text=str(h)
    trpr=t.rows[0]._tr.get_or_add_trPr();repeat=OxmlElement('w:tblHeader');repeat.set(qn('w:val'),'true');trpr.append(repeat)
    for row in rows:
        cells=t.add_row().cells
        for cell,value in zip(cells,row):cell.text=str(value)
    for i,row in enumerate(t.rows):
        for j,cell in enumerate(row.cells):
            cell.width=Mm(widths[j]);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
            cp=cell._tc.get_or_add_tcPr()
            if i==0:
                bd=OxmlElement('w:tcBorders');e=OxmlElement('w:bottom');e.set(qn('w:val'),'single');e.set(qn('w:sz'),'6');e.set(qn('w:color'),'000000');bd.append(e);cp.append(bd)
            for p in cell.paragraphs:
                p.paragraph_format.space_after=Pt(0);p.paragraph_format.space_before=Pt(0)
                p.paragraph_format.line_spacing=1.12;p.paragraph_format.keep_with_next=False
                if j==0:p.alignment=WD_ALIGN_PARAGRAPH.CENTER
                for r in p.runs:font(r,size,i==0)
    return t

def add_figure(cid,record):
    p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.page_break_before=True;p.paragraph_format.keep_with_next=True
    p.paragraph_format.space_after=Pt(4)
    inline=p.add_run().add_picture(record['png'],width=Mm(183))
    inline._inline.docPr.set('name',f'{cid} complete case figure')
    inline._inline.docPr.set('descr',f'{cid} evidence figure with actual structure and complete predefined dynamic predictions; caption specifies regions and limitations.')
    caption=Path(record['caption']).read_text(encoding='utf8').strip()
    caption=re.sub(r'^#+\s*','',caption)
    caption=re.sub(r'^图C\d+[：:\s]*','',caption)
    caption=re.sub(r'^C\d+\s+','',caption)
    caption=caption.replace('图4：','').replace('[28]','')
    cp=para('图'+cid+' '+caption,size=10,style='Caption')
    cp.paragraph_format.line_spacing=1.1;cp.paragraph_format.keep_together=True
    cp.paragraph_format.keep_with_next=False
    return caption

def scope_rows(cid,record):
    if cid=='C01':return [['RMSF','全部六个CDR','45位'],['Hχ','同类型pMHC界面','23位'],['接触几何','两态各3条轨迹','6条200 ns']]
    if cid in ['C02','C03','C04']:
        n={'C02':(13,11),'C03':(10,5),'C04':(9,7)}[cid]
        return [['RMSF','完整肽段',str(n[0])+'位'],['Hχ','两态χ有定义的肽位置',str(n[1])+'位'],['既有MD','两态各3条轨迹','6条200 ns']]
    if cid=='C31':return [['RMSF','四组成的完整KIX','4 × 87位'],['高活动集合','各状态最高20%','每态18位'],['专门MD','每组成1条轨迹','4条100 ns']]
    if cid=='C32':return [['Hχ','泛素全部有效位置','66位'],['逐角核对','χ1至χ4全部有效角','151角'],['角度曲线','Q2与Q41','3角 × 2位']]
    p=Path(record['audit'])
    audit=read(p) if p.exists() else {}
    a=audit.get('RMSF_valid_counts_by_contrast',audit.get('RMSF_rows','见图注'))
    b=audit.get('Hchi_valid_counts_by_contrast',audit.get('Hchi_rows','见图注'))
    def count_label(n):
        if isinstance(n,list):
            if len(n)==1:return str(n[0])+'位'
            if len(set(n))==1:return str(n[0])+'位 × '+str(len(n))+'对照'
            return '/'.join(str(x) for x in n)+'位'
        return str(n)+'位' if isinstance(n,int) else str(n)
    scope={'interface':'自然pMHC界面','cdr':'六个完整CDR','peptide':'完整肽段'}.get(audit.get('RMSF_scope'),'完整预设主链区域')
    if cid=='C17':scope='所有未替换共同位置'
    hchi_scope='未替换共同位置' if cid=='C17' else '同类型完整肽段' if cid=='C08' else '同类型完整界面'
    return [['RMSF',scope,count_label(a)],
            ['Hχ',hchi_scope,count_label(b)],
            ['已有MD','独立轨迹均值','每态10条' if cid=='C17' else '每态3条']]

def updated_discovery_summary(record):
    directory=Path(record['png']).parent
    summaries=[]
    for metric,name,unit in [('rmsf','主链','Å'),('hchi','侧链','nats')]:
        path=directory/f'panel_d_{metric}_all_contrasts.csv'
        rows=list(csv.DictReader(path.open(encoding='utf-8-sig')))
        observed=[float(r[metric+'_md_delta']) for r in rows]
        predicted=[float(r[metric+'_full_delta']) for r in rows]
        n=len(rows);assert n
        mae=sum(abs(a-b) for a,b in zip(observed,predicted))/n
        zero=sum(abs(a) for a in observed)/n
        gain=1-mae/zero if zero>1e-12 else None
        wording=('误差比全部预测不变低'+f'{gain*100:.1f}%' if gain is not None and gain>=0 else
                 '误差比全部预测不变高'+f'{-gain*100:.1f}%' if gain is not None else '参考变化近于零，未计算相对改善')
        summaries.append(f'{name}共{n}个配对位置，变化误差为{mae:.3f} {unit}，{wording}')
    return '更新的完整范围核对 '+ '；'.join(summaries)+'。这里合并图d列出的全部对照；位置数量表示比较范围，不是独立实验样本量。早期小集合中较高的排序相关，应与这些完整范围的幅度结果一起判断。'

def build():
    bundle=read(ROOT/'ATLAS_FIGURE_MAP.json')
    assert sha(bundle['main_manuscript'])==bundle['main_manuscript_sha256']
    narratives=read(ROOT/'evidence/narratives.json')
    cases=sorted(narratives['cases'],key=lambda c:c['catalog_case_id'])
    catalog=read(ROOT/'catalog/candidate_inventory.json')
    bycat={c['case_id']:c for c in catalog['cases']}
    assessments=read(ROOT/'evidence/case_assessment_table.json')
    p=doc.add_paragraph(style='Title');font(p.add_run('FlexPPI潜在案例图册'),18,True,'黑体')
    para('生物背景  模型证据  机制解释与适用边界',size=12)
    para('本图册汇集30个TCR–pMHC候选系列，以及KIX和泛素两个已有实例。MEL5亲和成熟与SB27/SB47的识别足迹差异提供较完整的双读数证据；KRAS G12D主要支持局部主链响应。其余系列既包括特定区域的正向结果，也包括信号较弱、模型版本间不一致、构建体差异及明确失败。候选的价值在于说明模型能够回答哪一种动态问题，并据此判断它是否适合进入论文。')
    para((ROOT/'evidence/reader_terms_paragraph.md').read_text(encoding='utf8').strip())
    para('TCR识别由肽与主要组织相容性分子共同形成的pMHC表面。人体的这类呈递分子称为HLA。受体的六条互补决定区（CDR）构成主要识别环，α、β两条链各含CDR1、CDR2和CDR3三条环。比较相同抗原上的不同TCR时，pMHC常有天然的一一残基对应；两个受体自身的同编号或同氨基酸位置，则还需要序列和结构关系支持。图中的链角色依据各结构单独识别，不以固定PDB链字母代替生物学身份。')
    heading('读图与证据判断',level=2)
    para('TCR图采用同一已有联合模型，其训练记录为ATLAS 1,866条与DynaRepo 564条，未使用ITrace。此前部分候选另用一个2,431条训练记录的历史版本复测，两者训练清单仅相差一条抗原单链记录；这类比较用于检查结论是否依赖已有训练版本，不能视为随机重复。KIX与泛素使用其对应的留出预测。模型曲线展示绝对值；变化核对中的每一点是残基位置，不是独立实验样本。')
    para('ITrace的已有MD通常为每状态三次独立重复，CatchBond为十次；KIX每状态一条专门轨迹。MD均值及重复差异提供模拟条件下的动态佐证。变化排序相关只衡量位置的顺序，变化误差则衡量幅度；“预测不变”对照将所有变化设为零，用于判断模型是否比忽略变化提供更多信息。未替换、同类型位置更适合检验环境响应；化学类型改变、缺失对应及无χ角位置分别标记。完整链与全部定义对照的数据均随图归档。')
    para('图C01至C17包含已有生物学复核及负例，图C18至C30为尚待完成机制核实的发现线索，C31和C32为非TCR实例。主对照与系列中的额外结构分别列明，不能把不同亲本、不同受体或不同HLA的结构串成一条工程路线。历史“SB27&R7”对应的有效对照是SB27/SB47及肽Q7，未列为独立受体案例。')
    heading('案例总览',before=True)
    para('表1汇总各例已有证据或可检验的问题，以及首先需要解决的缺口。“可检验”表示具有可用结构与预测的待验证问题，不能视为已证明的模型能力。编号可用于Word导航和图号检索。',size=11)
    t=table(assessments['headers'],assessments['rows'],[12,40,27,51,53],10)
    for i,row in enumerate(t.rows[1:]):
        p=row.cells[0].paragraphs[0];p.clear();internal_link(p,assessments['rows'][i][0],assessments['rows'][i][0])
    case_records=[]
    for i,c in enumerate(cases,1):
        cid=c['catalog_case_id'];record=bundle['figures'][cid]
        h=heading(cid+' '+c['title'],before=True);bookmark(h,cid,i)
        states=bycat.get(cid,{}).get('pdbs',c['pdbs'])
        para('结构条目 '+ '、'.join(s.upper() for s in states)+'。证据阶段 '+c['status']+'。',size=10.5)
        table(['读数或数据','图内分析范围','规模'],scope_rows(cid,record),[32,109,42],10)
        roles=c.get('paragraph_roles',['生物背景与对照','模型实际证据','机制与局限'])
        for j,text in enumerate(c['paragraphs']):
            # Paragraphs remain prose. Bold leads organize the distinct evidential roles.
            lead=roles[j] if j<len(roles) else '具体限制与未解决问题'
            p=para(text,lead=lead)
            if j==0:
                font(p.add_run(' ['+','.join(str(n) for n in bundle['case_references'][cid])+']'),10)
        if cid=='C02':para('同一研究的CA5结构4JRX作为第三受体背景保留；本例两态定量对照仅为2AK4与4JRY，不将CA5解释为SB27或SB47的工程衍生物。',size=10.5)
        scope_clarifications={
            'C06':'本例正文的7位侧链统计针对整条肽；图b、c另展示43位完整CDR主链和12位pMHC界面侧链，两者回答不同区域的响应问题。',
            'C12':'本例正文的23位主链及16位侧链统计均针对CDR3；图b、c分别展示47位完整CDR主链和14位pMHC界面侧链，两组数量对应不同分析范围。',
            'C15':'本例正文的22位统计针对两条CDR3；图b展开全部45位CDR，并保留三个抗原的所有列明对照。',
        }
        if cid in scope_clarifications:para(scope_clarifications[cid],size=10.5)
        if cid in [f'C{k:02}' for k in range(18,31)]:
            para('图中采用更新的完整肽段与同类型pMHC界面范围；早期少量跨受体同类型位置的筛选数值仅用于保留线索。整条链采用全局序列配准，所有空缺及化学替换均保留于源表，因此不能将早期小集合统计直接当作本图完整区域的成绩。',size=10.5)
            para(updated_discovery_summary(record),size=10.5)
        caption=add_figure(cid,record)
        case_records.append({'case_id':cid,'title':c['title'],'paragraphs':c['paragraphs'],'caption':caption,
                             'figure_sha256':record['sha256']['png'],'references':bundle['case_references'][cid]})
    heading('跨来源辅助检查',before=True)
    para('以下辅助比较涉及不同模拟来源或不同结合伙伴，作为相应案例的来源敏感性记录保留。它们不是新增的单变量生物学实验，也没有另行建立机制链。')
    for row in catalog['auxiliary_cross_source_checks']:
        text=row['name'].replace(' versus ','与')+'。'
        if 'rmsf_rho' in row:text+='RMSF变化排序相关为'+f"{row['rmsf_rho']:.3f}"+'。'
        if 'Hchi_rho' in row:text+='Hχ变化排序相关为'+f"{row['Hchi_rho']:.3f}"+'。'
        text+='跨来源的力场、采样和观测窗口不同，这些数值不能直接归因于TCR结合伙伴本身。'
        para(text)
    heading('参考文献与结构条目',before=True)
    for ref in bundle['references']:
        suffix=''
        if ref.get('verification'):suffix=' '+ref['verification']+'。'
        text=f"[{ref['number']}] {ref['title']}. "
        if ref.get('doi'):text+='doi:'+ref['doi']+'. '
        if ref.get('url'):text+=ref['url']+'.'
        para(text+suffix,size=10)
    heading('发现性结构配对索引',before=True)
    para('本索引保留此前数值检索涉及的343个结构对，便于定位尚未形成独立生物学叙述的原始线索。同一系列内可以出现多个对照；结构对数量不是独立案例或验证体系数量。完整分类、序列、条目来源与数值在对应源表中保留。')
    discovery=read(ROOT/'catalog/all_discovery_pair_index.json')
    rows=[]
    for i in range(0,len(discovery),3):
        row=[]
        for j in range(3):
            k=i+j
            row += [str(k+1),discovery[k]['a'].upper()+' / '+discovery[k]['b'].upper()] if k<len(discovery) else ['','']
        rows.append(row)
    table(['序号','结构对','序号','结构对','序号','结构对'],rows,[12,49,12,49,12,49],10)
    foot=sec.footer.paragraphs[0];foot.alignment=WD_ALIGN_PARAGRAPH.CENTER
    f=OxmlElement('w:fldSimple');f.set(qn('w:instr'),'PAGE');foot._p.append(f)
    doc.core_properties.title='FlexPPI潜在案例图册'
    doc.core_properties.subject='32个案例的背景 模型证据 机制与局限'
    doc.core_properties.author='';doc.core_properties.last_modified_by=''
    doc.save(OUT)
    with zipfile.ZipFile(OUT) as z:
        assert z.testzip() is None
        from lxml import etree
        xml=etree.fromstring(z.read('word/document.xml'))
        ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
        assert len(xml.xpath('.//a:blip',namespaces=ns))==32
    assert sha(bundle['main_manuscript'])==bundle['main_manuscript_sha256']
    audit={'status':'built_pending_render_and_review','document':str(OUT),'document_sha256':sha(OUT),
           'case_count':32,'figure_count':32,'reference_count':len(bundle['references']),
           'table_count':len(doc.tables),'main_manuscript_unchanged':True,'case_records':case_records,
           'source_manifest_sha256':sha(ROOT/'ATLAS_FIGURE_MAP.json'),
           'new_training':0,'new_inference':0,'new_MD':0}
    write(QA/'integration_audit.json',audit)
    (QA/'document_content.txt').write_text('\n'.join(p.text for p in doc.paragraphs),encoding='utf8')
    print(json.dumps({k:audit[k] for k in ['status','document','document_sha256','case_count','figure_count','table_count']},ensure_ascii=False))

if __name__=='__main__':build()
