"""Add completed ANM/GNM and paired alternate-label results to Chinese revision."""
from pathlib import Path
from copy import deepcopy
import json,hashlib,zipfile
from docx import Document
from docx.shared import Inches,Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
import pandas as pd
root=Path(__file__).resolve().parent
source=root.parent/'paper_v16_20261010_中文修订_已完成分析_三线表.docx'
output=root.parent/'paper_v17_20261010_中文修订_基线与标签敏感性_三线表.docx'
doc=Document(source)
reference_anchor=next(p for p in doc.paragraphs if p.text.startswith('[29]'))
for reference in [
 '[30] Ahmet Bakan, Lidio M. Meireles, and Ivet Bahar. ProDy: protein dynamics inferred from theory and experiments. Bioinformatics, 27(11):1575–1577, 2011.',
 '[31] She Zhang, James M. Krieger, Yan Zhang, Cihan Kaya, Burak Kaynak, Karolina Mikulska-Ruminska, Pemra Doruker, Hongchun Li, and Ivet Bahar. ProDy 2.0: increased scale and scope after 10 years of protein dynamics modelling with Python. Bioinformatics, 37(20):3657–3659, 2021. doi:10.1093/bioinformatics/btab187.'
]:
    new_reference=doc.add_paragraph(reference)
    reference_anchor._p.addnext(new_reference._p)
    reference_anchor=new_reference

def replace(p,text):
    for child in list(p._p):
        if child.tag==qn('w:r') and not child.findall('.//'+qn('w:commentReference')):p._p.remove(child)
    p.add_run(text)
def after(p,text,style='Normal'):
    new=doc.add_paragraph(text,style);p._p.addnext(new._p);return new
def table(caption,columns,rows,widths):
    p=doc.add_paragraph(caption,'Caption');p.paragraph_format.keep_with_next=True
    t=doc.add_table(rows=1,cols=len(columns));t.autofit=False
    for col,w in zip(t.columns,widths):col.width=Inches(w)
    for cell,label in zip(t.rows[0].cells,columns):cell.text=label
    for values in rows:
        for cell,value in zip(t.add_row().cells,values):cell.text=str(value)
    for i,row in enumerate(t.rows):
        row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
        if i==0:row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
        for cell,w in zip(row.cells,widths):
            cell.width=Inches(w)
            for p in cell.paragraphs:
                p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.space_before=Pt(4);p.paragraph_format.space_after=Pt(4)
                p.paragraph_format.first_line_indent=Pt(0)
                p.paragraph_format.left_indent=Pt(0);p.paragraph_format.right_indent=Pt(0)
                for run in p.runs:
                    run.font.size=Pt(10);run.font.bold=i==0;run.font.name='Times New Roman'
                    run._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'宋体')
    return t

anchor=next(p for p in doc.paragraphs if p.text.startswith('改以生物体系等权汇总后'))
replace(anchor,anchor.text+'ANM/GNM对照见补充材料D。')
p=next(p for p in doc.paragraphs if p.text.startswith('监督标签还受轨迹长度'))
replace(p,p.text.replace('统一轨迹处理、增加独立重复，以及联合学习与独立单目标训练的比较仍需补充。',
 '按来源预先抽样的共同时间窗口敏感性分析，在83条可比较记录上统一使用100 ns、200 ps间隔、逐链对齐和第一条轨迹，完整模型相对静态对照的两项MAE收益在三个来源中仍保持，但绝对误差与完整和链内模型的差距随标签处理改变（补充材料E）。该分析不统一力场或温度，也不等同于重新训练的全队列比较；增加独立重复及联合与独立单目标训练的比较仍需补充。'))

doc.add_page_break()
doc.add_paragraph('D 弹性网络基线','Heading 1')
doc.add_paragraph('使用ProDy 2.6.1[30,31]对冻结输入的完整多肽组成构建ANM和GNM，保留原始Cα顺序；两种网络的距离截断均为15 Å，弹簧常数为1，使用最低的20个非零模态，极短结构以可用模态数为上限。相对柔性分数为模态平方波动总和的平方根。全部3,063条结构均计算成功，无失败记录；有效计分位置为961,321个。每项分数先在记录内计算，再对记录等权平均，复合物中的链在同一记录内汇总。弹性网络输出不经过训练标签的幅值标定。标准化剖面误差指每条记录的参考与预测分别减去自身均值、除以自身标准差后，再计算逐位置平均绝对差，因而不带物理单位；最高10%位点命中率用于描述活动区域筛选。')
raw=pd.read_csv(root/'baselines/prody_common_evaluation/record_macro_metrics.csv')
raw=raw[raw.dataset=='all'].set_index('method')
names={'flexppi':'完整FlexPPI','intrachain':'链内FlexPPI','anm':'ANM','gnm':'GNM'}
rows=[]
for method in ['flexppi','intrachain','anm','gnm']:
    r=raw.loc[method];rows.append([names[method],f'{r.spearman:.3f}',f'{r.pearson:.3f}',f'{r.z_mae:.3f}',f'{r.top10_precision*100:.1f}'])
table('表S5：完整队列的弹性网络对照。各方法使用完全相同的有效RMSF位置；相关和标准化误差为记录等权均值，最高10%位点命中率单位为%。',
 ['方法','Spearman ρ','Pearson r','标准化MAE','前10%命中率'],rows,[1.50,1.15,1.10,1.25,1.30])
doc.add_paragraph('ANM、GNM的整体Spearman相关分别为0.715和0.663，完整模型为0.796；标准化MAE分别为0.450、0.487和0.371，最高10%位点命中率分别为56.7%、57.7%和65.1%。在ATLAS中，ANM和完整模型的排序相关相近，分别为0.785与0.788；DynaRepo分别为0.635与0.774，ITrace分别为0.453与0.938。由此可见，几何学习相对弹性网络的收益具有来源依赖；这些指标评价当前MD标签下的柔性剖面，不能把它们解释为对所有正常模分析用途的优劣判断。')

doc.add_page_break()
doc.add_paragraph('E 共同时间窗口的标签敏感性','Heading 1')
doc.add_paragraph('在查看重算标签和评价变化之前，按来源内样本标识排序，以固定随机种子20261010从ATLAS、DynaRepo、ITrace各无放回抽取30条记录，冻结90条名单。逐条读取原始轨迹的全部时间戳及文件哈希后发现，7条DynaRepo记录的保存时长约1 ns。依据时间戳可用性，在计算替代标签前将至少10 ns和时间戳单调作为可比较条件；这些短轨迹按该条件排除而不替换，最终保留30、23、30条，共83条。其余203条可用轨迹均支持0至100 ns的共同窗口；因ITrace保存间隔为200 ps，统一采用200 ps的实际观测网格，共501帧，不插值。对相同时间点的重复帧固定取首次保存帧；没有倒序时间戳。排除名单、每帧时间戳及共同方案均在计算替代标签前保存。')
doc.add_paragraph('主敏感性方案对每条记录只使用第一条轨迹，使窗口长度、采样间隔、每条记录的轨迹数及估计方法相同。各链全部可解析Cα逐帧Kabsch拟合至共同窗口首帧，计算围绕对齐均值的RMSF；Hχ使用与主分析相同的完整χ原子定义、三个角度分区、联合状态和自然对数。Ala、Gly保持无侧链监督，Pro按原冻结定义保留。完整、链内、静态、类型均值及BackFlip2预测均保持冻结，不重新训练，分别在原标签和新标签的完全相同有效交集计分。83条记录均重算成功；所有残基对应完整，RMSF及Hχ有效位置分别为38,599和33,493，旧新掩码无位置增减，未发现超过5 Å的χ共价键距离。另报告ATLAS和ITrace三条可用轨迹的共同窗口均值，DynaRepo保持一条轨迹，以区分重复数变化的影响。')
summary=pd.read_csv(root/'label_sensitivity/prediction_effect_record_mean_sample_sd.csv')
summary=summary[summary.variant=='uniform_rep1'].set_index(['dataset','metric'])
rows=[]
for metric,label in [('rmsf','RMSF'),('chi_entropy','Hχ')]:
    for ds,source_label,n in [('atlas','ATLAS',30),('dynarepo','DynaRepo',23),('itrace','ITrace',30)]:
        r=summary.loc[(ds,metric)]
        rows.append([label,source_label,n,f'{r.old_mae_mean:.3f}',f'{r.new_mae_mean:.3f}',f'{r.old_spearman_mean:.3f}',f'{r.new_spearman_mean:.3f}'])
table('表S6：83条可比较记录中完整模型在原标签与共同窗口第一条轨迹标签上的表现。误差单位为RMSF的Å和Hχ的nats；ρ为记录内Spearman相关的记录等权均值。此处为抽样子集结果，不能替代全部3,063条记录的主分数。',
 ['指标','来源','条数','原MAE','共同MAE','原ρ','共同ρ'],rows,[.70,1.05,.55,1.0,1.0,1.0,1.0])
doc.add_paragraph('统一使用第一条轨迹后，三个来源的RMSF与Hχ MAE均提高，记录内排序相关均降低（表S6）。完整模型相对静态对照的两项平均误差仍均较低；但较小的模型差距会变化，例如抽样ATLAS的RMSF原标签上链内模型误差较低，共同第一条轨迹方案下完整模型较低。因此，时间窗和重复数会影响绝对分数与小幅模型差异，不能把整体优势表述为每条记录或每个对照均稳定。')
doc.add_paragraph('使用三条可用轨迹均值后，ATLAS的RMSF/Hχ MAE为0.723 Å/0.172 nats，接近其原标签结果；ITrace为0.158 Å/0.110 nats，其RMSF误差相对原标签反而降低。DynaRepo因仅有一条本地轨迹，与主敏感性方案相同。该结果显示，标签变化并非所有来源中同方向的误差偏移。统一估计器、窗口和采样不能消除力场、温度、溶剂、蛋白组成或真实动态的差异；该子集分析也没有重新训练模型，更不能据此证明跨来源模拟条件已完全一致。逐记录标签变化、所有对照的配对分数和掩码审计随复现材料提供。')

for t in doc.tables:
    prop=t._tbl.tblPr
    for e in list(prop):
        if e.tag in {qn('w:tblBorders'),qn('w:shd')}:prop.remove(e)
    for i,row in enumerate(t.rows):
        for cell in row.cells:
            cp=cell._tc.get_or_add_tcPr()
            for e in list(cp):
                if e.tag in {qn('w:tcBorders'),qn('w:shd')}:cp.remove(e)
            borders=OxmlElement('w:tcBorders')
            for side in ['top','left','bottom','right','insideH','insideV']:
                outer=(i==0 and side=='top') or (i==len(t.rows)-1 and side=='bottom');head=i==0 and side=='bottom'
                e=OxmlElement('w:'+side);e.set(qn('w:val'),'single' if outer or head else 'nil')
                if outer or head:e.set(qn('w:sz'),str(12 if outer else 6));e.set(qn('w:color'),'000000');e.set(qn('w:space'),'0')
                borders.append(e)
            cp.append(borders)
ref_source=next(p for p in doc.paragraphs if p.text.startswith('[28]'))
for number in [29,30,31]:
    ref_target=next(p for p in doc.paragraphs if p.text.startswith(f'[{number}]'))
    if ref_target._p.pPr is not None:ref_target._p.remove(ref_target._p.pPr)
    if ref_source._p.pPr is not None:ref_target._p.insert(0,deepcopy(ref_source._p.pPr))
    if ref_source.runs and ref_source.runs[0]._r.rPr is not None:
        for run in ref_target.runs:
            if run._r.rPr is not None:run._r.remove(run._r.rPr)
            run._r.insert(0,deepcopy(ref_source.runs[0]._r.rPr))
doc.save(output)
with zipfile.ZipFile(source) as a,zipfile.ZipFile(output) as b:
    assert a.read('word/comments.xml')==b.read('word/comments.xml')
    assert all(a.read(name)==b.read(name) for name in a.namelist() if name.startswith('word/media/'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit=dict(source=str(source),output=str(output),output_sha256=sha(output),tables=len(doc.tables),
   completed_added=['ANM GNM all3063','common window labels 83 paired control comparisons'],
   input_audits={'prody':sha(root/'baselines/prody_common_evaluation/scoring_audit.json'),
    'labels':sha(root/'label_sensitivity/verification_audit.json'),'controls':sha(root/'label_sensitivity/controls_analysis_audit.json')},
   pending=['single-target 10fits','family10fits','new2430 external neural','RMSFnet','PLANET','KIX12 new trajectories'])
(root/'baseline_sensitivity_integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(audit,ensure_ascii=False))
