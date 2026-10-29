"""Integrate completed frozen scientific outputs into the Chinese manuscript."""
from pathlib import Path
from copy import deepcopy
import json, hashlib, zipfile
from lxml import etree
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT = Path(__file__).resolve().parent
RESULT = ROOT / 'completed_results_20261012'
SOURCE = ROOT.parent / 'paper_v18_20261010_中文修订_复用KIX单轨迹_三线表.docx'
OUTPUT = ROOT.parent / 'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'
doc = Document(SOURCE)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources = {}

def read(p):
    sources[str(p)] = sha(p)
    return pd.read_csv(p)

def replace(p, text):
    for child in list(p._p):
        if child.tag == qn('w:r') and not child.findall('.//' + qn('w:commentReference')):
            p._p.remove(child)
    p.add_run(text)

def find(start):
    items = [p for p in doc.paragraphs if p.text.startswith(start)]
    assert len(items) == 1, (start, len(items))
    return items[0]

def rep(start, text):
    replace(find(start), text)

def after(anchor, text, style='Normal'):
    p = doc.add_paragraph(text, style)
    anchor._p.addnext(p._p)
    return p

def heading(text):
    p = doc.add_paragraph(text, 'Heading 1')
    p.paragraph_format.page_break_before = True
    return p

def table(caption, columns, values, widths):
    p = doc.add_paragraph(caption, 'Caption')
    p.paragraph_format.keep_with_next = True
    t = doc.add_table(rows=1, cols=len(columns))
    t.autofit = False
    for col, w in zip(t.columns, widths):
        col.width = Inches(w)
    for cell, value in zip(t.rows[0].cells, columns):
        cell.text = value
    for values_row in values:
        for cell, value in zip(t.add_row().cells, values_row):
            cell.text = str(value)
    for i, row in enumerate(t.rows):
        row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
        if i == 0:
            row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
        for cell, w in zip(row.cells, widths):
            cell.width = Inches(w)
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_before = Pt(4)
                p.paragraph_format.space_after = Pt(4)
                p.paragraph_format.first_line_indent = Pt(0)
                p.paragraph_format.left_indent = Pt(0)
                p.paragraph_format.right_indent = Pt(0)
                for r in p.runs:
                    r.font.size = Pt(10)
                    r.font.bold = i == 0
                    r.font.name = 'Times New Roman'
                    r._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'), '宋体')
    return t

single = read(RESULT / 'summary/single_target_overall_record_macro.csv').set_index(['variant', 'metric'])
family = read(RESULT / 'summary/family_overall_record_macro.csv').set_index(['variant', 'metric'])
family_source = read(RESULT / 'summary/family_source_record_macro.csv').set_index(['variant', 'dataset', 'metric'])
external = read(RESULT / 'summary/external_2430_record_macro.csv').set_index(['variant', 'metric'])
baseline = read(ROOT / 'baselines/additional_original_oof/record_macro_metrics.csv').set_index(['method', 'dataset'])
baseline_family = read(ROOT / 'baselines/additional_family_oof/record_macro_metrics.csv').set_index(['method', 'dataset'])
assert len(single) == 4 and len(external) == 8
old_external_rows = [[c.text for c in row.cells] for row in doc.tables[1].rows]

abstract = find('蛋白质的局部动态既包括')
replace(abstract, abstract.text.replace('五折评价中', '原体系留出的五折评价中').replace(
    '去除相关训练体系后的外部验证进一步支持侧链状态预测的收益，主链幅值预测仍受到数据分布变化的影响。',
    '全组成链家族隔离后成绩下降，但完整模型相对链内模型仍保留平均收益；严格外部验证支持侧链状态预测的小幅收益，主链幅值误差未优于静态对照。联合训练相对单目标的平均误差收益较小且依赖指标。'))
related = find('早期的结构动态分析主要')
replace(related, related.text.replace('这些方法的输入需求和输出目标各有侧重',
    'PLANET-MD进一步从序列和可选结构表示预测RMSF、广义动态相关和结构异质性分布，其异质性分布与本文联合离散χ角状态熵具有不同定义[31]。这些方法的输入需求和输出目标各有侧重'))
overview = find('我们从三个互补层面评价')
replace(overview, overview.text.replace('并以BackFlip2发布模型检验RMSF预测的适用范围。',
    '并以BackFlip2、PLANET-MD及RMSF-net发布模型检验RMSF预测的适用范围；另以单目标对照及全组成链家族隔离评价检验联合训练收益和泛化范围。'))

split = find('序列相似性与重复结构可能')
after(split, '另构建全组成链家族隔离方案，对每条输入链执行关闭低复杂度屏蔽的MMseqs2序列搜索，采用至少30%一致性及70%双向覆盖，并将所有长度链的精确序列重复和同一来源生物体系关系共同连接。只要任意组成链关联，整个复合物所属记录就通过连通关系归为同组。3,063条训练记录形成1,662个连通家族，按记录数确定性平衡为633、608、608、607、607条的五折；同一连通家族及生物体系不跨折，标签与预测不参与分配。32条外部记录仅用于识别关联关系，不进入训练或测试记录。完整与链内模型均在新划分上从头训练40轮、种子7，以对应留出折预测计分。全部244条ITrace属于同一连通家族，仅在第0折测试，该折训练不含ITrace，因而同时涉及家族留出和来源外推，不能为该来源构造五折标准差。')
search = find('外部验证使用ThermoPCD')
replace(search, search.text.replace('使用MMseqs2搜索所有长度至少20个残基的链',
    '采用上述全组成链家族关系，使用关闭低复杂度屏蔽的MMseqs2搜索所有长度至少20个残基的链'))
rep('严格方案移除所有与外部队列关联', '严格外部方案删除与固定32条外部记录关联的633条训练记录，保留2,430条。完整与链内模型使用全链家族第0折在这2,430条记录上从头训练的单模型，种子7、第40轮权重；静态特征与类型均值对照在同一保留集合上拟合。外部链与保留训练集合在声明的搜索和精确序列条件下无命中。32条外部记录构成一个关联家族，按记录进行描述性汇总，不视作32个独立蛋白家族。普通原训练范围五模型平均、旧关联组删除632条的2,431条训练方案保留为补充敏感性；这些方案的关联规则、训练范围或模型数量不同，不用于单因素归因，也不按预测成绩选择报告。')
method_baseline = find('BackFlip2采用作者发布')
after(method_baseline, 'PLANET-MD[31]使用公开full及sequence-only权重，对每条链独立推断后合并对应记录。full输入序列和结构特征，采用权重对应的encoded结构编码；sequence-only不输入结构特征。保持发布模型默认300 K及float32精度，原RMSF输出nm乘以10转换为Å，不用测试标签拟合倍率。其结构异质性分布不等同于Hχ，因此仅评价RMSF。RMSF-net[10]采用PDB-only发布权重，将完整蛋白结构转换成4 Å分辨率、1.5 Å网格的原子序数加权高斯模拟密度；该适配未使用真实冷冻电镜密度，也未与原Chimera操作逐体素验证。RMSF-net输出保持相对分数，只评价排序及记录内标准化剖面误差，不赋予Å单位。PLANET、RMSF-net、ANM及GNM全部共同覆盖3,063条记录和961,321个RMSF有效位置；发布权重的训练语料及交叠不同于本文拟合对照，不解释为相同训练数据的架构比较。')
rep('表1：五折体系留出评价', '表1：原五折体系留出评价及发布模型对照。各项先逐记录计算，再对记录等权平均。MAE越低越好，ρ为Spearman相关；RMSF单位为Å，Hχ为nats。静态特征与类型均值使用对应折训练记录拟合；BackFlip2、PLANET及RMSF-net保留原发布权重，训练语料不同。无表示不输出该目标，不适用表示相对分数不具有Å幅值。')
for name, method in [('PLANET full 发布权重', 'planet_full'), ('PLANET seq 发布权重', 'planet_seq'), ('RMSF-net PDB-only', 'rmsf_net_pdb')]:
    b = baseline.loc[(method, 'all')]
    row = deepcopy(doc.tables[0].rows[-1]._tr)
    doc.tables[0]._tbl.append(row)
    for c, value in zip(doc.tables[0].rows[-1].cells, [name, f'{b.mae:.3f}' if pd.notna(b.mae) else '不适用', f'{b.spearman:.3f}', '无', '无']):
        c.text = value
    for c in doc.tables[0].rows[-1].cells:
        for p in c.paragraphs:
            for r in p.runs:
                r.font.size = Pt(10)
                r.font.name = 'Times New Roman'
                r._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'), '宋体')
anchor = find('完整模型相对仅链内模型的两项')
anchor = after(anchor, 'PLANET full的整体RMSF MAE为0.568 Å、ρ为0.828，优于完整FlexPPI的0.577 Å和0.796；但RMSE为1.186 Å，高于FlexPPI的0.792 Å，平均绝对误差与大误差尾部呈现不同结果。PLANET sequence-only的MAE为0.714 Å、ρ为0.771，RMSF-net PDB-only的ρ为0.658。BackFlip2仍具有较低总体RMSF MAE和较高排序相关。PLANET full在ATLAS的MAE/ρ为0.398 Å/0.876，在DynaRepo为0.704 Å/0.752，在ITrace为1.412 Å/0.732，来源差异明显（补充材料I）。未知训练交叠、独立链输入、默认温度和标签处理可能共同影响比较，不能只由来源差异推断单一机制。')
anchor = after(anchor, '单目标对照保持相同输入、编码器、留出划分和40轮训练预算，仅分别预测RMSF或Hχ。联合模型的MAE由单目标0.582降至0.577 Å、由0.169降至0.166 nats，平均收益约0.8%和1.8%；RMSF排序相关由0.800变为0.796，Hχ由0.938升至0.941。部分来源单目标误差更低，RMSF仅49.1%的记录在联合模型中MAE较低。因此，两类读数的互补用途与联合训练的性能收益应分开评价（补充材料G）。')
after(anchor, '全组成链家族隔离后，完整模型RMSF/Hχ MAE为0.616 Å/0.173 nats，链内模型为0.666 Å/0.180 nats，完整空间环境仍有平均收益；但相对于原体系留出，两项目标成绩下降。全部ITrace只在第0折测试，完整模型RMSF排序相关由原体系留出约0.938降至0.680；该变化同时涉及家族隔离和来源外推，不能全部归因于某一泄漏因素，也不能把原成绩解释为严格家族泛化（补充材料H）。')

for i, variant in enumerate(['full', 'intrachain', 'static_sequence_complex', 'aa_only'], start=1):
    values = [f'{external.loc[(variant, metric), column]:.3f}' for metric, column in [('rmsf', 'mae'), ('rmsf', 'spearman'), ('chi_entropy', 'mae'), ('chi_entropy', 'spearman')]]
    for c, v in zip(doc.tables[1].rows[i].cells[1:], values):
        c.text = v
rep('表2：去同源条件下', '表2：全组成链关联排除后的ThermoPCD外部验证。删除633条相关记录、保留2,430条训练记录；完整与链内为同一保留集合上的第0折单模型，静态特征与类型均值在同集合拟合。BackFlip2保留原发布权重，在相同外部记录和掩码上计分。各项为记录等权均值，MAE单位为Å或nats，ρ为Spearman相关；无表示不输出该目标。')
rep('在严格外部验证中，完整模型相对', '在全组成链关联排除后的严格外部验证中，完整模型在Hχ误差及两项目标的残基排序上优于链内和静态对照，但RMSF幅值误差略高。完整、链内和静态模型的RMSF MAE分别为0.474、0.471和0.471 Å，Hχ MAE为0.208、0.229和0.211 nats（表2，图4）。')
rep('具体而言，完整模型的RMSF MAE', '完整模型相对静态对照的Hχ误差降幅约1.5%，RMSF平均误差则高约0.003 Å；即使19条记录的RMSF误差较低，全队列平均仍略高。因此，外部结果支持部分动态读数的迁移和排序能力，尚不支持主链幅值或所有目标上的统一优势。BackFlip2在相同外部记录上的RMSF误差较高，与其主评价中较好的整体成绩形成来源依赖的结果。32条外部记录属于一个关联家族，上述均值和胜出记录数仅作描述，不进行独立家族统计推断。')
figure = RESULT / 'figures/fig4_strict_external_paired.png'
caption = find('图4：移除相关训练记录后')
previous = caption._p.getprevious()
while previous is not None and not previous.findall('.//' + qn('w:drawing')):
    previous = previous.getprevious()
assert previous is not None
inline = doc.add_paragraph()
inline.paragraph_format.keep_with_next = True
inline.add_run().add_picture(str(figure), width=Inches(6.3))
previous.getparent().replace(previous, inline._p)
caption.paragraph_format.keep_with_next = False
replace(caption, (RESULT / 'figures/fig4_strict_external_paired_caption_zh.md').read_text(encoding='utf8').replace('各为种子7的单模型', '各为第0折、种子7、第40轮的单模型'))

discussion = find('FlexPPI的主要价值在于')
replace(discussion, discussion.text.replace('BackFlip2在主评价中的总体RMSF成绩更好', 'BackFlip2在主评价中的总体RMSF成绩更好，PLANET full的MAE和排序也优于FlexPPI而RMSE较高'))
rep('完整复合物输入在两项任务中', '完整复合物输入在原体系留出及全组成链家族划分中均提供平均预测收益，但完整与链内模型分别训练，比较的是两种信息范围的建模方案；严格外部的RMSF误差并未优于链内或静态对照。联合目标对照显示平均MAE收益较小，RMSF排序及部分来源的误差并未改善。伙伴几何的平均价值、未见家族泛化与结合前后响应是需要分别检验的问题，不能由较好的绝对柔性图相互替代。KIX状态差值也显示活动区定位不保证响应方向准确。')
rep('数据划分也决定结论的适用范围', '数据划分决定结论的适用范围。原体系留出允许ITrace组成链家族跨折出现；全组成链连通家族隔离后的成绩下降，尤其ITrace主链排序由约0.938降至0.680。ITrace仅一个依赖家族、在单折测试且该折训练不含该来源，因此同时改变了家族和来源可见性。严格ThermoPCD外部排除减少序列关联，但其32条记录仍只有一个关联家族，不足以覆盖广泛独立家族。发布权重的训练交叠未经逐条审计，其分数不能作为相同家族划分下从头训练的结果；扩大独立家族和来源覆盖仍是进一步验证的方向。')
p = find('监督标签还受轨迹长度')
replace(p, p.text.replace('联合学习相对于独立单目标模型的贡献仍需单独检验。', '单目标对照已显示联合训练只有较小、依赖指标和来源的平均收益；两类动态读数的互补性不等同于普遍多任务协同。'))
conclusion = find('FlexPPI从静态主链结构共同估计')
replace(conclusion, conclusion.text.replace('体系留出和去同源外部评价表明，几何表示能够提供静态结构摘要之外的预测信息，但主链幅值预测仍有改进空间。', '原体系留出、全链家族隔离和严格外部评价共同显示，完整几何可以提供附加预测信息，但收益依赖目标、来源和验证条件；联合训练相对单目标的改善较小，主链幅值与结合响应仍有改进空间。'))

for start in ['图S1：数据组成', '图S10：外部验证', '图S11：普通五模型']:
    p = find(start)
    replace(p, p.text + ' 本图保留旧关联规则的历史敏感性结果；旧严格方案删除632条、保留2,431条，不代表主文删除633条、保留2,430条的新方案。')
for p in doc.paragraphs:
    if '[30,31]' in p.text:
        replace(p, p.text.replace('[30,31]', '[3,30]'))
reference = find('[31] She Zhang')
replace(reference, '[31] Samuel Sledzieski and Sonya M. Hanson. PLANET-MD: Ultra-fast Proteome-scale Prediction of Allosteric Networks in Proteins. bioRxiv preprint, 2026, 2025.06.12.659353, version 2. doi:10.1101/2025.06.12.659353.')
reference_style = find('[30] Ahmet Bakan')
if reference_style.runs[0]._r.rPr is not None:
    for run in reference.runs:
        if run._r.rPr is not None:
            run._r.remove(run._r.rPr)
        run._r.insert(0, deepcopy(reference_style.runs[0]._r.rPr))
for row in doc.tables[1].rows[1:5]:
    for col, cell in enumerate(row.cells[1:], start=1):
        sample_p = doc.tables[1].rows[-1].cells[col].paragraphs[0]
        for p in cell.paragraphs:
            if p._p.pPr is not None:
                p._p.remove(p._p.pPr)
            if sample_p._p.pPr is not None:
                p._p.insert(0, deepcopy(sample_p._p.pPr))
            for run in p.runs:
                run.font.size = Pt(10)
                run.font.name = 'Times New Roman'
                run._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'), '宋体')

heading('G 联合与单目标训练对照')
doc.add_paragraph('在原五折留出划分上分别从头训练仅RMSF及仅Hχ模型，各五个，共十次训练。两组均采用种子7、第40轮权重，与联合模型保持几何输入、编码器、来源与体系平衡抽样、优化器和训练轮数一致。每条折外预测对应的输入、标签、测试折和有效残基掩码逐位核对相同；RMSF与Hχ分别覆盖3,063条记录中的961,321和824,988个位置。各项先逐记录计分，再对记录等权平均；本对照使用固定种子，不把五折当成五次独立初始化。')
rows = []
for metric, label, arm in [('rmsf', 'RMSF', 'rmsf_only'), ('chi_entropy', 'Hχ', 'chi_entropy_only')]:
    for model, v in [('单目标', arm), ('联合目标', 'full')]:
        r = single.loc[(v, metric)]
        rows.append([label, model, f'{r.mae:.3f}', f'{r.rmse:.3f}', f'{r.spearman:.3f}'])
table('表S9：单目标与联合目标的记录等权结果。RMSF误差单位为Å，Hχ为nats，ρ为Spearman相关。各目标在两种模型间使用相同有效位置；本表不混用等权折均值与折间标准差。', ['目标', '训练方案', 'MAE', 'RMSE', 'ρ'], rows, [.8, 1.55, 1.25, 1.25, 1.45])
doc.add_paragraph('联合模型的RMSF MAE降低0.004621 Å、Hχ MAE降低0.002989 nats，约为单目标误差的0.8%与1.8%；但RMSF RMSE略升，排序相关略低。RMSF有1,503/3,063条记录、Hχ有1,863/3,063条记录的联合模型MAE较低。DynaRepo的RMSF及ITrace的Hχ单目标平均误差更低。上述观察为描述性平均收益，不证明每个来源或残基均获得改善，也不提供独立种子重复的显著性结论。')

heading('H 全组成链家族隔离评价')
doc.add_paragraph('全部1,662个训练连通家族均不跨折，完整与链内模型共完成十次40轮训练，每条记录仅使用对应留出折的单模型预测。五折大小及来源组成不同，主表采用3,063条记录等权均值，折均值与标准差另保存在源表，不将二者混配。共同有效残基与原体系留出一致。表S10列出整体和来源分层结果；ITrace仅在一个折中测试，不报告该来源五折标准差。')
rows = []
for dataset, name in [('all', '整体'), ('atlas', 'ATLAS'), ('dynarepo', 'DynaRepo'), ('itrace', 'ITrace')]:
    for v, label in [('full', '完整'), ('intrachain', '链内')]:
        r = [family.loc[(v, m)] if dataset == 'all' else family_source.loc[(v, dataset, m)] for m in ['rmsf', 'chi_entropy']]
        rows.append([name, label, f'{r[0].mae:.3f}', f'{r[0].spearman:.3f}', f'{r[1].mae:.3f}', f'{r[1].spearman:.3f}'])
table('表S10：全组成链家族隔离的记录等权结果。RMSF MAE为Å，Hχ MAE为nats；ρ为记录内Spearman相关的均值。ITrace共244条、属于同一依赖家族，仅第0折测试。', ['来源', '模型', 'RMSF MAE', 'RMSF ρ', 'Hχ MAE', 'Hχ ρ'], rows, [1.10, .7, 1.2, 1.1, 1.15, 1.05])
doc.add_paragraph('在家族方案中，完整模型相对链内模型的RMSF/Hχ MAE降低约7.5%/3.8%，两项排序也较好；相对于原体系留出，完整模型RMSF MAE由0.577升至0.616 Å，Hχ由0.166升至0.173 nats，排序相关分别由0.796/0.941降至0.760/0.936。ITrace的RMSF/Hχ MAE为0.378 Å/0.182 nats，RMSF排序为0.680。共享链形成的大连通组使来源分布改变，因此这些变化描述更严格外推条件，不能单独量化某一关联边或泄漏因素的作用。')

heading('I 发布模型与相对柔性基线')
doc.add_paragraph('PLANET full和sequence-only各完成4,973条链，RMSF-net完成3,063条完整结构输入；ANM/GNM也完整覆盖同一队列。各方法没有因失败或误差筛除记录。七方法共同计分为3,063条记录、961,321个RMSF位置；先汇总同一记录中的有效链，再对记录等权平均。PLANET、RMSF-net、ANM和GNM的预测在原划分与家族方案中不改变，后者仅替换完整与链内模型的留出预测，不意味着发布模型已按本文家族划分重训。')
rows = []
names = {'flexppi': '完整FlexPPI', 'intrachain': '链内FlexPPI', 'planet_full': 'PLANET full', 'planet_seq': 'PLANET seq', 'rmsf_net_pdb': 'RMSF-net PDB', 'anm': 'ANM', 'gnm': 'GNM'}
for m in ['flexppi', 'intrachain', 'planet_full', 'planet_seq', 'rmsf_net_pdb', 'anm', 'gnm']:
    a, b = baseline.loc[(m, 'all')], baseline_family.loc[(m, 'all')]
    rows.append([names[m], f'{a.mae:.3f}' if pd.notna(a.mae) else '不适用', f'{a.spearman:.3f}', f'{b.mae:.3f}' if pd.notna(b.mae) else '不适用', f'{b.spearman:.3f}', f'{a.z_mae:.3f}'])
table('表S11：共同RMSF位置上的发布模型及相对基线。误差单位为Å，ρ为Spearman；标准化MAE分别将每条记录的参考与预测减自身均值并除自身标准差，单位为无量纲。末列对应原体系留出；家族完整/链内标准化MAE为0.401/0.434，其余方法不变。', ['方法', '原MAE', '原ρ', '家族MAE', '家族ρ', '标准化MAE'], rows, [1.5, .9, .8, 1.0, .8, 1.3])
doc.add_paragraph('PLANET full在ATLAS的MAE/ρ为0.398 Å/0.876，DynaRepo为0.704 Å/0.752，ITrace为1.412 Å/0.732；sequence-only分别为0.499/0.826、0.782/0.668、2.171/0.719。PLANET full在ITrace的预测减参考平均偏差为+1.381 Å，sequence-only为+2.150 Å。完整来源分层和逐记录结果随复现材料提供；ATLAS占全部记录约63%，整体成绩不能替代来源评价。PLANET使用原发布默认温度、权重及独立链输入，训练交叠未逐条核实。RMSF-net是高斯模拟密度适配，未取得真实实验密度或原Chimera逐体素一致性证明，因此不能把这里的分数当作原冷冻电镜密度任务的直接复现。')

heading('J 旧严格外部协议敏感性')
doc.add_paragraph('原关联规则删除632条训练记录、保留2,431条；全组成链搜索新增关联后，主文删除633条、保留2,430条。两者外部队列保持32条不变，均属于一个关联家族。旧结果仅作协议敏感性保留（表S12及图S1、S10、S11）；新主文表2、图4采用2,430条方案。两种方案的关联规则、训练集合和训练运行同时变化，不能用差值单独归因一条新增关联或去同源效应。普通原训练范围五模型集成另见表S2，亦不与两个严格单模型方案混称。')
table('表S12：旧关联规则的2,431条训练方案，仅作历史敏感性。RMSF MAE为Å，Hχ MAE为nats；每项为相同32条外部记录的等权均值。主文采用更新的2,430条方案。', old_external_rows[0], old_external_rows[1:], [1.5, 1.25, 1.15, 1.25, 1.15])

for t in doc.tables:
    for e in list(t._tbl.tblPr):
        if e.tag in {qn('w:tblBorders'), qn('w:shd')}:
            t._tbl.tblPr.remove(e)
    for i, row in enumerate(t.rows):
        for cell in row.cells:
            cp = cell._tc.get_or_add_tcPr()
            for e in list(cp):
                if e.tag in {qn('w:tcBorders'), qn('w:shd')}:
                    cp.remove(e)
            borders = OxmlElement('w:tcBorders')
            for side in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
                outer = (i == 0 and side == 'top') or (i == len(t.rows)-1 and side == 'bottom')
                head = i == 0 and side == 'bottom'
                e = OxmlElement('w:' + side)
                e.set(qn('w:val'), 'single' if outer or head else 'nil')
                if outer or head:
                    e.set(qn('w:sz'), str(12 if outer else 6)); e.set(qn('w:color'), '000000'); e.set(qn('w:space'), '0')
                borders.append(e)
            cp.append(borders)
for p in doc.paragraphs:
    if '对每条输入链执行关闭低复杂度屏蔽的MMseqs2序列搜索' in p.text:
        replace(p, p.text.replace('对每条输入链执行关闭低复杂度屏蔽的MMseqs2序列搜索', '对长度至少20个残基的输入链执行关闭低复杂度屏蔽的MMseqs2序列搜索，短链按精确序列连接'))
    if p.text.startswith('外部验证主要检验模型能否迁移'):
        replace(p, '外部验证主要检验模型能否迁移到另一轨迹来源，同时减少相关训练体系造成的乐观评价。主分析删除与固定外部队列关联的633条记录，使用保留2,430条记录上训练的全链家族第0折完整与链内模型，并在同集合拟合静态与类型对照。关联组规模与保留集合见方法及补充材料H；图S1保留旧关联规则，仅作历史敏感性。外部链在声明的序列搜索和精确重复条件下不命中保留训练集合。表2和图4报告新的严格单模型方案。')
doc.save(OUTPUT)
with zipfile.ZipFile(SOURCE) as old, zipfile.ZipFile(OUTPUT) as new:
    original_comments = old.read('word/comments.xml')
    saved_comments = new.read('word/comments.xml')
    assert etree.tostring(etree.fromstring(original_comments), method='c14n') == etree.tostring(etree.fromstring(saved_comments), method='c14n')
    saved_items = [(info, original_comments if info.filename == 'word/comments.xml' else new.read(info.filename)) for info in new.infolist()]
temp = OUTPUT.with_suffix('.preserved.tmp')
with zipfile.ZipFile(temp, 'w') as package:
    for info, content in saved_items:
        package.writestr(info, content)
temp.replace(OUTPUT)
with zipfile.ZipFile(SOURCE) as old, zipfile.ZipFile(OUTPUT) as new:
    assert old.read('word/comments.xml') == new.read('word/comments.xml')
    assert all(old.read(n) == new.read(n) for n in old.namelist() if n.startswith('word/media/') and not n.endswith('/'))
assert len(doc.tables) == 14
audit = {'document':str(OUTPUT), 'sha256':sha(OUTPUT), 'tables':14, 'source_sha256':sources,
         'comments_preserved':True, 'old_images_preserved':True, 'external_figure_sha256':sha(figure),
         'new_primary_external_training_records':2430, 'new_primary_external_purged':633,
         'server_training_fits_completed':20, 'md_repeats_per_state':1,
         'required_remaining':'Independent manuscript content audit and full Word render/visual QA'}
(RESULT / 'integration_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf8')
print(json.dumps({'document':str(OUTPUT), 'sha256':sha(OUTPUT), 'tables':14}, ensure_ascii=False))
