"""Integrate the verified historical single trajectories without new MD results."""
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import zipfile
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[2]
SOURCE = ROOT.parent / 'paper_v17_20261010_中文修订_基线与标签敏感性_三线表.docx'
OUTPUT = ROOT.parent / 'paper_v18_20261010_中文修订_复用KIX单轨迹_三线表.docx'
CASE = PROJECT / 'experiments/_case_validation_sources/data'
QUALITY = PROJECT / 'artifacts/kix_md_20260914/pilot_results_20260916/quality'
STATE_ORDER = ['apo', 'myb', 'mll', 'both']
STATE_NAME = {'apo': '游离KIX', 'myb': 'cMyb', 'mll': 'MLL', 'both': '双伙伴'}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

doc = Document(SOURCE)

def replace(p, text):
    # Keep comment anchors and the original paragraph properties.
    for child in list(p._p):
        if child.tag == qn('w:r') and not child.findall('.//' + qn('w:commentReference')):
            p._p.remove(child)
    p.add_run(text)

def replace_start(start, text):
    matches = [p for p in doc.paragraphs if p.text.startswith(start)]
    assert len(matches) == 1, (start, len(matches))
    replace(matches[0], text)

def make_table(caption, columns, rows, widths):
    p = doc.add_paragraph(caption, 'Caption')
    p.paragraph_format.keep_with_next = True
    t = doc.add_table(rows=1, cols=len(columns))
    t.autofit = False
    for col, width in zip(t.columns, widths):
        col.width = Inches(width)
    for cell, label in zip(t.rows[0].cells, columns):
        cell.text = label
    for values in rows:
        for cell, value in zip(t.add_row().cells, values):
            cell.text = str(value)
    for i, row in enumerate(t.rows):
        row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
        if i == 0:
            row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
        for cell, width in zip(row.cells, widths):
            cell.width = Inches(width)
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_before = Pt(4)
                p.paragraph_format.space_after = Pt(4)
                p.paragraph_format.first_line_indent = Pt(0)
                p.paragraph_format.left_indent = Pt(0)
                p.paragraph_format.right_indent = Pt(0)
                for run in p.runs:
                    run.font.size = Pt(10)
                    run.font.bold = i == 0
                    run.font.name = 'Times New Roman'
                    run._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'), '宋体')
    return t

case_metrics = pd.read_csv(CASE / 'kix_fold1_case_metrics.csv').set_index('state')
residues = pd.read_csv(CASE / 'kix_fold1_all_residues.csv')
assert set(residues.metric) == {'rmsf'}
assert set(residues.model_fold) == {1}
assert set(residues.md_independent_repeats) == {1}
assert all(len(residues[residues.state == s]) == 87 for s in STATE_ORDER)
l12 = residues[residues.author_residue_number.between(614, 621)]
assert all(len(l12[l12.state == s]) == 8 for s in STATE_ORDER)
l12_mean = l12.groupby('state')[['md_value', 'model_value']].mean()

abstract = next(p for p in doc.paragraphs if p.text.startswith('蛋白质的局部动态既包括'))
replace(abstract, abstract.text.replace('持续高波动的L12连接环', '四种状态中均相对高波动的L12连接环'))
method = next(p for p in doc.paragraphs if p.text.startswith('KIX模拟使用OpenMM'))
replace(method, method.text.replace('每种状态实际完成一条生产轨迹', '本研究分析每种状态各一条已完成的生产轨迹').replace(
    '单条轨迹也不能用于估计独立模拟重复之间的不确定性。',
    '单条轨迹也不能用于估计独立模拟重复之间的不确定性。历史轨迹的完整性、采样诊断及结合状态差值见补充材料F。'))
result = next(p for p in doc.paragraphs if p.text.startswith('在未结合伙伴、结合cMyb'))
replace(result, result.text.replace('L12在MD中持续位于较高柔性区域', '在四种状态的完整100 ns轨迹中，L12均位于较高柔性区域') +
    '活动区域的定位不等同于结合状态变化的准确预测：L12的双伙伴减游离状态RMSF均值在MD中为+1.017 Å，在同一留出模型中为−1.589 Å；MLL减cMyb分别为+0.014 Å和−1.263 Å。后者的MD差值在前后50 ns中还会改变符号。因此，本案例支持活动连接环的筛选，不能据此证明模型正确恢复了结合诱导的柔性变化或变构响应（补充材料F）。')
discussion = next(p for p in doc.paragraphs if p.text.startswith('监督标签还受轨迹长度'))
replace(discussion, discussion.text.replace(
    '增加独立重复及联合与独立单目标训练的比较仍需补充。',
    '这些标签仍依赖有限采样，联合学习相对于独立单目标模型的贡献仍需单独检验。'))
case_discussion = next(p for p in doc.paragraphs if p.text.startswith('KIX实例展示从完整结构定位活动连接环'))
replace(case_discussion, case_discussion.text.replace(
    '结构实例用于解释模型输出，不能替代系统性外部评价；',
    'KIX每种状态仅有一条100 ns轨迹，前后半程存在采样差异，且状态差值的预测出现方向冲突；结构实例用于解释模型输出，不能替代系统性外部评价；'))

heading = doc.add_paragraph('F KIX单轨迹核验与结合状态差值', 'Heading 1')
heading.paragraph_format.page_break_before = True
doc.add_paragraph('本节分析四种KIX组成各一条已完成的MD轨迹。每种状态含10 ns带约束平衡和100 ns无约束生产；保存的生产轨迹各有100个连续1 ns分段，每段1,000帧，帧间隔1 ps，合计400,000个蛋白构象。逐段轨迹、蛋白拓扑及逐残基指标均通过文件完整性核验。全帧质量分析重算的完整100 ns RMSF和Hχ与原分析逐位一致，已核查的共价键几何未见断链或周期边界引起的巨跳。所有状态均按KIX自身87个Cα对齐生产首帧，RMSF围绕各自对齐后轨迹均值计算。模型使用与正文图5一致的固定第1折留出权重，结构输入和指标有效位置保持一致。')

rows = []
for state in STATE_ORDER:
    row = case_metrics.loc[state]
    mean = l12_mean.loc[state]
    rows.append([STATE_NAME[state], f'{row.mae:.3f}', f'{row.spearman:.3f}',
                 f'{int(row.top_overlap)}/18', f'{mean.md_value:.3f}', f'{mean.model_value:.3f}'])
make_table('表S7：四种KIX组成的单轨迹RMSF评价。MAE在全部87个KIX残基上计算；ρ为Spearman相关，命中数为最高20%位点集合的重合数。L12采用614–621八个位点的算术平均，MAE与均值单位均为Å。每行只有一条独立MD轨迹，不报告重复间误差条。',
           ['状态', '整体MAE', '整体ρ', '高20%命中', 'L12参考均值', 'L12预测均值'],
           rows, [.80, .80, .70, 1.10, 1.40, 1.50])

trace = pd.read_csv(QUALITY / 'trajectory_trace_binned.csv')
assert set(trace.n_frames) == {100}
q = trace.groupby(['state', 'metric'])['mean'].mean()
rmsd = [q.loc[(s, 'KIX_CA_RMSD_angstrom')] for s in STATE_ORDER]
thermo = pd.read_csv(QUALITY / 'thermodynamics_all_samples.csv')
temperature = thermo.groupby('state').temperature_kelvin.mean()
windows = pd.read_csv(QUALITY / 'residue_window_metrics.csv')
half = windows[(windows.window_kind == 'disjoint_half') & windows.author_residue_number.between(614, 621)]
h = half.groupby(['state', 'start_ns_exclusive', 'end_ns_inclusive']).rmsf_angstrom.mean()
half_pairs = '、'.join(f'{h.loc[(s, 0., 50.)]:.3f}/{h.loc[(s, 50., 100.)]:.3f}' for s in STATE_ORDER)
contacts = [q.loc[('myb', 'partner_A_native_residue_contact_fraction')],
            q.loc[('mll', 'partner_C_native_residue_contact_fraction')],
            q.loc[('both', 'partner_A_native_residue_contact_fraction')],
            q.loc[('both', 'partner_C_native_residue_contact_fraction')]]
doc.add_paragraph('采样诊断使用完整100 ns，另比较前后两个不重叠的50 ns时间窗口。四状态按游离、cMyb、MLL、双伙伴顺序，KIX Cα相对生产首帧的时间平均RMSD分别为' +
    '、'.join(f'{v:.3f}' for v in rmsd) + ' Å；这些值反映相对起始构象的偏移，不等同于RMSF。各状态10,000个热力学采样点的平均温度为' +
    f'{temperature.min():.2f}–{temperature.max():.2f} K。L12在前50 ns/后50 ns不重叠时间窗口内计算的RMSF均值依次为' + half_pairs +
    ' Å，显示尤其cMyb状态存在明显的时段依赖。初始界面接触定义为制备结构中KIX与伙伴间距离不超过4.5 Å的重原子对；残基对只在至少一个原始原子对仍满足该阈值时计为保持。cMyb单伙伴、MLL单伙伴，以及双伙伴中的cMyb、MLL界面，原始残基接触的平均保持比例分别为' +
    '、'.join(f'{v*100:.1f}%' for v in contacts) +
    '。这一统计不计入重新形成的接触，不能由接触比例下降直接断言解离；但也不能将这些轨迹概括为界面始终保持初始状态或已经充分收敛。')

delta_paragraph = doc.add_paragraph('结合状态差值统一定义为表中前一状态的L12均值减后一状态的L12均值；模型和MD分别计算，不用状态内百分位替代原单位幅度。全部六个两两差值均保留（表S8）。例如双伙伴减游离的MD差值为+1.017 Å，模型为−1.589 Å；MLL减cMyb的完整轨迹差值为+0.014 Å，模型为−1.263 Å，且MD在前后50 ns分别为−0.402和+0.511 Å。高残基排序相关与响应差值的方向不一致可以同时存在，因而四状态结果支持活动区域的定位，尚不足以支持结合响应幅值或变构机制的预测。')
delta_paragraph.paragraph_format.page_break_before = True
pairs = [('myb', 'apo'), ('mll', 'apo'), ('both', 'apo'), ('mll', 'myb'), ('both', 'myb'), ('both', 'mll')]
rows = []
delta_records = []
for first, second in pairs:
    delta = l12_mean.loc[first] - l12_mean.loc[second]
    error = abs(delta.model_value - delta.md_value)
    rows.append([STATE_NAME[first] + '−' + STATE_NAME[second], f'{delta.md_value:+.3f}',
                 f'{delta.model_value:+.3f}', f'{error:.3f}'])
    delta_records.append({'first': first, 'second': second, 'md_delta': float(delta.md_value),
                          'model_delta': float(delta.model_value), 'absolute_error': float(error)})
make_table('表S8：L12八个位点平均RMSF的全部结合状态差值，单位为Å。参考差值来自每状态一条完整100 ns轨迹，预测差值来自同一个固定留出模型；绝对偏差为预测差值与参考差值之差的绝对值。',
           ['状态差值', 'MD参考差值', '模型预测差值', '绝对偏差'], rows, [1.70, 1.55, 1.65, 1.40])
doc.add_paragraph('本案例共有四种条件，各条件只有一次独立启动；四种条件不是同一条件的四次独立重复。前后半程和累积窗口是同一轨迹的时间子集，尤其累积窗口彼此重叠，不能据此构造独立重复的置信区间或显著性检验，也不作为收敛证明。MD参考是特定制备、力场和采样时长下的模拟统计量。结合状态变化的方向与幅值、界面重排以及慢运动采样仍是本案例的限制，结果不用于推断真实结合亲和力、转录功能或热力学变构效应。')

for t in doc.tables:
    prop = t._tbl.tblPr
    for e in list(prop):
        if e.tag in {qn('w:tblBorders'), qn('w:shd')}:
            prop.remove(e)
    for i, row in enumerate(t.rows):
        for cell in row.cells:
            cp = cell._tc.get_or_add_tcPr()
            for e in list(cp):
                if e.tag in {qn('w:tcBorders'), qn('w:shd')}:
                    cp.remove(e)
            borders = OxmlElement('w:tcBorders')
            for side in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
                outer = (i == 0 and side == 'top') or (i == len(t.rows) - 1 and side == 'bottom')
                head = i == 0 and side == 'bottom'
                e = OxmlElement('w:' + side)
                e.set(qn('w:val'), 'single' if outer or head else 'nil')
                if outer or head:
                    e.set(qn('w:sz'), str(12 if outer else 6))
                    e.set(qn('w:color'), '000000')
                    e.set(qn('w:space'), '0')
                borders.append(e)
            cp.append(borders)

doc.save(OUTPUT)
with zipfile.ZipFile(SOURCE) as old, zipfile.ZipFile(OUTPUT) as new:
    assert old.read('word/comments.xml') == new.read('word/comments.xml')
    assert all(old.read(n) == new.read(n) for n in old.namelist() if n.startswith('word/media/'))
assert len(doc.tables) == 10
audit = {'source': str(SOURCE), 'output': str(OUTPUT), 'output_sha256': sha(OUTPUT),
         'new_md_results_used': False, 'repeats_per_state': 1, 'production_ns_per_state': 100,
         'model_fold': 1, 'tables': 10, 'l12_deltas': delta_records,
         'source_sha256': {str(p): sha(p) for p in [CASE / 'kix_fold1_case_metrics.csv',
             CASE / 'kix_fold1_all_residues.csv', QUALITY / 'summary.json',
             QUALITY / 'trajectory_trace_binned.csv', QUALITY / 'residue_window_metrics.csv',
             QUALITY / 'thermodynamics_all_samples.csv']},
         'comments_preserved': True, 'original_images_preserved': True,
         'render_and_independent_source_audit_pending': True}
(ROOT / 'kix_single_integration_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf8')
print(json.dumps({k: audit[k] for k in ['output', 'output_sha256', 'tables', 'new_md_results_used']}, ensure_ascii=False))
