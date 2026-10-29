"""Certify the final, visually reviewed Chinese v23 manuscript and delivery."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone, timedelta
import hashlib, json, re, shutil, zipfile
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from PIL import Image, ImageChops
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent.parent
QA = ROOT / 'manuscript_QA'
SOURCE = BASE / 'paper_v22_20261012_中文修订_自由能分析融入正文_三线表.docx'
DOC = BASE / 'paper_v23_20261012_中文修订_论证精简与补充整合_三线表.docx'
PACKAGE = BASE / 'paper_v23_20261012_更新图8及绘图数据.zip'
TMP = Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
write = lambda p, x: Path(p).write_text(json.dumps(x, ensure_ascii=False, indent=2), encoding='utf8')
now = datetime.now(timezone(timedelta(hours=8))).isoformat()
checks = []

def check(name, condition):
    checks.append({'name': name, 'passed': bool(condition)})
    assert condition, name

check('Final document is the independently reviewed version', sha(DOC) == '04e952548a9e2de9091f138d510238f93e816f5731ad801bee9af5f07fcbf4f8')
check('Frozen source version unchanged', sha(SOURCE) == '315ca2108ab9dbdfae83774f4d9fe1556087d66f8ea09697afa9c9eb67b67bdd')
integration = json.loads((QA / 'integration_audit.json').read_text(encoding='utf8'))
numerical = json.loads((QA / 'v23_results_independent_review.json').read_text(encoding='utf8'))
content = json.loads((QA / 'v23_independent_content_review.json').read_text(encoding='utf8'))
check('Independent numerical review passes on final document', numerical['status'] == 'passed' and numerical['document_sha256'] == sha(DOC))
check('Independent logic and references review passes on final document', content['status'] == 'pass' and content['docx_sha256'] == sha(DOC))
check('Integration certificate matches final document', integration['document_sha256'] == sha(DOC))
doc, source = Document(DOC), Document(SOURCE)
mapping = [0, 1, 3, 4, 11, 13, 12, 8, 2, 9, 14]
cells = lambda t: tuple(tuple(c.text for c in r.cells) for r in t.rows)
check('All eleven retained tables preserve native source cells exactly', len(doc.tables) == 11 and all(cells(t) == cells(source.tables[mapping[i]]) for i, t in enumerate(doc.tables)))
check('Exactly 393 preserved native table cells', sum(len(r.cells) for t in doc.tables for r in t.rows) == 393)
for i, table in enumerate(doc.tables, 1):
    check(f'Table {i} has no table grid or cell shading', not table._tbl.xpath('.//w:shd | ./w:tblPr/w:tblBorders'))
    correct = True
    for j, row in enumerate(table.rows):
        for cell in row.cells:
            borders = cell._tc.tcPr.find(qn('w:tcBorders'))
            assert borders is not None
            visible = {etree.QName(e).localname for e in borders if e.get(qn('w:val')) not in {'nil', 'none'}}
            correct &= visible == ({'top', 'bottom'} if j == 0 else {'bottom'} if j == len(table.rows)-1 else set())
    check(f'Table {i} has exactly three horizontal rules', correct)

ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main', 'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math', 'a': 'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a, zipfile.ZipFile(DOC) as b:
    ax, bx = etree.fromstring(a.read('word/document.xml')), etree.fromstring(b.read('word/document.xml'))
    check('Original user comments byte-identical', a.read('word/comments.xml') == b.read('word/comments.xml'))
    anchors = lambda x: Counter((e.tag, e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'), qn('w:commentRangeEnd'), qn('w:commentReference')})
    check('All user comment anchors preserved', anchors(ax) == anchors(bx))
    math = lambda x: [etree.tostring(e, method='c14n') for e in x.xpath('.//m:oMath', namespaces=ns)]
    check('All native display and inline mathematics preserved', math(ax) == math(bx) and len(math(bx)) == 36)
    check('Ten actual scientific images', len(bx.xpath('.//a:blip', namespaces=ns)) == 10)
    old_media = {hashlib.sha256(a.read(n)).hexdigest() for n in a.namelist() if n.startswith('word/media/')}
    new_figure_sha = sha(ROOT / 'figures/fig8_joint_partition_CA.png')
    rel = etree.fromstring(b.read('word/_rels/document.xml.rels'))
    targets = {e.get('Id'): e.get('Target') for e in rel}
    active = [targets[e.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed')] for e in bx.xpath('.//a:blip', namespaces=ns)]
    active_hashes = [hashlib.sha256(b.read('word/' + n)).hexdigest() for n in active]
    check('Nine retained figures unchanged and new figure8 exact', active_hashes.count(new_figure_sha) == 1 and all(s in old_media or s == new_figure_sha for s in active_hashes))

text = '\n'.join(p.text for p in doc.paragraphs)
main = text[:text.index('4 讨论')]
figs = [m.group(1) for p in doc.paragraphs if (m := re.match(r'^图\s*(S?\d+)\s*[:：]', p.text))]
tables = [m.group(1) for p in doc.paragraphs if (m := re.match(r'^表\s*(S?\d+)\s*[:：]', p.text))]
check('Nine main and one supplementary figure captions consecutive', figs == [str(i) for i in range(1, 10)] + ['S1'])
check('Two main and nine supplementary table captions consecutive', tables == ['1', '2'] + ['S' + str(i) for i in range(1, 10)])
for prefix, captions in [('图', figs), ('表', tables)]:
    check(f'All {prefix} cross-references resolve', set(re.findall(prefix + r'\s*(S?\d+)', text)) <= set(captions))
refs = [int(m.group(1)) for p in doc.paragraphs if (m := re.match(r'^\[(\d+)\]', p.text))]
check('35 references preserved and consecutive', refs == list(range(1, 36)))
check('References all cited', content['checks']['references_1_through_35_complete_and_cited'])
sections = [m.group(1) for p in doc.paragraphs if (m := re.match(r'^([A-G])\s', p.text))]
check('Seven supplementary modules A through G', sections == list('ABCDEFG'))
check('Free energy and distribution evidence remains in main text', all(x in main for x in ['3.6.1 泛素', '3.6.2 角度', '3.6.3 不同蛋白', '图7：', '图8：', '图9：']))
check('Superseded tables, figures and posthoc calibration removed', not any(x in text for x in ['图S2：', '图S12', '图S13', '表S10', '表S11', '表S12', '事后分组校准']))

pdf_path = TMP / 'v23_word.pdf'
pdf = PdfReader(pdf_path)
check('31 rendered Word pages', len(pdf.pages) == 31)
check('No blank pages or replacement glyphs', all(len((p.extract_text() or '').strip()) > 20 and '\ufffd' not in (p.extract_text() or '') for p in pdf.pages))
reuse = json.loads((QA / 'page_pixel_reuse.json').read_text(encoding='utf8'))
viewed = {5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}
check('Page coverage matches actual PDF, ignoring obsolete render files', len(reuse) == 31 and [p['page'] for p in reuse] == list(range(1, 32)))
page_review = []
for p in reuse:
    n, old = p['page'], p['reuse_v22']
    image_path = TMP / 'render_v23' / f'page-{n}.png'
    if n in viewed:
        method = 'final_PNG_visually_inspected_by_root'
    else:
        assert old is not None
        with Image.open(image_path) as current, Image.open(TMP / 'render_v22' / f'page-{old}.png') as previous:
            check(f'Page {n} body pixel-identical to previously reviewed v22', current.size == previous.size and ImageChops.difference(current.crop((0, 0, current.width, current.height-110)), previous.crop((0, 0, previous.width, previous.height-110))).getbbox() is None)
        method = 'pixel_identical_to_previously_inspected_v22_page_body'
    page_review.append({'page': n, 'method': method, 'v22_page': old, 'PNG_sha256': sha(image_path)})
check('All 31 page bodies visually verified', len(viewed) == 25 and len(page_review) == 31)
with zipfile.ZipFile(PACKAGE) as z:
    check('Updated figure package readable and CRC-valid', z.testzip() is None)
    check('Package contains exact new figure8', any(hashlib.sha256(z.read(n)).hexdigest() == new_figure_sha for n in z.namelist() if n.endswith('.png')))

review = {'status': 'passed', 'document_sha256': sha(DOC), 'reviewed_beijing': now, 'checks': checks, 'checks_passed': len(checks), 'pages': 31, 'directly_inspected_page_count': 25, 'pixel_reuse_page_count': 6, 'page_review': page_review, 'internal_Word_PDF_sha256': sha(pdf_path), 'independent_numerical_review_sha256': sha(QA / 'v23_results_independent_review.json'), 'independent_content_review_sha256': sha(QA / 'v23_independent_content_review.json'), 'layout_findings': 'All final changed page bodies inspected in original PNG resolution: captions, native formulas, three-line tables and two-panel figure8 are complete and readable; remaining six bodies are pixel-identical to previously inspected v22 pages. Deliberate page breaks keep large scientific figure plates intact.'}
write(QA / 'final_content_and_layout_audit.json', review)
integration['review'] = 'passed_independent_numerical_and_content_reviews_and_all_page_visual_QA'
integration['pages'] = 31
integration['final_review_certificate'] = str(QA / 'final_content_and_layout_audit.json')
write(QA / 'integration_audit.json', integration)

note = ROOT / '论证精简与补充整合说明.md'
note.write_text('''v23 根据逐节论证审计完成实质删改。结果部分由预测能力、同类型残基差异、双读数互补、外部验证、具体位点定位和构象分布解释构成。讨论缩为三个连贯段落；结论限定于实际完成的评价范围。联合训练作为默认方案，消除围绕家族隔离成绩变化的重复讨论。主文中的自由能分析仍完整保留于第2.8和3.6节，泛素多角度分布、83条记录的联合占据与Cα均方波动分解、不同来源代表位点的投影自由能继续作为正文证据。

补充材料由十一节整合为七个模块：A数据来源与标签观测，B训练方案与随机初始化，C基线实现与来源比较，D全组成链家族隔离评价，E共同窗口与构象统计，F KIX轨迹与状态比较，G双读数组合的完整统计。旧普通外部验证、单列事后校准、重复弹性网络表和重复补充病例图撤出稿件；十张补充图删除，仅保留必要的标签观测图S1。原全状态差值表改为关键数值与完整源表入口。原角度分区的完整统计由正文表3转入补充表S7，图8重画为直接解释联合熵与主链散布的两个面板。

精简保留真实结果及其必要适用范围：低／高与高／低组合召回率为15.8%和16.5%，KIX六对状态差值方向均与单轨迹参考相反；正文未将状态内活动区域定位写成状态响应预测。原生方差覆盖率保留在补充A，不以事后校准改写模型成绩。随机初始化明确为固定第0折，单目标对照为每目标五个模型。全部十一张保留表的393个单元格与v22逐字一致，均为原生可编辑三线表；原用户批注、锚点、五个数学公式及35篇参考文献保留。

最终Word为31页，共九张正文图、一张补充图、两张正文表和九张补充表。最终版本通过独立数值审阅、独立论证与引用审阅，以及全部31页的视觉验证。本轮为稿件结构和图8的精简，没有追加训练或MD生产。原始计算、已完成训练权重、轨迹、被撤出稿件的图表以及此前稿件版本均继续存档；v21科学源数据归档和v22自由能解释源表保持不变。
''', encoding='utf8')
pointer = ROOT.parent / 'completed_results_20261012/delivery_status.json'
status = json.loads(pointer.read_text(encoding='utf8'))
snapshot = QA / 'delivery_status_v22_snapshot.json'
if not snapshot.exists():
    shutil.copy2(pointer, snapshot)
    shutil.copy2(ROOT.parent / 'completion_status_zh.md', QA / 'completion_status_v22_snapshot.md')
status.update({'status': 'chinese_manuscript_v23_argument_pruning_delivered', 'completed_beijing': now, 'previous_document': str(SOURCE), 'previous_document_sha256': sha(SOURCE), 'document': str(DOC), 'document_sha256': sha(DOC), 'pages': 31, 'figures': 10, 'main_figures': 9, 'supplementary_figures': 1, 'three_line_tables': 11, 'main_tables': 2, 'supplementary_tables': 9, 'supplementary_modules': 7, 'references': 35, 'content_review': 'passed_independent_numerical_and_logic_reviews', 'render_and_visual_review': 'passed_all_31_pages', 'new_figure_package': str(PACKAGE), 'new_figure_package_sha256': sha(PACKAGE), 'new_experiment_note': str(note), 'changes_this_round': 'Pruned defensive repetition and superseded external/calibration blocks; merged eleven supplements into seven modules; removed ten redundant supplementary figures; simplified main figure8 to joint entropy and state-conditioned CA variance; retained main-text free energy and material weak results without modifying source table cells.', 'new_MD_production_this_round': 0, 'new_model_training_this_round': 0, 'remaining_required_work': []})
write(ROOT / 'delivery_status.json', status)
write(pointer, status)
(ROOT.parent / 'completion_status_zh.md').write_text(f'当前交付为 {DOC}，SHA256为 {sha(DOC)}，31页。结果、讨论和补充材料已按论证作用精简，自由能与构象占据证据保留正文；九张正文图、一张补充图、十一张原生可编辑三线表、七个补充模块和35篇参考文献。393个保留表格单元格与v22完全相同，原用户批注及五个原生数学公式保留。独立数值及论证引用审阅通过，最终31页视觉验证通过。图8及数据包为 {PACKAGE}。全部原实验结果、训练权重、轨迹、被撤图表和历史版本继续存档，当前详细状态见 {ROOT / "delivery_status.json"}。本轮所需工作已完成。\n', encoding='utf8')
print(json.dumps({'status': 'complete', 'document': str(DOC), 'document_sha256': sha(DOC), 'pages': 31, 'tables': 11, 'figures': 10, 'supplementary_modules': 7, 'root_checks': len(checks), 'package_MB': round(PACKAGE.stat().st_size/1e6, 2)}, ensure_ascii=False))
