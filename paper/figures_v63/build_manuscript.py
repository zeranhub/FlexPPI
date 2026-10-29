"""Insert the selected Figure 2 and update its scientific descriptions in v61."""
from copy import deepcopy
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import json
import sys

from lxml import etree as ET
from PIL import Image

sys.stdout.reconfigure(encoding='utf-8')
ROOT = Path(__file__).resolve().parent
BASE = Path('E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v61_20261023_中文修订_图1输入输出标注修正.docx')
OUTPUT = ROOT.parent / 'paper_v63_20261023_中文修订_图2紧凑构图与文字同步.docx'
FIGURE = ROOT.parent / 'figure2_compact_layouts_20261023/fig2_compact_merged.png'
SUPPLEMENT = ROOT.parent / 'figure2_C_selected_20261023/fig2_source_AP_companion.png'
NS = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'wp': 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
    'pic': 'http://schemas.openxmlformats.org/drawingml/2006/picture',
    'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'rel': 'http://schemas.openxmlformats.org/package/2006/relationships',
    'w14': 'http://schemas.microsoft.com/office/word/2010/wordml',
}
sha = lambda data: hashlib.sha256(data).hexdigest()

CAPTION = (
    '图2：FlexPPI的数值预测与高动态位点识别。Full、Intra、Static和AA分别为完整模型、链内模型、静态特征回归和训练折氨基酸类型均值对照。'
    '（a）Full预测与MD参考的逐位置对应关系。蓝色RMSF使用下、左轴（Å），深赭色Hχ使用上、右轴（nats），分别按各自坐标读取。'
    '二维颜色包含3,063条记录全部961,321个RMSF及824,988个Hχ有效位置的网格概率质量，各记录等权、记录内位置等权；双色标尺共用对数质量标度。'
    '实线为参考分箱内的预测均值，同色带为1,000次来源分层重抽记录所得的95%区间，灰点线为理想一致关系。RMSF两轴在1 Å以内线性、以外对数，Hχ两轴为线性。MAE和Spearman相关ρ为逐记录结果的等权平均。'
    '（b）RMSF、总体Hχ及同类型Hχ三个任务的记录级AP分布。每组从左到右为Full、Intra、Static、AA，每种方法均保留3,063条记录；点为记录，小提琴表示分布，内粗条为第25–75百分位，白短划为中位数。'
    '（c）与b中三个任务依次对应的平均PR曲线。RMSF和总体Hχ以记录内MD第90百分位及以上为阳性，保留并列值；同类型Hχ在记录内各氨基酸类型中采用同一规则，42,108个有效记录–类型组先在记录内对类型等权，再对记录等权。'
    '灰横线为实际阳性率，依次为10.26%、10.29%和14.63%，末端保留全部入选状态。三个纵轴采用各自的连续分段线性刻度，均完整覆盖0–1，具体映射见补充材料C。'
    'AP按Full、Intra、Static、AA顺序列出，由各记录的原始完整阈值曲线计算后平均，不从显示曲线的视觉面积计算。'
)

DISPLAY_METHOD = (
    '图2b给出各任务全部记录的AP分布，四分位范围描述记录间差异。图2c完整保留0–1的Precision范围，采用连续分段线性纵轴以展开曲线集中的区间：'
    '以0–1对应的显示高度为100%，RMSF的0–0.5和0.5–1分别占25%和75%；总体Hχ的0–0.5、0.5–0.9和0.9–1分别占10%、35%和55%；'
    '同类型Hχ的0–0.45、0.45–0.75和0.75–1分别占12%、78%和10%。各轴独立标注真实数值，未省略任何数值区间；'
    '坐标变换仅用于显示，AP仍按原始阈值曲线逐记录计算，再等权平均。'
)

SOURCE_CAPTION = (
    '图S2：同一五折体系留出评价的来源分层平均AP。三组列分别为RMSF、总体Hχ及同类型Hχ，每组按ATLAS、DynaRepo、ITrace排列，n为记录数。'
    '各来源内先逐记录计分，再对记录等权平均；阳性定义及同类型汇总方式与图2相同。Full、Intra、Static、AA的含义同图2。热图统一使用0–1色标，深色表示较高AP。'
)


def text(p):
    return ''.join(p.xpath('.//w:t/text() | .//m:t/text()', namespaces=NS))


def replace_text(p, old, new):
    matches = [t for t in p.findall('.//w:t', NS) if old in (t.text or '')]
    assert len(matches) == 1, (old, len(matches))
    matches[0].text = matches[0].text.replace(old, new)


def set_text(p, value):
    """Keep original paragraph formatting while replacing simple text only."""
    assert not p.findall('.//m:oMath', NS)
    prototype = p.find('w:r', NS)
    rpr = deepcopy(prototype.find('w:rPr', NS)) if prototype is not None else None
    for child in list(p):
        if child.tag != '{%s}pPr' % NS['w']:
            p.remove(child)
    run = ET.SubElement(p, '{%s}r' % NS['w'])
    if rpr is not None:
        run.append(rpr)
    node = ET.SubElement(run, '{%s}t' % NS['w'])
    node.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    node.text = value


def clone_text(prototype, value):
    p = deepcopy(prototype)
    p.attrib.pop('{%s}paraId' % NS['w14'], None)
    set_text(p, value)
    return p


def set_picture(p, width, height, name, description, identifier=None, relationship=None):
    for extent in p.xpath('.//wp:extent | .//pic:spPr/a:xfrm/a:ext', namespaces=NS):
        extent.set('cx', str(width))
        extent.set('cy', str(height))
    docpr = p.find('.//wp:docPr', NS)
    docpr.set('name', name)
    docpr.set('descr', description)
    if identifier is not None:
        docpr.set('id', str(identifier))
        picpr = p.find('.//pic:cNvPr', NS)
        picpr.set('id', str(identifier))
        picpr.set('name', name)
    if relationship is not None:
        p.find('.//a:blip', NS).set('{%s}embed' % NS['r'], relationship)


def main():
    assert sha(BASE.read_bytes()) == '666ecc56bb3159e2146d8a2a22cec5db09069697c28a9d89c12d41b80750593c'
    assert sha(FIGURE.read_bytes()) == '2f449b964fa36f9a829b5c96b618e219229e75c18f59e0a9cf2f485aee54e6db'
    assert Image.open(FIGURE).size == (3242, 2657)
    assert Image.open(SUPPLEMENT).size == (3242, 797)
    with ZipFile(BASE) as source:
        original = {name: source.read(name) for name in source.namelist()}
        doc = ET.fromstring(original['word/document.xml'])
        rels = ET.fromstring(original['word/_rels/document.xml.rels'])
        body = doc.find('w:body', NS)
        ps = body.findall('w:p', NS)
        assert len(ps) == 178
        assert text(ps[75]).startswith('图2：')
        original_math = [ET.tostring(m) for m in doc.findall('.//m:oMath', NS)]
        assert len(original_math) == 36
        original_paragraphs = [ET.tostring(p) for p in ps]
        replace_text(ps[67], '预测与MD参考的逐位置对应关系见图2a、b。',
                     '图2a以各自坐标展示完整模型两类读数与MD参考的逐位置对应关系。')
        replace_text(ps[68], '（图2c、d）', '（图2b、c）')
        replace_text(ps[68], '完整模型在两个目标上的AP均高于三项对照，',
                     '完整模型在两项任务上的记录级AP中位数分别为0.755和0.939，高于链内模型的0.724和0.934及静态特征回归的0.669和0.919，但各方法的记录间分布仍有重叠（图2b）。完整模型在两个目标上的记录等权平均AP均高于三项对照，')
        set_text(ps[75], CAPTION)
        replace_text(ps[77], '（图2f）', '（图S2）')
        replace_text(ps[82], '（图2e）', '（图2b、c）')
        replace_text(ps[82], '该评价覆盖3,063条记录中的42,108个记录–类型组，',
                     '完整模型的记录级AP中位数为0.617，亦高于链内模型的0.604和静态特征回归的0.548（图2b）。该评价覆盖3,063条记录中的42,108个记录–类型组，')

        width = 6035040
        figure_height = round(width * 150 / 183)
        set_picture(ps[74], width, figure_height,
                    'Figure 2 prediction agreement and high-dynamics screening',
                    'a: Full RMSF and side-chain entropy predictions with separate paired axes; b: record-level AP distributions; c: task-specific precision-recall curves with independently labelled continuous piecewise-linear precision scales')

        # Add the unchanged source-stratified evidence to supplement C.
        anchor = ps[162]
        body.insert(body.index(anchor), clone_text(ps[159], DISPLAY_METHOD))
        picture = deepcopy(ps[74])
        picture.attrib.pop('{%s}paraId' % NS['w14'], None)
        max_id = max(int(n.get('id')) for n in doc.findall('.//wp:docPr', NS))
        new_rid = 'rId53'
        assert rels.find('rel:Relationship[@Id="%s"]' % new_rid, NS) is None
        set_picture(picture, width, round(width * 45 / 183),
                    'Supplementary Figure S2 source-stratified average precision',
                    'Mean AP for four methods across ATLAS, DynaRepo and ITrace in three tasks, with a common zero-to-one colour scale',
                    identifier=max_id + 1, relationship=new_rid)
        body.insert(body.index(anchor), picture)
        body.insert(body.index(anchor), clone_text(ps[151], SOURCE_CAPTION))
        ET.SubElement(rels, '{%s}Relationship' % NS['rel'],
                      Id=new_rid, Type=NS['r'] + '/image', Target='media/figureS2_source_AP.png')

        assert [ET.tostring(m) for m in doc.findall('.//m:oMath', NS)] == original_math
        edited_indices = [i for i, p in enumerate(ps) if ET.tostring(p) != original_paragraphs[i]]
        assert edited_indices == [67, 68, 74, 75, 77, 82], edited_indices
        assert len(body.findall('w:p', NS)) == 181
        assert len(doc.findall('.//wp:inline', NS)) == 7
        all_text = text(doc)
        for obsolete in ['图2a、b', '图2c、d', '图2e', '图2f']:
            assert obsolete not in all_text, obsolete
        assert all_text.count('图S2：') == 1
        assert all_text.count('（图S2）') == 1
        changed = {
            'word/document.xml': ET.tostring(doc, encoding='UTF-8', xml_declaration=True, standalone=True),
            'word/_rels/document.xml.rels': ET.tostring(rels, encoding='UTF-8', xml_declaration=True, standalone=True),
            'word/media/image2.png': FIGURE.read_bytes(),
        }
        added = {'word/media/figureS2_source_AP.png': SUPPLEMENT.read_bytes()}
        # PNG content type already exists; preserve all other package members.
        assert b'Extension="png"' in original['[Content_Types].xml']
        with ZipFile(OUTPUT, 'w') as destination:
            for info in source.infolist():
                destination.writestr(info, changed.get(info.filename, original[info.filename]))
            for name, data in added.items():
                destination.writestr(name, data, compress_type=ZIP_DEFLATED)
    with ZipFile(OUTPUT) as result:
        differences = [name for name, data in original.items() if result.read(name) != data]
        assert set(differences) == set(changed)
        assert result.read('word/media/image2.png') == FIGURE.read_bytes()
        assert result.read('word/media/figureS2_source_AP.png') == SUPPLEMENT.read_bytes()
    qa = ROOT / 'document_QA'
    qa.mkdir(exist_ok=True)
    audit = {
        'status': 'structural_checks_passed_render_pending',
        'baseline': str(BASE), 'baseline_sha256': sha(BASE.read_bytes()),
        'output': str(OUTPUT), 'document_sha256': sha(OUTPUT.read_bytes()),
        'modified_members': differences, 'added_members': list(added),
        'edited_original_paragraph_indices': edited_indices,
        'native_math_objects_preserved_byte_identically': 36,
        'all_other_existing_members_byte_identical': True,
        'figure2_width_mm': width / 36000, 'figure2_height_mm': figure_height / 36000,
        'figure2_sha256': sha(FIGURE.read_bytes()),
        'figureS2_sha256': sha(SUPPLEMENT.read_bytes()),
        'main_figure_panels': ['a', 'b', 'c'], 'supplementary_figure': 'S2',
        'new_figure_caption': CAPTION, 'new_display_methods': DISPLAY_METHOD,
        'new_source_caption': SOURCE_CAPTION,
    }
    (qa / 'document_edit_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    print(OUTPUT)
    print('Preserved all 36 native equations and all unrelated package members.')


if __name__ == '__main__':
    main()
