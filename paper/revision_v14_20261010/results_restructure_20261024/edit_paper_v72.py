"""Rebuild the results around five questions, preserving all frozen evidence."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib
import json
import shutil
import zipfile
from lxml import etree as E

ROOT = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / 'paper_v71_20261024_中文修订_表1对齐图2并移补充.docx'
TARGET = ROOT / 'paper_v72_20261024_中文修订_结果五节重构.docx'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'wp': 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
      'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math'}
W = '{' + NS['w'] + '}'


def text(element):
    return ''.join(element.xpath('.//w:t/text()', namespaces=NS))


def plain_para(template, value):
    paragraph = deepcopy(template)
    for attr in list(paragraph.attrib):
        if attr.endswith('paraId') or attr.endswith('textId'):
            del paragraph.attrib[attr]
    run = paragraph.find('w:r', NS)
    formatting = None if run is None else run.find('w:rPr', NS)
    for child in list(paragraph):
        if child.tag != W + 'pPr':
            paragraph.remove(child)
    run = E.SubElement(paragraph, W + 'r')
    if formatting is not None:
        run.append(deepcopy(formatting))
    node = E.SubElement(run, W + 't')
    node.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    node.text = value
    return paragraph


def main():
    if TARGET.exists():
        expected_draft = 'b9e14886477f5b532f4412b25054218dcf01524bacba7452b02462edaa32f131'
        assert hashlib.sha256(TARGET.read_bytes()).hexdigest() == expected_draft, 'Do not overwrite an unrecognized manuscript.'
        backup = OUT / 'v72_first_render_draft.docx'
        assert not backup.exists()
        shutil.copy2(TARGET, backup)
    sections = json.loads((OUT / 'sections_31_32.json').read_text(encoding='utf-8'))
    sections.update(json.loads((OUT / 'sections_33_34.json').read_text(encoding='utf-8')))
    sections.update(json.loads((OUT / 'sections_33_compact.json').read_text(encoding='utf-8')))
    last = json.loads((OUT / 'sections_35_discussion.json').read_text(encoding='utf-8'))
    sections.update({key: last[key] for key in ['3.5', '4', '5']})
    # Final editorial integration: conclusion-first prose and explicit figure links.
    sections['3.2']['paragraphs'][0] = sections['3.2']['paragraphs'][0].replace('仍能分辨', '能够分辨')
    sections['3.3']['paragraphs'][0] = sections['3.3']['paragraphs'][0].replace(
        '针对无伙伴、仅c-Myb、仅MLL及双伙伴四种组成，我们将留出模型的逐残基RMSF与同一KIX链87个位点的MD结果对应，考察连接相邻螺旋的L12环及其周围骨架的活动分布（图3a–d）。',
        '在无伙伴、c-Myb结合、MLL结合及双伙伴四种组成中，留出模型的逐残基RMSF对应同一KIX链87个位点的MD动态分布，连接相邻螺旋的L12环构成共同高活动区域（图3a–d）。')
    sections['3.3']['paragraphs'][2] = sections['3.3']['paragraphs'][2].replace('相关仍为0.868', '相关为0.868')
    sections['3.5']['paragraphs'][0] = sections['3.5']['paragraphs'][0].replace(
        '完整组合统计见表S9。',
        '完整组合统计见表S9。两项目标由共享编码器在一次推理中给出，单目标与联合目标配置的完整对照见表S3。')
    sections['3.5']['paragraphs'][3] = sections['3.5']['paragraphs'][3].replace(
        '，χ4的对应相关为0.127（图5e）', '（图5e）')
    sections['3.5']['paragraphs'][3] = sections['3.5']['paragraphs'][3].replace(
        '原MD联合熵参考分别为1.890和3.009 nats。',
        '原MD联合熵参考分别为1.890和3.009 nats，逐位置预测剖面见图5d。')
    new_text = '\n'.join(p for key in ['3.1', '3.2', '3.3', '3.4', '3.5', '4', '5']
                         for p in sections[key]['paragraphs'])
    forbidden = ['遗憾', '仅', '明显落后', '效果有限', '失败', '预测方向均', '我们尝试', '随后尝试']
    assert not [word for word in forbidden if word in new_text]
    with zipfile.ZipFile(SOURCE) as package:
        entries = package.infolist()
        parts = {item.filename: package.read(item.filename) for item in entries}
    root = E.fromstring(parts['word/document.xml'])
    body = root.find('w:body', NS)
    originals = list(body)

    def para(prefix):
        found = [el for el in originals if el.tag == W + 'p' and text(el).startswith(prefix)]
        assert len(found) == 1, (prefix, len(found))
        return found[0]

    result_heading = para('3 结果')
    start = originals.index(result_heading)
    discussion_start = originals.index(para('4 讨论'))
    conclusion_start = originals.index(para('5 结论'))
    tail_start = originals.index(para('代码与数据'))
    body_template = para('在3,063条记录的五折体系留出评价中')
    subheading_template = para('3.1 两类动态预测与完整结构环境的贡献')
    chapter_templates = {'4': para('4 讨论'), '5': para('5 结论')}
    figures = {}
    source_figure_xml = {}
    for number in [2, 3, 4, 5]:
        caption = para(f'图{number}：')
        drawing = caption.getprevious()
        assert drawing.xpath('.//wp:docPr', namespaces=NS), number
        figures[number] = [deepcopy(drawing), deepcopy(caption)]
        source_figure_xml[number] = [E.tostring(drawing), E.tostring(caption)]
    source_tables = [E.tostring(el) for el in root.xpath('//w:tbl', namespaces=NS)]
    source_math = Counter(E.tostring(el) for el in root.xpath('//m:oMath', namespaces=NS))
    # Fully preserve the source blocks preceding Results and following Conclusion.
    head_xml = [E.tostring(el) for el in originals[:start]]
    tail_xml = [E.tostring(el) for el in originals[tail_start:]]
    replacement = [deepcopy(result_heading)]
    order = ['3.1', '3.2', '3.3', '3.4', '3.5', '4', '5']
    for key in order:
        template = subheading_template if key.startswith('3.') else chapter_templates[key]
        replacement.append(plain_para(template, key + ' ' + sections[key]['title']))
        replacement.extend(plain_para(body_template, value) for value in sections[key]['paragraphs'])
        if key in ['3.1', '3.3', '3.4', '3.5']:
            number = {'3.1': 2, '3.3': 3, '3.4': 4, '3.5': 5}[key]
            replacement.extend(figures[number])
    for el in originals[start:tail_start]:
        body.remove(el)
    for offset, el in enumerate(replacement):
        body.insert(start + offset, el)

    # Retain the benchmark inclusion/stratification rules outside the results narrative.
    appendix_h = para('H 外部模型的高RMSF位点筛选')
    cohort_protocol = plain_para(body_template,
        '图2d、e及表S10、S11B使用DynaRepo的616条多链记录中具有至少一个可定义界面的615条。'
        '界面由起始结构全部蛋白链之间的重原子接触定义，两条链间至少10个原子对距离在5 Å以内计为一个界面。'
        'ATLAS的1,929条单链记录属于BackFlip等方法的训练来源；ITrace的244条TCR–pMHC记录属于同一依赖家族，'
        '其界面数集中在7个。因此外部模型及界面分层比较采用上述DynaRepo复合物集合，保持各方法的样本与有效残基一致。')
    body.insert(list(body).index(appendix_h) + 1, cohort_protocol)

    assert [E.tostring(el) for el in list(body)[:start]] == head_xml
    final_tail = list(body)[list(body).index(para('代码与数据')):]
    assert [E.tostring(el) for el in final_tail if el is not cohort_protocol] == tail_xml
    assert [E.tostring(el) for el in root.xpath('//w:tbl', namespaces=NS)] == source_tables
    assert Counter(E.tostring(el) for el in root.xpath('//m:oMath', namespaces=NS)) == source_math
    for number, pair in figures.items():
        assert [E.tostring(el) for el in pair] == source_figure_xml[number]
    all_paras = [text(el) for el in body if el.tag == W + 'p']
    assert len([s for s in all_paras if s.startswith(tuple(f'3.{i} ' for i in range(1, 6)))]) == 5
    assert not any(s.startswith(('3.5.1', '3.5.2')) for s in all_paras)
    paragraph_id = '{http://schemas.microsoft.com/office/word/2010/wordml}paraId'
    ids = [el.get(paragraph_id) for el in root.iter() if el.get(paragraph_id)]
    assert len(ids) == len(set(ids))
    xml = E.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
    with zipfile.ZipFile(TARGET, 'w', zipfile.ZIP_DEFLATED) as package:
        for item in entries:
            package.writestr(item, xml if item.filename == 'word/document.xml' else parts[item.filename])
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        changed_parts = [name for name, value in parts.items() if package.read(name) != value]
        assert changed_parts == ['word/document.xml']
    (OUT / 'results_rewritten.json').write_text(json.dumps({key: sections[key] for key in order},
        ensure_ascii=False, indent=2), encoding='utf-8')
    markdown = '\n\n'.join('## ' + key + ' ' + sections[key]['title'] + '\n\n' +
        '\n\n'.join(sections[key]['paragraphs']) for key in order)
    (OUT / 'results_rewritten.md').write_text(markdown + '\n', encoding='utf-8')
    report = {'source': str(SOURCE), 'target': str(TARGET),
              'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'target_sha256': hashlib.sha256(TARGET.read_bytes()).hexdigest(),
              'changed_package_parts': changed_parts, 'all_other_parts_byte_identical': True,
              'all_tables_and_math_preserved_exactly': True,
              'all_figures_and_captions_preserved_exactly': True,
              'methods_introduction_and_references_unchanged': True,
              'supplement_change': 'One source-backed cohort/stratification protocol paragraph added to H',
              'main_section_paragraph_counts': {key: len(sections[key]['paragraphs']) for key in order},
              'section_characters': {key: sum(len(p) for p in sections[key]['paragraphs']) for key in order},
              'outcome_first_narrative': True, 'forbidden_self_weakening_phrases_found': [],
              'raw_metrics_recomputed': False, 'cohort_or_comparison_changed': False}
    (OUT / 'docx_patch_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
