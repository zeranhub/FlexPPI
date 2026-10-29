"""Align the main comparison with Fig. 2 and move revised Table 1 to S11."""
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import zipfile
from lxml import etree as E

ROOT = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / 'paper_v70_20261024_中文修订_图2紧凑终版.docx'
TARGET = ROOT / 'paper_v71_20261024_中文修订_表1对齐图2并移补充.docx'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
W = '{' + NS['w'] + '}'


def text(element):
    return ''.join(element.xpath('.//w:t/text()', namespaces=NS))


def set_text(element, value):
    """Keep paragraph properties while replacing its complete textual content."""
    old_run = element.find('w:r', NS)
    old_format = None if old_run is None else old_run.find('w:rPr', NS)
    for child in list(element):
        if child.tag != W + 'pPr':
            element.remove(child)
    run = E.SubElement(element, W + 'r')
    if old_format is not None:
        run.append(deepcopy(old_format))
    node = E.SubElement(run, W + 't')
    node.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    node.text = value


def replace_in_paragraph(element, old, new):
    value = text(element)
    assert value.count(old) == 1, (old, value)
    set_text(element, value.replace(old, new))


def add_para(template, value):
    result = deepcopy(template)
    for attr in list(result.attrib):
        if attr.endswith('paraId') or attr.endswith('textId'):
            del result.attrib[attr]
    set_text(result, value)
    return result


def main():
    assert not TARGET.exists(), 'Keep prior document versions.'
    external = json.loads((OUT / 'external_prediction_source.json').read_text(encoding='utf-8'))
    with zipfile.ZipFile(SOURCE) as package:
        entries = package.infolist()
        parts = {item.filename: package.read(item.filename) for item in entries}
    root = E.fromstring(parts['word/document.xml'])
    body = root.find('w:body', NS)
    originals = list(body)
    before = {id(el): E.tostring(el) for el in originals}

    def para(prefix):
        found = [el for el in originals if el.tag == W + 'p' and text(el).startswith(prefix)]
        assert len(found) == 1, (prefix, len(found))
        return found[0]

    old_caption = para('表1：')
    old_table = body[list(body).index(old_caption) + 1]
    assert old_table.tag == W + 'tbl'
    old_note = para('注：Full为完整结构模型')
    old_rows = old_table.findall('w:tr', NS)
    # Preserve the full original methods paragraph, including its OMML formula.
    legacy_backflip = deepcopy(para('BackFlip-2采用2.1版本'))
    first = legacy_backflip.xpath('.//w:t', namespaces=NS)[0]
    first.text = '在全队列扩展比较中，' + first.text.replace('BackFlip-2采用', '逐链BackFlip-2采用', 1)
    legacy_other = deepcopy(para('PLANET-MD[11]比较'))
    legacy_performance = deepcopy(para('在已有RMSF方法中'))
    for copied in [legacy_backflip, legacy_other, legacy_performance]:
        for attr in list(copied.attrib):
            if attr.endswith('paraId') or attr.endswith('textId'):
                del copied.attrib[attr]
    legacy_source = add_para(para('在ATLAS、DynaRepo和ITrace三个来源中'),
        '在全队列的来源分层扩展比较中，FlexPPI相对逐链BackFlip-2.1的ATLAS RMSF误差较高，'
        'DynaRepo和ITrace幅值误差较低，DynaRepo中的残基排序相关接近。'
        '逐链BackFlip-2.1在全队列的RMSF MAE和Spearman相关分别为0.537 Å和0.848；'
        '在ThermoPCD的32条记录中分别为0.749 Å和0.593。上述逐链结果与表S11B的assembly结果分别报告。')

    set_text(para('BackFlip-2采用2.1版本'),
        '外部RMSF比较采用BackFlip-2.1、BackFlip-1.0、DynaProt-M、FastProtFlex、ANM和GNM，'
        '与完整FlexPPI共同在615条DynaRepo复合物的相同有效残基上评价。BackFlip两个版本均采用assembly模式，'
        '一次输入全部链的主链坐标；DynaProt-M与FastProtFlex使用作者发布的权重及默认推理设置。'
        'DynaProt-M按作者实现逐链推理后合并。'
        'ANM与GNM由完整输入的Cα构建弹性网络。')
    set_text(para('PLANET-MD[11]比较'),
        '具有绝对单位输出的预测统一换算为Å后计算幅值误差。FastProtFlex、ANM和GNM输出相对柔性分数，'
        '仅比较排序、筛选及记录内标准化误差，不报告Å单位MAE。已有方法直接使用作者提供的模型；'
        '静态特征和类型均值对照使用本文对应训练集合拟合。弹性网络参数及全队列额外基线的实现与结果见补充材料C。')
    metric_para = para('评价先在每条记录的有效残基上计算')
    set_text(metric_para, text(metric_para) +
        '标准化MAE将每条记录的参考与预测分别减去自身均值、除以自身标准差，再计算平均绝对差，为无量纲。')
    replace_in_paragraph(para('在3,063条记录的五折体系留出评价中'), '（表1）', '（表S11A）')
    replace_in_paragraph(para('完整结构编码相对仅链内编码'), '（表1）', '（表S11A）')
    set_text(para('在已有RMSF方法中'),
        '外部模型比较采用图2d、e的七种方法及相同615条DynaRepo复合物。'
        '数值预测结果见表S11B，PR与ROC筛选指标见表S10；全队列额外方法的扩展比较见补充材料C。')
    set_text(para('在ATLAS、DynaRepo和ITrace三个来源中'),
        '在ATLAS、DynaRepo和ITrace三个来源中，FlexPPI的两项平均误差均低于静态特征模型。')
    interface = para('为检验结构编码在多链环境中的作用')
    set_text(interface, text(interface).replace('BackFlip-2', 'BackFlip-2.1（assembly）'))
    replace_in_paragraph(para('该比较只使用DynaRepo复合物记录'),
        '排除后剩余616条记录，其中615条具有至少一个可定义的界面。',
        '因此比较限定于DynaRepo的616条多链记录，其中615条具有至少一个可定义的界面。')
    set_text(para('各对照方法按其官方实现处理同一份完整输入'),
        '各对照方法按其官方实现处理同一份完整输入。BackFlip-2.1与BackFlip-1.0均采用assembly模式，'
        '一次输入全部链的主链坐标；逐链模式不使用链间坐标，未用于图2d、e。'
        'FastProtFlex、ANM和GNM输出相对柔性分数，其绝对RMSF MAE不适用，'
        '但可参与排序、筛选及记录内标准化误差比较。')
    sidechain = para('现有柔性预测方法均不输出侧链构象熵')
    suffix = text(sidechain).split('这使联合预测成为', 1)[1]
    set_text(sidechain,
        '本文主比较中的BackFlip-2.1、BackFlip-1.0、DynaProt-M和FastProtFlex均不提供Hχ输出；'
        'ANM和GNM给出主链相对柔性分数，因此侧链任务采用Full、Intra、Static和AA进行内部对照。'
        '这使联合预测成为' + suffix)
    replace_in_paragraph(para('在ThermoPCD的32条外部记录中'), '（表1B）', '（表S11C）')

    # Preserve older expanded comparisons in C, clearly separated from Fig. 2.
    c_heading = para('C 基线实现与来源比较')
    c_intro = add_para(para('评价先在每条记录的有效残基上计算'),
        '本节保留3,063条记录、961,321个RMSF有效位置上的额外基线比较，'
        '包括逐链BackFlip-2.1、PLANET-MD和RMSF-net。该全队列扩展评价与图2d、e及表S11B中'
        '615条DynaRepo复合物的外部方法比较分别报告。')
    position = list(body).index(c_heading) + 1
    for el in [c_intro, legacy_backflip, legacy_other, legacy_performance, legacy_source]:
        body.insert(position, el)
        position += 1
    replace_in_paragraph(para('表S4：'),
        '表S4：共同RMSF位置上的RMSF预测方法与相对柔性基线。',
        '表S4：全队列共同RMSF位置上的额外方法与相对柔性基线扩展比较。')

    # Derive A/C directly from the original table; only B uses the frozen 615 cohort.
    def row_cells(row):
        return [text(cell) for cell in row.findall('w:tc', NS)]
    original_data = {'A': {}, 'C': {}}
    group = None
    for row in old_rows[2:]:
        cells = row_cells(row)
        if len(cells) == 1:
            group = 'A' if cells[0].startswith('A') else 'C'
        elif cells[0] in ['氨基酸类型均值', '静态特征回归', 'FlexPPI (Intra)', 'FlexPPI (Full)']:
            original_data[group][cells[0]] = cells[1:]
    control_order = [('FlexPPI (Full)', 'FlexPPI (Full)'), ('FlexPPI (Intra)', 'FlexPPI (Intra)'),
                     ('静态特征回归', 'Static'), ('氨基酸类型均值', 'AA')]
    body_rows = []
    expected_rows = []
    group_template = old_rows[2]
    data_template = old_rows[3]
    for key, heading in [('A', 'A  五折体系留出：内部对照（3,063条记录）'),
                         ('B', 'B  DynaRepo复合物：外部方法比较（615条记录）'),
                         ('C', 'C  ThermoPCD：内部对照（32条记录）')]:
        group_row = deepcopy(group_template)
        set_text(group_row.find('w:tc/w:p', NS), heading)
        body_rows.append(group_row)
        expected_rows.append([heading])
        if key != 'B':
            records = [[label, *original_data[key][old]] for old, label in control_order]
        else:
            names = ['Full', 'backflip-2.1 (assembly)', 'backflip-1.0 (assembly)',
                     'DynaProt-M', 'FastProtFlex', 'ANM', 'GNM']
            by_method = {item['method']: item for item in external['methods']}
            assert set(by_method) == set(names)
            records = []
            for method in names:
                item = by_method[method]
                absolute = item['absolute_angstrom']
                assert (item['mae'] is not None) == absolute
                label = 'FlexPPI (Full)' if method == 'Full' else method.replace('backflip', 'BackFlip')
                records.append([label, f"{item['mae']:.3f}" if absolute else '—',
                                f"{item['spearman']:.3f}", '—', '—'])
        for values in records:
            row = deepcopy(data_template)
            for cell, value in zip(row.findall('w:tc', NS), values):
                set_text(cell.find('w:p', NS), value)
                for rp in cell.xpath('.//w:rPr', namespaces=NS):
                    for bold in list(rp.findall('w:b', NS)):
                        rp.remove(bold)
                    if values[0] == 'FlexPPI (Full)':
                        E.SubElement(rp, W + 'b')
            body_rows.append(row)
            expected_rows.append(values)
    new_table = deepcopy(old_table)
    for row in new_table.findall('w:tr', NS)[2:]:
        new_table.remove(row)
    for row in body_rows:
        new_table.append(row)
    # Keep a conventional three-rule table with a final bottom border.
    for cell in body_rows[-1].findall('w:tc', NS):
        borders = cell.find('w:tcPr/w:tcBorders', NS)
        edge = borders.find('w:bottom', NS)
        edge.set(W + 'val', 'single')
        edge.set(W + 'sz', '12')
        edge.set(W + 'color', '000000')
    # Keep the table within one supplement page; no artificial fixed row heights.
    for row in new_table.findall('w:tr', NS):
        for para_el in row.xpath('./w:tc/w:p', namespaces=NS):
            props = para_el.find('w:pPr', NS)
            if props is not None:
                for node in props.findall('w:keepNext', NS):
                    props.remove(node)

    new_heading = add_para(para('H 外部模型的高RMSF位点筛选'), 'I 数值预测与对照结果')
    new_caption = add_para(old_caption,
        '表S11：与图2模型口径一致的数值预测结果。A为3,063条记录的五折体系留出内部比较，'
        'B为图2d、e相同615条DynaRepo复合物上的七种方法，C为ThermoPCD的32条外部记录。'
        '各项先逐记录计算，再对记录等权平均；MAE越低越好，Spearman相关ρ越高越好。')
    new_note = add_para(old_note,
        '注：Full为完整结构模型，Intra为链内模型，Static为静态特征回归，AA为训练集氨基酸类型均值对照。'
        'B中BackFlip两个版本均采用assembly模式，DynaProt-M按作者实现逐链推理后合并；'
        'FastProtFlex、ANM和GNM输出相对分数，Å单位MAE不适用。'
        'B只比较RMSF，Hχ未在该子集中汇总；“—”表示不适用或未评价。其AP与AUROC见表S10。'
        'C中四种方法均使用排除关联记录后的2,430条训练记录。')
    for el in [old_caption, old_table, old_note]:
        body.remove(el)
    # Let the remaining 3.1 results fill the space left by the removed table.
    # The figure and its caption remain together after the PR/ROC results.
    figure_nodes = body.xpath('./w:p[.//*[local-name()="docPr" and @id="2"]]', namespaces=NS)
    assert len(figure_nodes) == 1
    figure_group = [figure_nodes[0], para('图2：')]
    for el in figure_group:
        body.remove(el)
    after_results = list(body).index(para('在相同615条DynaRepo复合物上')) + 1
    for el in figure_group:
        body.insert(after_results, el)
        after_results += 1
    last = body.find('w:sectPr', NS)
    insert_at = len(body) if last is None else list(body).index(last)
    for el in [new_heading, new_caption, new_table, new_note]:
        body.insert(insert_at, el)
        insert_at += 1

    actual_rows = [row_cells(row) for row in new_table.findall('w:tr', NS)[2:]]
    assert actual_rows == expected_rows
    assert len(new_table.findall('w:tr', NS)) == 20
    all_text = ''.join(root.xpath('//w:t/text()', namespaces=NS))
    assert '表1' not in all_text.replace(' ', '')
    assert all(f'表S11{key}' in all_text for key in 'ABC')
    assert E.tostring(para('表S10：')) == before[id(para('表S10：'))]
    changed_original = [i for i, el in enumerate(originals)
                        if el.getparent() is body and E.tostring(el) != before[id(el)]]
    xml = E.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
    with zipfile.ZipFile(TARGET, 'w', zipfile.ZIP_DEFLATED) as package:
        for item in entries:
            package.writestr(item, xml if item.filename == 'word/document.xml' else parts[item.filename])
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        changed_parts = [name for name, data in parts.items() if package.read(name) != data]
        assert changed_parts == ['word/document.xml']
    report = {'source': str(SOURCE), 'target': str(TARGET),
              'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'target_sha256': hashlib.sha256(TARGET.read_bytes()).hexdigest(),
              'changed_package_parts': changed_parts, 'all_other_parts_byte_identical': True,
              'figure2_and_all_other_media_byte_identical': True,
              'modified_original_block_indices': changed_original,
              'old_main_table_removed': True, 'new_supplement_table': 'S11 in Appendix I',
              'S11_rows': expected_rows, 'S10_preserved': True,
              'legacy_extended_results_preserved_in_appendix_C': True}
    (OUT / 'docx_patch_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'S11_rows'}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
