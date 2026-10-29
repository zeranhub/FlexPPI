"""Shorten Figure 2 caption and consolidate low-information supplement tables."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import csv
import hashlib
import io
import json
import re
import zipfile
from lxml import etree as E

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'paper_v78_20261024_中文修订_方法精简与顺序调整.docx'
TARGET = ROOT.parent / 'paper_v79_20261024_中文修订_图2图注与补充表精简.docx'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math'}
W = '{' + NS['w'] + '}'


def text(el):
    return ''.join(el.xpath('.//w:t/text()', namespaces=NS))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def paragraph(template, value):
    el = deepcopy(template)
    properties = el.find('w:r/w:rPr', NS)
    for child in list(el):
        if child.tag != W+'pPr':
            el.remove(child)
    run = E.SubElement(el, W+'r')
    if properties is not None:
        run.append(deepcopy(properties))
    node = E.SubElement(run, W+'t')
    node.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    node.text = value
    return el


def replace_range(el, start, end, new):
    offset = 0
    inserted = False
    for node in el.xpath('.//w:t', namespaces=NS):
        value = node.text or ''
        left, right = offset, offset + len(value)
        if right > start and left < end:
            prefix = value[:max(0, start-left)]
            suffix = value[max(0, end-left):] if right > end else ''
            node.text = prefix + (new if not inserted else '') + suffix
            inserted = True
        offset = right
    assert inserted


MAPPING = {2: '1A', 3: '1B', 4: '1C', 5: '2', 6: '3', 7: '4', 8: '5', 10: '6A', 11: '6B'}
REF = re.compile(r'(?P<prefix>表|[、及与])S(?P<n>1[01]|[1-9])(?P<panel>[AB]?)')


def renumber_refs(el):
    matches = list(REF.finditer(text(el)))
    if not matches:
        return el
    result = deepcopy(el)
    for match in reversed(matches):
        n = int(match['n'])
        assert n in MAPPING, text(el)
        new = match['prefix']+'S'+MAPPING[n]+match['panel']
        replace_range(result, match.start(), match.end(), new)
    return result


def rows(table):
    return [[text(cell) for cell in row.findall('w:tc', NS)] for row in table.findall('w:tr', NS)]


def keep_with_next(el, value=True):
    pp = el.find('w:pPr', NS)
    if pp is None:
        pp = E.SubElement(el, W+'pPr')
    k = pp.find('w:keepNext', NS)
    if k is None:
        k = E.SubElement(pp, W+'keepNext')
    k.set(W+'val', '1' if value else '0')


def merged_table(template, groups, widths):
    """Combine small tables as labelled panels, retaining every source cell's content."""
    table = deepcopy(template)
    for row in table.findall('w:tr', NS):
        table.remove(row)
    grid = table.find('w:tblGrid', NS)
    for col in list(grid):
        grid.remove(col)
    for width in widths:
        col = E.SubElement(grid, W+'gridCol')
        col.set(W+'w', str(width))
    tp = table.find('w:tblPr', NS)
    tw = tp.find('w:tblW', NS)
    if tw is None:
        tw = E.SubElement(tp, W+'tblW')
    tw.set(W+'w', str(sum(widths)))
    tw.set(W+'type', 'dxa')
    layout = tp.find('w:tblLayout', NS)
    if layout is None:
        layout = E.SubElement(tp, W+'tblLayout')
    layout.set(W+'type', 'fixed')
    header_cell = template.find('w:tr/w:tc', NS)
    for title, source, spans in groups:
        panel = E.Element(W+'tr')
        cell = deepcopy(header_cell)
        for child in list(cell):
            if child.tag != W+'tcPr':
                cell.remove(child)
        cp = cell.find('w:tcPr', NS)
        for tag in ['gridSpan', 'vMerge']:
            for node in cp.findall('w:'+tag, NS):
                cp.remove(node)
        span = E.SubElement(cp, W+'gridSpan'); span.set(W+'val', str(len(widths)))
        cw = cp.find('w:tcW', NS)
        if cw is None:
            cw = E.SubElement(cp, W+'tcW')
        cw.set(W+'w', str(sum(widths))); cw.set(W+'type', 'dxa')
        shade = cp.find('w:shd', NS)
        if shade is None:
            shade = E.SubElement(cp, W+'shd')
        shade.set(W+'fill', 'E7EDF4'); shade.set(W+'val', 'clear')
        cell.append(paragraph(header_cell.find('w:p', NS), title))
        panel.append(cell); table.append(panel)
        for original_row in source.findall('w:tr', NS):
            row = deepcopy(original_row)
            rp = row.find('w:trPr', NS)
            if rp is not None:
                for node in rp.findall('w:tblHeader', NS)+rp.findall('w:trHeight', NS):
                    rp.remove(node)
            offset = 0
            cells = row.findall('w:tc', NS)
            assert len(cells) == len(spans)
            for source_cell, count in zip(cells, spans):
                properties = source_cell.find('w:tcPr', NS)
                if properties is None:
                    properties = E.Element(W+'tcPr'); source_cell.insert(0, properties)
                for node in properties.findall('w:gridSpan', NS):
                    properties.remove(node)
                if count > 1:
                    span = E.SubElement(properties, W+'gridSpan'); span.set(W+'val', str(count))
                cw = properties.find('w:tcW', NS)
                if cw is None:
                    cw = E.SubElement(properties, W+'tcW')
                cw.set(W+'w', str(sum(widths[offset:offset+count]))); cw.set(W+'type', 'dxa')
                offset += count
            assert offset == len(widths)
            table.append(row)
    all_rows = table.findall('w:tr', NS)
    for i, row in enumerate(all_rows):
        rp = row.find('w:trPr', NS)
        if rp is None:
            rp = E.Element(W+'trPr'); row.insert(0, rp)
        if rp.find('w:cantSplit', NS) is None:
            E.SubElement(rp, W+'cantSplit')
        for p in row.xpath('.//w:p', namespaces=NS):
            keep_with_next(p, i < len(all_rows)-1)
    return table


def main():
    assert sha(SOURCE) == 'd0facaa397c0f0853b8fdf5c38fb102a89b756d8fcfda4edbfdb9b3431580743'
    assert not TARGET.exists()
    with zipfile.ZipFile(SOURCE) as package:
        entries = package.infolist()
        parts = {entry.filename: package.read(entry.filename) for entry in entries}
    root = E.fromstring(parts['word/document.xml'])
    body = root.find('w:body', NS)
    original = list(body)

    def find(prefix):
        matches = [el for el in original if el.tag == W+'p' and text(el).startswith(prefix)]
        assert len(matches) == 1, (prefix, len(matches))
        return matches[0]

    old_captions = {n: find(f'表S{n}：') for n in range(1, 12)}
    tables = {n: original[original.index(cap)+1] for n, cap in old_captions.items()}
    assert all(el.tag == W+'tbl' for el in tables.values())
    saved_math = Counter(E.tostring(el) for el in root.xpath('//m:oMath', namespaces=NS))
    saved_drawings = [E.tostring(el) for el in root.xpath('//w:drawing', namespaces=NS)]
    formula = next(el for el in original if el.xpath('.//m:oMath', namespaces=NS)
                   and any('(5)' in ''.join(el.itertext()) for _ in [0]))
    formula_bytes = E.tostring(formula)

    archive = ROOT / 'source_data'
    archive.mkdir(exist_ok=True)
    manifest = {}
    for n, table in tables.items():
        (archive / f'original_table_S{n}.xml').write_bytes(E.tostring(table, xml_declaration=True, encoding='UTF-8'))
        csv_name = {1: 'label_observation_design.csv', 9: 'common_window_label_validation.csv'}.get(n, f'original_table_S{n}.csv')
        csv_path = archive / csv_name
        with csv_path.open('w', encoding='utf-8-sig', newline='') as stream:
            csv.writer(stream).writerows(rows(table))
        manifest[f'S{n}'] = {'caption': text(old_captions[n]), 'csv': csv_name,
                             'sha256': sha(csv_path), 'rows': rows(table)}

    old_fig2 = find('图2：FlexPPI的数值预测与高动态位点识别')
    fig2_caption = ('图2：FlexPPI的柔性预测与高动态位点识别。Full、Intra、Static、AA分别为完整、链内、静态特征及训练集类型均值模型。'
                    '（a）预测与MD对应。蓝色RMSF读下、左轴（Å），赭色Hχ读上、右轴（nats）；共用对数概率质量色标，实线和色带为分箱预测均值及95%区间，灰点线为理想一致关系。RMSF轴在1 Å以内线性、以外对数。'
                    '（b）RMSF、总体及同类型Hχ的记录级AP分布；点为记录，粗条为四分位范围，白划为中位数。'
                    '（c）同类型Hχ PR，采用覆盖0–1的连续分段线性纵轴。'
                    '（d）615条DynaRepo复合物的外部模型比较，左PR、右ROC；c、d灰横线为阳性率，ROC灰对角线为随机参考。'
                    '（e）按界面数分组的RMSF Spearman相关、AP及标准化MAE；柱为记录均值，误差线为95%家族整簇重采样区间。计分与作图方法见2.6节及补充材料C。')
    assert len(fig2_caption) < len(text(old_fig2))*0.55

    robust = merged_table(tables[2], [
        ('A 随机初始化：原体系划分第0折，训练40轮', tables[2], [1]*6),
        ('B 单目标与联合目标：五折记录等权评价', tables[3], [2,1,1,1,1]),
        ('C 全组成链家族隔离：3,063条记录', tables[4], [1]*6)
    ], [1800, 1500, 1600, 1600, 1600, 1600])
    joint = merged_table(tables[10], [
        ('A 同类型细箱熵评价：按可定义χ角数分层', tables[10], [1]*7),
        ('B 同类型相对双读数组合：全部3,063条记录', tables[11], [2,1,1,1,1,1])
    ], [1500,1000,1450,1450,1450,1300,1550])
    for n, merged in [(2,robust), (3,robust), (4,robust), (10,joint), (11,joint)]:
        available = Counter(tuple(row) for row in rows(merged))
        required = Counter(tuple(row) for row in rows(tables[n]))
        assert not (required-available), n

    b_head = find('B 训练与划分验证')
    c_head = find('C 预测性能与方法比较')
    b_intro = ('表S1汇总训练与划分稳健性。A固定原体系第0折，以种子7、11、23分别训练Full与Intra；'
               'B的单目标模型每目标训练五次，与联合模型使用相同输入、编码器、40轮预算及961,321个RMSF和824,988个Hχ有效位置；'
               'C按2.5节的1,662个连通家族进行五折评价。ITrace的244条记录同属一个依赖家族，均在第0折测试。')
    b_note = ('联合模型在一次推理中提供两类读数，RMSF与Hχ MAE相对单目标分别降低0.8%和1.8%；'
              '家族隔离评价中，Full相对Intra的两项MAE分别降低7.5%和3.8%。各设定及全部来源结果见表S1。')
    b_cap = ('表S1：训练与划分稳健性。A为固定划分的随机初始化，B为单／联合目标，C为全组成链家族隔离；'
             '各分块采用其注明的训练与评价集合。RMSF误差为Å，Hχ误差为nats，ρ为记录等权Spearman相关。')
    new_b = [b_head, paragraph(find('固定原体系划分的第0折训练'), b_intro),
             paragraph(old_captions[2], b_cap), robust,
             paragraph(find('联合模型的RMSF与Hχ平均'), b_note)]
    joint_cap = ('表S6：联合柔性读数验证。A为83条记录的同类型细箱熵评价，按χ角数分层，ρ为去同类型均值后的记录等权Spearman相关，高减低为按Full Hχ预测划分半组的真实熵差（nats/角）。'
                 'B为全队列双读数组合；位置数为原始计数，比例和召回率为记录等权值（%），RMSF、Hχ MAE分别为Å、nats。A各层类型组成不同；B未分类位置保留在比例分母中。')

    a_intro = find('ATLAS以重复轨迹统计建立标签')
    new_a = text(a_intro).replace('各来源的观测设计、有效Hχ位置与零Hχ数量见表S1。',
                                '各来源观测设计、有效Hχ位置及零Hχ数量保存在配套数据label_observation_design.csv中。')
    e_window_cap = old_captions[9]
    e_window_note = find('表S9并列报告固定预测')
    window_note = ('固定预测在原汇总标签与第一条共同窗口标签上的核查保存在配套数据common_window_label_validation.csv中，完整保留各来源的MAE和相关。'
                   '共同窗口下Full两项MAE均低于Static。采用三条共同窗口轨迹均值标签时，ATLAS的RMSF／Hχ MAE为0.723 Å／0.172 nats，ITrace为0.158 Å／0.110 nats，DynaRepo采用一条轨迹。')
    pair_methods = find('两类读数的配对分析使用全部3,063')
    pair_note = find('表S11提供各参考类别的完整计数')
    angle_note = find('共同窗口中的33,493个位点有351个')

    replacements = {
        old_fig2: paragraph(old_fig2, fig2_caption),
        a_intro: paragraph(a_intro, new_a),
        e_window_note: paragraph(e_window_note, window_note),
        pair_methods: None,
        pair_note: None,
        old_captions[10]: paragraph(old_captions[10], joint_cap),
        tables[10]: joint,
        old_captions[5]: paragraph(old_captions[5], text(old_captions[5]).replace('表S5：', '表S2：', 1)),
        old_captions[6]: paragraph(old_captions[6], text(old_captions[6]).replace('表S6：', '表S3：', 1)),
        old_captions[7]: paragraph(old_captions[7], text(old_captions[7]).replace('表S7：', '表S4：', 1)),
        old_captions[8]: paragraph(old_captions[8], text(old_captions[8]).replace('表S8：', '表S5：', 1))
    }
    remove = {old_captions[1], tables[1], old_captions[9], tables[9], old_captions[11], tables[11]}
    pair_intro = renumber_refs(pair_methods)
    pair_details = renumber_refs(pair_note)
    graph_note = ('图2a密度包含3,063条记录的961,321个RMSF及824,988个Hχ有效位置，记录等权、记录内位置等权；'
                  '两套色标共用对数概率质量标度。分箱均值的95%色带由1,000次来源分层重抽记录得到。RMSF两轴在1 Å以内线性、以外对数，Hχ两轴线性。'
                  '图2b各方法均含3,063条记录，内粗条为第25–75百分位，白划为中位数。'
                  '图2c含42,108个有效记录–类型组，先类型后记录等权，实际平均阳性率为14.63%；图2d七种方法共用260,340个有效残基，阳性率为10.14%。'
                  '图2e按1个、2–3个、至少4个界面分组，误差线为95%家族整簇重采样区间，记录数与家族数见表S4注。')
    figure_s1 = find('图S1：内部对照的来源分层')
    insert_before = {
        figure_s1: [],
        old_captions[10]: [pair_intro, pair_details],
    }
    # Insert graph-statistics prose before the supplementary figure's image, preserving image-caption adjacency.
    s1_image = original[original.index(figure_s1)-1]
    assert s1_image.xpath('.//w:drawing', namespaces=NS)
    insert_before[s1_image] = [paragraph(find('图2c及图S1b、c的Precision轴'), graph_note)]
    intro = find('将侧链多角联合类别c的占比记为pc')
    intro_copy = deepcopy(intro); keep_with_next(intro_copy)
    replacements[intro] = intro_copy

    bs, be = original.index(b_head), original.index(c_head)
    rebuilt = []
    for i, el in enumerate(original):
        if i == bs:
            rebuilt.extend(new_b)
            continue
        if bs < i < be or el in remove:
            continue
        rebuilt.extend(insert_before.get(el, []))
        if el in replacements:
            item = replacements[el]
            if item is not None:
                rebuilt.append(item)
        else:
            rebuilt.append(renumber_refs(el) if el.tag == W+'p' else el)
    for child in list(body):
        body.remove(child)
    for child in rebuilt:
        body.append(child)

    assert E.tostring(formula) == formula_bytes
    assert Counter(E.tostring(el) for el in root.xpath('//m:oMath', namespaces=NS)) == saved_math
    assert [E.tostring(el) for el in root.xpath('//w:drawing', namespaces=NS)] == saved_drawings
    final_tables = body.findall('w:tbl', NS)
    assert len(final_tables) == 6
    captions = [text(el) for el in body if re.match(r'^表S\d+：', text(el))]
    assert [re.match(r'^表S(\d+)：', t).group(1) for t in captions] == [str(n) for n in range(1,7)]
    for n, expected_index in [(5,1), (6,2), (7,3), (8,4)]:
        assert rows(tables[n]) == rows(final_tables[expected_index])
    for token in ['70.050', '15.8', '16.5', '0.790', '0.792', '0.568', '0.828', '1.186', '0.792']:
        assert token in '\n'.join(text(el) for el in body), token
    assert all(n <= 6 for n in [int(m.group('n')) for el in body for m in REF.finditer(text(el))])
    xml = E.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
    with zipfile.ZipFile(TARGET, 'w', zipfile.ZIP_DEFLATED) as package:
        for entry in entries:
            package.writestr(entry, xml if entry.filename == 'word/document.xml' else parts[entry.filename])
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        assert [name for name, value in parts.items() if package.read(name) != value] == ['word/document.xml']
    manifest['common_window_label_validation.csv_notes'] = text(e_window_note)
    (archive/'table_archive_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    (archive/'README.md').write_text('原补充表的完整数值与OOXML备份。旧S1和S9改为配套数据，其余表格数值全部保留并合并／改号。\n'
                                    '映射：旧S2/S3/S4→新S1A/B/C；旧S5→新S2；旧S6→新S3；旧S7→新S4；旧S8→新S5；旧S10/S11→新S6A/B。\n',encoding='utf-8')
    audit = {'source':str(SOURCE),'target':str(TARGET),'target_sha256':sha(TARGET),
             'figure_2_caption_characters':{'before':len(text(old_fig2)),'after':len(fig2_caption)},
             'supplement_tables':{'before':11,'after':6,'mapping':MAPPING,'archived_as_data':[1,9]},
             'retained_table_cells_preserved_exactly':True,
             'all_original_table_data_archived':True,
             'all_native_math_and_images_preserved_exactly':True,
             'equation_5_original_paragraph_preserved_exactly':True,
             'equation_intro_kept_with_equation':True,
             'archive':str(archive),'visual_qa':'pending'}
    (ROOT/'v79_patch_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(audit,ensure_ascii=False))


if __name__ == '__main__':
    main()
