"""Limit main-text 2.7 to MD parameters and preserve case methods in supplement D."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib
import json
import re
import zipfile
from lxml import etree as E

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'paper_v76_20261024_中文修订_方法去除实例概述.docx'
TARGET = ROOT.parent / 'paper_v77_20261024_中文修订_2.7仅保留MD参数.docx'
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
        if child.tag != W + 'pPr':
            el.remove(child)
    run = E.SubElement(el, W + 'r')
    if properties is not None:
        run.append(deepcopy(properties))
    node = E.SubElement(run, W + 't')
    node.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    node.text = value
    return el


def replace_phrase(el, old, new):
    """Replace one phrase while retaining all original run properties."""
    el = deepcopy(el)
    assert text(el).count(old) == 1
    start = text(el).index(old)
    end = start + len(old)
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
    return el


def main():
    assert sha(SOURCE) == '394f667b47cad4478eb062a25157bab59d6a2eef4a1b810df9accf425b4ba3ab'
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

    h27 = find('2.7 结构体系构建与轨迹分析')
    h28 = find('2.8 联合预测与构象占据的验证')
    hd = find('D 复合物实例补充验证')
    he = find('E 联合柔性与构象占据验证')
    start, end = original.index(h27), original.index(h28)
    ds, de = original.index(hd), original.index(he)
    old_methods = original[start+1:end]
    assert len(old_methods) == 6
    kix, md, tcr, contacts, ub, sasa = old_methods
    sd = original[ds:de]
    assert len(sd) == 6
    result_tcr = find('主链与侧链柔性的联合解析能够将TCR')
    fig4 = find('图4：MEL5体系中受体亲和成熟')
    saved_math = Counter(E.tostring(el) for el in root.xpath('//m:oMath', namespaces=NS))
    saved_tables = [E.tostring(el) for el in root.xpath('//w:tbl', namespaces=NS)]
    saved_drawings = [E.tostring(el) for el in root.xpath('//w:drawing', namespaces=NS)]

    md_text = ('本研究新增的全原子MD模拟使用OpenMM 8.4、CHARMM36力场及其配套TIP3P水模型。'
               '以pH 7.0的构建选项加入氢原子，组氨酸采用HID形式，链端带电且未封端。'
               '各体系置于立方体周期性水盒，初始蛋白与盒边至少相隔1.0 nm，加入0.1 mol/L NaCl并电荷中和。'
               '长程静电使用粒子网格Ewald方法（PME），非键作用截断为1.2 nm、开关距离为1.0 nm；'
               '约束所有涉及氢原子的共价键，并使用刚性水分子，积分步长为2 fs。'
               '经两阶段L-BFGS能量最小化后，在298 K下采用Langevin恒温器，阻尼系数为5 ps⁻¹，采用恒温恒容（NVT）模拟。'
               '先运行10 ns蛋白重原子位置约束平衡，约束势为U=k d²、k=4,184 kJ mol⁻¹ nm⁻²，'
               '再解除位置约束运行100 ns生产模拟，每1 ps保存蛋白坐标。体系构建、已有轨迹来源及分析规则见补充材料D。')
    new_methods = [paragraph(h27, '2.7 分子动力学模拟参数'), paragraph(md, md_text)]
    for name in ['KIX', 'cMyb', 'MLL', 'MEL5', 'TCR', '泛素', 'UBA', '2AGH', 'ITrace', 'A003R']:
        assert name not in md_text
    for setting in ['OpenMM 8.4', 'CHARMM36', 'TIP3P', 'pH 7.0', 'HID', '1.0 nm', '0.1 mol/L',
                    'PME', '1.2 nm', '2 fs', 'L-BFGS', '298 K', '5 ps⁻¹', 'NVT', '10 ns',
                    'U=k d²', '4,184 kJ mol⁻¹ nm⁻²', '100 ns', '1 ps']:
        assert setting in md_text

    kix_analysis = ('KIX体系在完整三链结构中补齐缺失重原子后派生各状态，五个组氨酸均采用HID形式。'
                    '各状态按2.7节进行模拟；每种状态分析一条100 ns生产轨迹，以全部87个Cα对齐生产首帧后按式(1)计算RMSF；'
                    '另核查逐帧共价键几何，并比较前后半段的指标。留出模型在相同位置的定位与幅值见表S8。')
    ub_validation = text(sd[5])
    prefix = '泛素–UBA的结构环境分析采用与图3一致的A003R完整双链起始坐标。'
    assert ub_validation.startswith(prefix)
    ub_validation = ub_validation[len(prefix):].replace('邻域接触和计算设置见配套源表及2.7节。',
                                                    '侧链面积及邻域接触数据见配套源表。')
    new_d = [paragraph(hd, 'D 复合物体系与补充验证'), deepcopy(kix), paragraph(sd[1], kix_analysis),
             deepcopy(sd[2]), deepcopy(sd[3]), deepcopy(sd[4]), deepcopy(tcr), deepcopy(contacts),
             deepcopy(ub), deepcopy(sasa), paragraph(sd[5], ub_validation)]
    for el in [kix, tcr, contacts, ub, sasa]:
        assert any(text(el) == text(item) for item in new_d)

    new_result_tcr = replace_phrase(result_tcr, '模型设置见方法2.5及2.7节。',
                                   '模型设置见2.5节，结构输入与对照定义见补充材料D。')
    new_fig4 = replace_phrase(fig4, '状态定义见方法。', '状态定义见补充材料D。')
    replacement = {result_tcr: new_result_tcr, fig4: new_fig4}
    rebuilt = []
    for i, el in enumerate(original):
        if i == start:
            rebuilt.extend(new_methods)
        elif start < i < end:
            continue
        elif i == ds:
            rebuilt.extend(new_d)
        elif ds < i < de:
            continue
        else:
            rebuilt.append(replacement.get(el, el))
    for child in list(body):
        body.remove(child)
    for child in rebuilt:
        body.append(child)

    assert Counter(E.tostring(el) for el in root.xpath('//m:oMath', namespaces=NS)) == saved_math
    assert [E.tostring(el) for el in root.xpath('//w:tbl', namespaces=NS)] == saved_tables
    assert [E.tostring(el) for el in root.xpath('//w:drawing', namespaces=NS)] == saved_drawings
    excluded = set(original[start:end] + original[ds:de] + [result_tcr, fig4])
    assert all(any(el is item for item in rebuilt) for el in original if el not in excluded)
    full = '\n'.join(text(el) for el in body)
    assert len(re.findall(r'(?<!\d)2\.7(?!\d)', full)) == 2
    for name in ['2AGH', 'A00AQ', '3HG1', '4JFF', '4JFD', '4JFE', '2QHO', 'A003R',
                 'NE2', '150°', '5,000', '1,000', '27类', 'dot_solvent=1', 'dot_density=4']:
        assert name in '\n'.join(text(el) for el in new_d)
    xml = E.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
    with zipfile.ZipFile(TARGET, 'w', zipfile.ZIP_DEFLATED) as package:
        for entry in entries:
            package.writestr(entry, xml if entry.filename == 'word/document.xml' else parts[entry.filename])
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        assert [name for name, value in parts.items() if package.read(name) != value] == ['word/document.xml']
    audit = {'source': str(SOURCE), 'target': str(TARGET), 'target_sha256': sha(TARGET),
             'main_2_7': [text(el) for el in new_methods],
             'supplement_D': [text(el) for el in new_d],
             'cross_references_updated': {'3.4': '2.5 and supplement D', 'Figure 4': 'supplement D',
                                          'supplement D SASA': 'local method and source tables'},
             'only_parts_changed': ['word/document.xml'],
             'all_native_math_tables_and_figures_preserved_exactly': True,
             'all_other_body_elements_preserved_exactly': True,
             'case_methods_preserved_in_supplement': True, 'visual_qa': 'pending'}
    (ROOT / 'v77_patch_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in audit.items() if k not in ['main_2_7', 'supplement_D']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
