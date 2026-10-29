"""Shorten comparison and occupancy methods; retain native equation 5 and reorder MD."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib
import json
import re
import zipfile
from lxml import etree as E

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'paper_v77_20261024_中文修订_2.7仅保留MD参数.docx'
TARGET = ROOT.parent / 'paper_v78_20261024_中文修订_方法精简与顺序调整.docx'
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
    el = deepcopy(el)
    assert text(el).count(old) == 1, old
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
    assert sha(SOURCE) == '50d11b1873c4b54163d412e71e49114dfbd5d02f668c19b4d5f30f61237ef0e3'
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

    h26 = find('2.6 对照方法与统计指标')
    hmd = find('2.7 分子动力学模拟参数')
    hjoint = find('2.8 联合预测与构象占据的验证')
    h3 = find('3 结果')
    a, b, c, d = map(original.index, [h26, hmd, hjoint, h3])
    old26, oldmd, oldjoint = original[a+1:b], original[b+1:c], original[c+1:d]
    assert len(old26) == 6 and len(oldmd) == 1 and len(oldjoint) == 8
    equation = oldjoint[1]
    assert len(equation.xpath('.//m:oMath', namespaces=NS)) == 1
    equation_bytes = E.tostring(equation)
    saved_math = Counter(E.tostring(el) for el in root.xpath('//m:oMath', namespaces=NS))
    saved_tables = [E.tostring(el) for el in root.xpath('//w:tbl', namespaces=NS)]
    saved_drawings = [E.tostring(el) for el in root.xpath('//w:drawing', namespaces=NS)]

    concise26 = [
        '内部以氨基酸类型均值（AA）、47项静态特征回归（Static）及链内模型（Intra）为对照，AA和Static使用对应训练集合拟合。外部RMSF比较采用BackFlip-2.1、BackFlip-1.0、DynaProt-M、FastProtFlex、ANM和GNM，与Full共同在615条DynaRepo复合物的相同有效残基上计分。BackFlip采用完整多链assembly输入，DynaProt-M逐链推理后合并；配置与实现见补充材料C。',
        'MAE和RMSE评价幅值误差，Spearman与Pearson相关评价预测与参考的一致性；各指标先逐记录计算，再对记录等权平均。标准化MAE将每条记录的预测与参考分别减去自身均值、除以自身标准差后计算平均绝对差。FastProtFlex、ANM和GNM输出相对分数，不报告Å单位MAE；同类型与体系等权汇总规则见补充材料C。',
        '高动态阳性取各记录MD参考的第90百分位及以上，保留并列参考；以AP和AUROC评价筛选。Hχ另在同记录同类型至少4个位点内重定阳性，先对类型、再对记录等权平均。PR、ROC由逐记录曲线等权汇总；AP与AUROC由各记录原始分数计算后平均，显示曲线面积不替代计分。并列分数、曲线构建及置信区间见补充材料C。'
    ]
    concise_joint_intro = '将侧链多角联合类别c的占比记为pc，最大占比记为pmax，以最常见类别为零点定义无量纲相对经验投影自由能，式(2)可写为'
    concise_joint = [
        '式(5)连接侧链标签与多角联合状态占据。另以MD逐角轨迹的低相对自由能角区总宽度与细箱熵检验预测所识别的取向范围[3,24,25]。按5°箱统计概率密度，令F*(q)=−ln[p̂(q)/p̂max]，累加F*≤1的全部箱宽（包括分离角区），并计算Shannon熵。两指标均对四个箱格原点及全部可定义χ角取平均，单位分别为°/角和nats/角；具体定义见补充材料E。',
        '双读数互补性在3,063条记录的824,988个共同有效位置上评价，采用记录内Spearman相关及同类型高低组合召回。构象占据验证使用83条记录，统一取0–100 ns、200 ps间隔的501帧；主分析取每记录第一条轨迹，并以可用重复指标均值复核。保留同记录同类型至少4个位点的33,142个共同有效位置，对轨迹指标及Full Hχ、Intra Hχ、Full RMSF留出预测分别减去记录内同类型均值，再计算相关并对记录等权平均。高低组按预测排序确定，使用MD指标评价组差。样本筛选、分组与统计细节见补充材料E。'
    ]
    new_main = [h26] + [paragraph(old26[0], t) for t in concise26]
    new_main += [paragraph(hjoint, '2.7 联合预测与构象占据验证'),
                 paragraph(oldjoint[0], concise_joint_intro), equation]
    new_main += [paragraph(oldjoint[2], t) for t in concise_joint]
    new_main += [paragraph(hmd, '2.8 分子动力学模拟参数'), oldmd[0]]

    old26_length = sum(len(text(el)) for el in old26)
    new26_length = sum(map(len, concise26))
    oldjoint_length = sum(len(text(el)) for el in oldjoint)
    newjoint_length = len(concise_joint_intro) + sum(map(len, concise_joint))
    assert new26_length < old26_length * 0.70
    assert newjoint_length < oldjoint_length * 0.70

    hc = find('C 预测性能与方法比较')
    he = find('E 联合柔性与构象占据验证')
    c_config = [
        text(old26[0]),
        '外部比较中，BackFlip两个版本均在assembly模式下一次输入全部链的主链坐标；DynaProt-M与FastProtFlex使用作者发布的权重及默认推理设置，DynaProt-M按作者实现逐链推理后合并。ANM与GNM由完整输入的Cα构建弹性网络。具有绝对单位的输出统一换算为Å后计算幅值误差；FastProtFlex、ANM和GNM输出相对分数，用于排序、筛选及记录内标准化误差评价。',
        '除MAE、RMSE及Spearman/Pearson相关外，另计算预测减参考的平均偏差。相关要求至少4个有效位置且预测与参考均有变化。体系等权汇总先在同一来源的生物体系内平均相关记录，再对体系等权；同类型误差先在每条记录计算，再对包含该类型有效位置的记录等权平均。AA在同类型内为常数，采用MAE评价。固定划分的随机初始化对照分别报告各模型的留出成绩。',
        '内部四种方法在主评价的相同有效位置计分，并列预测作为同一阈值组。AP按召回率增量乘对应精确率求和；AUROC为所有阳性–阴性位点对中阳性分数更高的比例，并列计半分。PR在共同召回率处取首次达到该召回率的阈值组精确率，再对记录等权平均。外部ROC保留原始阈值顺序及垂直段，在共同假阳性率处逐段线性插值后对记录等权平均；外部PR、ROC与界面分层共用615条DynaRepo复合物及相同阳性、并列规则。',
        'Hχ同类型筛选要求同记录同类型至少4个位点，阳性重新取该组参考第90百分位及以上，先对有效类型、再对记录等权平均；单类记录或类型组保留覆盖统计而不计分。来源分层AP在各来源内对记录等权平均。内部比较的95%区间按来源分层重抽记录10,000次。预测–参考密度及分箱均值对记录等权、记录内位置等权。逐残基曲线按作者编号和实际链边界排列，无效位置断开。'
    ]
    c_insert = [paragraph(find('表S5汇总与图2一致'), t) for t in c_config]

    e_cohort = find('种子20261010在来源内各抽30条记录')
    e_cohort_new = text(e_cohort) + 'Full与Intra采用对应留出折的原始预测；共同位置要求原RMSF与Hχ标签、轨迹角度及Cα均有效，且预测值有限。提取角度时保留化学名称校正及逐残基映射。'
    e_width = ('逐角以5°等宽箱统计MD概率密度p̂(q)，定义F*(q)=−ln[p̂(q)/p̂max]，最高密度位置为零。'
               '低相对自由能角区总宽度为全部F*≤1角度箱的宽度之和，等价于密度至少为最高值1/e的角区总长度，分离的合格箱均计入。'
               '对0°、1.25°、2.5°、3.75°四个箱格原点分别计算后平均，再对该残基全部可定义角平均，单位为°/角；角度按圆周处理，主阈值固定为1。'
               '平均细箱角度占据熵先对每个可定义χ角计算5°箱Shannon熵，再依次对相同四原点及全部角平均，单位为nats/角。'
               '两指标使用完全相同的位置和观测窗口。记录内同类型至少4个共同有效位置时计分；各轨迹指标及Full Hχ、Intra Hχ、Full RMSF预测分别去除记录内同类型均值，在相同位置计算Spearman相关后对记录等权平均。')
    e_groups = text(oldjoint[-1]).replace('位置筛选仅使用各预测读数：', '构象占据的高低组筛选仅使用各预测读数：', 1)
    e_pair = text(oldjoint[3])
    e_group_results = find('共同窗口中的33,493个位点有351个')
    e_pair_results = find('正文图5a的两类读数配对分析采用824,988')
    e_pair_results_new = text(e_pair_results).replace('正文图5a的两类读数配对分析采用824,988个共同有效位置。', '', 1)

    fig2 = find('图2：FlexPPI的数值预测与高动态位点识别')
    fig5 = find('图5：联合动态读数与MD构象占据的对应')
    dkix = find('KIX体系在完整三链结构中补齐')
    dub = find('泛素–UBA分析使用DynaRepo记录A003R')
    replacements = {
        e_cohort: paragraph(e_cohort, e_cohort_new),
        e_group_results: replace_phrase(e_group_results, '按2.8节处理', '按本节规则处理'),
        e_pair_results: paragraph(e_pair_results, e_pair_results_new),
        fig2: replace_phrase(fig2, '阈值及平均方法见2.6节', '阈值与计分见2.6节，曲线构建见补充材料C'),
        fig5: replace_phrase(fig5, '2.8节指标', '2.7节指标'),
        dkix: replace_phrase(dkix, '按2.7节进行模拟', '按2.8节进行模拟'),
        dub: replace_phrase(dub, '按2.8节对四个箱格原点取平均', '按补充材料E对四个箱格原点取平均')
    }
    inserts_after = {hc: c_insert,
                     e_cohort: [paragraph(e_cohort, e_width), paragraph(e_cohort, e_groups)]}
    inserts_before = {e_pair_results: [paragraph(e_pair_results, e_pair)]}
    rebuilt = []
    for i, el in enumerate(original):
        if i == a:
            rebuilt.extend(new_main)
            continue
        if a < i < d:
            continue
        rebuilt.extend(inserts_before.get(el, []))
        rebuilt.append(replacements.get(el, el))
        rebuilt.extend(inserts_after.get(el, []))
    for child in list(body):
        body.remove(child)
    for child in rebuilt:
        body.append(child)

    assert E.tostring(equation) == equation_bytes
    assert Counter(E.tostring(el) for el in root.xpath('//m:oMath', namespaces=NS)) == saved_math
    assert [E.tostring(el) for el in root.xpath('//w:tbl', namespaces=NS)] == saved_tables
    assert [E.tostring(el) for el in root.xpath('//w:drawing', namespaces=NS)] == saved_drawings
    excluded = set(original[a:d] + list(replacements))
    for el in original:
        if el not in excluded:
            assert any(el is item for item in rebuilt)
    assert oldmd[0] in rebuilt
    full = '\n'.join(text(el) for el in body)
    sections = [text(el) for el in body if re.match(r'2\.[678] ', text(el))]
    assert sections == ['2.6 对照方法与统计指标', '2.7 联合预测与构象占据验证', '2.8 分子动力学模拟参数']
    refs = [{'body_index': i, 'text': text(el)} for i, el in enumerate(body)
            if re.search(r'(?<!\d)2\.[678](?!\d)', text(el))]
    for old_token in ['120轮', '0.07', '31个', '50', 'L2正则为1', '63', '203条', '33,493',
                      '33,142', '824,988', '25与第75', '1.25°', '2.5°', '3.75°', '20261013',
                      '10,000次', '1/e', '奇数类型', '输入残基顺序']:
        assert old_token in full, old_token
    xml = E.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
    with zipfile.ZipFile(TARGET, 'w', zipfile.ZIP_DEFLATED) as package:
        for entry in entries:
            package.writestr(entry, xml if entry.filename == 'word/document.xml' else parts[entry.filename])
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        assert [name for name, value in parts.items() if package.read(name) != value] == ['word/document.xml']
    audit = {'source': str(SOURCE), 'target': str(TARGET), 'target_sha256': sha(TARGET),
             'section_order': sections,
             'main_text_characters': {'2.6_before': old26_length, '2.6_after': new26_length,
                                      'joint_before': oldjoint_length, 'joint_after': newjoint_length},
             'new_2_6': concise26, 'new_2_7': [concise_joint_intro] + concise_joint,
             'equation_5_paragraph_preserved_exactly': True,
             'all_native_math_tables_and_figures_preserved_exactly': True,
             'MD_parameter_paragraph_preserved_exactly': True,
             'technical_details_preserved_in_supplement_C_E': True,
             'only_parts_changed': ['word/document.xml'], 'cross_references': refs,
             'visual_qa': 'pending'}
    (ROOT / 'v78_patch_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in audit.items() if k not in ['new_2_6', 'new_2_7', 'cross_references']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
