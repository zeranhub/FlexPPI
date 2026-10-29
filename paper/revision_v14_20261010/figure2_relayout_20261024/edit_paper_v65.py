"""Patch only Figure 2, its related prose/references, and new Figure S3.

Use a ZIP/OOXML patch to preserve every untouched package part byte-for-byte.
"""
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import re
import zipfile
from lxml import etree as ET
from PIL import Image

ROOT = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT = ROOT / 'figure2_relayout_20261024'
SOURCE = ROOT / 'paper_v64_20261024_中文修订_图2加多链界面子图.docx'
TARGET = ROOT / 'paper_v65_20261024_中文修订_图2取舍与外部PR_ROC正文版.docx'
NS = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
      'pic':'http://schemas.openxmlformats.org/drawingml/2006/picture'}
W = '{'+NS['w']+'}'
REL = 'http://schemas.openxmlformats.org/package/2006/relationships'


def text(el):
    return ''.join(el.xpath('.//w:t/text()', namespaces=NS))


def set_text(el, value):
    nodes=el.xpath('.//w:t', namespaces=NS)
    if not nodes:
        r=ET.SubElement(el,W+'r');nodes=[ET.SubElement(r,W+'t')]
    nodes[0].text=value
    nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for t in nodes[1:]: t.text=''


def replace_text(el, old, new):
    value=text(el)
    assert value.count(old)==1, (old,value[:100])
    set_text(el,value.replace(old,new))


def append_text(el, value):
    nodes=el.xpath('.//w:t',namespaces=NS)
    nodes[-1].text=(nodes[-1].text or '')+value


def resize_drawing(el, cx, image):
    with Image.open(image) as im: w,h=im.size
    cy=round(cx*h/w)
    for ext in el.xpath('.//wp:extent | .//a:xfrm/a:ext',namespaces=NS):
        ext.set('cx',str(cx));ext.set('cy',str(cy))
    return cy


def main():
    assert SOURCE.is_file()
    assert (OUT/'fig2_relayout.png').is_file() and (OUT/'figS3_internal_PR.png').is_file()
    with zipfile.ZipFile(SOURCE) as z:
        entries={i.filename:z.read(i.filename) for i in z.infolist()}
        infos=z.infolist()
    root=ET.fromstring(entries['word/document.xml'])
    body=root.find('w:body',NS)
    original=list(body)
    before=[ET.tostring(x) for x in original]
    assert text(original[76]).startswith('图2：')
    assert text(original[173]).startswith('图S2：')
    assert text(original[174]).startswith('D 全组成链')
    # Local method/prose edits are tied to the moved panels and added curves.
    replace_text(original[47], '四种方法在相同有效位置计分', '内部四种方法在主评价的相同有效位置计分')
    replace_text(original[47], '95%区间按来源分层重抽记录10,000次', '内部比较的95%区间按来源分层重抽记录10,000次')
    append_text(original[47], '外部PR与ROC比较使用与界面分层相同的615条DynaRepo复合物，沿用相同的逐记录阳性定义与并列分数组。ROC保留原始阈值顺序及垂直段，在共同假阳性率处逐段线性插值，再对记录等权平均；图例AP与AUROC均由各记录原始分数计算后平均，不以显示曲线积分替代。')
    replace_text(original[68], '（图2b、c）', '（图2b及图S3）')
    replace_text(original[79], '（图2d）', '（图2f）')
    replace_text(original[171], '图2c完整保留0–1的Precision范围', '图2c与图S3完整保留0–1的Precision范围')
    append_text(original[171], '同类型Hχ保留在图2c，RMSF与总体Hχ曲线移至图S3；图2d、e的外部模型曲线使用常规线性坐标。')
    # Keep the complete original descriptions of panels a and b.
    old_caption=text(original[76])
    head=old_caption.split('（c）',1)[0]
    caption=(head+
        '（c）同类型Hχ的平均PR曲线。以记录内各氨基酸类型的MD第90百分位及以上为阳性，保留并列值；42,108个有效记录–类型组先在记录内对类型等权，再对3,063条记录等权。灰横线为实际阳性率14.63%，末端保留全部入选状态；纵轴为完整覆盖0–1的连续分段线性刻度，具体映射见补充材料C。AP按Full、Intra、Static、AA顺序列出，由各记录的原始阈值曲线计算后平均。'
        '（d,e）与f相同的615条DynaRepo复合物上，FlexPPI与六个外部方法的高RMSF位点PR和ROC曲线，七种方法共用260,340个有效残基。阳性为各记录MD参考RMSF第90百分位及以上的位置；PR按共同召回率首达阈值组取精确率，ROC按共同假阳性率在原始阈值折线上插值，再对记录等权平均。图例报告逐记录AP或AUROC的均值；d中灰横线为实际平均阳性率10.14%，e中灰对角线为随机参考。BackFlip两版本均采用完整结构输入的assembly模式。'
        '（f）按模拟结构中蛋白–蛋白界面数分组的DynaRepo复合物结果：三条柱依次为1个界面（n=414）、2–3个（n=144）和3个以上（n=57），柱为记录均值，误差线为95%家族聚类重抽样区间，三幅依次报告Spearman相关、高RMSF位点AP和标准化MAE。')
    set_text(original[76],caption)
    extra=deepcopy(original[78])
    set_text(extra,'在相同615条DynaRepo复合物上，完整模型对高RMSF位点的记录等权平均AP为0.659，AUROC为0.917；FastProtFlex的对应值为0.625和0.903，BackFlip-2.1为0.599和0.891（图2d、e）。两类曲线给出不同阈值下的筛选表现，其他方法的平均AP及AUROC同时列于图例。该结果对应复合物子集，与全队列的幅值估计和内部对照分别评价。')
    # Avoid duplicate Word paragraph identifiers in the cloned paragraph.
    for attr in list(extra.attrib):
        if attr.rsplit('}',1)[-1] in ('paraId','textId'): del extra.attrib[attr]
    original[80].addnext(extra)
    # Main figure replaces its existing image relationship, preserving all others.
    main_drawing=original[75]
    cx=int(main_drawing.xpath('.//wp:extent',namespaces=NS)[0].get('cx'))
    main_cy=resize_drawing(main_drawing,cx,OUT/'fig2_relayout.png')
    main_pr=main_drawing.xpath('.//wp:docPr',namespaces=NS)[0]
    main_pr.set('descr','a prediction/reference correspondence; b original record AP distributions; c within-type Hchi PR; d external RMSF PR; e external RMSF ROC; f unchanged three interface-stratified metrics')
    # New supplementary figure follows S2 and precedes the next appendix section.
    rels=ET.fromstring(entries['word/_rels/document.xml.rels'])
    ids=[int(x.get('Id')[3:]) for x in rels if re.fullmatch(r'rId\d+',x.get('Id',''))]
    rid='rId'+str(max(ids)+1)
    new_rel=ET.SubElement(rels,'{'+REL+'}Relationship')
    new_rel.set('Id',rid);new_rel.set('Type',NS['r']+'/image')
    new_rel.set('Target','media/figureS3_internal_PR.png')
    s3=deepcopy(original[172])
    for attr in list(s3.attrib):
        if attr.rsplit('}',1)[-1] in ('paraId','textId'): del s3.attrib[attr]
    s3.xpath('.//a:blip',namespaces=NS)[0].set('{'+NS['r']+'}embed',rid)
    drawing_ids=[int(x.get('id')) for x in root.xpath('.//wp:docPr',namespaces=NS)]
    pr=s3.xpath('.//wp:docPr',namespaces=NS)[0]
    pr.set('id',str(max(drawing_ids)+1));pr.set('name','Supplementary Figure S3 internal RMSF and overall Hchi PR')
    pr.set('descr','Two original internal-baseline PR curves moved from Figure 2c without changes to their data or precision-axis mappings.')
    pic=s3.xpath('.//pic:cNvPr',namespaces=NS)
    if pic: pic[0].set('name','figureS3_internal_PR.png')
    s3_cy=resize_drawing(s3,cx,OUT/'figS3_internal_PR.png')
    s3_caption=deepcopy(original[173])
    set_text(s3_caption,'图S3：全队列内部对照的RMSF与总体Hχ高动态位点PR曲线。（a）RMSF；（b）总体Hχ。Full、Intra、Static、AA与图2相同；以各记录MD参考第90百分位及以上为阳性，保留真实值与预测分数的并列组，各曲线在共同召回率处取首次达到该召回率的阈值组精确率，再对3,063条记录等权平均。图例或上方数值为逐记录原始阈值计算后平均的AP；灰横线为实际平均阳性率，分别为10.26%和10.29%，末端保留全部入选状态。纵轴沿用连续分段线性刻度并完整覆盖0–1，映射见补充材料C；显示曲线面积不代替平均AP。同类型Hχ曲线见图2c。')
    original[173].addnext(s3);s3.addnext(s3_caption)
    modified_original=[i for i,x in enumerate(original) if ET.tostring(x)!=before[i]]
    assert set(modified_original)=={47,68,75,76,79,171},modified_original
    replacement={
        'word/document.xml':ET.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),
        'word/_rels/document.xml.rels':ET.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True),
        'word/media/image2.png':(OUT/'fig2_relayout.png').read_bytes(),
        'word/media/figureS3_internal_PR.png':(OUT/'figS3_internal_PR.png').read_bytes()}
    with zipfile.ZipFile(TARGET,'w',zipfile.ZIP_DEFLATED) as z:
        for info in infos:
            z.writestr(info,replacement.get(info.filename,entries[info.filename]))
        z.writestr('word/media/figureS3_internal_PR.png',replacement['word/media/figureS3_internal_PR.png'])
    with zipfile.ZipFile(TARGET) as z:
        assert z.testzip() is None
        changed=[name for name in entries if z.read(name)!=entries[name]]
        assert set(changed)=={'word/document.xml','word/_rels/document.xml.rels','word/media/image2.png'},changed
        all_text=''.join(ET.fromstring(z.read('word/document.xml')).xpath('//w:t/text()',namespaces=NS))
        assert '图S3：' in all_text and '（图2f）' in all_text and '（图2d、e）' in all_text
    report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'target_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),
            'modified_original_body_indices':modified_original,'inserted_paragraphs':3,
            'changed_package_parts':changed,'new_package_parts':['word/media/figureS3_internal_PR.png'],
            'all_other_package_parts_byte_identical':True,'main_figure_extent_emu':[cx,main_cy],
            'supp_figure_extent_emu':[cx,s3_cy], 'main_caption_characters':len(caption)}
    (OUT/'docx_patch_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__': main()
