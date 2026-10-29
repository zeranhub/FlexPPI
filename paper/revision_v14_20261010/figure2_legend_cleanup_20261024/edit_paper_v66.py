"""Apply the user's legend annotations and move external metrics to Table S10.

Patch the existing OOXML package. All unrelated package parts remain identical.
"""
from pathlib import Path
from copy import deepcopy
import csv
import hashlib
import json
import zipfile
from lxml import etree as ET
from PIL import Image

ROOT=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT=ROOT/'figure2_legend_cleanup_20261024'
SOURCE=ROOT/'paper_v65_20261024_中文修订_图2取舍与外部PR_ROC正文版.docx'
TARGET=ROOT/'paper_v66_20261024_中文修订_图2图例内置与指标移补充.docx'
METRICS=ROOT/'figure2_external_curves_20261024'/'summary_metrics.csv'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
W='{'+NS['w']+'}'
METHODS=[('Full','FlexPPI (Full)'),('backflip-2.1 (assembly)','BackFlip-2.1 (assembly)'),
    ('backflip-1.0 (assembly)','BackFlip-1.0 (assembly)'),('DynaProt-M','DynaProt-M'),
    ('FastProtFlex','FastProtFlex'),('ANM','ANM'),('GNM','GNM')]


def text(el): return ''.join(el.xpath('.//w:t/text()',namespaces=NS))

def set_text(el,value):
    nodes=el.xpath('.//w:t',namespaces=NS)
    assert nodes
    nodes[0].text=value
    nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for node in nodes[1:]: node.text=''


def replace(el,old,new):
    old_text=text(el);assert old_text.count(old)==1,(old,old_text[:120])
    set_text(el,old_text.replace(old,new))


def clean_identifiers(el):
    for node in el.iter():
        for key in list(node.attrib):
            if key.rsplit('}',1)[-1] in ('paraId','textId'):del node.attrib[key]


def clone_paragraph(source,value):
    result=deepcopy(source);clean_identifiers(result);set_text(result,value);return result


def table_from_template(template,values):
    # Match the manuscript's existing three-line table typography and borders.
    table=deepcopy(template);clean_identifiers(table)
    grid=table.find('w:tblGrid',NS)
    widths=[4608,2160,2160]
    for col in list(grid):grid.remove(col)
    for width in widths:
        node=ET.SubElement(grid,W+'gridCol');node.set(W+'w',str(width))
    rows=table.findall('w:tr',NS)
    assert len(rows)==8
    content=[['方法','AP ↑','AUROC ↑']]+values
    for ri,(row,record) in enumerate(zip(rows,content)):
        cells=row.findall('w:tc',NS);assert len(cells)==4
        row.remove(cells[3]);cells=cells[:3]
        for ci,(cell,value,width) in enumerate(zip(cells,record,widths)):
            cell.find('w:tcPr/w:tcW',NS).set(W+'w',str(width))
            set_text(cell,value)
            # Preserve the frozen template's font and paragraph alignment.
    assert len(table.findall('w:tblGrid/w:gridCol',NS))==3
    assert all(len(r.findall('w:tc',NS))==3 for r in rows)
    return table


def main():
    assert (OUT/'fig2_legend_cleanup.png').is_file()
    with zipfile.ZipFile(SOURCE) as z:
        infos=z.infolist();entries={i.filename:z.read(i.filename) for i in infos}
    doc=ET.fromstring(entries['word/document.xml']);body=doc.find('w:body',NS)
    original=list(body);before=[ET.tostring(x) for x in original]
    def paragraph(prefix):
        matches=[x for x in original if ET.QName(x).localname=='p' and text(x).startswith(prefix)]
        assert len(matches)==1,(prefix,len(matches));return matches[0]
    caption=paragraph('图2：')
    replace(caption,'AP按Full、Intra、Static、AA顺序列出，由各记录的原始阈值曲线计算后平均。',
        '内部对照方法的图例置于c内，供b、c共同读取。')
    replace(caption,'图例报告逐记录AP或AUROC的均值；',
        '七种方法共用e内的图例，逐记录AP与AUROC的均值见补充表S10；')
    replace(caption,'三条柱依次为1个界面（n=414）、2–3个（n=144）和3个以上（n=57），柱为记录均值，',
        '三条柱依次为1个界面（414条记录、196个家族）、2–3个界面（144条记录、58个家族）和至少4个界面（57条记录、34个家族）；同一家族可出现在不同组。柱为记录均值，')
    method=paragraph('高动态位点筛选')
    replace(method,'图例AP与AUROC均由各记录原始分数计算后平均','AP与AUROC均由各记录原始分数计算后平均')
    result=paragraph('在相同615条DynaRepo复合物上')
    replace(result,'其他方法的平均AP及AUROC同时列于图例。','各方法的平均AP及AUROC汇总于补充表S10。')
    # The main figure is the only existing image to change.
    image_el=next(x for x in original if x.xpath('.//wp:docPr[@id="2"]',namespaces=NS))
    extent=image_el.xpath('.//wp:extent',namespaces=NS)[0]
    cx=int(extent.get('cx'))
    with Image.open(OUT/'fig2_legend_cleanup.png') as im:
        cy=round(cx*im.height/im.width)
    for node in image_el.xpath('.//wp:extent | .//a:xfrm/a:ext',namespaces=NS):
        node.set('cx',str(cx));node.set('cy',str(cy))
    image_el.xpath('.//wp:docPr',namespaces=NS)[0].set('descr',
        'a prediction/reference correspondence; b record AP distributions; c within-type Hchi PR and internal method legend; d external RMSF PR; e external RMSF ROC with shared external method legend; f original interface-stratified metrics. Numerical external AP/AUROC are in Supplementary Table S10.')
    # Append the numerical table after S9 so all existing appendix/table labels stay.
    with METRICS.open(encoding='utf-8-sig',newline='') as f:raw=list(csv.DictReader(f))
    lookup={(r['method'],r['metric']):float(r['mean']) for r in raw if r['scope']=='all'}
    values=[[name,f'{lookup[(method,"AP")]:.3f}',f'{lookup[(method,"AUROC")]:.3f}'] for method,name in METHODS]
    header=clone_paragraph(paragraph('G 双读数组合'),'H 外部模型的高RMSF位点筛选')
    note=clone_paragraph(paragraph('正文图5a的两类读数配对分析'),
        '表S10汇总图2d、e中七种方法的高RMSF位点筛选结果。阳性为各记录MD参考RMSF第90百分位及以上的位置，并列值全部保留。AP与AUROC均在各记录的相同有效残基上计算后对记录等权平均；AP依据原始阈值组计算，显示的平均曲线面积不替代这些数值。')
    tab_caption=clone_paragraph(paragraph('表S4：'),
        '表S10：DynaRepo复合物高RMSF位点筛选的平均AP与AUROC。七种方法共同覆盖615条记录、261个序列家族及260,340个有效残基，两项指标均越高越好。BackFlip两个版本均使用完整结构输入的assembly模式。')
    # Keep the new caption with the new table.
    ppr=tab_caption.find('w:pPr',NS)
    keep=ppr.find('w:keepNext',NS)
    if keep is None: keep=ET.SubElement(ppr,W+'keepNext')
    keep.set(W+'val','1')
    template=next(x for x in original if ET.QName(x).localname=='tbl' and 'PLANET full' in text(x))
    table=table_from_template(template,values)
    assert ET.QName(body[-1]).localname=='sectPr'
    for node in [header,note,tab_caption,table]:body.insert(len(body)-1,node)
    modified=[i for i,node in enumerate(original) if ET.tostring(node)!=before[i]]
    assert set(modified)=={47,75,76,81},modified
    changes={'word/document.xml':ET.tostring(doc,xml_declaration=True,encoding='UTF-8',standalone=True),
             'word/media/image2.png':(OUT/'fig2_legend_cleanup.png').read_bytes()}
    with zipfile.ZipFile(TARGET,'w',zipfile.ZIP_DEFLATED) as z:
        for info in infos:z.writestr(info,changes.get(info.filename,entries[info.filename]))
    with zipfile.ZipFile(TARGET) as z:
        assert z.testzip() is None
        changed=[name for name,value in entries.items() if z.read(name)!=value]
        assert set(changed)==set(changes)
        check=ET.fromstring(z.read('word/document.xml'))
        all_text=''.join(check.xpath('//w:t/text()',namespaces=NS))
        assert all_text.count('表S10：')==1
        assert 'AP按Full、Intra、Static、AA顺序列出' not in all_text
        assert '其他方法的平均AP及AUROC同时列于图例' not in all_text
    with (OUT/'table_S10_source.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(['method','AP','AUROC']);writer.writerows(values)
    report={'source':str(SOURCE),'target':str(TARGET),
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'target_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),
        'modified_original_body_indices':modified,'changed_package_parts':changed,
        'all_other_package_parts_byte_identical':True,'new_body_nodes':4,
        'supplementary_table':'S10','table_rows':8,'table_columns':3,
        'metrics_source':str(METRICS),'figure_extent_emu':[cx,cy]}
    (OUT/'docx_patch_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':main()

