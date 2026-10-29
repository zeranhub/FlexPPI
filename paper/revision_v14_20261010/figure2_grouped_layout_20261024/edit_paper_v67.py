"""Group external PR/ROC, restore b legend, and keep the caption compact."""
from pathlib import Path
import hashlib,json,zipfile
from lxml import etree as ET
from PIL import Image

ROOT=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT=ROOT/'figure2_grouped_layout_20261024'
SOURCE=ROOT/'paper_v66_20261024_中文修订_图2图例内置与指标移补充.docx'
TARGET=ROOT/'paper_v67_20261024_中文修订_图2组合子图与图例重排.docx'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main'}


def text(el):return ''.join(el.xpath('.//w:t/text()',namespaces=NS))

def set_text(el,value):
    nodes=el.xpath('.//w:t',namespaces=NS);assert nodes
    nodes[0].text=value
    nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for x in nodes[1:]:x.text=''


def replace(el,old,new):
    source=text(el);assert source.count(old)==1,(old,source[:100])
    set_text(el,source.replace(old,new))


def main():
    assert (OUT/'fig2_grouped.png').is_file()
    with zipfile.ZipFile(SOURCE) as z:
        infos=z.infolist();parts={i.filename:z.read(i.filename) for i in infos}
    root=ET.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
    original=list(body);before=[ET.tostring(x) for x in original]
    def para(prefix):
        matches=[x for x in original if ET.QName(x).localname=='p' and text(x).startswith(prefix)]
        assert len(matches)==1,(prefix,len(matches));return matches[0]
    caption=para('图2：');head=text(caption).split('（c）',1)[0]
    new_caption=(head+
        '（c）同类型Hχ PR曲线。42,108个有效记录–类型组先在记录内对类型等权，再对3,063条记录等权；阳性定义同2.6节，灰横线为实际阳性率14.63%。纵轴沿用完整覆盖0–1的连续分段线性刻度（补充材料C），保留全部入选末端。'
        '（d）615条DynaRepo复合物的外部模型RMSF筛选比较，左为PR、右为ROC，七种方法共用260,340个有效残基，曲线逐记录等权平均。PR灰横线为阳性率10.14%，ROC灰对角线为随机参考；阈值及平均方法见2.6节，AP与AUROC见表S10。BackFlip两个版本均采用assembly模式。'
        '（e）按界面数分组的DynaRepo复合物结果，三幅依次为Spearman相关、高RMSF位点AP和标准化MAE。1个、2–3个和至少4个界面的记录数分别为414、144和57；柱为记录均值，误差线为95%家族整簇重抽样区间，家族数见补充材料H。')
    set_text(caption,new_caption)
    replace(para('为检验结构编码在多链环境中的作用'), '（图2f）','（图2e）')
    replace(para('在相同615条DynaRepo复合物上'), '（图2d、e）','（图2d）')
    replace(para('图2b给出各任务全部记录的AP分布'), '图2d、e的外部模型曲线','图2d的外部模型PR和ROC曲线')
    note=para('表S10汇总图2d、e中七种方法')
    replace(note,'图2d、e','图2d')
    set_text(note,text(note)+'同一队列按1个、2–3个和至少4个界面分组时，分别有414条记录（196个家族）、144条记录（58个家族）和57条记录（34个家族）；同一家族可出现在不同组。')
    replace(para('表S10：'),'BackFlip两个版本均使用完整结构输入的assembly模式。',
        'BackFlip两个版本均使用完整结构输入的assembly模式，曲线对应图2d。')
    drawing=next(x for x in original if x.xpath('.//wp:docPr[@id="2"]',namespaces=NS))
    cx=int(drawing.xpath('.//wp:extent',namespaces=NS)[0].get('cx'))
    with Image.open(OUT/'fig2_grouped.png') as im:cy=round(cx*im.height/im.width)
    for el in drawing.xpath('.//wp:extent | .//a:xfrm/a:ext',namespaces=NS):
        el.set('cx',str(cx));el.set('cy',str(cy))
    drawing.xpath('.//wp:docPr',namespaces=NS)[0].set('descr',
        'a prediction agreement; b record AP distributions with its own color legend; c within-type Hchi PR with internal legend; d paired external RMSF PR and ROC as one subfigure; e original interface-stratified three metrics with centered shared legend. AP/AUROC values are in Table S10.')
    modified=[i for i,x in enumerate(original) if ET.tostring(x)!=before[i]]
    changed={'word/document.xml':ET.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),
             'word/media/image2.png':(OUT/'fig2_grouped.png').read_bytes()}
    with zipfile.ZipFile(TARGET,'w',zipfile.ZIP_DEFLATED) as z:
        for info in infos:z.writestr(info,changed.get(info.filename,parts[info.filename]))
    with zipfile.ZipFile(TARGET) as z:
        assert z.testzip() is None
        actual=[name for name,data in parts.items() if z.read(name)!=data]
        assert set(actual)==set(changed)
        all_text=''.join(ET.fromstring(z.read('word/document.xml')).xpath('//w:t/text()',namespaces=NS))
        assert '图2f' not in all_text and '图2d、e' not in all_text
        assert '图2d、e' not in text(note)
    report={'source':str(SOURCE),'target':str(TARGET),
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'target_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),
        'changed_package_parts':actual,'all_other_parts_byte_identical':True,
        'modified_original_body_indices':modified,'main_caption_characters':len(new_caption),
        'figure_extent_emu':[cx,cy],
        'subfigures':['a','b','c','d:PR+ROC','e:three interface metrics'],
        'table_S10_preserved':True}
    (OUT/'docx_patch_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
