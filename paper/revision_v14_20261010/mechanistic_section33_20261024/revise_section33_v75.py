"""Add evidence-backed biological context and explicitly labelled local hypotheses."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import csv
import hashlib
import json
import shutil
import zipfile
from lxml import etree as E

ROOT=Path(__file__).resolve().parent
PAPER=ROOT.parent
SOURCE=PAPER/'paper_v74_20261024_中文修订_补充材料精简合并.docx'
TARGET=PAPER/'paper_v75_20261024_中文修订_3.3背景与机制推断.docx'
ASSETS=Path(r'E:\Project\FlexPPI\paper\figures_v59\03_KIX与泛素局部动态验证\source_archive')
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
W='{'+NS['w']+'}'

def txt(el):return ''.join(el.xpath('.//w:t/text()',namespaces=NS))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def plain(template,value):
    el=deepcopy(template)
    for key in list(el.attrib):
        if key.endswith(('paraId','textId')):del el.attrib[key]
    rp=el.find('w:r/w:rPr',NS)
    for child in list(el):
        if child.tag!=W+'pPr':el.remove(child)
    run=E.SubElement(el,W+'r')
    if rp is not None:run.append(deepcopy(rp))
    node=E.SubElement(run,W+'t')
    node.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    node.text=value
    return el

def main():
    assert sha(SOURCE)=='ccdd897c7128867076b61b3f1ca00e6db60aeffe06e50904c4c945db6eb81680'
    if TARGET.exists():
        assert sha(TARGET)=='8b32a97c44e04236a3afb77aa4a45c786b180a48c71fba41e3acf2e4aa620b36'
        backup=ROOT/'v75_first_render_draft.docx'
        assert not backup.exists()
        shutil.copy2(TARGET,backup)
    section=json.loads((ROOT/'section33_v75.json').read_text(encoding='utf-8'))
    with zipfile.ZipFile(SOURCE) as package:
        entries=package.infolist()
        parts={i.filename:package.read(i.filename) for i in entries}
    root=E.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
    old=list(body)
    source_tables=[E.tostring(el) for el in root.xpath('//w:tbl',namespaces=NS)]
    source_math=Counter(E.tostring(el) for el in root.xpath('//m:oMath',namespaces=NS))
    source_images=Counter(E.tostring(el) for el in root.xpath('//w:drawing',namespaces=NS))
    start=next(i for i,el in enumerate(old) if txt(el).startswith('3.3 '))
    end=next(i for i,el in enumerate(old) if txt(el).startswith('3.4 '))
    heading,template=old[start:start+2]
    caption=next(el for el in old[start:end] if txt(el).startswith('图3：'))
    drawing=caption.getprevious()
    assert drawing.xpath('.//w:drawing',namespaces=NS)
    cap_value=txt(caption).replace('图3：静态结构的局部动态验证。','图3：复合物主链活动与侧链构象约束。')
    cap_value+='Q41/Q62的起始结构接触与溶剂可及面积见补充材料D。'
    new=[plain(heading,'3.3 '+section['title'])]
    # The combined figure sits between the two case analyses, so expanding the
    # biological context does not strand a short text fragment before a full figure.
    for i,value in enumerate(section['paragraphs']):
        new.append(plain(template,value))
        if i==1:new.extend([deepcopy(drawing),plain(caption,cap_value)])
    for el in old[start:end]:body.remove(el)
    for offset,el in enumerate(new):body.insert(start+offset,el)
    # The next full-page figure is interleaved with its case discussion so the
    # expanded preceding section does not leave two short paragraphs alone.
    h34=next(el for el in body if txt(el).startswith('3.4 '))
    h35=next(el for el in body if txt(el).startswith('3.5 '))
    case34=list(body)[list(body).index(h34)+1:list(body).index(h35)]
    cap4=next(el for el in case34 if txt(el).startswith('图4：'))
    image4=cap4.getprevious()
    paragraphs34=[el for el in case34 if el.tag==W+'p' and txt(el) and el is not cap4]
    assert len(paragraphs34)==4
    body.remove(image4);body.remove(cap4)
    image_slot=list(body).index(paragraphs34[1])+1
    body.insert(image_slot,image4);body.insert(image_slot+1,cap4)
    sasa_old=next(el for el in body if txt(el).startswith('泛素侧链溶剂可及面积由'))
    sasa_method='泛素侧链的起始结构溶剂可及面积采用PyMOL 3.0.0的get_area计算，以完整复合物的935个重原子为环境，设dot_solvent=1、探针半径1.4 Å、dot_density=4，累加各残基除N、Cα、C、O以外的侧链重原子面积，作为预测后的结构环境解释量。局部邻域记录目标侧链与其他残基重原子的最近距离不超过5 Å的位置。侧链取向叠加使用实际MD代表帧，经局部主链对齐展示；逐残基预测、联合状态频率、接触与面积源表及代表帧来源随配套数据提供。'
    body.replace(sasa_old,plain(sasa_old,sasa_method))
    appendix_e=next(el for el in body if txt(el).startswith('E 联合柔性与构象占据验证'))
    appendix_note='泛素–UBA的结构环境分析采用与图3一致的A003R完整双链起始坐标。Q41与Q62的侧链溶剂可及面积分别为0.0和100.67 Å²，描述起始结构的暴露环境；邻域接触和计算设置见配套源表及2.7节。五个连续100 ns窗口中，Q41与Q62的Hχ范围分别为1.829–1.977与2.488–2.774 nats，各窗口均为Q62较高；这些窗口来自同一条500 ns轨迹。模型与MD的逐位置读数、27类联合状态占据及窗口统计随源数据提供。'
    body.insert(list(body).index(appendix_e),plain(template,appendix_note))
    # All previously frozen data, figures and native equations remain intact.
    assert [E.tostring(el) for el in root.xpath('//w:tbl',namespaces=NS)]==source_tables
    assert Counter(E.tostring(el) for el in root.xpath('//m:oMath',namespaces=NS))==source_math
    assert Counter(E.tostring(el) for el in root.xpath('//w:drawing',namespaces=NS))==source_images
    new_text='\n'.join(section['paragraphs'])
    assert not any(s in new_text for s in ['遗憾','仅','明显落后','效果有限','证明变构','重现结合诱导'])
    for el in old:
        if el in old[start:end] or el is sasa_old or el.tag!=W+'p':continue
        assert any(E.tostring(el)==E.tostring(candidate) for candidate in body)
    id_name='{http://schemas.microsoft.com/office/word/2010/wordml}paraId'
    ids=[el.get(id_name) for el in root.iter() if el.get(id_name)]
    assert len(ids)==len(set(ids))
    # Archive the existing evidence, not newly simulated or refitted measurements.
    evidence=ROOT/'source_data';evidence.mkdir(exist_ok=True)
    source_paths=[ASSETS/'ubiquitin/reference_environment_contacts.csv',
                  ASSETS/'ubiquitin/reference_solvent_accessibility.csv',
                  ASSETS/'ubiquitin/render_from_bundle.py',
                  ASSETS/'data/ubiquitin_Q41_Q62_validation_summary.csv',
                  ASSETS/'data/A003R_Q41_Q62_label_windows.csv']
    inputs=[]
    for source in source_paths:
        target=evidence/source.name;shutil.copy2(source,target)
        assert sha(source)==sha(target)
        inputs.append({'source':str(source),'archive':str(target),'sha256':sha(source)})
    sasa=list(csv.DictReader(source_paths[1].open(encoding='utf-8-sig')))
    assert [round(float(r['sidechain_SASA_A2']),2) for r in sasa]==[0.0,100.67]
    window=list(csv.DictReader(source_paths[-1].open(encoding='utf-8-sig')))
    (evidence/'README.md').write_text('# 3.3结构机制证据\n\n本目录复用既有结构接触、溶剂可及面积、MD窗口及留出预测源表；未运行新模拟或新模型。溶剂可及面积对应完整双链起始结构，五个连续100 ns窗口来自同一500 ns轨迹。正文的邻域扰动描述为待检验的机制假说。\n\n背景沿用原文献[22] De Guzman等2006年JMB（https://pubmed.ncbi.nlm.nih.gov/16253272/）与[23] Kozlov等2007年JBC（https://pubmed.ncbi.nlm.nih.gov/17897937/）。\n',encoding='utf-8')
    xml=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
    with zipfile.ZipFile(TARGET,'w',zipfile.ZIP_DEFLATED) as package:
        for entry in entries:package.writestr(entry,xml if entry.filename=='word/document.xml' else parts[entry.filename])
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        changed=[name for name,value in parts.items() if package.read(name)!=value]
        assert changed==['word/document.xml']
    audit={'source':str(SOURCE),'target':str(TARGET),'target_sha256':sha(TARGET),
           'changed_package_parts':changed,'main_revision':'3.3 biological context and explicitly labelled hypotheses',
           'supporting_revisions':['2.7 existing SASA method expanded','Figure 3 caption and placement','Figure 4 placement for pagination only','Supplement D structural provenance'],
           'all_tables_native_math_and_images_preserved_exactly':True,
           'new_model_predictions_or_simulations':False,'source_evidence':inputs,
           'literature':[{'citation':22,'url':'https://pubmed.ncbi.nlm.nih.gov/16253272/','role':'published cooperative-binding background'},
                         {'citation':23,'url':'https://pubmed.ncbi.nlm.nih.gov/17897937/','role':'published EDD UBA recognition background'}],
           'causal_scope':'L12 adaptation and neighborhood perturbations explicitly framed as hypotheses',
           'visual_qa':'pending'}
    (ROOT/'v75_patch_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k not in ['source_evidence','literature']},ensure_ascii=False))

if __name__=='__main__':main()
