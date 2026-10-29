"""Keep section 2.7 procedural and remove duplicated case-study overviews."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib
import json
import zipfile
from lxml import etree as E

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'paper_v75_20261024_中文修订_3.3背景与机制推断.docx'
TARGET=ROOT.parent/'paper_v76_20261024_中文修订_方法去除实例概述.docx'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
W='{'+NS['w']+'}'

def txt(el):return ''.join(el.xpath('.//w:t/text()',namespaces=NS))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def plain(template,value):
    el=deepcopy(template)
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
    assert sha(SOURCE)=='74b1f7903a601038b1e80a4baa2f74c40b0774c12a990d390ccdead8ca9bcbf1'
    assert not TARGET.exists()
    with zipfile.ZipFile(SOURCE) as package:
        entries=package.infolist();parts={entry.filename:package.read(entry.filename) for entry in entries}
    root=E.fromstring(parts['word/document.xml']);body=root.find('w:body',NS)
    originals=list(body)
    def find(prefix):
        matches=[el for el in originals if el.tag==W+'p' and txt(el).startswith(prefix)]
        assert len(matches)==1,(prefix,len(matches))
        return matches[0]
    title=find('2.7 结构实例与分子动力学验证')
    kix=find('结构实例考察主链活动区域定位')
    tcr=find('TCR实例采用ITrace')
    ub=find('侧链实例使用泛素与EDD')
    start=originals.index(title);end=originals.index(find('2.8 联合预测与构象占据的验证'))
    head=[E.tostring(el) for el in originals[:start]]
    tail=[E.tostring(el) for el in originals[end:]]
    native_math=Counter(E.tostring(el) for el in root.xpath('//m:oMath',namespaces=NS))
    source_tables=[E.tostring(el) for el in root.xpath('//w:tbl',namespaces=NS)]
    removed_intro='结构实例考察主链活动区域定位、受体亲和成熟后的跨界面动态响应，以及同类型残基侧链朝向多样性的区分。为检验静态结构预测是否能够定位明确的主链活动区域，我们对CBP的KIX结构域开展专门分子动力学模拟。'
    assert txt(kix).startswith(removed_intro)
    revised={title:'2.7 结构体系构建与轨迹分析',
             kix:txt(kix)[len(removed_intro):],
             tcr:txt(tcr).replace('TCR实例采用ITrace中的','TCR–pMHC分析采用ITrace中的',1),
             ub:'泛素–UBA分析使用DynaRepo记录A003R，其起始复合物对应PDB 2QHO[23]；输入完整两链结构，采用第0折留出预测及主评价标签。前5,000帧分为五个连续1,000帧窗口，联合状态熵按窗口计算后取均值；状态占据率与累积覆盖合并全部帧。谷氨酰胺三个χ角各分为三个等宽区间，共27类联合状态。按各位点联合状态占据率选取最高三类的实际MD代表帧，以该位点9 Å内的局部主链对齐。联合Hχ预测与MD单角细箱熵按χ1–χ4分别配对；在每个χ角层内，分别减去同一种氨基酸在预测值与MD指标中的平均值，再计算Spearman相关。逐角比较在泛素六个谷氨酰胺中按留出模型Hχ预测选取最低和最高两个位置；全部逐角指标统一取0–100 ns、200 ps间隔的501帧，并按2.8节对四个箱格原点取平均。'}
    changes=[]
    for original,value in revised.items():
        assert not original.xpath('.//m:oMath|.//w:drawing',namespaces=NS)
        changes.append({'body_index':originals.index(original),'before':txt(original),'after':value})
        body.replace(original,plain(original,value))
    assert [E.tostring(el) for el in list(body)[:start]]==head
    assert [E.tostring(el) for el in list(body)[end:]]==tail
    assert Counter(E.tostring(el) for el in root.xpath('//m:oMath',namespaces=NS))==native_math
    assert [E.tostring(el) for el in root.xpath('//w:tbl',namespaces=NS)]==source_tables
    for el in list(body)[start:end]:
        if el.tag==W+'p':assert '实例' not in txt(el)
    xml=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
    with zipfile.ZipFile(TARGET,'w',zipfile.ZIP_DEFLATED) as package:
        for entry in entries:package.writestr(entry,xml if entry.filename=='word/document.xml' else parts[entry.filename])
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        assert [name for name,value in parts.items() if package.read(name)!=value]==['word/document.xml']
    audit={'source':str(SOURCE),'target':str(TARGET),'target_sha256':sha(TARGET),
           'only_section_2_7_changed':True,'changed_paragraphs':changes,
           'all_results_supplement_tables_figures_and_native_math_preserved_exactly':True,
           'removed_content_role':['duplicated biological-purpose overview','figure walk-through text already described in captions/results'],
           'retained_content_role':['PDB and record identifiers','input chain selection and conditions','simulation parameters','held-out model assignments','trajectory windows','mapping and scoring rules','representative-frame and site selection','SASA and contact algorithms'],
           'visual_qa':'pending'}
    (ROOT/'v76_patch_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k!='changed_paragraphs'},ensure_ascii=False))

if __name__=='__main__':main()
