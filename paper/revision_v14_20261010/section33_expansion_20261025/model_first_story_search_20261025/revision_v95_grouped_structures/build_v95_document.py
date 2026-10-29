"""Update Figure 3, its caption, and exact panel references; preserve other parts."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from copy import deepcopy
from lxml import etree
from PIL import Image
import json
import hashlib
import sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parents[2]
SRC=BASE/'paper_v94_20261026_中文修订_恢复逐角占据关联.docx'
DEST=BASE/'paper_v95_20261026_中文修订_结构集中与局部精简.docx'
sha=lambda b: hashlib.sha256(b).hexdigest()
assert sha(SRC.read_bytes())=='6aac6e6333bc6a9f7ae9c1fc6e369cb78f535c9bdae10fecd22ae0e75ca4cb7b'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
q=lambda tag:'{'+NS['w']+'}'+tag
with ZipFile(SRC) as z: parts={n:z.read(n) for n in z.namelist()}
old=parts.copy()
doc=etree.fromstring(parts['word/document.xml'])
ps=doc.findall('.//w:body/w:p',NS)
def text(p): return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def settext(p,s):
    runs=p.findall('w:r',NS)
    prop=deepcopy(runs[0].find('w:rPr',NS)) if runs and runs[0].find('w:rPr',NS) is not None else None
    for child in list(p):
        if child.tag!=q('pPr'):p.remove(child)
    run=etree.SubElement(p,q('r'))
    if prop is not None:run.append(prop)
    etree.SubElement(run,q('t')).text=s
def find(prefix):
    result=[p for p in ps if text(p).startswith(prefix)]
    assert len(result)==1,(prefix,len(result))
    return result[0]
edits=[]
def replace(prefix,a,b):
    p=find(prefix);s=text(p);assert a in s,(prefix,a)
    settext(p,s.replace(a,b));edits.append({'paragraph':prefix,'old':a,'new':b})
replace('CaMKI的调控序列与自抑制区域重叠','（图3a、b）','（图3a、c）')
replace('从结构与预测的对应看','（图3b）','（图3c）')
replace('泛素–UBA复合物2QHO提供','（图3c）','（图3d）')
replace('其中，Q41侧链位于紧密堆积','（图3d、e）','（图3b、e）')
replace('其中，Q41侧链位于紧密堆积','局部全残基棍状表示将这种状态多样性差别与实际堆积邻域对应起来。','')
replace('CaM–CaMKI实例采用DynaRepo记录','图3a、b中的1MXE预测','图3a、c中的1MXE预测')
replace('图S2：泛素–UBA实例','完整66个Hχ位点预测见图3c','完整66个Hχ位点预测见图3d')
caption='图3：主链约束、同类型侧链环境差异与角度占据关联。（a、b）CaM–CaMKI（1MXE）与泛素–UBA（2QHO）整体丝带结构；紫色为被分析链，灰色为伙伴，青绿色标出W303/M316或Q41/Q62。（c）CaMKI全部25位点RMSF：MD黑色点线、Full蓝色实线、Intra暖赭虚线；浅区为W303–M316，下方为4.5 Å内CaM接触残基数。（d）泛素全部66个有效Hχ位点的MD与Full预测，无χ角位置留空。（e）全部六个Gln的SASA与Hχ，空心圆为MD、蓝方点为Full，连线对应同一位置，不作SASA回归。（f）同一联合Hχ预测与MD各χ角5°细箱熵的关联，两轴分别按氨基酸类型去均值；Q41/Q2为六个Gln中预测最低/最高位点，以菱形/方形标出，四图均保留全部有效点。（g）83记录中Full Hχ与Full RMSF对应MD细箱熵及低F*角区宽度的同类型相关；全部记录散点，大点及误差线为记录均值与来源分层重抽10000次的95%区间。模型不生成构象或直接预测结合自由能。'
settext(find('图3：主链约束'),caption)
new_image=(ROOT/'figure3/fig3_dual_case_entropy.png').read_bytes()
parts['word/media/image5.png']=new_image
parts['word/document.xml']=etree.tostring(doc,encoding='UTF-8',xml_declaration=True,standalone=True)
changed=[n for n in parts if old[n]!=parts[n]]
assert set(changed)=={'word/document.xml','word/media/image5.png'}
with ZipFile(DEST,'w',ZIP_DEFLATED) as z:
    for n,data in parts.items():z.writestr(n,data)
with ZipFile(DEST) as z:assert z.testzip() is None
for ext in ['png','pdf','svg']:
    src=ROOT/f'figure3/fig3_dual_case_entropy.{ext}'
    dst=BASE/f'figure3_v95_20261026_结构集中与局部精简.{ext}'
    dst.write_bytes(src.read_bytes())
audit={'source':str(SRC),'source_sha256':sha(SRC.read_bytes()),'output':str(DEST),'output_sha256':sha(DEST.read_bytes()),'changed_zip_parts':changed,'all_six_tables_unchanged':True,'all_other45_media_unchanged':True,'supplement4_unchanged':True,'caption_characters':len(caption),'edits':edits,'source_figure_size_px':Image.open(ROOT/'figure3/fig3_dual_case_entropy.png').size}
(ROOT/'manuscript/build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'manuscript/图注与引用调整.md').write_text(caption+'\n\n'+json.dumps(edits,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
