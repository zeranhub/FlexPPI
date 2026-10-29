from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from copy import deepcopy
from lxml import etree
from PIL import Image
import json,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parents[2]
SRC=BASE/'paper_v95_20261026_中文修订_结构集中与局部精简.docx'
DEST=BASE/'paper_v96_20261026_中文修订_图3视觉层级重排.docx'
sha=lambda b:hashlib.sha256(b).hexdigest()
assert sha(SRC.read_bytes())=='070ff7277549ebf2290ced20e86bb9bf056db9dc0cddfe72217a877b42a75067'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
q=lambda tag:'{'+NS['w']+'}'+tag
with ZipFile(SRC) as z:parts={n:z.read(n) for n in z.namelist()}
before=parts.copy();doc=etree.fromstring(parts['word/document.xml']);ps=doc.findall('.//w:body/w:p',NS)
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def settext(p,s):
    runs=p.findall('w:r',NS);prop=deepcopy(runs[0].find('w:rPr',NS)) if runs and runs[0].find('w:rPr',NS) is not None else None
    for child in list(p):
        if child.tag!=q('pPr'):p.remove(child)
    run=etree.SubElement(p,q('r'))
    if prop is not None:run.append(prop)
    etree.SubElement(run,q('t')).text=s
def find(prefix):
    results=[p for p in ps if text(p).startswith(prefix)];assert len(results)==1,(prefix,len(results));return results[0]
edits=[]
def replace(prefix,old,new):
    p=find(prefix);s=text(p);assert old in s,(prefix,old);settext(p,s.replace(old,new));edits.append(dict(paragraph=prefix,old=old,new=new))
replace('CaMKI的调控序列与自抑制区域重叠','（图3a、c）','（图3a、b）')
replace('从结构与预测的对应看','（图3c）','（图3b）')
replace('泛素–UBA复合物2QHO提供','（图3d）','（图3c）')
replace('泛素–UBA复合物2QHO提供','（图3e）','（图3d）')
replace('其中，Q41侧链位于紧密堆积','（图3b、e）','（图3a、d）')
replace('将同一个联合Hχ预测进一步与MD逐χ角','（图3f）','（图3e）')
replace('控制记录内氨基酸类型后','（图3g）','（图3f）')
replace('CaM–CaMKI实例采用DynaRepo记录','图3a、c中的1MXE预测','图3a、b中的1MXE预测')
replace('图S2：泛素–UBA实例','完整66个Hχ位点预测见图3d','完整66个Hχ位点预测见图3c')
replace('Q41/Q2按联合Hχ预测极值选出','逐χ角的模型–MD关联已恢复至图3f','逐χ角的模型–MD关联见图3e')
caption='图3：主链约束、同类型侧链环境差异与角度占据关联。（a）CaM–CaMKI（1MXE，左）与泛素–UBA（2QHO，右）整体结构；紫色为被分析链，灰色为伙伴，青绿色标出W303/M316或Q41/Q62。（b）CaMKI全部25个位点的RMSF，浅区为W303–M316核心；黑色空心点线为MD，蓝色实线为Full，暖赭虚线为Intra，该MD/Full配色也用于c、d。（c）泛素全部66个有效Hχ位点的MD与Full预测，无χ角位置留空。（d）全部六个Gln的SASA与Hχ；空心圆为MD，蓝色实心圆为Full，连线对应同一位置，仅标注Q41/Q62，不作SASA回归。（e）同一联合Hχ预测与MD各χ角5°细箱熵的关联，两轴分别按氨基酸类型去均值；Q41/Q2为六个Gln中预测最低/最高位点，以菱形/方形标出，保留全部有效点。（f）83条记录中Full Hχ与Full RMSF对应MD细箱熵及低F*角区宽度的同类型相关；浅点保留全部记录，大点及误差线为记录均值与来源分层重抽10000次的95%区间。模型不生成构象或直接预测结合自由能。'
settext(find('图3：主链约束'),caption)
im=Image.open(ROOT/'figure3/fig3_dual_case_entropy.png');cx=164*36000;cy=round(cx*im.height/im.width)
nodes=doc.xpath('.//w:p[.//a:blip[@r:embed="rId13"]]',namespaces=NS);assert len(nodes)==1
for node in nodes[0].xpath('.//wp:extent|.//a:xfrm/a:ext',namespaces=NS):node.set('cx',str(cx));node.set('cy',str(cy))
parts['word/media/image5.png']=(ROOT/'figure3/fig3_dual_case_entropy.png').read_bytes()
parts['word/document.xml']=etree.tostring(doc,encoding='UTF-8',xml_declaration=True,standalone=True)
changed=[n for n in parts if parts[n]!=before[n]];assert set(changed)=={'word/document.xml','word/media/image5.png'}
with ZipFile(DEST,'w',ZIP_DEFLATED) as z:
    for n,data in parts.items():z.writestr(n,data)
with ZipFile(DEST) as z:assert z.testzip() is None
for ext in ['png','pdf','svg']:
    (BASE/f'figure3_v96_20261026_视觉层级重排.{ext}').write_bytes((ROOT/f'figure3/fig3_dual_case_entropy.{ext}').read_bytes())
audit=dict(source=str(SRC),source_sha256=sha(SRC.read_bytes()),output=str(DEST),output_sha256=sha(DEST.read_bytes()),changed_zip_parts=changed,all_other45_media_unchanged=True,caption_characters=len(caption),image_width_mm=cx/36000,image_height_mm=cy/36000,edits=edits)
(ROOT/'manuscript/build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'manuscript/图注与引用调整.md').write_text(caption+'\n\n'+json.dumps(edits,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
