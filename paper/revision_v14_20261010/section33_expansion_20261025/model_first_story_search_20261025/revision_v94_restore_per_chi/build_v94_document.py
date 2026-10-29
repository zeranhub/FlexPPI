from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from copy import deepcopy
from lxml import etree
from PIL import Image
import json,hashlib
R=Path(__file__).resolve().parent;BASE=R.parents[2];O=R/'manuscript';O.mkdir(exist_ok=True)
SRC=BASE/'paper_v93_20261026_中文修订_双案例与构象占据整合.docx'
DEST=BASE/'paper_v94_20261026_中文修订_恢复逐角占据关联.docx'
sha=lambda b:hashlib.sha256(b).hexdigest()
assert sha(SRC.read_bytes())=='6c603bfff03d033610d228854c05c0cb65bce4ceee3e9f9dae2bf8b3433aab44'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
q=lambda t:'{'+NS['w']+'}'+t
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def settext(p,s):
    runs=p.findall('w:r',NS);prop=deepcopy(runs[0].find('w:rPr',NS)) if runs and runs[0].find('w:rPr',NS) is not None else None
    for child in list(p):
        if child.tag!=q('pPr'):p.remove(child)
    run=etree.SubElement(p,q('r'))
    if prop is not None:run.append(prop)
    etree.SubElement(run,q('t')).text=s
with ZipFile(SRC) as z:parts={i.filename:z.read(i.filename) for i in z.infolist()}
original=parts.copy();doc=etree.fromstring(parts['word/document.xml']);body=doc.find('w:body',NS)
before=[etree.tostring(t,method='c14n') for t in body.findall('w:tbl',NS)]
ps=body.findall('w:p',NS)
def find(prefix):
    result=[p for p in ps if text(p).startswith(prefix)];assert len(result)==1,(prefix,len(result));return result[0]
caption='图3：主链约束、同类型侧链环境差异与角度占据关联。（a）1MXE中CaM–CaMKI丝带结构及W303/M316，紫色为被分析链，灰色为伙伴。（b）全部25位点RMSF：MD黑色点线、Full蓝色实线、Intra暖赭虚线；浅区为W303–M316，下方为4.5 Å内CaM接触残基数。（c）泛素全部66个有效Hχ位点的MD与Full预测，无χ角位置留空。（d）泛素–UBA整体与Q41/Q62的5 Å完整残基邻域棍状放大，目标为青绿色，两局部同视角同比例。（e）全部六个Gln的SASA与Hχ，空心圆为MD、蓝方点为Full，连线对应同一位置，不作SASA回归。（f）同一联合Hχ预测与MD各χ角5°细箱熵的关联，两轴分别按氨基酸类型去均值；Q41/Q2为六个Gln中预测最低/最高位点，以菱形/方形标出，四图均保留全部有效点。（g）83记录中Full Hχ与Full RMSF对应MD细箱熵及低F*角区宽度的同类型相关；全部记录散点，大点及误差线为记录均值与来源分层重抽10000次的95%区间。模型不生成构象或直接预测结合自由能。'
settext(find('图3：主链约束'),caption)
paragraph='将同一个联合Hχ预测进一步与MD逐χ角的5°细箱熵对应，并对两轴分别去除氨基酸类型均值，可检验联合输出是否保留了单角层面的占据差异（图3f）。χ1–χ3呈中等正相关（ρ=0.50、0.45、0.52，n=66、52、23），χ4的样本较少且关联较弱（ρ=0.13，n=10）。四幅图共享同一个联合Hχ预测，并不表示模型分别预测四个χ角。标出的Q41和Q2为六个Gln中按联合Hχ预测选出的最低和最高位点，与前述Q41/Q62的结构环境对照不同；这一对应将联合Hχ预测与实际单角占据联系起来，支持其作为侧链构象多样性摘要的解释。'
clone=deepcopy(find('泛素–UBA复合物2QHO提供'));settext(clone,paragraph)
head=find('3.4 侧链预测');body.insert(list(body).index(head),clone)
p=find('控制记录内氨基酸类型后');old=text(p);assert '（图3f）' in old;settext(p,old.replace('（图3f）','（图3g）'))
p=find('进一步仅按预测值在每条记录');old=text(p);assert '（图3g）' in old
new=old.replace('（图3g）','（图S4c）').replace('泛素逐角分布的局部核查见图S4c、d','泛素逐角MD分布的局部核查见图S4d')
settext(p,new)
s4='图S4：柔性读数与角度占据的补充验证。（a）3063条记录中MD参考与模型预测的RMSF–Hχ同类型相关经验累积分布，每条曲线保留全部3063个顶点；它描述两种读数的相关结构，不是条件独立或预测准确性的证据。（b）83条记录中，Full Hχ、Intra Hχ及Full RMSF与MD细箱熵和低F*角区宽度的同类型相关，全部83记录成对显示；大点为均值，误差线为来源分层重抽记录10000次的95%区间，种子20261026。（c）在同类型内仅按Full Hχ预测划分的低/高两半，其实际MD细箱熵和低F*角区宽度，保留全部83条记录的配对值及同样定义的均值区间。（d）六个Gln中按联合Hχ预测最低的Q41与最高的Q2所对应的MD逐角F*曲线，预测分别为2.000和2.990 nats；窗口为0–100 ns、501帧、原点0°的5°箱。全部432角度箱保留，211个零计数箱断线，不加伪计数或平滑；虚线F*=1。曲线由MD占据统计得到，不是模型输出，也不换算为结合自由能。'
settext(find('图S4：双柔性'),s4)
p=find('Q41/Q2是按六个Gln预测极端选出');settext(p,'Q41/Q2按联合Hχ预测极值选出，与主图Q41/Q62的环境对照具有不同选择依据；逐χ角的模型–MD关联已恢复至图3f。b、c采用原体系划分的折外预测，d采用原泛素预测与已有MD。残基、角度箱及帧不作为独立生物重复。')
images=[]
for part,path in [('word/media/image5.png',R/'figure3/fig3_dual_case_entropy.png'),('word/media/image40.png',R/'supplement4/figS4_joint_readouts_angular_occupancy.png')]:
    rid='rId13' if 'image5.' in part else 'rId51';nodes=body.xpath('.//w:p[.//a:blip[@r:embed="'+rid+'"]]',namespaces=NS);assert len(nodes)==1
    im=Image.open(path);cx=164*36000;cy=round(cx*im.height/im.width)
    for node in nodes[0].xpath('.//wp:extent|.//a:xfrm/a:ext',namespaces=NS):node.set('cx',str(cx));node.set('cy',str(cy))
    parts[part]=path.read_bytes();images.append(dict(part=part,png=str(path),sha256=sha(parts[part]),width_mm=164,height_mm=cy/36000))
assert before==[etree.tostring(t,method='c14n') for t in body.findall('w:tbl',NS)]
parts['word/document.xml']=etree.tostring(doc,encoding='UTF-8',xml_declaration=True,standalone=True)
changed=[name for name,data in parts.items() if original[name]!=data]
assert set(changed)=={'word/document.xml','word/media/image5.png','word/media/image40.png'}
with ZipFile(DEST,'w',ZIP_DEFLATED) as z:
    for name,data in parts.items():z.writestr(name,data)
with ZipFile(DEST) as z:assert z.testzip() is None
audit=dict(source=str(SRC),source_sha256=sha(SRC.read_bytes()),output=str(DEST),output_sha256=sha(DEST.read_bytes()),changed_zip_parts=changed,all_six_tables_unchanged=True,all_other44_media_unchanged=True,images=images,restored_main_panel='f, per-chi 2x2',previous_selection_plot='Supplementary S4c',added_text=paragraph,main_caption_chars=len(caption))
(O/'build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'图注与新增段落.md').write_text(caption+'\n\n'+paragraph+'\n\n'+s4,encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
