"""Reorganize frozen supplementary evidence and remap all manuscript references."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib
import json
import re
import shutil
import zipfile
from lxml import etree as E
from PIL import Image

ROOT = Path(__file__).resolve().parent
PAPER = ROOT.parent
SOURCE = PAPER / 'paper_v73_20261024_中文修订_AA至Full分层贡献.docx'
TARGET = PAPER / 'paper_v74_20261024_中文修订_补充材料精简合并.docx'
NS = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
W = '{'+NS['w']+'}'
REL = '{http://schemas.openxmlformats.org/package/2006/relationships}'
TABLE_MAP = {'1':'1','2':'2','3':'3','4':'6','5':'4','6':'9','7':'10',
             '8':'8','9':'11','10':'7','11A':'5A','11B':'7','11C':'5B','11':'5'}

def txt(el):
    return ''.join(el.xpath('.//w:t/text()', namespaces=NS))

def setprop(parent, tag, **attrs):
    node = parent.find('w:'+tag, NS)
    if node is None:
        node = E.SubElement(parent, W+tag)
    for k,v in attrs.items():
        node.set(W+k,str(v))
    return node

def ppr(el):
    node = el.find('w:pPr', NS)
    if node is None:
        node = E.Element(W+'pPr')
        el.insert(0,node)
    return node

def plain(template, value):
    el=deepcopy(template)
    for key in list(el.attrib):
        if key.endswith(('paraId','textId')):
            del el.attrib[key]
    fmt=el.find('w:r/w:rPr', NS)
    for child in list(el):
        if child.tag != W+'pPr':
            el.remove(child)
    run=E.SubElement(el,W+'r')
    if fmt is not None:
        run.append(deepcopy(fmt))
    node=E.SubElement(run,W+'t')
    node.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    node.text=value
    return el

def replace_spans(el, regex, callback):
    nodes=el.xpath('.//w:t',namespaces=NS)
    offsets=[]; total=0
    for node in nodes:
        value=node.text or ''
        offsets.append((total,total+len(value),node))
        total+=len(value)
    value=''.join(n.text or '' for n in nodes)
    for match in reversed(list(regex.finditer(value))):
        a,b=match.span(); replacement=callback(match)
        affected=[(x,y,n) for x,y,n in offsets if x<b and y>a]
        if not affected:continue
        first=affected[0]; last=affected[-1]
        if first[2] is last[2]:
            n=first[2]; old=n.text or ''
            n.text=old[:a-first[0]]+replacement+old[b-first[0]:]
        else:
            n=first[2]; n.text=(n.text or '')[:a-first[0]]+replacement
            for _,_,n in affected[1:-1]:n.text=''
            n=last[2]; n.text=(n.text or '')[b-last[0]:]

def remap(el):
    replace_spans(el,re.compile(r'图\s*S([23])(?!\d)'),
                  lambda m:'图S1a' if m[1]=='2' else '图S1b、c')
    replace_spans(el,re.compile(r'S(11[A-C]?|10|[1-9])(?![0-9A-Za-z])'),
                  lambda m:'S'+TABLE_MAP[m[1]])
    for pattern,value in [(r'表S7及S7','表S7'),(r'表S7、S7','表S7')]:
        replace_spans(el,re.compile(pattern),lambda m,v=value:v)
    letters={'D':'B','E':'E','F':'D','G':'E','H':'C','I':'C'}
    replace_spans(el,re.compile(r'补充材料\s*([D-I])'),lambda m:'补充材料'+letters[m[1]])
    return el

def rows(table):
    return [[txt(c) for c in r.findall('w:tc',NS)] for r in table.findall('w:tr',NS)]

def main():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='e19aac6b5b5c828d00afb351d0f3af9407f716186e5c0bded7c48ec4906dad0d'
    if TARGET.exists():
        assert hashlib.sha256(TARGET.read_bytes()).hexdigest()=='7195bf31068735b6fa8c9932268df7665a292a78946478b62b7151ff2f4aa1f8'
        backup=ROOT/'v74_first_render_draft.docx'
        assert not backup.exists()
        shutil.copy2(TARGET,backup)
    merged_image=ROOT/'figS1_merged.png'
    assert merged_image.exists()
    with zipfile.ZipFile(SOURCE) as package:
        entries=package.infolist()
        parts={item.filename:package.read(item.filename) for item in entries}
    root=E.fromstring(parts['word/document.xml'])
    body=root.find('w:body',NS)
    old=list(body)
    supplement_start=next(i for i,el in enumerate(old) if txt(el)=='补充材料')
    assert supplement_start==139 and len(old)==200
    all_math=Counter(E.tostring(el) for el in root.xpath('//m:oMath',namespaces=NS))
    source_main_images=[E.tostring(el) for el in old[:supplement_start] if el.xpath('.//w:drawing',namespaces=NS)]
    # Preserve the original supplementary diagnostics and full implementation text.
    archive=ROOT/'source_data'
    archive.mkdir(exist_ok=True)
    (archive/'supplement_v73_original.xml').write_bytes(E.tostring(E.Element('original')))
    snapshot={'source':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'supplement_original_text':[txt(el) for el in old[supplement_start:] if el.tag==W+'p'],
              'tables_original':{str(i):rows(old[i]) for i in [143,150,152,162,172,176,179,184,189,194,197]}}
    (archive/'supplement_v73_complete_source.json').write_text(json.dumps(snapshot,ensure_ascii=False,indent=2),encoding='utf-8')
    original_tree=E.Element('originalSupplement')
    for el in old[supplement_start:]:original_tree.append(deepcopy(el))
    (archive/'supplement_v73_original.xml').write_bytes(E.tostring(original_tree,xml_declaration=True,encoding='UTF-8'))
    rels=E.fromstring(parts['word/_rels/document.xml.rels'])
    rel_targets={el.get('Id'):el.get('Target') for el in rels}
    diagnostic_rid=old[144].xpath('.//a:blip/@r:embed',namespaces=NS)[0]
    diagnostic_part='word/'+rel_targets[diagnostic_rid]
    (archive/'label_observation_diagnostics.png').write_bytes(parts[diagnostic_part])
    (archive/'README.md').write_text(
        '# 配套源数据说明\n\n标签观测诊断图已从论文补充图中移出，保存在label_observation_diagnostics.png。'
        '该图沿用原图S1，数据与图形均未修改。\n\n'
        'supplement_v73_complete_source.json及supplement_v73_original.xml保存精简前的全部补充文字和表格，'
        '其中包括完整实现设置、来源分层汇总、KIX状态差值和双读数组合统计。'
        '原始比较的样本与推理口径继续分别报告，未按结果优劣筛选。\n',encoding='utf-8')

    body_template=old[191]
    caption_template=old[142]
    heading_template=old[140]

    def para(value, caption=False):
        el=plain(caption_template if caption else body_template,value)
        pr=ppr(el)
        setprop(pr,'pageBreakBefore',val=0)
        setprop(pr,'keepNext',val=1 if caption else 0)
        setprop(pr,'keepLines',val=1 if caption else 0)
        setprop(pr,'spacing',before=80 if caption else 0,after=60 if caption else 100,
                line=260 if caption else 300,lineRule='auto')
        return el

    def heading(value):
        el=plain(heading_template,value)
        pr=ppr(el);setprop(pr,'pageBreakBefore',val=0);setprop(pr,'keepNext',val=1)
        setprop(pr,'spacing',before=180,after=100)
        return el

    def original(i):
        el=remap(deepcopy(old[i]))
        if el.tag==W+'p':
            pr=ppr(el);setprop(pr,'pageBreakBefore',val=0)
            if txt(el).startswith(('表S','图S')):
                setprop(pr,'spacing',before=80,after=60,line=260,lineRule='auto')
            else:setprop(pr,'spacing',after=100,line=300,lineRule='auto')
        return el

    # Split the old dual-task table into the unchanged A and C data panels.
    old_table=old[197]
    source_rows=old_table.findall('w:tr',NS)
    table5=deepcopy(old_table)
    for r in table5.findall('w:tr',NS):table5.remove(r)
    for i in [0,1,2,3,4,5,6,15,16,17,18,19]:table5.append(deepcopy(source_rows[i]))
    panel_b=table5.findall('w:tr',NS)[7].find('w:tc',NS)
    cell_template=panel_b.find('w:p',NS)
    for child in list(panel_b):
        if child.tag != W+'tcPr':panel_b.remove(child)
    panel_b.append(plain(cell_template,'B  ThermoPCD：外部验证（32条记录）'))
    # Merge the same seven 615-complex rows without recomputing displayed values.
    old_ap=rows(old[194])
    old_values=rows(old[197])
    b_values=old_values[8:15]
    assert [r[0] for r in b_values]==[r[0] for r in old_ap[1:]]
    joined=[['方法','MAE (Å) ↓','ρ ↑','AP ↑','AUROC ↑']]
    joined += [[b[0],b[1],b[2],ap[1],ap[2]] for b,ap in zip(b_values,old_ap[1:])]
    table7=deepcopy(old_table)
    for r in table7.findall('w:tr',NS):table7.remove(r)
    row_template=source_rows[3]
    for ri,values in enumerate(joined):
        row=deepcopy(row_template)
        trpr=row.find('w:trPr',NS)
        if trpr is None:
            trpr=E.Element(W+'trPr');row.insert(0,trpr)
        setprop(trpr,'cantSplit')
        if ri==0:setprop(trpr,'tblHeader')
        for cell,value in zip(row.findall('w:tc',NS),values):
            template=cell.find('w:p',NS)
            for child in list(cell):
                if child.tag!=W+'tcPr':cell.remove(child)
            cell.append(plain(template,value))
            for run in cell.xpath('.//w:r',namespaces=NS):
                rp=run.find('w:rPr',NS)
                if rp is None:
                    rp=E.Element(W+'rPr');run.insert(0,rp)
                setprop(rp,'b',val=1 if ri in (0,1) else 0)
            setprop(ppr(cell.find('w:p',NS)),'keepNext',val=1 if ri==0 else 0)
            cp=cell.find('w:tcPr',NS)
            borders=setprop(cp,'tcBorders')
            setprop(borders,'top',val='single' if ri==0 else 'nil',sz=12,color='000000')
            setprop(borders,'bottom',val='single' if ri in (0,len(joined)-1) else 'nil',sz=6 if ri==0 else 12,color='000000')
        table7.append(row)

    sequence=[deepcopy(old[139]),heading('A 数据与标签设计'),
        para('ATLAS以重复轨迹统计建立标签，DynaRepo使用连续窗口，ITrace使用完整重复轨迹；各来源的观测设计、有效Hχ位置与零Hχ数量见表S1。观测标准差描述不同重复或窗口的读数差异。单次观测的样本方差以零占位；多次观测读数相同及观测次数未保存的零方差记录在配套诊断中分别标记。零Hχ表示参与汇总的各窗口或重复内出现一种联合状态。'),
        original(142),original(143),
        para('逐记录PDB来源注释及原API响应随数据提供[26]。标签观测方差与模型残差尺度按2.4节用于训练；标签均值、观测标准差、角数归一化熵及观测次数诊断图保存在配套源数据中。'),
        heading('B 训练与划分验证'),original(148),original(149),original(150),original(151),original(152),
        para('联合模型的RMSF与Hχ平均绝对误差相对单目标分别降低约0.8%和1.8%。默认共享模型在一次推理中提供两类读数，各训练目标的MAE、RMSE和排序相关完整列于表S3。'),
        original(170),original(171),original(172),
        heading('C 预测性能与方法比较'),
        para('表S5汇总与图2一致的双任务内部对照和ThermoPCD外部验证；表S6为全队列的额外RMSF方法比较，表S7为相同615条DynaRepo复合物的外部模型比较。各组分别采用其注明的样本、有效位置及推理配置。'),
        para('表S5：双任务数值预测结果。A为3,063条记录的五折体系留出内部比较，B为ThermoPCD的32条外部记录；各项先逐记录计算，再对记录等权平均。MAE越低越好，Spearman相关ρ越高越好。',True),table5,
        para('注：Full为完整结构模型，Intra为链内模型，Static为静态特征回归，AA为训练集氨基酸类型均值对照。B中四种方法均使用排除关联记录后的2,430条训练记录。'),
        para('全队列扩展比较覆盖3,063条记录、961,321个RMSF有效位置。逐链BackFlip-2.1采用代码提交887c6a7；PLANET-MD使用full及sequence-only配置[11]，RMSF-net使用PDB-only配置[9]。其逐链结果与表S7的复合物输入结果分别报告。'),
        original(156),
        para('PLANET-MD各链独立预测后合并，采用300 K与float32精度，并将nm换算为Å。RMSF-net将结构转换成4 Å分辨率、1.5 Å网格的原子序数加权高斯密度，评价相对分数。ANM与GNM采用ProDy 2.6.1，由完整输入的Cα构建，距离截断15 Å、弹簧常数1，保留最低20个非零模态，对平方波动之和开方[4,27]。已有方法使用作者权重，内部对照使用相应训练集合拟合。'),
        original(161),original(162),
        para('逐链BackFlip-2.1在全队列的RMSF MAE和Spearman相关为0.537 Å和0.848，在ThermoPCD的32条记录中为0.749 Å和0.593。PLANET-MD full与FlexPPI的RMSE为1.186 Å和0.792 Å。PLANET-MD两种配置的来源分层、最高10%位点命中率及逐记录结果与完整实现说明随源数据保留。'),
        para('图2d、e及表S7采用DynaRepo的616条多链记录中具有至少一个可定义界面的615条。界面由起始结构全部蛋白链间的重原子接触定义，两条链间至少10个原子对距离在5 Å以内计为一个界面。ATLAS的1,929条单链记录属于BackFlip等方法的训练来源；ITrace的244条TCR–pMHC记录属于同一依赖家族，界面数集中在7个。外部模型与界面分层比较使用上述DynaRepo集合，保持样本和有效残基一致。'),
        para('表S7：615条DynaRepo复合物的RMSF数值预测与高波动位点筛选。七种方法共同覆盖261个序列家族及260,340个有效残基；各项按记录等权平均。MAE越低越好，ρ、AP和AUROC越高越好；图2d、e采用相同样本与方法。',True),table7,
        para('注：BackFlip两个版本采用完整结构的assembly输入；DynaProt-M按作者实现逐链推理后合并。FastProtFlex、ANM和GNM输出相对分数，Å单位MAE以“—”表示不适用。高RMSF阳性取各记录参考第90百分位及以上，保留并列值；AP按原始阈值组计算。1个、2–3个和至少4个界面的组分别含414条（196个家族）、144条（58个家族）和57条（34个家族），同一家族可跨组。'),
    ]
    # Reuse the source figure container, linking a new PNG without touching old assets.
    image_para=deepcopy(old[165])
    image_pr=ppr(image_para);setprop(image_pr,'pageBreakBefore',val=0);setprop(image_pr,'keepNext',val=1)
    rid='rIdSupplementMergedV74'
    E.SubElement(rels,REL+'Relationship',Id=rid,Type='http://schemas.openxmlformats.org/officeDocument/2006/relationships/image',Target='media/figS1_merged_v74.png')
    for node in image_para.xpath('.//a:blip',namespaces=NS):node.set('{'+NS['r']+'}embed',rid)
    im=Image.open(merged_image);ratio=im.height/im.width
    extent=image_para.find('.//wp:extent',NS);cx=int(extent.get('cx'));cy=round(cx*ratio)
    extent.set('cy',str(cy))
    for node in image_para.xpath('.//a:xfrm/a:ext',namespaces=NS):node.set('cx',str(cx));node.set('cy',str(cy))
    for node in image_para.xpath('.//wp:docPr',namespaces=NS):node.set('name','Supplementary Figure S1');node.set('descr','来源分层平均AP与内部高动态位点PR曲线')
    sequence += [image_para,
        para('图S1：内部对照的来源分层与高动态位点筛选。（a）RMSF、总体Hχ及同类型Hχ的来源分层平均AP，每组依次为ATLAS、DynaRepo、ITrace；n为记录数，热图统一使用0–1色标。（b、c）全队列RMSF与总体Hχ的PR曲线。以各记录MD参考第90百分位及以上为阳性，保留真实值与预测分数的并列组；各曲线在共同召回率处取首次达到该召回率的阈值组精确率，再对3,063条记录等权平均。图中数值为逐记录原始阈值计算后平均的AP；灰横线为实际平均阳性率，分别为10.26%和10.29%，末端保留全部入选状态。Full、Intra、Static、AA及同类型汇总方式与图2一致。',True),
        para('图2c及图S1b、c的Precision轴采用完整覆盖0–1的连续分段线性显示：RMSF的0–0.5与0.5–1分别占显示高度25%与75%；总体Hχ的0–0.5、0.5–0.9与0.9–1分别占10%、35%与55%；同类型Hχ的0–0.45、0.45–0.75与0.75–1分别占12%、78%与10%。各轴标注真实数值，显示曲线面积不代替平均AP。图2b四分位范围描述记录间差异；图2d的外部PR与ROC采用常规线性坐标。'),
        heading('D 复合物实例补充验证'),original(182),original(183),original(184),
        para('表S8报告KIX各状态内的定位与幅值。状态间比较的六对L12差值、前后半程统计及相应MD与预测数值保存在病例源表。'),
        heading('E 联合柔性与构象占据验证'),original(174),original(175),original(176),
        para('表S9并列报告固定预测在原汇总标签及第一条共同窗口标签上的表现；完整模型相对静态对照的两项平均误差均较低。采用三条共同窗口轨迹的均值标签时，ATLAS的RMSF／Hχ MAE为0.723 Å／0.172 nats，ITrace为0.158 Å／0.110 nats，DynaRepo采用一条轨迹。'),
        original(178),original(179),original(180),original(187),original(188),original(189),deepcopy(old[199])]

    # Caption following an image ends the keep-with-next chain.
    merged_caption=sequence[sequence.index(image_para)+1]
    setprop(ppr(merged_caption),'keepNext',val=0)
    # Keep each short table and its caption together, avoiding orphaned headers
    # and single trailing rows on otherwise unrelated pages.
    for element in sequence:
        if element.tag!=W+'tbl':continue
        table_rows=element.findall('w:tr',NS)
        for ri,row in enumerate(table_rows):
            for paragraph in row.xpath('.//w:p',namespaces=NS):
                setprop(ppr(paragraph),'keepNext',val=0 if ri==len(table_rows)-1 else 1)
                setprop(ppr(paragraph),'keepLines',val=1)
    kix_intro=next(el for el in sequence if txt(el).startswith('四种KIX组成各使用'))
    setprop(ppr(kix_intro),'keepNext',val=1)
    # Preserve all original data in unchanged tables and both split/merged tables.
    for source_index in [143,150,152,162,172,176,179,184,189]:
        source_values=rows(old[source_index])
        assert any(el.tag==W+'tbl' and rows(el)==source_values for el in sequence),source_index
    assert rows(table5)[3:7]==old_values[3:7]
    assert rows(table5)[8:12]==old_values[16:20]
    assert rows(table7)==joined
    # Update references only in the pre-existing manuscript, then replace supplement.
    before_text=[txt(el) for el in old[:supplement_start]]
    for el in old[:supplement_start]:remap(el)
    reference_changes=[{'body_index':i,'old':before_text[i],'new':txt(el)} for i,el in enumerate(old[:supplement_start]) if before_text[i]!=txt(el)]
    for el in old[supplement_start:]:body.remove(el)
    for el in sequence:body.append(el)
    assert Counter(E.tostring(el) for el in root.xpath('//m:oMath',namespaces=NS))==all_math
    # Main figure XML may change caption references; image containers are preserved.
    new_main_images=[E.tostring(el) for el in list(body)[:supplement_start] if el.xpath('.//w:drawing',namespaces=NS)]
    assert new_main_images==source_main_images
    assert len(root.xpath('//w:tbl',namespaces=NS))==11
    captions=[txt(el) for el in body if el.tag==W+'p' and re.match(r'表S\d+：',txt(el))]
    assert [int(re.match(r'表S(\d+)：',x)[1]) for x in captions]==list(range(1,12))
    assert [txt(el) for el in body if el.tag==W+'p' and re.match(r'[A-I] ',txt(el))]==[
        'A 数据与标签设计','B 训练与划分验证','C 预测性能与方法比较','D 复合物实例补充验证','E 联合柔性与构象占据验证']
    assert not any(re.search(r'补充材料\s*[F-I]|图S[23](?!\d)',txt(el)) for el in body)
    modified={'word/document.xml':E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),
              'word/_rels/document.xml.rels':E.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True)}
    with zipfile.ZipFile(TARGET,'w',zipfile.ZIP_DEFLATED) as package:
        for entry in entries:package.writestr(entry,modified.get(entry.filename,parts[entry.filename]))
        package.writestr('word/media/figS1_merged_v74.png',merged_image.read_bytes())
    with zipfile.ZipFile(TARGET) as package:assert package.testzip() is None
    audit={'source':str(SOURCE),'target':str(TARGET),'target_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),
           'sections':5,'tables':11,'supplementary_figures':1,'table_mapping':TABLE_MAP,
           'all_original_table_measurements_retained':True,'all_native_math_preserved':True,
           'main_figure_containers_preserved_exactly':True,'source_diagnostic_archived':str(archive),
           'reference_changes':reference_changes,'figure_source':str(merged_image),
           'visual_qa':'pending'}
    (ROOT/'v74_patch_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k!='reference_changes'},ensure_ascii=False))

if __name__=='__main__':main()
