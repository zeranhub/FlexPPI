"""Add the approved bound-state antibody narrative, preserving all existing media."""
from pathlib import Path
from copy import deepcopy
from zipfile import ZipFile, ZIP_DEFLATED
from hashlib import sha256
from lxml import etree
import csv, json

ROOT=Path(__file__).resolve().parent
DATA=ROOT.parent
PAPER=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
SOURCE=PAPER/'paper_v128_20261028_中文修订_单点突变与残基动态叙事.docx'
TARGET=PAPER/'paper_v129_20261029_中文修订_抗原识别与抗体动态.docx'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
q=lambda n:'{'+NS['w']+'}'+n
with ZipFile(SOURCE) as z:
    items=z.infolist();parts={x.filename:z.read(x.filename) for x in items}
assert sha256(SOURCE.read_bytes()).hexdigest()=='1dd44957215bebde0ff2ff71cfa44ce2b8398d089da1f71597396b35993174cb'
doc=etree.fromstring(parts['word/document.xml']);body=doc.find('w:body',NS)
ps=body.findall('w:p',NS)
text=lambda p:''.join(p.xpath('.//w:t/text()',namespaces=NS))
def find(prefix):
    nodes=[p for p in ps if text(p).startswith(prefix)]
    assert len(nodes)==1,(prefix,len(nodes))
    return nodes[0]
def paragraph(value,template):
    p=etree.Element(q('p'))
    pp=template.find('w:pPr',NS)
    if pp is not None:p.append(deepcopy(pp))
    r=etree.SubElement(p,q('r'));rp=template.find('w:r/w:rPr',NS)
    if rp is not None:r.append(deepcopy(rp))
    etree.SubElement(r,q('t')).text=value
    return p
changes=[]
def replace(prefix,value):
    old=find(prefix);new=paragraph(value,old);body.replace(old,new)
    changes.append({'operation':'replace','before':text(old),'after':value})
    return new
def insert_after(anchor,value,template=None):
    p=paragraph(value,template if template is not None else anchor)
    anchor.addnext(p)
    changes.append({'operation':'insert_after','anchor':text(anchor),'after':value})
    return p

TITLE='3.4 抗原识别中TCR与抗体的残基动态'
replace('3.4 TCR亲和成熟与肽替换中的残基动态响应',TITLE)
replace('这些层面的贡献将在后续分析中落实到具体结构位置。',
        '这些层面的贡献将在后续分析中落实到具体结构位置。3.3节结合整体RMSF–Hχ分布、CaM–CaMKI与泛素–UBA两个实例及MD角度占据，依次检验双尺度柔性预测的准确性、结构环境对应及侧链构象含义；3.4节结合MEL5亲和成熟、肽替换和B5209B–Robo1抗体–抗原复合物，将两类读数进一步联系到识别环及跨链接触组织。')

METHOD=(
    'B5209B–Robo1抗体–抗原实例采用DynaRepo记录A00LZ，其起始复合物来源为PDB 6A77[37]。'
    '预测使用全组成链家族划分第0折的Full与Intra第40轮权重，该记录及其关联家族位于测试集合。'
    '在同一结合态输入上核对全部六个CDR，并结合原有轨迹的逐帧接触与侧链状态占据解释局部读数；输入、编号、参考范围及分析方法见补充材料H。'
)
insert_after(find('排除与外部集合关联的633条训练记录后'),METHOD,find('外部评价采用ThermoPCD的32条'))

P1=(
    'B5209B Fab–Robo1 Ig5复合物进一步展示了抗体识别环中主链与侧链约束的差异（PDB：6A77）[37]。'
    '全组成链家族留出预测在完整H-CDR3的15个主链位置及12个有效侧链位置上的MAE分别为0.116 Å和0.160 nats。'
    'R101的预测RMSF为1.149 Å，接近MD参考的1.123 Å，链内模型为1.640 Å；其预测与参考Hχ分别为3.431和3.419 nats。'
    '起始结构中，R101主链羰基O与轻链Y49羟基O相距3.21 Å，侧链则朝向溶剂。'
    'MD中，R101主链与轻链接触的帧占比为90.3%，侧链仅为2.3%（≤4 Å）。'
    '这一对应表明，主链参与链间接触可以与同一残基侧链保留较多取向并存，残基层面的接触标签不能代替对两类动态自由度的区分。'
)
P2=(
    '同一H-CDR3中的R96提供了同类型侧链参照。R96位于密集的重链内部环境，起始结构侧链溶剂可及面积为16.20 Å²，R101为189.13 Å²；'
    'R96的预测与MD参考Hχ分别为1.969和1.801 nats，均低于R101。'
    '相同化学类型在同一识别环中的差异，将侧链状态读数联系到局部堆积与暴露环境。'
    '结合R101的主链接触，这些对应将抗体识别环的主链波动与侧链多样性联系到不同的局部环境和原子接触方式，为分析结合态识别环的接触组织提供依据。'
    '完整六个识别环及同类型对照的结果见补充材料H。'
)
anchor=insert_after(find('重复MD将这些结合态差异进一步联系到未替换位置参与的接触重组。'),P1)
keep_lines=anchor.find('w:pPr/w:keepLines',NS)
assert keep_lines is not None
keep_lines.set(q('val'),'1')
insert_after(anchor,P2)

DISCUSSION=(
    'B5209B–Robo1抗体实例中，R101的主链接触与高侧链多样性并存，同型R96的内部堆积则对应较低Hχ，'
    '进一步将双尺度读数联系到不同的局部约束。'
)
insert_after(find('MEL5家族留出实例'),DISCUSSION)

REFERENCE=(
    '[37] T. Yamashita, E. Mizohata, S. Nagatoishi, T. Watanabe, M. Nakakido, H. Iwanari, Y. Mochizuki, T. Nakayama, '
    'Y. Kado, Y. Yokota, H. Matsumura, T. Kawamura, T. Kodama, T. Hamakubo, T. Inoue, H. Fujitani, and K. Tsumoto. '
    'Affinity Improvement of a Cancer-Targeted Antibody through Alanine-Induced Adjustment of Antigen-Antibody Interface. '
    'Structure, 27(3):519–527.e5, 2019. doi:10.1016/j.str.2018.11.002.'
)
assert not any('10.1016/j.str.2018.11.002' in text(p).lower() for p in ps)
insert_after(find('[36]'),REFERENCE)

SUPPLEMENT=[
    'H B5209B–Robo1识别环动态与接触核查',
    '输入与参考。A00LZ使用冻结预处理结构的全部519个残基，含Robo1 Ig5链A 91位、Fab轻链B 211位和重链C 217位。'
    '重轻链按ANARCII编号[36]界定全部六个CDR，正文使用与原结构对应的作者残基编号；H-CDR3为重链95–109共15位，具有12个有效Hχ位置。'
    '结构中的组氨酸和异亮氨酸原子名按既有规范化版本处理，坐标和7778个原子的对应保持不变。'
    'Full与Intra采用2.5节第0折、种子7、第40轮检查点，在同一519位结合态结构上独立推理，使用同一有效掩码计分。',
    'MD参考与采样范围。参考来自本地可用的一条1000帧XTC，按2.2节逐链对齐重算RMSF，按原固定χ角分箱定义重算联合Hχ。'
    'H-CDR3 RMSF与冻结标签的最大差为4.2×10⁻⁷ Å，有效Hχ的最大差为1.1×10⁻⁷ nats。'
    '本地仅有一条轨迹，其头部时标记录为0–999 ps，与来源元数据的三条100 ns模拟描述不能直接对等；未取得原始保存帧映射，因此不据此宣称已核实物理采样时长、三份独立重复或收敛。'
    '四个连续250帧块仅用于检查本条轨迹内结果的稳定性，不作为独立重复，也不作显著性检验。',
    '接触与暴露。对完整H-CDR3的15位与全部519位蛋白残基保留逐帧最小重原子距离，采用周期盒最小镜像约定。'
    '主口径以最近距离≤4 Å记为接触，并保留4.5和5 Å结果；链内序列间隔不超过两位的残基不计入非局部接触汇总。'
    '主链定义为N、Cα、C、O及OXT，其余非氢原子定义为侧链；分别统计各组与轻链或抗原的接触帧占比，同时接触取同一帧两项条件的交集。'
    '静态侧链SASA采用1.4 Å探针、960个球面采样点，并以240点作敏感性核查，面积均在完整复合物环境中计算。'
    'R101主链与轻链的4 Å接触占比在四个连续块为95.2%、88.4%、87.2%、90.4%，侧链为0%、9.2%、0%、0%；'
    '接触占比描述采样帧中的几何邻接，不等同接触寿命、氢键占据或结合能。',
    '完整区域结果。全部六个CDR均保留，未按模型表现删除识别环。RMSF MAE按轻链CDR1、CDR2、CDR3及重链CDR1、CDR2、CDR3的顺序，'
    'Full为0.272、0.080、0.182、0.302、0.119、0.116 Å，Intra为0.129、0.047、0.058、0.422、0.146、0.290 Å。'
    'H-CDR3的Hχ MAE分别为0.160和0.255 nats；完整519位的RMSF MAE分别为0.236和0.390 Å，抗体428位为0.198和0.215 Å。'
    '三个轻链CDR的RMSF误差均以Intra更小，不能将H-CDR3的改善概括为全部识别环一致改善。'
    'R96/R101是H-CDR3中的全部两个Arg，Intra的Hχ为1.782/3.358 nats，也恢复了同型差异；该对照用于解释局部环境对应，不单独证明Full的跨链新增优势。'
    'Full未完整恢复H-CDR3三个Tyr的同型排序，因此不将局部绝对误差优势解释为全部同型侧链环境差异均获恢复。',
    '解释范围与复现。该实例是WT Fab–抗原结合态的结构—动态对应，不是完整膜BCR，也不提供结合前后或WT/P103A差值的动态验证。'
    'Hχ表征固定粗分箱下侧链取向占据多样性，不等同总结合熵；Pro的零Hχ不表示其主链不动。'
    '本文据R101的原子接触和R96/R101同型读数讨论局部约束，不据此认定能量热点、突变敏感性或亲和力贡献。'
    '完整519位与六环源表、全部H-CDR3—伙伴距离、主链/侧链占据、侧链联合状态、不同拟合参考及输入和检查点SHA256均随复现材料提供。'
    '独立核查以手工最小镜像法复算10帧全部15×519个距离（最大差8.4×10⁻⁶ Å），另对六个关注位点整条轨迹复算主链/侧链接触及联合熵，均与原统计一致。'
]
last=body.findall('w:p',NS)[-1]
last=insert_after(last,SUPPLEMENT[0],find('G TCR结合模型对照'))
for value in SUPPLEMENT[1:]:
    last=insert_after(last,value,find('全部分数、输入及训练覆盖在配套CSV'))

parts['word/document.xml']=etree.tostring(doc,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
    for item in items:z.writestr(item,parts[item.filename])
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
    assert b.testzip() is None
    changed=[n for n in a.namelist() if a.read(n)!=b.read(n)]
    assert changed==['word/document.xml']
    assert a.namelist()==b.namelist()
    original=etree.fromstring(a.read('word/document.xml'))
    for xpath in ['//w:tbl','//w:fldSimple','//w:instrText','//w:commentRangeStart','//w:commentRangeEnd','//w:sectPr','//w:drawing']:
        before=original.xpath(xpath,namespaces=NS);after=doc.xpath(xpath,namespaces=NS)
        assert len(before)==len(after),xpath
        assert all(etree.tostring(x)==etree.tostring(y) for x,y in zip(before,after)),xpath
    media=[n for n in a.namelist() if n.startswith('word/media/')]
    assert len(media)==50 and all(a.read(n)==b.read(n) for n in media)
    newtexts=[text(p) for p in body.findall('w:p',NS)]
    for i in range(79,87):
        assert text(ps[i]) in newtexts
    assert text(find('图4：MEL5体系')) in newtexts
    assert text(find('5 结论')) in newtexts
    assert text(ps[2]) in newtexts
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
        'target_sha256':sha256(TARGET.read_bytes()).hexdigest(),'title':TITLE,
        'changed_original_package_parts':changed,'original_media_preserved':len(media),
        'new_images_inserted':0,'existing_figure4_caption_unchanged':True,
        'all_original_TCR_result_paragraphs_unchanged':True,
        'abstract_and_conclusion_unchanged':True,'original_tables_fields_comments_sections_drawings_unchanged':True,
        'changes':changes,'render_validation':'pending'}
(ROOT/'revision_changes.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'新增叙事与方法_中文.md').write_text(TITLE+'\n\n'+P1+'\n\n'+P2+'\n\n讨论\n\n'+DISCUSSION+'\n\n方法2.5新增\n\n'+METHOD+'\n\n'+REFERENCE+'\n\n'+'\n\n'.join(SUPPLEMENT),encoding='utf-8')
print(json.dumps({k:report[k] for k in ['target','title','new_images_inserted','original_media_preserved','target_sha256']},ensure_ascii=False,indent=2))
