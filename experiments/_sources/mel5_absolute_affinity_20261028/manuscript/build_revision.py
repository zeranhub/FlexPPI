"""Produce a new Chinese revision, preserving all original media and tables."""
from pathlib import Path
from copy import deepcopy
from zipfile import ZipFile, ZIP_DEFLATED
from hashlib import sha256
import csv
import json
import math
from lxml import etree
from docx import Document
from docx.shared import Cm

ROOT = Path(__file__).resolve().parent
EXPERIMENT = ROOT.parent
PAPER = Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010")
SOURCE = PAPER / "paper_v125_20261028_中文修订_MEL5核心优势叙述.docx"
TARGET = PAPER / "paper_v126_20261028_中文修订_绝对亲和力对照与肽突变验证.docx"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
      "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
      "pic": "http://schemas.openxmlformats.org/drawingml/2006/picture"}
qn = lambda name: "{" + NS[name.split(":")[0]] + "}" + name.split(":")[1]
with ZipFile(SOURCE) as archive:
    items = archive.infolist()
    parts = {item.filename: archive.read(item.filename) for item in items}
document = etree.fromstring(parts["word/document.xml"])
body = document.find("w:body", NS)
paragraphs = body.findall("w:p", NS)
text = lambda p: "".join(p.xpath(".//w:t/text()", namespaces=NS))
assert len(paragraphs) == 187
original_text = [text(p) for p in paragraphs]
assert original_text[83].startswith("受体亲和成熟与肽替换还呈现不同的区域响应")
with (EXPERIMENT / "affinity/absolute_affinity_results.csv").open(encoding="utf-8-sig", newline="") as handle:
    affinity = list(csv.DictReader(handle))
methods = list(dict.fromkeys(row["model"] for row in affinity))
assert len(methods) == 2 and len(affinity) == 8, "Both prespecified absolute models must be complete"
prodigy = {row["pdb"]: row for row in affinity if row["model"].startswith("PRODIGY")}
gnn = {row["pdb"]: row for row in affinity if row["model"].startswith("ProAffinity")}
assert set(prodigy) == set(gnn) == {"3hg1", "4jff", "4jfd", "4jfe"}
num = lambda row, key: float(row[key])
gnn_p = [num(gnn[pdb], "predicted_pKD") for pdb in ["3hg1", "4jff", "4jfd", "4jfe"]]

def fresh_paragraph(value, template=83, page_break=False):
    p = etree.Element(qn("w:p"))
    ppr = paragraphs[template].find("w:pPr", NS)
    ppr = deepcopy(ppr) if ppr is not None else etree.Element(qn("w:pPr"))
    for child in list(ppr):
        if child.tag in [qn("w:sectPr"), qn("w:pageBreakBefore")]:
            ppr.remove(child)
    if page_break:
        etree.SubElement(ppr, qn("w:pageBreakBefore"))
    p.append(ppr)
    run = etree.SubElement(p, qn("w:r"))
    old_rpr = paragraphs[template].find("w:r/w:rPr", NS)
    if old_rpr is not None:
        run.append(deepcopy(old_rpr))
    t = etree.SubElement(run, qn("w:t"))
    t.text = value
    return p

changes = []
def replace(index, value):
    old = paragraphs[index]
    new = fresh_paragraph(value, index)
    body.replace(old, new)
    changes.append({"source_index": index, "before": original_text[index], "after": value})
    return new

replace(78, "3.4 TCR亲和成熟与肽替换中的残基动态响应")
affinity_paragraph = (
    "肽单点替换提供了直接接触相近而亲和力显著不同的对照。在成熟受体背景中，G4A与I7A使SPR解离常数由0.6 nM分别增至36与31 μM，亲和力下降约60,000与51,700倍；原研究将不同热力学特征联系到溶剂重排[28]。我们对四种晶体复合物分别预测绝对亲和力，将TCR α、β链作为一方，完整pMHC作为另一方。PRODIGY对ELA、G4A和I7A的预测解离常数分别为0.203、0.200和0.113 nM，未恢复实验中的下降；作者简化推理版ProAffinity-GNN的对应pKD分别为"
    f"{gnn_p[1]:.3f}、{gnn_p[2]:.3f}和{gnn_p[3]:.3f}，同样未再现实验效应量；完整四态结果、输入处理及适用范围见图S5和补充材料F[31,32]。"
    "这一比较检验亲和力数值的恢复，残基柔性读数则用于进一步定位不同结合态中的动态响应。"
)
anchor = replace(83, affinity_paragraph)
dynamic_paragraph = (
    "肽替换的动态对应具有状态与区域选择性。I7A对照中，完整45个CDR位置的RMSF变化有37位与三重复MD均值同向，变化排序Spearman相关为0.498；CDR2α的全部六位方向一致、相关为0.943，MD区域均值由1.153升至1.427 Å，预测由0.878升至0.989 Å。CDR3β的MD均值由0.926升至1.194 Å，三条I7A重复的区域均值均高于三条ELA重复，但预测由0.648升至0.671 Å，低估了变化幅度。G4A的完整CDR变化排序相关为−0.095、方向一致29/45位；并非两种肽替换的全部响应均被恢复。完整肽、识别环和同型界面侧链的逐位置与重复核查见图S6及配套源表。"
)
contact_paragraph = (
    "重复MD将这些结合态差异进一步联系到未替换位置参与的接触重组。在≤4 Å的最小重原子距离口径下，I7A中P7–Thrβ96的平均接触占据由36.0%降至1.6%，未替换P8与Glyβ97、Glyβ99的占据分别由1.7%和15.3%升至45.6%和39.2%；三项在两态三重复范围间均分离，并在4.5 Å阈值下保持相同方向。G4A中，P8–HLA Arg97的平均占据由1.0%升至37.8%，重复间幅度差异较大（图S7）。这些结果为肽背景相关的局部约束和接触重分配提供轨迹证据；结合态RMSF与Hχ不等同于总结合熵，以上对应尚不计量亲和力巨变的因果贡献。"
)
for value in [dynamic_paragraph, contact_paragraph]:
    new = fresh_paragraph(value)
    anchor.addnext(new)
    anchor = new
    changes.append({"inserted_after_source_index": 83, "text": value})
replace(89, original_text[89] +
        " 肽替换的完整核查显示，I7A的识别环主链变化与MD具有区域性对应，G4A的恢复则较弱。亲和力评分、残基动态与接触占据分别描述效应量、响应位置和轨迹中的局部组织，三者的结合为机制研究提供可检验线索；蛋白结合态柔性与原研究强调的溶剂热力学应分别解释。")

references = [
    "[31] L. C. Xue, J. P. G. L. M. Rodrigues, P. L. Kastritis, A. M. J. J. Bonvin, and A. Vangone. PRODIGY: a web server for predicting the binding affinity of protein–protein complexes. Bioinformatics, 32(23):3676–3678, 2016. doi:10.1093/bioinformatics/btw514.",
    "[32] Z. Zhou, Y. Yin, H. Han, Y. Jia, J. H. Koh, A. W.-K. Kong, and Y. Mu. ProAffinity-GNN: A Novel Approach to Structure-Based Protein–Protein Binding Affinity Prediction via a Curated Data Set and Graph Neural Networks. Journal of Chemical Information and Modeling, 64(23):8796–8808, 2024. doi:10.1021/acs.jcim.4c01850."
]
for value in references:
    paragraphs[127].addprevious(fresh_paragraph(value, 126))

supplement = []
def append_text(value, template=165, page_break=False):
    p = fresh_paragraph(value, template, page_break)
    body.insert(len(body) - 1, p)
    supplement.append(value)
    return p

append_text("F 绝对亲和力对照与肽替换动态验证", 156, True)
append_text(
    "绝对亲和力分别输入RCSB原始3HG1、4JFF、4JFD和4JFE晶体，经同一流程保留首个模型、标准蛋白残基及一个替代构象。相互作用两方固定为D、E链的TCR与A、B、C链的HLA、β2M及十肽，全部已解析重原子保留。PRODIGY 2.4.0与FreeSASA 2.2.1按作者实现计算，温度为与SPR一致的25 °C[31]。pKD=−log10(KD/M)，越高表示结合越强。各状态独立预测绝对亲和力，不使用突变ΔΔG专用模型或本例实验值校准。"
)
append_text(
    "ProAffinity-GNN使用作者公开简化版（提交0dd07bb）及原始ESM2-650M、GNN权重，严格加载且未重训[32]；该版本按已解析结构提取序列，作者说明精度可能受影响。原程序pMHC特征循环错误复用末条肽的输出，修正一个赋值变量output1为output2后核对全部图节点数与边索引。PDBQT采用OpenBabel加氢及AutoDockTools的标准Cys/Met硫SA规则，晶体重原子身份和坐标不变；其准备工具不同于作者建议的ADFR。可获得的训练PDB清单未含四个目标编号，但尚未确认同源家族隔离。版本、权重哈希、原始代码与修正记录随源数据提供，结果限定于该回顾性实现。"
)
append_text("表S7：四态绝对亲和力对照。KD统一以nM表示，ProAffinity-GNN由作者输出的负十进制亲和力标度pKa换算。结果对应规定的结构准备和公开权重，不表示独立体系上的总体性能。", 162)

table = etree.Element(qn("w:tbl"))
tblpr = etree.SubElement(table, qn("w:tblPr"))
etree.SubElement(tblpr, qn("w:tblW"), {qn("w:w"): "9500", qn("w:type"): "dxa"})
etree.SubElement(tblpr, qn("w:tblLayout"), {qn("w:type"): "fixed"})
borders = etree.SubElement(tblpr, qn("w:tblBorders"))
for side in ["top", "left", "bottom", "right", "insideH", "insideV"]:
    etree.SubElement(borders, qn("w:" + side), {qn("w:val"): "single", qn("w:sz"): "4", qn("w:color"): "D9D9D9"})
grid = etree.SubElement(table, qn("w:tblGrid"))
widths = [2600, 1300, 1850, 1850, 1900]
for width in widths:
    etree.SubElement(grid, qn("w:gridCol"), {qn("w:w"): str(width)})
rows = [["状态", "PDB", "SPR KD (nM)", "PRODIGY KD (nM)", "ProAffinity KD (nM)"]]
for pdb, state in [("3hg1", "MEL5 / ELA"), ("4jff", "α24β17 / ELA"),
                   ("4jfd", "α24β17 / G4A"), ("4jfe", "α24β17 / I7A")]:
    rows.append([state, pdb.upper(), f"{num(prodigy[pdb], 'experimental_KD_M') * 1e9:.4g}",
                 f"{num(prodigy[pdb], 'predicted_KD_M') * 1e9:.4g}",
                 f"{num(gnn[pdb], 'predicted_KD_M') * 1e9:.4g}"])
for row_index, values in enumerate(rows):
    tr = etree.SubElement(table, qn("w:tr"))
    for width, value in zip(widths, values):
        tc = etree.SubElement(tr, qn("w:tc"))
        tcpr = etree.SubElement(tc, qn("w:tcPr"))
        etree.SubElement(tcpr, qn("w:tcW"), {qn("w:w"): str(width), qn("w:type"): "dxa"})
        etree.SubElement(tcpr, qn("w:vAlign"), {qn("w:val"): "center"})
        if row_index == 0:
            etree.SubElement(tcpr, qn("w:shd"), {qn("w:fill"): "E8EEF4"})
        p = etree.SubElement(tc, qn("w:p"))
        pp = etree.SubElement(p, qn("w:pPr"))
        etree.SubElement(pp, qn("w:ind"), {qn("w:left"): "0", qn("w:right"): "0", qn("w:firstLine"): "0", qn("w:firstLineChars"): "0"})
        etree.SubElement(pp, qn("w:jc"), {qn("w:val"): "center"})
        etree.SubElement(pp, qn("w:snapToGrid"), {qn("w:val"): "0"})
        etree.SubElement(pp, qn("w:spacing"), {qn("w:before"): "70", qn("w:after"): "70", qn("w:line"): "240", qn("w:lineRule"): "auto"})
        rr = etree.SubElement(p, qn("w:r"))
        rp = etree.SubElement(rr, qn("w:rPr"))
        etree.SubElement(rp, qn("w:rFonts"), {qn("w:ascii"): "Arial", qn("w:hAnsi"): "Arial", qn("w:eastAsia"): "宋体"})
        etree.SubElement(rp, qn("w:sz"), {qn("w:val"): "18"})
        if row_index == 0:
            etree.SubElement(rp, qn("w:b"))
        etree.SubElement(rr, qn("w:t")).text = value
body.insert(len(body) - 1, table)

relationships = etree.fromstring(parts["word/_rels/document.xml.rels"])
REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
rid_number = max(int(x.get("Id")[3:]) for x in relationships
                 if x.get("Id", "").startswith("rId") and x.get("Id")[3:].isdigit()) + 1
docpr_number = max(int(x.get("id")) for x in document.xpath("//wp:docPr", namespaces=NS)) + 1
new_media = {}
def append_figure(filename, page_break=False, width_cm=17.0):
    global rid_number, docpr_number
    path = EXPERIMENT / "figures" / filename
    temp = Document()
    temp.add_picture(str(path), width=Cm(width_cm))
    p = etree.fromstring(etree.tostring(temp.paragraphs[-1]._p))
    pp = etree.SubElement(p, qn("w:pPr")) if p.find("w:pPr", NS) is None else p.find("w:pPr", NS)
    etree.SubElement(pp, qn("w:jc"), {qn("w:val"): "center"})
    etree.SubElement(pp, qn("w:spacing"), {qn("w:before"): "0", qn("w:after"): "40", qn("w:line"): "240", qn("w:lineRule"): "auto"})
    etree.SubElement(pp, qn("w:keepNext"))
    if page_break:
        etree.SubElement(pp, qn("w:pageBreakBefore"))
    p.remove(pp)
    p.insert(0, pp)
    rid = f"rId{rid_number}"
    rid_number += 1
    image_name = "v126_" + filename
    for blip in p.xpath(".//a:blip", namespaces=NS):
        blip.set(qn("r:embed"), rid)
    for docpr in p.xpath(".//wp:docPr", namespaces=NS):
        docpr.set("id", str(docpr_number))
        docpr.set("name", image_name)
        docpr_number += 1
    etree.SubElement(relationships, "{" + REL_NS + "}Relationship", {
        "Id": rid, "Type": NS["r"] + "/image", "Target": "media/" + image_name})
    new_media["word/media/" + image_name] = path.read_bytes()
    body.insert(len(body) - 1, p)

append_figure("figS5_absolute_affinity.png", width_cm=16.9)
append_text("图S5：四态绝对亲和力与SPR实验比较。pKD=−log10(KD/M)，较高值表示较强结合；同一状态的实验和各模型结果分别显示。每个状态使用一套晶体结构，各方法使用公开实现及权重。没有用四个相关状态计算独立体系上的性能显著性，也没有将结构推理次数作为实验重复。完整原始输出与全精度数值见源数据。", 174)

append_figure("figS6_peptide_dynamic_validation.png", page_break=True, width_cm=17.0)
append_text("图S6：肽替换动态响应的完整核查。（a、b）G4A和I7A全部十个肽位置的RMSF变化，MD误差线为两态独立均值差的SE，即√[(SD²前+SD²后)/3]。（c、d）全部45个CDR位置的预测–MD变化；（e、f）全部30个pMHC界面位置中同型且χ角有效的23、22位Hχ变化。ρ为位置变化的Spearman相关，direction为同向位置数；残基不作为独立生物重复。（g、h）完整CDR3β的三条独立轨迹区域均值，黑点及横线为重复值与均值，蓝方块为预测。全部差值为肽替换态减成熟ELA态。", 174)

append_text("肽接触占据核查", 156, True)
append_text(
    "对4JFF、4JFD、4JFE各三条已有轨迹保留0–200 ns全部保存帧，不运行新MD。候选为三态准备起始结构中距完整十肽任一重原子7 Å以内的非肽蛋白残基并集，共88位，其中HLA 63位、TCRα 15位、TCRβ 10位；全部候选在三态中氨基酸类型一致。对完整十肽与88邻居的全部880种组合，采用每帧周期盒最小镜像计算最小重原子距离，并以≤4 Å和≤4.5 Å统计占据。先按每条独立轨迹统计，再对三重复等权汇总；占据率不是接触寿命或结合能。图S7的局部对照为完整筛查后的探索性核查，不作显著性或因果检验。"
)
append_figure("figS7_peptide_contact_redistribution.png", width_cm=17.0)
append_text("图S7：肽替换的接触占据重分配。（a、b）G4A、I7A与ELA的全部880对肽–邻居平均占据，≤4 Å口径，零接触组合保留。（c–f）完整筛查后核查的四对残基，显示每态全部三条独立轨迹的占据及均值，蓝圆与橙三角分别为≤4 Å和≤4.5 Å。TCR采用MEL5共同编号；Thrβ96、Glyβ97、Glyβ99对应成熟输入E97、E98、E100。没有跨状态配对模拟，三个重复之间不连线。", 174)
append_text(
    "动态分析复用图4家族留出预测，四态12条已有轨迹的标签均值与方差在3×10⁻⁶内复现。完整源表保留全部十肽、45个CDR、30个pMHC界面位置及不一致结果。处理后TCRα端部存在差异：4JFD较4JFF少一个末端附近Phe，4JFE首两位排列不同，因此属于既有结构状态比较。processed轨迹仅含蛋白，不能验证动态水桥或溶剂熵；结合态RMSF与Hχ不换算为总结合熵或结合自由能。"
)

parts["word/document.xml"] = etree.tostring(document, xml_declaration=True, encoding="UTF-8", standalone=True)
parts["word/_rels/document.xml.rels"] = etree.tostring(relationships, xml_declaration=True, encoding="UTF-8", standalone=True)
content_types = etree.fromstring(parts["[Content_Types].xml"])
CTNS = "http://schemas.openxmlformats.org/package/2006/content-types"
if not any(x.get("Extension") == "png" for x in content_types):
    etree.SubElement(content_types, "{" + CTNS + "}Default", {"Extension": "png", "ContentType": "image/png"})
    parts["[Content_Types].xml"] = etree.tostring(content_types, xml_declaration=True, encoding="UTF-8", standalone=True)
parts.update(new_media)
with ZipFile(TARGET, "w", compression=ZIP_DEFLATED, compresslevel=6) as archive:
    for item in items:
        archive.writestr(item, parts[item.filename])
    for name, data in new_media.items():
        archive.writestr(name, data)
with ZipFile(SOURCE) as a, ZipFile(TARGET) as b:
    assert b.testzip() is None
    original_media = [name for name in a.namelist() if name.startswith("word/media/")]
    assert all(a.read(name) == b.read(name) for name in original_media)
    changed = [name for name in a.namelist() if a.read(name) != b.read(name)]
    assert set(changed) <= {"word/document.xml", "word/_rels/document.xml.rels", "[Content_Types].xml"}
    old_root = etree.fromstring(a.read("word/document.xml"))
    for query in ["//w:tbl", "//w:fldSimple", "//w:instrText", "//w:commentRangeStart", "//w:commentRangeEnd", "//w:sectPr"]:
        old_nodes = old_root.xpath(query, namespaces=NS)
        new_nodes = document.xpath(query, namespaces=NS)
        assert len(new_nodes) >= len(old_nodes)
        assert all(etree.tostring(x) == etree.tostring(y) for x, y in zip(old_nodes, new_nodes)), query
report = {"source": str(SOURCE), "target": str(TARGET),
          "source_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
          "target_sha256": sha256(TARGET.read_bytes()).hexdigest(),
          "changed_original_parts": changed, "new_media": list(new_media),
          "all_original_media_and_tables_fields_comments_and_sections_preserved": True,
          "affinity_methods": methods, "changed_paragraphs": changes,
          "references_added": references, "supplement_text": supplement,
          "render_validation": "pending"}
(ROOT / "revision_changes.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
(ROOT / "3.4新增正文.md").write_text("\n\n".join([affinity_paragraph, dynamic_paragraph, contact_paragraph]) + "\n", encoding="utf-8")
print(json.dumps({"target": str(TARGET), "changed_original_parts": changed,
                  "new_figures": len(new_media), "original_media_preserved": len(original_media)}, ensure_ascii=False))
