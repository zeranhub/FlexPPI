"""Create a manuscript integration proposal; never write or modify DOCX files."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import re
from lxml import etree

BASE = Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010")
SOURCE = BASE / "paper_v42_20261017_中文修订_蛋白背景与实例动机_三线表.docx"
OUT = BASE / "experiment_completion/figure_merge_v43_20261018/text_review"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
with ZipFile(SOURCE) as z:
    doc = etree.fromstring(z.read("word/document.xml"))
paragraphs = ["".join(p.xpath(".//w:t/text()", namespaces=NS)) for p in doc.xpath("//w:body/w:p", namespaces=NS)]


def find(prefix):
    values = [p for p in paragraphs if p.startswith(prefix)]
    if len(values) != 1:
        raise ValueError((prefix, len(values)))
    return values[0]


def swap(prefix, pairs):
    old = find(prefix)
    new = old
    for a, b in pairs:
        if a not in new:
            raise ValueError((prefix, a))
        new = new.replace(a, b)
    return {"find_paragraph_start": prefix, "old": old, "new": new}


replacements = [
    swap("侧链实例使用泛素", [("图5显示", "图3f显示"), ("图7a", "图5d"), ("图7b", "图5e"), ("图7c", "图5f")]),
    swap("完整结构编码相对仅链内编码", [("图2进一步", "图2e、f进一步")]),
    swap("同类型残基比较检验", [("（图3）", "（图2c、d）")]),
    swap("转录共激活因子CBP", [("图5a–d", "图3a–d")]),
    swap("四种KIX组成中", [("图5c、d", "图3c、d")]),
    swap("KRAS是一种", [("图6a", "图4a")]),
    swap("在全部10个肽位置上", [("图6b", "图4b")]),
    swap("去除发生替换的G6/D6后", [("图6b", "图4b")]),
    swap("正文图4采用", [("正文图4采用", "正文图5a的两类读数配对分析采用")]),
]
type_hchi_old = find("氨基酸类型均值对照的整体Hχ相关已达")
type_hchi_new = "氨基酸类型均值对照的整体Hχ相关已达0.888。完整模型在18种具有侧链标签的氨基酸中均取得低于静态特征模型的类型内平均误差，其中16种低于链内模型（图2d）；与类型均值对照相比，17种的误差更低。这些比较支持模型在固定氨基酸类型后提高逐位置侧链状态多样性的估计精度。"
replacements.append({"find_paragraph_start": "氨基酸类型均值对照的整体Hχ相关已达", "old": type_hchi_old, "new": type_hchi_new})

intro = "主链波动与侧链联合状态多样性描述同一残基的不同动态特征。全体折外预测中，RMSF与Hχ原值的参考／预测Spearman相关为0.233／0.260；减去记录内同类型均值后为0.391／0.517（图5a）。按记录内同类型分位形成低／低、低／高、高／低及高／高四类组合，模型的记录平均召回率分别为54.6%、15.8%、16.5%和61.4%（表S9）。我们将两项预测进一步对应到MD连续角度占据及主占据取向范围，检验侧链读数能否为这些构象特征提供比主链波动更合适的位置排序。"

population1 = "在83条记录的33,142个同类型有效位置上，完整模型Hχ预测与MD低相对自由能角区总宽度的记录平均相关为0.400，链内Hχ为0.391，完整模型RMSF为0.239。完整侧链输出相对主链输出的配对相关差为0.161，95%重采样区间为0.144–0.179（图5b）。各输出采用相同位置和记录内同类型控制，对照不同读数与主占据取向范围的联系。"
population2 = "按侧链预测在各记录同类型内形成高、低半组，比较两组真实低自由能角区总宽度。F*≤1时，低组与高组分别覆盖35.41和43.17°/角，差值为7.76°/角，95%区间为7.26–8.26，83条记录均为正值。链内侧链输出与完整模型主链输出确定的组差分别为7.54和5.09°/角（图5c）。位置组由原静态结构预测确定，实际MD宽度用于评价这些位置，支持侧链输出筛选主占据取向覆盖更广的同类型位点。"
population3 = "对同一组位置，完整模型Hχ、链内Hχ及完整模型RMSF与真实平均细箱角度占据熵的记录平均相关分别为0.510、0.493和0.321；完整侧链输出相对主链输出的配对相关差为0.189，95%区间为0.172–0.206（图5b）。仅按完整侧链预测确定的高、低半组真实细箱熵差为0.293 nats/角，95%区间为0.277–0.307，83条记录均为正值；链内侧链与完整主链输出的组差分别为0.281和0.198 nats/角（表S7及配套数据）。低自由能角区总宽度检验主占据范围，细箱熵检验整体占据多样性；两项结果支持侧链输出比主链输出更好地区分同类型位置的实际占据差异。"

ub1 = find("泛素是由76个氨基酸").replace("图5e、g", "图3e、g").replace("图7a", "图5d")
ub2 = find("Q41与Q62均为具有三个侧链角").replace("图5g", "图3g").replace("图5f、h", "图3f、h")
ub3 = find("在泛素各χ角层内控制氨基酸类型后").replace("图7b", "图5e").replace("图7c", "图5f")

caption5 = "图5：联合动态读数与侧链构象占据的对应。（a）3,063条记录中，RMSF与Hχ在减去记录内同类型均值后的Spearman相关分布；参考与预测均按记录等权汇总。（b）83条记录、33,142个共同有效位置上，完整模型Hχ相对其RMSF输出与MD指标的配对相关差；两项MD指标为平均细箱角度占据熵和低相对自由能角区总宽度。每点为一条记录，均值与95%区间采用按来源分层的配对重抽样。（c）分别按完整Hχ、链内Hχ和完整RMSF在记录内同类型中选择高、低半组，其真实MD低自由能角区总宽度之差，单位为°/角。（d）泛素全部66个Hχ有效残基的联合熵预测与MD参考，仅标出Q2与Q41；无可定义χ角的Gly/Ala不绘点，连线在此断开。（e）同一个模型联合Hχ读数与χ1–χ4各层MD单角细箱熵的关联，n依次为66、52、23和10；每层按氨基酸类型分别减去预测与MD指标均值。（f）预测选出的Q41与Q2的χ1–χ3实际MD自由能投影F*=−ln(p̂/p̂max)，虚线为F*=1。自由能曲线与角区宽度均由轨迹占据计算。泛素曲线使用0–100 ns、200 ps间隔的501帧；d的联合熵参考仍为原五个连续窗口均值。"

proposal = {
    "status": "proposal_ready_no_DOCX_edit",
    "source_document": str(SOURCE),
    "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "main_figure_map": {"1": "原图1原样", "2": "旧图2+3的六块能力证据", "3": "原图5原样", "4": "原图6 KRAS原样", "5": "原图4a+83条桥梁+原图7", "S1": "原图S1原样"},
    "final_panel_maps": {
        "Figure2": {"a": "RMSF逐记录MAE ECDF", "b": "Hχ逐记录MAE ECDF", "c": "RMSF类型内相对Intra/Static MAE降幅", "d": "Hχ类型内相对Intra/Static MAE降幅", "e": "RMSF按链数的Intra减Full MAE", "f": "Hχ按链数的Intra减Full MAE"},
        "Figure5": {"a": "旧图4a两类读数的类型内相关", "b": "83条记录中Full Hχ减Full RMSF对MD占据指标的配对ρ差", "c": "三个读数各自选位后的MD低F角宽高低组差", "d": "旧图7a泛素全部66位", "e": "旧图7b联合Hχ对各χ层MD占据", "f": "旧图7c Q41/Q2的χ1–χ3实际MD曲线"}
    },
    "heading_replacements": {"3.1 两类动态预测与完整结构环境的贡献": "不变", "3.2 同类型残基的局部动态差异": "不变", "3.4 KIX结构域中主链活动区域的定位": "3.3 KIX结构域中主链活动区域的定位", "3.5 KRAS G12D单点抗原变化的主链动态响应": "3.4 KRAS G12D单点抗原变化的主链动态响应", "3.3 主链与侧链读数的互补组合": "从旧位置移除，内容压缩到新3.5引入", "3.6 侧链动态预测的构象对应与群体验证": "3.5 联合读数的构象对应与群体验证", "3.6.2 侧链预测对低自由能取向范围与角度占据的验证": "3.5.1 侧链预测对低自由能取向范围与角度占据的验证", "3.6.1 泛素同类型侧链多样性的预测与轨迹对照": "3.5.2 泛素同类型侧链多样性的预测与轨迹对照"},
    "paragraph_replacements_outside_new35": replacements,
    "remove_main_text_paragraph": {"old": find("图3同时展示模型在不同动态幅度下的误差"), "reason": "原幅度诊断散点不再活动图，不保留不存在的图引。该诊断的图及源表仍随原实验材料归档；不虚构新图S2。", "optional_append_to_Supplement_C": "动态幅度分层的配套数据中，RMSF较大的记录通常具有较大的绝对误差；按参考均值归一化后，该趋势在ATLAS和DynaRepo中减弱，在ITrace中仍较明显。"},
    "new_section_35_ordered_paragraphs": [
        {"style": "section_heading", "text": "3.5 联合读数的构象对应与群体验证"},
        {"style": "body", "text": intro},
        {"style": "subsection_heading", "text": "3.5.1 侧链预测对低自由能取向范围与角度占据的验证"},
        {"style": "body", "text": population1},
        {"style": "body", "text": population2},
        {"style": "body", "text": population3},
        {"style": "subsection_heading", "text": "3.5.2 泛素同类型侧链多样性的预测与轨迹对照"},
        {"style": "body", "text": ub1},
        {"style": "body", "text": ub2},
        {"style": "body", "text": ub3}
    ],
    "figure5_caption_proposal_must_match_final_plot": caption5,
    "legacy_caption_mapping": {"图2：五折折外预测及输入链数分层。": "由Figure2 worker的新六panel图注替换", "图3：动态幅度与预测误差的关系及同类残基内的结构信息收益。": "删除旧活动图与旧图注；相关数据归档", "图4：主链波动与侧链状态多样性的互补性及组合恢复。": "删除旧活动图与旧图注，四类完整结果保留表S9", "图5：静态结构的局部动态验证。": "仅将开头图5改为图3；其余图注原样", "图6：KRAS G12D单点替换后的肽段主链动态响应预测。": "仅将开头图6改为图4；其余图注原样", "图7：联合熵预测、逐角占据关联与MD自由能。": "由Figure5 worker的a–f新图注替换，旧d/e/f细则保留在新图注及既有方法"},
    "scientific_guards_for_editor_only": [
        "图5a为两项读数之间的相关，不是预测对参考的准确性，不从预测ρ高于参考推导互补信息保真。",
        "图5b显示对物理占据指标的相关差，不能称作控制RMSF后的条件独立贡献或独立信息增量。",
        "图5c各高低组只由各自模型预测选择，MD仅用于评价选位，不颠倒选择方向。",
        "模型不输出自由能曲线或低F角宽；f的曲线和b/c的MD指标均由真实轨迹计算。",
        "图5e使用同一个逐残基联合Hχ读数，不是χ1、χ2、χ3、χ4各有一个预测头。",
        "83为记录数而非249个独立重复；主分析为每记录第一条共同窗口，不能混available_replicas结果。",
        "图2c/d只有Intra和Static类型内收益，不声称其绘出了AA对照；17/18 AA改善事实保留正文与原源表。",
        "表S9四类完整召回54.6/15.8/16.5/61.4%均保留；不把组合存在、相关不为1或低召回写成联合筛位准确性的强证据。"
    ],
    "must_preserve_structurally": {"native_tables": 10, "native_OMML_equations": 36, "references": 28, "comment_parts_and_anchors": "原样", "unchanged_visible_media": ["word/media/image1.png", "word/media/image5.png", "word/media/image41.png", "word/media/image7.png"], "final_visible_figures": 6},
    "numbering_checks": ["最终仅图1–5及图S1活动图，正文无图6/图7引用。", "结果顺序为3.1、3.2、3.3KIX、3.4KRAS、3.5及其3.5.1群体/3.5.2泛素。", "补充G原正文图4改为正文图5a，不能将整张新图5都说成824,988位置。", "方法2.8不变，方法2.7仅改泛素图引用。"]
}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "manuscript_integration_proposal.json").write_text(json.dumps(proposal, ensure_ascii=False, indent=2), encoding="utf-8")
print(OUT / "manuscript_integration_proposal.json")
