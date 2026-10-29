"""Validate exported ablation artifacts and assemble the Chinese handoff report."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET

import fitz
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/information_ablation_figures/verified"
AUDIT = ROOT / "output/information_ablation_analysis/verified_audit"
FIGURES = {
    "fig_neural_ablation_main": "主图：四组神经网络消融与配对差值",
    "fig_neural_ablation_distributions": "补充图：所有蛋白的配对误差分布",
    "fig_neural_ablation_heterogeneity": "补充图：数据来源与五折方向一致性",
    "fig_neural_ablation_dependence": "补充图：相关蛋白整组重采样的区间敏感性",
    "fig_ablation_all_controls": "补充图：全部九种方法的比较背景",
}
NAMES = {"flexppi": "完整模型", "nn_sequence": "仅序列", "nn_geometry": "隐藏残基身份", "nn_intrachain": "仅链内"}

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def link(path, label):
    return f"[{label}]({Path(path).resolve().as_posix()})"

def main():
    exports = []
    for name in FIGURES:
        paths = [OUT / f"{name}.{ext}" for ext in ("pdf", "svg", "png")]
        assert all(p.exists() and p.stat().st_size > 1000 for p in paths)
        with fitz.open(paths[0]) as document:
            assert len(document) == 1
            page = document[0]
            width_mm = page.rect.width * 25.4 / 72
            height_mm = page.rect.height * 25.4 / 72
            assert abs(width_mm - 183) < .01 and height_mm < 247
            spans = [s for b in page.get_text("dict")["blocks"] if "lines" in b
                     for line in b["lines"] for s in line["spans"]]
            minimum_font = min(s["size"] for s in spans)
            assert minimum_font >= 5
            clipping = [s["text"] for s in spans if s["bbox"][0] < -.1 or s["bbox"][1] < -.1
                        or s["bbox"][2] > page.rect.width+.1 or s["bbox"][3] > page.rect.height+.1]
            assert not clipping, (name, clipping)
            fonts = page.get_fonts()
            assert all(len(document.extract_font(f[0])[3]) > 0 for f in fonts)
            font_names = [f[3] for f in fonts]
        ET.parse(paths[1])
        with Image.open(paths[2]) as img:
            dpi = img.info["dpi"]
            assert min(dpi) > 299.9
            content_density = float((np.asarray(img.convert("RGB")).min(axis=2) < 245).mean())
            assert .015 < content_density < .5
            png_size = img.size
        exports.append(dict(figure=name,width_mm=width_mm,height_mm=height_mm,
                            min_pdf_font_pt=minimum_font,all_fonts_embedded=True,font_names=font_names,
                            clipped_text=clipping,png_dpi=dpi,png_pixels=png_size,
                            content_density=content_density,sha256={p.suffix:digest(p) for p in paths}))
    source=pd.read_csv(OUT/"source_data_all_proteins.csv")
    original=pd.read_csv(ROOT/"output/information_ablation_figures/source_data_protein_error.csv")
    pd.testing.assert_frame_equal(source,original,check_exact=False,rtol=1e-12,atol=1e-12)
    assert not source.duplicated(["protein_id","metric"]).any()
    statistics=pd.read_csv(OUT/"performance_protein_equal.csv")
    contrasts=pd.read_csv(OUT/"paired_ablation_effects.csv")
    dependence=pd.read_csv(OUT/"source_data_dependence_intervals.csv")
    assert len(statistics)==18 and len(contrasts)==6 and len(dependence)==6
    for row in contrasts.itertuples():
        block=source[source.metric.eq(row.metric)]
        assert abs((block[row.method]-block.flexppi).mean()-row.mean_delta)<1e-12
    independent=json.loads((AUDIT/"audit_independent.json").read_text(encoding="utf-8"))
    qa={"status":"PASS_WITH_INTERPRETATION_LIMITS","delivery_date":"2026-10-01",
        "export_checks":exports,"source_rows_retained":len(source),"source_rows_dropped":0,
        "source_copy_matches_original":True,"all_18_method_means_present":True,
        "paired_difference_direction_verified":True,
        "review":"All exported previews were visually inspected; full distribution tails retained; no overlapping labels in final main figure.",
        "style_validator_manual_resolution":{
            "AP-0":"False alarm from exact whitespace match for DIVERGING. Required colors are present verbatim.",
            "CL-1":"False alarm: the regex also matches capsize=2. Actual PDF minimum font is >=5.95 pt.",
            "CL-2":"Validator does not parse 183/25.4 expressions. All actual PDF page widths are exactly 183 mm."},
        "statistical_scope":"Descriptive effect estimation conditional on one fitted seed; marginal 95% intervals, no formal p-value tests.",
        "scientific_limits":["142 exact protein identities span folds; validation is not strictly protein-disjoint",
                             "Sequence-only capacity is 61.26% of full",
                             "Local ablation checkpoint provenance cannot be independently rechecked"]}
    (OUT/"export_QA.json").write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding="utf-8")
    qa_text="""# 图表质量核验

交付检查日期：2026-10-01。图表采用 Python / Matplotlib，统一宽度 183 mm；每张同时交付含嵌入字体的 PDF、可编辑文字的 SVG 与 300 dpi PNG。

Pass 0（样式与反模式）通过：使用明确语义的蓝／橙／紫／灰配色，移除多余边框，采用直接标注与图外图例，没有彩虹色或默认分类色。主图样本量为数千个蛋白，误差条由全部观测计算；补充累积分布保留了每一个蛋白及全部极端值。

Pass 1（导出与字体）通过：读取实际 PDF 验证五图均为 183 mm 宽，高度不超过 247 mm；嵌入字体可提取，实际最小文字约 5.95 pt，无文字边界超过页面。PNG 元数据为 300 dpi，SVG XML 可解析。skill 的字符串检查器误把 capsize=2 识别为小字号，并因 DIVERGING 常量的空格报缺失，已通过实际文件与代码核对排除；183/25.4 表达式也由实际 PDF 尺寸核验。

Pass 2（证据与数据）通过：全部 7,737 条蛋白×指标观测与原表对齐，全部九方法的 18 个均值保留，六个差值重新核对为“消融 − 完整模型”。主图描述平均误差及配对效应，分布图保留个体差异，来源／逐折图检查异质性，依赖敏感性图说明重采样单位对区间的影响。没有用挑选子集或修改数据增强效果。

Pass 3（实际视觉检查）通过：逐张打开 PNG 核对文字、数值、图例、坐标与颜色。已修复主图 Full wins 百分比与差值标签重叠；图注明确分布图的 ±0.01 线性区域、同一训练种子、逐折结果不是独立重复及来源分组的解释限制。总像素内容比例检查通过，但最终判断以实际观察为准。

核验结论为“图表通过，科学解释需保留已有局限”：这里的区间描述给定模型下测试样本变化，不代表多随机种子重复训练稳定性；142 个相同序列蛋白跨折，不能据此声称严格未见蛋白泛化；仅序列组还改变模型架构和容量。整组重采样扩大区间而保持点估计，不能修复数据划分或训练来源证据缺口。没有生成显著性星号或将边际区间解释为多重比较校正。
"""
    (OUT/"figure_QA.md").write_text(qa_text,encoding="utf-8")
    lines=["# 消融实验核验与图表交付", "", "交付日期：2026-10-01。已复核三组重新训练的消融结果，完整模型与各组的逐残基预测在记录、标签和有效掩码上对齐；原统计结果可以独立复算。建议以主图呈现四组神经网络的比较，并将其他四图作为补充证据。现有结果支持不同信息限制产生不同误差表现，尤其是隐藏残基身份后侧链熵误差增大；RMSF 的身份消融则没有稳定方向。", "", "| 文件 | 用途 | 导出 |", "|---|---|---|"]
    for name,label in FIGURES.items():
        lines.append(f"| {label} | {'建议正文使用' if name.endswith('main') else '补充材料'} | " + " · ".join(link(OUT/f"{name}.{ext}",ext.upper()) for ext in ("pdf","svg","png")) + " |")
    lines.extend(["", "三组各有 3,068 条记录，五折记录数为 613、613、614、614、614；RMSF 和侧链熵有效残基分别为 959,347 与 798,385。没有发现重复残基或缺失有效预测，各组文件哈希符合现有记录。全部九方法的 7,737 个蛋白／指标组合独立重建后，逐蛋白误差最大差异约 7.33×10⁻⁸；原均值、配对增益与区间差异均在数值精度范围内。RMSF 用 3,870 个精确序列蛋白身份汇总，侧链熵用 3,867 个身份，因为三个身份没有有效侧链熵标签；这是沿用共同有效掩码，并非人为筛除差结果。", "", "平均绝对误差（MAE）表示预测值与参考值的距离，越小越好。RMSF 是残基位置随时间摆动的尺度，单位 Å；侧链熵 Hχ 描述侧链转动状态的分散程度，使用自然对数，单位 nats。各蛋白先汇总链内误差，再平衡同一记录中的链拷贝与同一序列的不同记录，最后给每个蛋白相同权重，避免长蛋白靠残基数量主导结论。", "", "| 组别 | RMSF MAE / Å | Hχ MAE / nats | 参数量 |", "|---|---:|---:|---:|"])
    for method in ["flexppi","nn_geometry","nn_intrachain","nn_sequence"]:
        vals=statistics[statistics.method.eq(method)].set_index("metric")
        lines.append(f"| {NAMES[method]} | {vals.loc['rmsf','mean_mae']:.6f} | {vals.loc['chi_entropy','mean_mae']:.6f} | {'836,676' if method=='nn_sequence' else '1,365,812'} |")
    lines.extend(["", "主图差值统一定义为 ΔMAE = MAE(消融组) − MAE(完整模型)，正值表示完整模型更好。隐藏身份的 RMSF 差值为 −0.001390 Å，95% 蛋白重采样区间跨零，不能描述成移除身份后 RMSF 恶化。仅链内组的 RMSF / Hχ 差值分别为 +0.089964 Å / +0.009751 nats；仅序列组分别为 +0.141971 Å / +0.052445 nats。隐藏身份对 Hχ 的差值为 +0.443942 nats，但由于残基类型本身决定侧链转动角的数量，这个差距也包含类型依赖，不能全部解释成动态信息丢失。完整配对区间与胜出比例见 " + link(OUT/"paired_ablation_effects.csv","配对效应表") + "。", "", "神经网络三组消融的配置与完整模型仅在组别和输出目录上不同，均采用 seed=7、相同五折与 40 epoch 预算；完整模型五折训练日志达到 epoch 40。仅序列模型的注意力架构和容量同时改变，参数是完整模型的 61.26%，因此不能把这组差值当作几何信息的纯独立贡献。隐藏身份仍保留结构中与残基类型相关的线索，仅链内组则去掉跨链消息。现有隔离测试首次有一个测试排列错误地改变链内顺序，经最小修复后 14/14 通过；未改模型或训练结果。", "", "划分按记录、生物系统与现有 cluster_id 隔离，但有 142 个相同蛋白序列身份跨折，其中 97 个长度至少 20。归档划分与冻结导出哈希一致；各折训练／测试中 sequence_clusters 也有重叠。因此这些结果应称为记录层面的折外评估，不能表述成严格未见序列或同源家族泛化测试。本地没有消融 checkpoint、原预测 manifest 和服务端逐折复现日志，本次核验确认的是现有预测表及已记录的复现状态，没有重新执行消融 checkpoint 推理。", "", "为检查相关蛋白对区间的影响，把共享相同序列的生物分组连接成 2,035 个单位，整组重采样，每次保留抽中组内全部蛋白并以总差值除以总蛋白数，保持蛋白等权。最大单位包含 625 个蛋白，约占 16.2%；部分区间明显扩大，尤其 Hχ 的身份消融。五项正向效应的整组区间仍高于零，RMSF 身份效应仍跨零；这里是条件于现有模型的描述性敏感性分析，没有正式显著性检验，也不能修复序列重叠或代替多种子训练。", "", "| 指标 | 组别 | 平均 ΔMAE | 蛋白重采样 95% 区间 | 整组重采样 95% 区间 |", "|---|---|---:|---|---|"])
    for row in dependence.itertuples():
        metric="RMSF / Å" if row.metric=="rmsf" else "Hχ / nats"
        lines.append(f"| {metric} | {NAMES[row.method]} | {row.mean_delta:+.6f} | [{row.protein_ci_low:.6f}, {row.protein_ci_high:.6f}] | [{row.cluster_ci_low:.6f}, {row.cluster_ci_high:.6f}] |")
    lines.extend(["", "来源分组保留 ATLAS、DynaRepo、iTrace 及两种混合来源身份，未随意把跨来源蛋白归到某一个来源。隐藏身份对 RMSF 在三个主要来源上方向较小且不同，Hχ 则均有较大增幅。仅链内组的 RMSF 差值在多链组为 +0.165640 Å、单链组为 +0.014504 Å，但 ATLAS 全为单链而 iTrace 全为多链，来源与链数存在混杂，这只能描述为异质性，不能当作因果交互。逐折图对完整五折的全部有效残基计算 MAE，所以采用残基权重，与主图蛋白等权的数值不同；它显示方向一致性，不是五次独立重复训练。", "", "所有图均为 183 mm 宽，PDF 字体嵌入，SVG 保持文字可编辑，PNG 为 300 dpi。主图展示均值及 5,000 次配对蛋白重采样的边际 95% 区间，分布图不裁掉极端蛋白，来源与依赖补充图也保留全部对应观测。五图实际视觉与导出核验见 " + link(OUT/"figure_QA.md","质量报告") + "，详细元数据见 " + link(OUT/"statistics_and_reproducibility.json","统计与复现记录") + "。", "", "独立数据核验见 " + link(AUDIT/"audit_independent.md","数据核验报告") + "，实现和训练设计见 " + link(AUDIT/"design_audit.md","实验设计报告") + "，整组重采样方法见 " + link(AUDIT/"cluster_bootstrap_sensitivity.md","敏感性分析报告") + "。图的每个统计量可追溯到 " + link(OUT/"source_data_all_proteins.csv","全部蛋白源数据") + " 与各补充 source_data CSV。", "", "复现绘图：", "", "```powershell", "python -X utf8 scripts/information_sources/make_verified_ablation_figures.py", "python -X utf8 scripts/information_sources/finalize_ablation_deliverables.py", "```", "", "绘图依赖原结果表及本次核验产生的分组／逐折／整组重采样 CSV；首次计算整组敏感性可运行 output/information_ablation_analysis/verified_audit/cluster_bootstrap_sensitivity.py。随机种子与重采样次数已记录。", "", "建议图注：四组 FlexPPI 在相同有效残基上进行五折记录层面折外评估。a–b 为蛋白等权平均绝对误差，c–d 为消融相对完整模型的配对误差差值；正值表示完整模型误差更低。误差条为 5,000 次蛋白重采样的边际 95% 区间，Full wins 为完整模型在蛋白层面误差更低的比例。每组使用同一训练随机种子和 40 epoch 预算；仅序列组容量为完整模型的 61.26%，相同序列可能跨折，因此区间不代表多种子或严格未见蛋白泛化的稳定性。"])
    (OUT/"README_核验与图表说明.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    report_path=OUT/"statistics_and_reproducibility.json"
    report=json.loads(report_path.read_text(encoding="utf-8"))
    report.update(status="verified_visual_and_export_QA_with_interpretation_limits",delivery_date="2026-10-01",
                  final_QA="export_QA.json",n_proteins_spanning_folds=142,
                  connected_group_bootstrap={"n_units":2035,"largest_unit_proteins":625,
                                            "interval_file":"source_data_dependence_intervals.csv","resampling":"paired whole connected groups, retaining equal-protein estimator"})
    report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(f"Validated {len(exports)} figures; retained {len(source)} protein/metric rows; Chinese delivery report saved.")

if __name__=="__main__":
    main()
