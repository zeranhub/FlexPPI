"""Bind independent numerical and narrative review to the final v44 DOCX."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
from lxml import etree

ROOT = Path("E:/Project/FlexPPI/experiments/_sources/predictive_figure_v44_20261018")
DOC = ROOT.parent.parent / "paper_v44_20261018_中文修订_预测一致性与高动态识别_三线表.docx"
DOC_SHA = "4156791ba646d8137754386ea387b920aaaf1bdac041614414a17028915ffcab"
FIG = ROOT / "figure2/fig2_predictive_capability.png"
FIG_SHA = "ac1466e072178be049d29671693060c9c6eba7118f500f259312a2732872ebf0"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    package = json.loads((ROOT / "manuscript_QA/independent_package_audit.json").read_text(encoding="utf-8"))
    metrics = json.loads((ROOT / "text_review/independent_metrics_audit.json").read_text(encoding="utf-8"))
    proposal = json.loads((ROOT / "text_review/manuscript_proposal.json").read_text(encoding="utf-8"))
    summary = json.loads((ROOT / "analysis/summary.json").read_text(encoding="utf-8"))
    with ZipFile(DOC) as z:
        tree = etree.fromstring(z.read("word/document.xml"))
    paragraphs = ["".join(p.xpath(".//w:t/text()", namespaces=NS)) for p in tree.xpath("//w:body/w:p", namespaces=NS)]
    find = lambda prefix: next(p for p in paragraphs if p.startswith(prefix))
    method = find("高动态位点筛选以每条记录MD参考")
    result = find("进一步将连续预测用于高动态位点筛选")
    within = find("在同一记录、同一种氨基酸内重新定义高熵位置后")
    fig2 = find("图2：")
    new_texts = [x["new"] for x in proposal["replacements"]] + [x["text"] for x in proposal["insertions"]] + [fig2]
    rmsf, hchi = summary["metrics"]["rmsf"], summary["metrics"]["chi_entropy"]
    main_values = [f"{m[method_name]['mean_ap']:.3f}" for m in [rmsf, hchi] for method_name in ["Full", "Intra", "Static", "AA"]]
    main_values += [f"{rmsf['Full']['mean_auroc']:.3f}", f"{hchi['Full']['mean_auroc']:.3f}"]
    ci_values = [f"{v:.3f}" for m in [rmsf, hchi] for v in m["Full"]["ap_ci"]]
    type_values = [f"{summary['within_type'][m]['mean_ap']:.3f}" for m in ["Full", "Intra", "Static", "AA"]]
    type_values += [f"{summary['within_type']['Full']['mean_auroc']:.3f}"]
    checks = {
        "final_document_SHA_bound": sha(DOC) == DOC_SHA,
        "final_Figure2_SHA_bound": sha(FIG) == FIG_SHA,
        "independent_package_audit_passed_same_DOCX": package["status"] == "passed" and package["document_sha256"] == DOC_SHA,
        "independent_metrics_audit_passed_final_summary": metrics["status"] == "pass" and metrics["summary_sha256"] == sha(ROOT / "analysis/summary.json") and metrics["protocol_sha256"] == sha(ROOT / "analysis/PRESPECIFIED_PROTOCOL.md"),
        "all_new_reviewed_text_in_DOCX_once": all(paragraphs.count(s) == 1 for s in new_texts),
        "q90_reference_threshold_and_ties_explicit": "第90百分位及以上为阳性" in method and "并列参考全部纳入" in method and "并列预测作为同一阈值组" in method,
        "AP_step_and_AUC_pair_definitions_correct": "召回率增量乘对应精确率求和" in method and "所有阳性–阴性位点对" in method and "并列计半分" in method,
        "mean_AP_and_macro_PR_area_distinguished": "曲线面积不代替平均AP" in method and "逐记录等权平均" in method,
        "record_source_bootstrap_unit_and_interval_scope_correct": "两指标的95%区间按来源分层重抽记录10,000次" in method,
        "all_main_AP_AUC_values_match_computed_results": all(s in result for s in main_values),
        "Full_AP_CIs_match_computed_results": all(s in result for s in ci_values),
        "Full_AP_higher_than_all_three_controls_both_targets": all(m["Full"]["mean_ap"] > m[k]["mean_ap"] for m in [rmsf, hchi] for k in ["Intra", "Static", "AA"]),
        "within_type_auxiliary_task_separate_and_aggregation_defined": "Hχ另在同记录同类型至少4个位点内重定第90百分位" in method and "先对有效类型、再对记录等权平均" in method and "同一记录、同一种氨基酸内重新定义高熵位置" in within,
        "within_type_AP_AUC_and_coverage_match": all(s in within for s in type_values) and "3,063" in within and "42,108" in within,
        "AA_constant_score_not_mislabeled_random_AP": not any("随机AP" in s or "随机分类" in s for s in new_texts),
        "density_original_predictions_not_recalibrated": not any("校准后" in s or "校准后的" in s for s in new_texts),
        "density_record_weighting_and_predicted_binned_mean_explicit": "记录等权" in fig2 and "按MD参考值分箱的预测均值" in fig2,
        "density_CI_1000_record_bootstraps_correct": "1,000次来源分层重抽记录的95%区间" in fig2,
        "PR_complete_selection_endpoint_explicit": "全部位点入选的终点" in fig2,
        "no_accuracy_ACC_substitution_for_regression": not any("准确率" in s or "ACC" in s for s in new_texts),
        "new_type_MAE_panel_references_correct": "图2e、f" in find("同类型残基比较检验模型") and "图2f" in find("氨基酸类型均值对照的整体Hχ相关已达") and "图2e" in find("RMSF方面，完整模型在20种氨基酸类型"),
        "existing_RMSF_comparator_strength_retained": "低于FlexPPI的0.577" in find("在已有RMSF方法中，BackFlip"),
        "all_other_figures_cases_tables_and_background_preserved": package["checks"]["all_other_body_nodes_exactly_preserved"] and package["checks"]["six_images_same_order_and_Fig5_unchanged"],
    }
    report = {
        "status": "passed" if all(checks.values()) else "failed",
        "document": str(DOC), "document_sha256": sha(DOC),
        "figure2": str(FIG), "figure2_sha256": sha(FIG),
        "metrics_summary_sha256": sha(ROOT / "analysis/summary.json"),
        "independent_numeric_audit": str(ROOT / "text_review/independent_metrics_audit.json"),
        "checks": checks,
        "manual_review": {
            "Figure2_final_PNG_inspected": True,
            "numerical_agreement_panels": "All retained values displayed as record-equal density; RMSF high-reference underestimation is visible and is not described as calibrated accuracy.",
            "ranking_panels": "Full mean AP and AUROC match original OOF rankings; q90 prevalence and full-selection endpoints are retained; no best-performing record selection.",
            "type_panels": "Original type MAE gains retained, including negative Full vs Intra Hchi values and Ala/Gly no-angle status.",
            "within_type_narrative": "Independent within-record/type q90 AP demonstrates structure-dependent ranking beyond type prior without calling a constant-score baseline a random finite-sample classifier.",
            "scope": "No new training, MD, inference, prediction calibration, threshold optimization or test-set model selection.",
        },
        "page_layout_QA": "Handled separately by root/assigned page reviewers; this report does not claim page-level layout completion.",
    }
    out = ROOT / "manuscript_QA/independent_scientific_audit.json"
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "document_sha256": sha(DOC), "failed_checks": [k for k, v in checks.items() if not v]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
