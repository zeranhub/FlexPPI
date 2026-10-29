"""Adapt the previously audited export QA without inheriting old visual claims."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
source = (ROOT.parent / "qa_current.py").read_text(encoding="utf-8")


def replace_once(old, new):
    global source
    assert source.count(old) == 1, old[:100]
    source = source.replace(old, new)


replace_once('GROUPS = {"single_chain": 2193, "other_multichain": 626, "tcr_pmhc_itrace": 244}\nRESIDUES = {"single_chain": 497407, "other_multichain": 264246, "tcr_pmhc_itrace": 199668}',
    'BINARY_SCORE_AUDIT = json.loads((ROOT / "evaluation/scoring_audit.json").read_text(encoding="utf-8"))\nGROUPS = BINARY_SCORE_AUDIT["group_counts"]\nRESIDUES = BINARY_SCORE_AUDIT["group_residues"]')
start = source.index('VISUAL_REVIEWS = {')
end = source.index('\n\n\ndef sha256_bytes', start)
source = source[:start] + 'VISUAL_PATH = ROOT / "visual_reviews.json"\nVISUAL_REVIEWS = json.loads(VISUAL_PATH.read_text(encoding="utf-8")) if VISUAL_PATH.exists() else {}' + source[end:]
replace_once('        ROOT / "references/training_exposure_protocol_sources_sha256.json",',
    '        ROOT / "references/training_exposure_protocol_sources_sha256.json",\n        ROOT / "classification/membership_frozen.csv",\n        ROOT / "classification/membership_freeze_audit.json",\n        ROOT / "classification/policy_frozen.json",')
replace_once('        "training_exposure_audit_sha256": hashes[str(ROOT / "references/training_exposure_protocol_sources_sha256.json")],',
    '        "training_exposure_audit_sha256": hashes[str(ROOT / "references/training_exposure_protocol_sources_sha256.json")],\n        "membership_sha256": hashes[str(ROOT / "classification/membership_frozen.csv")],')
old = '''    tcr_rows = [row for row in points if row["group"] == "tcr_pmhc_itrace"]
    check("TCR_SINGLE_FAMILY_NO_INDEPENDENT_CI", all(int(row["n_families"]) == 1 and not math.isfinite(number(row["ci_low"])) and not math.isfinite(number(row["ci_high"])) for row in tcr_rows),
          "TCR-pMHC has one family component; current figure correctly uses descriptive points without independent-family intervals.")'''
new = '''    single_family_rows = [row for row in points if int(row["n_families"]) == 1]
    check("ANY_SINGLE_FAMILY_NO_INDEPENDENT_CI", all(not math.isfinite(number(row["ci_low"])) and not math.isfinite(number(row["ci_high"])) for row in single_family_rows),
          "Any single-family group must have no independent-family CI. TCR records remain a dependency family inside the larger complex group.")
    membership = read_csv_bytes(inputs[str(ROOT / "classification/membership_frozen.csv")])
    freeze = json.loads(inputs[str(ROOT / "classification/membership_freeze_audit.json")].decode("utf-8-sig"))
    check("BIOLOGICAL_MEMBERSHIP_FROZEN_ZERO_UNRESOLVED", len(membership) == len({row["sample_id"] for row in membership}) == 3063
          and dict(Counter(row["benchmark_group"] for row in membership)) == GROUPS
          and freeze["status"] == "frozen" and freeze["unresolved_records"] == 0
          and freeze["membership_sha256"] == hashes[str(ROOT / "classification/membership_frozen.csv")] == score_audit["membership_sha256"]
          and freeze["policy_sha256"] == hashes[str(ROOT / "classification/policy_frozen.json")],
          {"groups": GROUPS, "membership_sha256": score_audit["membership_sha256"], "unresolved_records": freeze["unresolved_records"]})
    check("PREDICTIONS_TARGETS_AND_SCORES_UNCHANGED", score_audit["predictions_recomputed"] is False
          and score_audit["targets_or_masks_changed"] is False and score_audit["precomputed_scores_numerically_unchanged"] is True,
          "Binary biological grouping is a metadata and aggregation update, not a model rerun.")
    check("ALL_RECORDS_INCLUDED_IN_ECDF", plot_audit.get("all_observations_in_ecdf") is True,
          "Plot runtime verifies each selected method uses the complete group, with no sampled or clipped observations.")'''
replace_once(old, new)
replace_once('f"数据检查：主图 {len(plotted)} 个完整方法/配置，三个互斥组为 2,193 / 626 / 244 条，总计 3,063 条；source.csv 共 {len(points)} 行。核对了全部方法×组×适用指标、有效残基覆盖、summary 数值和完成状态；相对分数方法没有被写成绝对 Å MAE。未将 incomplete 方法放入主图。"',
    'f"数据检查：主图 {len(plotted)} 个完整方法/配置，两组为 {GROUPS}，总计 3,063 条；source.csv 共 {len(points)} 行。逐条分类全部明确，核对全部方法×组×适用指标、有效残基覆盖和summary数值。相对分数未写成绝对 Å MAE。"')
source = source.replace('Tight bounding-box export slightly changes page width from nominal 183 mm.', 'Canvas-bound exports preserve the 183 mm physical width; actual PDF and PNG dimensions are measured.')
source = source.replace('名义 183 mm 经 tight crop 后的实际差异已记录，未隐藏。', '画布固定183 mm，并独立测量PDF与PNG实际尺寸。')
source = source.replace("Exported PDF/PNG widths are separately measured because bbox_inches='tight' changes the physical crop.", "Exported PDF/PNG widths are separately measured; the new script exports the fixed canvas boundary.")
replace_once('def safe_numeric(node):\n', 'SCRIPT_NUMERICS = {}\n\n\ndef safe_numeric(node):\n    if isinstance(node, ast.Name) and node.id in SCRIPT_NUMERICS:\n        return float(SCRIPT_NUMERICS[node.id])\n    if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name):\n        values = SCRIPT_NUMERICS.get(node.value.id)\n        if isinstance(values, (tuple, list)):\n            return float(values[int(safe_numeric(node.slice))])\n')
replace_once('    widths = []\n', '    SCRIPT_NUMERICS.update({name: item["value"] for name, item in assignments.items()})\n    widths = []\n')
(ROOT / "qa_current.py").write_text(source, encoding="utf-8")
print("Binary QA prepared; no legacy visual hashes reused.")
