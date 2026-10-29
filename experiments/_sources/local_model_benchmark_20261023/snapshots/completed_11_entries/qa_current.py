"""Audit current saved benchmark figures; never load residue-level score arrays.

Uses the installed Academic Figure Skill validator programmatically because its
file has no CLI entry point. Original findings are preserved, then supplemented
by AST, exported-file, and small-table checks. Visual observations are bound to
the exact two PNG hashes actually inspected; a future render is marked pending.
"""
import ast
from collections import Counter
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import re
import sys

from PIL import Image
import pdfplumber
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
VALIDATOR = Path("C:/Users/DELL/.codex/skills/academic-figure-skill/scripts/qa_validator.py")
GROUPS = {"single_chain": 2193, "other_multichain": 626, "tcr_pmhc_itrace": 244}
RESIDUES = {"single_chain": 497407, "other_multichain": 264246, "tcr_pmhc_itrace": 199668}
NOMINAL_WIDTH_MM = 183.0

# These observations were made by viewing the actual PNGs, not predicted from
# source code. The script deliberately does not apply them to changed images.
VISUAL_REVIEWS = {
    "benchmark_current.png": {
        "sha256": "501fde5cf189547e565acbaa98bde08382d065f95fa9170888f56437df814572",
        "reviewed_utc": "2026-10-23T16:34:54+00:00",
        "status": "inspected_current_hash_with_notes",
        "observations": [
            "Viewed the actual current 3-column by 3-row dot/interval figure. Eleven method/configuration rows are present in consistent order, including assembly and chains versions of both BackFlip-2.1 and BackFlip-1.",
            "Blue/ochre/teal semantic colors are visible. Direct method labels and marker shapes supplement color.",
            "No title, row-label, axis-label, or external-legend overlap/clipping was observed at the inspected image resolution.",
            "ANM/GNM absolute-MAE cells visibly say relative score; TCR-pMHC points are shown without independent-family confidence intervals.",
            "The longer BackFlip row labels fit without observed collision or clipping. Some intervals are smaller than their markers; this is not evidence that intervals are absent.",
            "This is an eleven-completed-entry progress snapshot. Partial FastProtFlex and other unfinished models are absent from the primary figure, rather than being compared on smaller surviving subsets.",
        ],
    },
    "benchmark_distributions_current.png": {
        "sha256": "c0e4de60093aee48d75e2155b0c068fdb119ae75d5877f559c364e3f4fbf26b3",
        "reviewed_utc": "2026-10-23T16:34:54+00:00",
        "status": "inspected_current_hash_with_notes",
        "observations": [
            "Viewed the actual current three-panel empirical-CDF figure. Titles, axes, and the external bottom legend are visible without observed clipping.",
            "The five displayed methods are FlexPPI, BackFlip-2.1 assembly, BackFlip-2.1 chains, DynaProt-M, and ANM. The five-entry, two-row external legend fits without observed clipping.",
            "The learned-baseline curves share an ochre family color but have solid/dashed/dotted styles visible in the actual figure. Direct legend labels identify each; ANM is teal with a distinct dash pattern.",
            "Single-chain curves partly overlap; no local zoom or nonlinear axis is used. No color-vision simulation or paper printout was inspected.",
        ],
    },
}


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256(path):
    return sha256_bytes(Path(path).read_bytes())


def read_csv_bytes(data):
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))


def true(value):
    return str(value).lower() == "true"


def number(value):
    return float(value) if str(value).strip() else float("nan")


def safe_numeric(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        value = safe_numeric(node.operand)
        return value if isinstance(node.op, ast.UAdd) else -value
    if isinstance(node, ast.BinOp):
        left, right = safe_numeric(node.left), safe_numeric(node.right)
        if isinstance(node.op, ast.Div):
            return left / right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
    raise ValueError("Not a literal arithmetic expression")


def imported_skill_findings(source):
    spec = importlib.util.spec_from_file_location("academic_figure_qa_validator", VALIDATOR)
    qa_validator = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = qa_validator
    spec.loader.exec_module(qa_validator)
    # Explicitly call the shipped check functions; simply executing that module
    # produces no report because the installed file has no __main__ entry point.
    findings = qa_validator.check_ap0_style_baseline(source)
    for check in (
        qa_validator.check_ap1_default_palette, qa_validator.check_ap2_jet_rainbow,
        qa_validator.check_ap3_four_sided_borders, qa_validator.check_ap4_legend_occlusion,
        qa_validator.check_ap5_low_res_export, qa_validator.check_ap6_missing_points,
        qa_validator.check_ap7_default_font, qa_validator.check_cl1_fontsize,
        qa_validator.check_cl2_dimensions, qa_validator.check_cl3_dpi,
        qa_validator.check_cl4_font_embedding, qa_validator.check_cl5_spine_linewidth,
        qa_validator.check_cl6_tick_direction, qa_validator.check_cl7_export_completeness,
    ):
        findings.append(check(source))
    return [asdict(finding) for finding in findings]


def pdf_font_records(reader):
    rows = []
    for page_index, page in enumerate(reader.pages, 1):
        fonts = page.get("/Resources", {}).get("/Font", {})
        fonts = fonts.get_object() if hasattr(fonts, "get_object") else fonts
        for resource_name, reference in fonts.items():
            font = reference.get_object()
            descendants = font.get("/DescendantFonts")
            effective = descendants[0].get_object() if descendants else font
            descriptor = effective.get("/FontDescriptor")
            descriptor = descriptor.get_object() if descriptor else {}
            streams = [name for name in ("/FontFile", "/FontFile2", "/FontFile3") if name in descriptor]
            rows.append({"page": page_index, "resource": str(resource_name),
                         "basefont": str(font.get("/BaseFont", effective.get("/BaseFont", ""))),
                         "subtype": str(font.get("/Subtype", "")), "descendant_subtype": str(effective.get("/Subtype", "")),
                         "embedded": bool(streams), "embedded_stream_keys": streams,
                         "embedded_bytes": sum(len(descriptor[name].get_object().get_data()) for name in streams)})
    return rows


def main():
    input_paths = [
        ROOT / "plot_benchmark.py", FIGURES / "benchmark_current_source.csv", FIGURES / "plot_audit.json",
        ROOT / "evaluation/scoring_audit.json", ROOT / "evaluation/group_summary.csv",
        ROOT / "evaluation/coverage_summary.csv",
    ]
    input_paths += [FIGURES / (stem + extension) for stem in ("benchmark_current", "benchmark_distributions_current") for extension in (".png", ".pdf", ".svg")]
    inputs = {str(path): path.read_bytes() for path in input_paths}
    hashes = {path: sha256_bytes(data) for path, data in inputs.items()}
    source = inputs[str(ROOT / "plot_benchmark.py")].decode("utf-8-sig")
    raw = imported_skill_findings(source)
    tree = ast.parse(source)
    assignments = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    try:
                        assignments[target.id] = {"line": node.lineno, "value": ast.literal_eval(node.value)}
                    except (ValueError, TypeError):
                        pass
    widths = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            for keyword in node.keywords:
                if keyword.arg == "figsize" and isinstance(keyword.value, (ast.Tuple, ast.List)):
                    widths.append({"source_line": node.lineno, "width_mm": safe_numeric(keyword.value.elts[0]) * 25.4,
                                   "expression": ast.unparse(keyword.value.elts[0])})
    custom_checks = []

    def check(check_id, passed, detail):
        custom_checks.append({"check_id": check_id, "pass": bool(passed), "detail": detail})

    expected_diverging = ["#2166AC", "#F7F7F7", "#B2182B"]
    check("AP0_AST_PALETTE", assignments.get("DIVERGING", {}).get("value") == expected_diverging and "CATEGORICAL" in assignments,
          {"DIVERGING": assignments.get("DIVERGING"), "CATEGORICAL": assignments.get("CATEGORICAL")})
    check("CL2_AST_NOMINAL_WIDTH", len(widths) == 2 and all(abs(row["width_mm"] - 183) < 1e-9 for row in widths), widths)
    explanations = [
        {"check_id": "AP-0", "raw_finding_preserved": True, "disposition": "resolved_literal_whitespace_false_positive",
         "explanation": "The shipped checker requires the exact substring 'DIVERGING = ['. Actual source uses 'DIVERGING   = ['. AST parsing verifies the complete named palette with the exact expected values; no source or raw finding was changed."},
        {"check_id": "CL-2", "raw_finding_preserved": True, "disposition": "resolved_insufficient_arithmetic_parsing",
         "explanation": "The shipped regular expressions do not evaluate figsize=(183/25.4,...). Safe AST arithmetic gives 183 mm for both source declarations. Exported PDF/PNG widths are separately measured because bbox_inches='tight' changes the physical crop."},
        {"check_id": "PDF_TEXT_SIZE", "disposition": "corrected_rotation_sensitive_bbox_size",
         "explanation": "pdfplumber char['size'] is a vertical bounding-box measure: for a 90-degree label it can equal the glyph advance width (e.g. the 7pt letter i is reported as 1.554). Font-size compliance uses actual PDF Tf values from pypdf visitor_text, multiplied by the text and graphics matrix scale. Raw pdfplumber minima remain recorded for transparency."},
    ]
    actual_bar_calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr in {"bar", "barh", "boxplot"}]
    if not actual_bar_calls:
        explanations.append({"check_id": "AP-6", "raw_finding_preserved": True, "disposition": "not_applicable_no_bar_or_box_calls",
                             "explanation": "Any lexical warning is caused by words such as 'bars' in prose. AST contains no bar/barh/boxplot call; the figures use points/intervals and empirical CDFs."})
    explanations.append({"check_id": "AP-4", "raw_finding_preserved": True, "disposition": "supplemented_by_actual_current_image_review",
                         "explanation": "The shipped legend detector only searches 30 characters after .legend(. Both calls use bbox_to_anchor later in their arguments. The hash-bound current PNG review confirms external bottom legends; future renders need a new visual review."})

    plot_audit = json.loads(inputs[str(FIGURES / "plot_audit.json")].decode("utf-8-sig"))
    score_audit = json.loads(inputs[str(ROOT / "evaluation/scoring_audit.json")].decode("utf-8-sig"))
    points = read_csv_bytes(inputs[str(FIGURES / "benchmark_current_source.csv")])
    summary = read_csv_bytes(inputs[str(ROOT / "evaluation/group_summary.csv")])
    coverage = read_csv_bytes(inputs[str(ROOT / "evaluation/coverage_summary.csv")])
    plotted = set(plot_audit["methods"])
    complete = set(score_audit["complete_methods"])
    point_methods = {row["method"] for row in points}
    check("PRIMARY_GROUP_COUNTS", plot_audit["records"] == score_audit["records"] == sum(GROUPS.values()) == 3063 and plot_audit["primary_groups"] == GROUPS,
          {"plot_groups": plot_audit["primary_groups"], "records": score_audit["records"]})
    check("GROUPING_CONFIRMED", score_audit.get("group_confirmed") is True, "User-confirmed mutually exclusive grouping is required.")
    check("ALL_COMPLETE_METHODS_PLOTTED", plotted == complete == point_methods,
          {"plotted": sorted(plotted), "complete_in_scoring_audit": sorted(complete), "source_csv": sorted(point_methods)})
    check("NO_INCOMPLETE_METHOD_IN_MAIN", plotted <= complete and plot_audit.get("unfinished_methods_in_primary_figure") is False and score_audit.get("no_primary_incomplete_methods") is True,
          {"incomplete_in_main": sorted(plotted - complete), "not_yet_complete_metadata_entries": sorted(set(score_audit["methods"]) - complete)})
    check("PROGRESS_SNAPSHOT_EXPLICIT", plot_audit.get("progress_snapshot") is True,
          "This report does not claim all planned models have finished.")
    expected_source_bindings = {
        "plot_script_sha256": hashes[str(ROOT / "plot_benchmark.py")],
        "summary_sha256": hashes[str(ROOT / "evaluation/group_summary.csv")],
        "scoring_audit_sha256": hashes[str(ROOT / "evaluation/scoring_audit.json")],
    }
    source_bindings = {}
    for field, expected_digest in expected_source_bindings.items():
        recorded_digest = plot_audit.get(field)
        matched = isinstance(recorded_digest, str) and recorded_digest == expected_digest
        detail = {"field": field, "recorded_sha256": recorded_digest,
                  "current_file_sha256": expected_digest, "matches": matched,
                  "disposition_if_missing_or_different": "regeneration_required"}
        source_bindings[field] = detail
        check("PLOT_AUDIT_SOURCE_BINDING:" + field, matched, detail)
    generated_value = plot_audit.get("generated_utc")
    generated_valid = False
    generated_error = None
    try:
        assert isinstance(generated_value, str) and generated_value
        generated_time = datetime.fromisoformat(generated_value.replace("Z", "+00:00"))
        assert generated_time.tzinfo is not None and generated_time.utcoffset().total_seconds() == 0
        assert (generated_time - datetime.now(timezone.utc)).total_seconds() <= 300
        generated_valid = True
    except (AssertionError, ValueError, TypeError) as error:
        generated_error = "Missing or invalid UTC generation timestamp: " + type(error).__name__
    generated_check = {"recorded_generated_utc": generated_value, "valid_UTC_timestamp": generated_valid,
                       "error": generated_error, "disposition_if_missing_or_invalid": "regeneration_required"}
    check("PLOT_AUDIT_GENERATED_UTC", generated_valid, generated_check)
    source_binding_valid = generated_valid and all(row["matches"] for row in source_bindings.values())
    ecdf_guard_tokens = ["selected_values.min()", "np.floor(selected_values.min()*10)/10", "ecdf_left=max(-1.,min(-.2", "xlim=(ecdf_left,1)"]
    ecdf_guard_in_code = all(token in source for token in ecdf_guard_tokens)
    check("ECDF_FULL_NEGATIVE_RANGE_GUARD_IN_CURRENT_CODE", ecdf_guard_in_code,
          "Read the updated code's data-driven left limit and its actual xlim use. This does not prove an older saved export used the updated code.")
    artifact_refresh = {
        "plot_script_mtime_ns": (ROOT / "plot_benchmark.py").stat().st_mtime_ns,
        "plot_script_newer_than_export": any((ROOT / "plot_benchmark.py").stat().st_mtime_ns > path.stat().st_mtime_ns for path in input_paths if path.suffix in {".png", ".pdf", ".svg"}),
        "plot_audit_records_new_ECDF_bounds": "ecdf_x_limits" in plot_audit,
        "saved_ECDF_bounds": plot_audit.get("ecdf_x_limits"),
        "source_hash_bindings": source_bindings,
        "generation_timestamp": generated_check,
        "source_hash_and_generation_bindings_valid": source_binding_valid,
        "interpretation": "Missing/mismatched source hashes, invalid generation timestamp, a newer plot script, or missing range audit require regeneration. Matching hashes bind the actual current script, group summary, and scoring audit to the render; old absent fields are not accepted as PASS. Current saved export hashes are separately inspected below.",
    }
    key = lambda row: (row["method"], row["group"], row["metric"])
    expected_keys = {(method, group, metric) for method in plotted for group in GROUPS for metric in ("spearman", "top10_ap", "mae")
                     if metric != "mae" or score_audit["methods"][method]["absolute_angstrom"]}
    actual_keys = [key(row) for row in points]
    check("SOURCE_METHOD_GROUP_METRIC_GRID", set(actual_keys) == expected_keys and len(actual_keys) == len(expected_keys),
          {"rows": len(points), "expected_rows": len(expected_keys), "missing": sorted(expected_keys - set(actual_keys)), "extra": sorted(set(actual_keys) - expected_keys)})
    check("SOURCE_FULL_RECORD_COUNTS", all(int(row["n_records"]) == GROUPS.get(row["group"], -1) and true(row["primary_complete_cohort"]) for row in points),
          "Every displayed mean uses all records of its stratum, rather than each method's advantageous surviving subset.")
    check("SOURCE_VALUES_FINITE", all(math.isfinite(number(row["mean"])) for row in points), "Every displayed point estimate is finite.")
    summary_by_key = {key(row): row for row in summary}
    source_matches = True
    for row in points:
        same = summary_by_key.get(key(row))
        source_matches &= same is not None
        if same:
            for column in ("mean", "ci_low", "ci_high"):
                a, b = number(row[column]), number(same[column])
                source_matches &= (math.isnan(a) and math.isnan(b)) or math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-12)
            source_matches &= row["n_records"] == same["n_records"] and row["n_families"] == same["n_families"]
    check("SOURCE_MATCHES_SMALL_SCORING_SUMMARY", source_matches, "Compared every displayed key/value/count to group_summary.csv; no residue-level or per-record score file was loaded.")
    coverage_lookup = {(row["method"], row["group"]): row for row in coverage}
    coverage_ok = True
    coverage_details = []
    for method in sorted(plotted):
        for group, count in GROUPS.items():
            row = coverage_lookup.get((method, group))
            passed = row is not None and int(row["records"]) == int(row["complete_records"]) == count and int(row["expected_residues"]) == int(row["covered_residues"]) == RESIDUES[group]
            coverage_ok &= passed
            coverage_details.append({"method": method, "group": group, "pass": passed})
    check("FULL_RESIDUE_COVERAGE_FOR_PLOTTED_METHODS", coverage_ok, coverage_details)
    tcr_rows = [row for row in points if row["group"] == "tcr_pmhc_itrace"]
    check("TCR_SINGLE_FAMILY_NO_INDEPENDENT_CI", all(int(row["n_families"]) == 1 and not math.isfinite(number(row["ci_low"])) and not math.isfinite(number(row["ci_high"])) for row in tcr_rows),
          "TCR-pMHC has one family component; current figure correctly uses descriptive points without independent-family intervals.")

    exports = []
    visual = []
    for stem in ("benchmark_current", "benchmark_distributions_current"):
        png = FIGURES / (stem + ".png")
        pdf = FIGURES / (stem + ".pdf")
        svg = FIGURES / (stem + ".svg")
        with Image.open(io.BytesIO(inputs[str(png)])) as image:
            dpi = image.info.get("dpi")
            dimensions = list(image.size)
            physical_mm = [dimensions[i] / dpi[i] * 25.4 for i in range(2)] if dpi else None
            nonblank = any(lo != hi for lo, hi in image.convert("RGB").getextrema())
        reader = PdfReader(io.BytesIO(inputs[str(pdf)]))
        fonts = pdf_font_records(reader)
        pages = []
        with pdfplumber.open(io.BytesIO(inputs[str(pdf)])) as document:
            for page_index, page in enumerate(document.pages, 1):
                raw_sizes = [float(char["size"]) for char in page.chars if str(char.get("text", "")).strip()]
                font_runs = []
                def record_font_run(text, cm, tm, font_dictionary, font_size):
                    if not text.strip():
                        return
                    # A PDF text glyph's vertical direction (0,1) is first
                    # transformed by tm, then cm. Pure 90-degree rotation has
                    # scale=1 and must not shrink the reported font size.
                    vx, vy = float(tm[2]), float(tm[3])
                    transformed_x = float(cm[0]) * vx + float(cm[2]) * vy
                    transformed_y = float(cm[1]) * vx + float(cm[3]) * vy
                    scale = math.hypot(transformed_x, transformed_y)
                    font_runs.append({"text": text.strip(), "declared_Tf_pt": float(font_size),
                                      "effective_size_pt": float(font_size) * scale,
                                      "graphics_text_scale": scale,
                                      "font": str(font_dictionary.get("/BaseFont", "")) if font_dictionary else ""})
                reader.pages[page_index - 1].extract_text(visitor_text=record_font_run)
                sizes = [run["effective_size_pt"] for run in font_runs]
                assert sizes, "No actual PDF font runs found"
                outside = [char for char in page.chars if char["x0"] < -.25 or char["x1"] > page.width + .25 or char["top"] < -.25 or char["bottom"] > page.height + .25]
                pages.append({"page": page_index, "width_pt": page.width, "height_pt": page.height,
                              "width_mm": page.width / 72 * 25.4, "height_mm": page.height / 72 * 25.4,
                              "text_characters": len(page.chars), "minimum_text_size_pt": min(sizes), "maximum_text_size_pt": max(sizes),
                              "font_size_basis": "pypdf actual PDF Tf text operator with combined graphics/text matrix scale",
                              "font_size_text_run_counts": dict(sorted(Counter(round(size, 3) for size in sizes).items())),
                              "minimum_pdfplumber_bbox_size_raw": min(raw_sizes),
                              "rotated_text_characters": sum(not char.get("upright", True) for char in page.chars),
                              "font_runs": font_runs,
                              "characters_outside_page": len(outside)})
        svg_text = inputs[str(svg)].decode("utf-8")
        check(stem + ":PNG_DPI_AND_DIMENSIONS", dpi is not None and all(abs(value - 300) < 1 for value in dpi) and min(dimensions) > 500 and nonblank,
              {"pixels": dimensions, "dpi": dpi, "physical_mm": physical_mm, "nonblank": nonblank})
        check(stem + ":PDF_EMBEDDED_FONTS", bool(fonts) and all(font["embedded"] and "Type3" not in font["subtype"] for font in fonts), fonts)
        check(stem + ":PDF_TEXT_SIZE_FLOOR", all(page["minimum_text_size_pt"] >= 5 - 1e-5 for page in pages),
              {"minimum_pt": min(page["minimum_text_size_pt"] for page in pages), "basis": "Actual PDF text extraction, not source regex."})
        check(stem + ":PDF_PAGE_AND_CROP", len(pages) == 1 and all(abs(page["width_mm"] - NOMINAL_WIDTH_MM) <= 3 and page["height_mm"] <= 247 and page["characters_outside_page"] == 0 for page in pages), pages)
        check(stem + ":PNG_PDF_PHYSICAL_AGREEMENT", physical_mm is not None and all(abs(physical_mm[i] - [pages[0]["width_mm"], pages[0]["height_mm"]][i]) < .5 for i in range(2)),
              "Allow 0.5 mm for different PDF/raster text-metric tight crops; record actual dimensions rather than assuming exactly 183 mm.")
        check(stem + ":SVG_EDITABLE_TEXT_PRESENT", "<text" in svg_text, {"text_elements": svg_text.count("<text")})
        exports.append({"stem": stem, "png": {"sha256": hashes[str(png)], "pixels": dimensions, "dpi": dpi, "physical_mm": physical_mm},
                        "pdf": {"sha256": hashes[str(pdf)], "pages": pages, "fonts": fonts}, "svg": {"sha256": hashes[str(svg)]}})
        review = VISUAL_REVIEWS.get(png.name)
        if review and review["sha256"] == hashes[str(png)]:
            visual.append({"file": png.name, **review})
        else:
            visual.append({"file": png.name, "sha256": hashes[str(png)], "status": "pending_actual_visual_review",
                           "reason": "Image hash differs from the PNG that was actually inspected. Automated export/data checks do not substitute for visual inspection."})
    check("INPUT_FILES_STABLE_DURING_QA", all(sha256(Path(path)) == digest for path, digest in hashes.items()),
          "Rerun after controller finishes rendering if any scoring/plot file changed during this snapshot audit.")
    resolved_ids = {"AP-0", "CL-2"}
    unresolved_failures = [row for row in raw if row["category"] == "FAIL" and row["check_id"] not in resolved_ids]
    automated_ok = all(row["pass"] for row in custom_checks) and not unresolved_failures
    visually_current = all(row["status"] == "inspected_current_hash_with_notes" for row in visual)
    needs_refresh = not source_binding_valid or artifact_refresh["plot_script_newer_than_export"] or not artifact_refresh["plot_audit_records_new_ECDF_bounds"]
    if not source_binding_valid:
        status = "REGENERATION_REQUIRED_PLOT_SOURCE_HASH_BINDING_INVALID"
    elif automated_ok and visually_current and needs_refresh:
        status = "SAVED_SNAPSHOT_QA_PASS_UPDATED_PLOT_REGENERATION_REQUIRED"
    elif automated_ok and visually_current:
        status = "PASS_CURRENT_SNAPSHOT_WITH_DOCUMENTED_VALIDATOR_LIMITATIONS"
    elif automated_ok:
        status = "AUTOMATED_CHECKS_PASS_VISUAL_REVIEW_PENDING"
    else:
        status = "FAIL_CHECKS_REQUIRE_REVIEW"
    report = {
        "status": status, "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Only the saved current snapshot and exact hashes below; not future model results or future renders.",
        "qa_script_sha256": sha256(__file__), "validator_path": str(VALIDATOR), "validator_sha256": sha256(VALIDATOR),
        "validator_invocation": "importlib import plus explicit check_ap0..check_cl7 function calls; no nonexistent CLI relied upon",
        "original_validator_findings": raw, "original_validator_counts": dict(Counter(row["category"] for row in raw)),
        "supplementary_explanations": explanations, "checks": custom_checks,
        "passed_custom_checks": sum(row["pass"] for row in custom_checks), "total_custom_checks": len(custom_checks),
        "plot_methods": plot_audit["methods"], "plot_groups": GROUPS, "source_csv_rows": len(points),
        "artifact_refresh": artifact_refresh,
        "exports": exports, "visual_reviews": visual, "input_sha256": hashes,
        "limitations": ["No future image is considered visually reviewed merely because this script runs successfully.",
                        "No color-vision simulation or printed-page inspection was performed.",
                        "ECDF all-observation inclusion was inspected in source code/audit assertions; this memory-light QA does not reload the large per-record or residue-level scoring files.",
                        "PNG is a 300 dpi preview; PDF/SVG are vector masters for line art. Tight bounding-box export slightly changes page width from nominal 183 mm."],
    }
    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "qa_raw_skill_findings.json").write_text(json.dumps(raw, indent=2, ensure_ascii=False), encoding="utf-8")
    (FIGURES / "qa_export_metrics.json").write_text(json.dumps(exports, indent=2, ensure_ascii=False), encoding="utf-8")
    (ROOT / "QA_REPORT.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    lines = ["# 当前横评图 QA", "", f"状态：`{status}`。本报告只对应下列已保存文件的哈希，不代表尚未运行完的模型或未来重绘图已经通过检查。", "",
             f"原 skill validator 已通过 import 显式执行，共 {len(raw)} 条原始 findings：{dict(Counter(row['category'] for row in raw))}。原文保留在 JSON 和 `figures/qa_raw_skill_findings.json`。补充文件/数据检查 {sum(row['pass'] for row in custom_checks)}/{len(custom_checks)} 通过。", "",
             "AP-0 的配色 FAIL 来自精确空白匹配：代码写成 `DIVERGING   = [`，检查器只接受 `DIVERGING = [`；AST 已验证变量及颜色值完整存在。CL-2 无法解析 `183/25.4`，AST 给出两个画布的名义宽度均为 183 mm；同时独立测量实际导出的 PDF/PNG，不用名义值代替实际尺寸。AP-4 的 30 字符 legend 搜索与 AP-6 对注释中 bars 的词法搜索限制也保留并解释。", "",
             "PDF 字号检查已修正旋转文字的解析误判：pdfplumber 的 size 对竖排坐标标签可能是字形宽度，7 pt 的 i 会被读成 1.554；本次直接读取 PDF 的 Tf 字号，并考虑文字与图形矩阵缩放，原始 bbox size 也留档。", "",
             f"绘图来源绑定：plot_audit 的 generated_utc 以及 plot_script_sha256 / summary_sha256 / scoring_audit_sha256 全部逐项校验；当前绑定状态为 {source_binding_valid}。任何旧字段缺失、哈希不符或UTC时间无效均要求重绘，不作为通过。", "",
             f"数据检查：主图 {len(plotted)} 个完整方法/配置，三个互斥组为 2,193 / 626 / 244 条，总计 3,063 条；source.csv 共 {len(points)} 行。核对了全部方法×组×适用指标、有效残基覆盖、summary 数值和完成状态；相对分数方法没有被写成绝对 Å MAE。未将 incomplete 方法放入主图。", "", "| 文件 | PNG像素 / DPI | PDF实际宽×高 mm | PDF最小字号 | 字体嵌入 |", "|---|---|---|---|---|"]
    for entry in exports:
        page = entry["pdf"]["pages"][0]
        lines.append(f"| {entry['stem']} | {entry['png']['pixels']} / {entry['png']['dpi']} | {page['width_mm']:.3f} × {page['height_mm']:.3f} | {page['minimum_text_size_pt']:.2f} pt | {all(font['embedded'] for font in entry['pdf']['fonts'])} |")
    lines += ["", "实际 PDF 均检查了字体文件流、文本字号和字符边界；名义 183 mm 经 tight crop 后的实际差异已记录，未隐藏。", "", "视觉检查与哈希："]
    if needs_refresh:
        reasons = []
        if not source_binding_valid:
            reasons.append("generated_utc或绘图脚本/summary/scoring audit的来源哈希绑定缺失、无效或不符")
        if artifact_refresh["plot_script_newer_than_export"]:
            reasons.append("当前绘图脚本比已导出图新")
        if not artifact_refresh["plot_audit_records_new_ECDF_bounds"]:
            reasons.append("plot_audit尚未记录ecdf_x_limits")
        lines += ["", "需要重新导出的原因：" + "；".join(reasons) + "。不能将这些已保存图视为当前计算与新ECDF保护已验证，需要重新绘制并按新哈希视觉复核。QA本身不修改或重新执行绘图脚本。"]
    for review in visual:
        lines += ["", f"- `{review['file']}`：`{review['status']}`，SHA-256 `{review['sha256']}`。"]
        lines.extend("  - " + observation for observation in review.get("observations", []))
    failed = [row for row in custom_checks if not row["pass"]]
    if failed:
        lines += ["", "未通过的补充检查："] + [f"- {row['check_id']}：{row['detail']}" for row in failed]
    lines += ["", "限制：未做色觉模拟或纸张打印检查；为控制内存没有读取大体积逐记录/逐残基评分文件。未来图像哈希变化时，本脚本自动把视觉复核改为 pending，不会沿用旧图的通过记录。"]
    (ROOT / "QA_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "custom_checks": f"{sum(row['pass'] for row in custom_checks)}/{len(custom_checks)}", "original_validator_counts": report["original_validator_counts"],
                      "plot_methods": len(plotted), "source_rows": len(points), "failed_checks": [row["check_id"] for row in failed]}, indent=2), flush=True)
    return 0 if automated_ok else 2


if __name__ == "__main__":
    sys.exit(main())
