from pathlib import Path
import sys, json, hashlib
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image, ImageChops

OLD = ROOT.parent / 'figure7_combined_prediction_v33_20261014'
OUT = ROOT / 'figures'
STEM = 'fig7_sequence_order_predictions_and_MD_curves'
OLDSTEM = 'fig7_combined_predictions_and_MD_curves'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj = lambda p: json.loads(Path(p).read_text(encoding='utf8'))
audit = readj(OUT / (STEM + '.audit.json'))
new = pd.read_csv(OUT / (STEM + '.points.csv'), float_precision='round_trip')
old = pd.read_csv(OLD / 'figures' / (OLDSTEM + '.points.csv'), float_precision='round_trip')
checks = []

def check(name, test):
    assert bool(test), name
    checks.append({'check': name, 'pass': True})

check('576_source_points_retained', len(new) == len(old) == 576)
check('all_panels_have_identical_point_counts', new.groupby('panel').size().to_dict() == old.groupby('panel').size().to_dict() == {'a_Gln': 12, 'a_Hchi': 132, 'b': 144, 'c': 144, 'd': 144})
for panel in ['a_Hchi', 'b', 'c', 'd']:
    a = new[new.panel == panel].reset_index(drop=True)
    b = old[old.panel == panel].reset_index(drop=True)
    check(panel + '_identical_source_rows', a.equals(b))
for (site, kind), group in new[new.panel == 'a_Gln'].groupby(['residue_id', 'kind']):
    original = old[(old.panel == 'a_Gln') & (old.residue_id == site) & (old.kind == kind)]
    check('Gln_' + str(site) + '_' + kind + '_x_and_error_unchanged', np.array_equal(group[['x', 'error']].to_numpy(), original[['x', 'error']].to_numpy()))
seq = pd.read_csv(OUT / (STEM + '.all6_sequence_order.csv'), float_precision='round_trip')
check('all6_sequence_order', seq.author_residue_id.tolist() == [2, 31, 40, 41, 49, 62])
check('Q40_retained_with_unchanged_error', abs(float(seq[seq.author_residue_id == 40].signed_error_nats.iloc[0]) - .3676173686981201) < 1e-14)
check('all6_signed_errors_preserved', np.allclose(seq.signed_error_nats, seq.predicted_Hchi_nats - seq.frozen_MD_joint_Hchi_nats, rtol=0, atol=1e-14))
source_hashes = {}
for p in sorted((ROOT / 'source').iterdir()):
    if p.is_file():
        check('unchanged_source_' + p.name, sha(p) == sha(OLD / 'source' / p.name))
        source_hashes[p.name] = sha(p)
check('v33_PNG_not_modified', sha(OLD / 'figures' / (OLDSTEM + '.png')) == '36682433a643991cae532f69ce949c0993f03f8f4fe55354a1885c13a7427fae')
check('xlim_not_expanded_to_reduce_error', [1.65, 3.15] == audit['a_right_xlim_unchanged_from_v33'])
check('original_metrics_unchanged', audit['primary_metrics_unchanged'] == readj(OLD / 'source/fig7_all_residues_Gln_and_MD_free_energy.statistics.json'))
oldim = Image.open(OLD / 'figures' / (OLDSTEM + '.png')).convert('RGB')
newim = Image.open(OUT / (STEM + '.png')).convert('RGB')
check('identical_export_pixel_dimensions', oldim.size == newim.size)
regions = {'a_left_complete_and_all_labels': (0, 0, 1450, 850), 'b_to_d_all_curves_labels_and_axes': (0, 850, newim.width, newim.height)}
for name, region in regions.items():
    check('pixel_identical_' + name, ImageChops.difference(oldim.crop(region), newim.crop(region)).getbbox() is None)
info = Image.open(OUT / (STEM + '.png')).info
physical = [newim.width / info['dpi'][0] * 25.4, newim.height / info['dpi'][1] * 25.4]
text = (OUT / 'plot_sequence_order_predictions_and_MD_curves.py').read_text(encoding='utf8')
check('no_signed_error_label_plotting', 'be.text' not in text and '预测−参考' not in text)
check('no_window_errorbars', 'errorbar(' not in text)
check('all6_display_order_is_sequence_order', audit['all6_order_basis'] == 'author_residue_id_ascending')
check('minimum_font_7point5', audit['minimum_font_pt'] >= 7.5)
check('all_text_within_canvas', audit['text_outside_canvas'] == [])
audit.update(status='complete', visual_review='passed: direct inspection of full PNG; paired points and all six labels clear; no overlap, clipping, or red error text; a-left and b-d pixel-identical to v33', configured_canvas_mm=[183., 135.], actual_export_png_pixels=list(newim.size), actual_export_png_dpi=list(info['dpi']), actual_export_dimensions_mm=physical, export_crop_policy='Inherited v33 rcParams savefig.bbox=tight; retained by request for identical unmodified panels.', independent_checks=checks, independent_check_count=len(checks), pixel_difference_bbox=ImageChops.difference(oldim, newim).getbbox(), unchanged_regions=regions, local_source_sha256=source_hashes, source_row_comparison='All 576 x and measured/predicted y values unchanged. Only a_Gln categorical row positions reorder six complete pairs.', signed_error_column_display='removed uniformly for all six residues', new_hypotheses_or_defensive_prose_added=False)
(OUT / (STEM + '.audit.json')).write_text(json.dumps(audit, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
(ROOT / 'figure_contract_zh.md').write_text('图7 v34小幅修改：全部66个有效残基继续检验总体预测；全部六个谷氨酰胺继续检验同一残基类型内的区别，按残基序号排列；b–d保持模型预先选出的六个位点预测极值Q41与Q2及其MD自由能曲线。删除所有六个位点的带符号误差文字列，不改任何数值、样本或Hχ轴范围。保留Q40和六位点MAE。版式、字体、颜色和下排继承v33。\n', encoding='utf8')
statistics = {
    'new_analysis_performed': False,
    'all_metrics_frozen_from_v33': audit['primary_metrics_unchanged'],
    'panel_a_left': {
        'n': 66,
        'n_definition': 'All Hchi-valid ubiquitin chain-A residues, including three valid zero-entropy prolines. Not 66 independent proteins.',
        'model': 'Frozen full-model test output, fold 0. No refit or new seed.',
        'reference': 'Mean original joint Hchi over five consecutive 1000-frame windows of the same existing MD trajectory. These are not independent MD repeats.',
        'MAE_definition': 'Mean absolute model-minus-reference joint-Hchi difference across all included residues, nats.',
        'rho_definition': 'Spearman rank correlation across all included residues.',
        'corrected_rho_definition': 'Previously frozen Spearman correlation after subtracting each training-set amino-acid-type mean from both model and reference.',
        'spread_or_CI': 'None displayed.',
        'test_or_multiple_correction': 'No new significance test; no multiple-comparison correction.'
    },
    'panel_a_right': {
        'n': 6,
        'n_definition': 'All six ubiquitin glutamines Q2, Q31, Q40, Q41, Q49, Q62; no omission.',
        'order': 'Ascending author residue number.',
        'display': 'MD reference hollow circles and full-model filled squares, joined within each residue.',
        'xlim': [1.65, 3.15],
        'spread_or_CI': 'None displayed.',
        'signed_error_storage': 'Preserved in all6_sequence_order.csv signed_error_nats and points.csv error; removed uniformly from plot text.'
    },
    'panels_b_d': {
        'n': 'Two residues x three angles x 72 five-degree bins = 432 original curve rows.',
        'frame_scope': 'Existing 0–100 ns, 200 ps interval, 501 frames, exactly as v33.',
        'selection': 'Frozen-model minimum Q41 and maximum Q2 among all six GLN, not selected by observed error.',
        'Fstar_definition': '-ln(bin_probability / maximum_bin_probability); zero-count bins remain undefined and unconnected.',
        'low_Fstar_width': 'Sum of 5-degree bins with finite Fstar <= 1; grid origin 0 degrees.',
        'spread_or_CI': 'None. No new MD run or statistical test.'
    },
    'source_traceability': {
        'source/all66_residue_summary.csv': 'All 66 joint-Hchi paired values and residue identifiers.',
        'source/ubiquitin_all66_chi_501frames.npz': 'Independent frozen paired-value equality check.',
        'source/histograms_501frames.csv': 'All 432 original Fstar curve rows.',
        'source/fig7_all_residues_Gln_and_MD_free_energy.statistics.json': 'Previously frozen overall and all-six metrics.'
    },
    'export': {'configured_canvas_mm': [183, 135], 'actual_png_pixels': list(newim.size), 'dpi': list(info['dpi']), 'physical_dimensions_mm': physical, 'crop': 'Original tight crop retained for v33 visual consistency; dimensions reported explicitly.'}
}
(OUT / (STEM + '.statistics.json')).write_text(json.dumps(statistics, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
report = f'''# 图7 v34 核验与统计说明

所有576条图中源数据保留：a左132条（66个位点×预测/参考），a右12条（全部六个GLN×预测/参考），b–d各144条。与v33逐值核验通过，Q40预测、参考及+0.3676173686981201 nats差值均保留。改动仅为a右按Q2、Q31、Q40、Q41、Q49、Q62排序，统一删除六个位点的误差文字列，配对点轴宽度由25 mm增为44 mm；Hχ范围仍为1.65–3.15，没有加时段误差条。

已完成{len(checks)}项独立核验。a左完整区域与b–d全部区域的PNG逐像素等同于v33；四份复制源文件SHA256一致，旧v33 PNG未修改。直接查看PNG，六个位点、配对点、图例和曲线均清晰，无文字碰撞或裁切；字体最小7.5 pt。保留橙色Q41和紫色Q2，并以圆/方和实/虚线提供非颜色区分。未新增机制解释或防御性文字。

声明画布为183×135 mm；沿用v33的tight裁边后，导出PNG为{newim.width}×{newim.height}像素、约300 dpi，实际外沿约{physical[0]:.2f}×{physical[1]:.2f} mm。PDF/SVG同步导出，SVG不含嵌入位图。沿用裁边是为了保持未改区域的位置和像素一致。

统计量、样本含义、参考时间范围、模型折号和逐文件溯源见同目录statistics.json。没有增加检验、显著性标记、误差条、训练或MD实验；五个连续时段没有作为五次独立实验处理。逐点误差仍在CSV中，图中六位点MAE仍为0.105 nats，包含Q40。
'''
(OUT / (STEM + '.QA_statistics_zh.md')).write_text(report, encoding='utf8')
print(json.dumps({'status': 'complete', 'checks': len(checks), 'pixel_difference_bbox': audit['pixel_difference_bbox'], 'PNG_pixels': list(newim.size), 'actual_mm': physical, 'root': str(ROOT)}, ensure_ascii=False))
