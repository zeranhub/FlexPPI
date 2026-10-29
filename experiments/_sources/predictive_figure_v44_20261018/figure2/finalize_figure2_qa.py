from pathlib import Path
import sys, json, hashlib
ROOT = Path(__file__).resolve().parent
EXP = ROOT.parent.parent
sys.path.insert(0, str(EXP / 'completed_results_20261012' / 'plot_runtime'))
import numpy as np
import pandas as pd
from PIL import Image
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
AN = ROOT.parent / 'analysis'
audit = json.loads((ROOT / 'data_audit.json').read_text(encoding='utf8'))
summary = json.loads((AN / 'summary.json').read_text(encoding='utf8'))
inputs = json.loads((AN / 'input_alignment_audit.json').read_text(encoding='utf8'))
assert inputs['status'] == 'passed' and summary['status'] == 'complete'
assert audit['svg_raster_image_count'] == 0 and not audit['canvas_text_outside']
assert audit['minimum_font_at_manuscript_width_pt'] >= 8
assert audit['source_summary_sha256'] == sha(AN / 'summary.json')
cap = (ROOT / 'caption_zh.md').read_text(encoding='utf8').strip()
assert len(cap) <= 450
endpoints = pd.read_csv(AN / 'all_selected_endpoint.csv')
assert len(endpoints) == 8
for r in endpoints.itertuples():
    assert r.recall == 1 and abs(r.precision - summary['metrics'][r.target][r.method]['prevalence_mean']) < 1e-12
old = pd.read_csv(EXP / 'figure_merge_v43_20261018/figure2/source_all_visible_points.csv', float_precision='round_trip')
old = old[old.panel.isin(['c', 'd'])].copy()
old['panel'] = old.panel.map({'c': 'e', 'd': 'f'})
new = pd.read_csv(ROOT / 'source_type_visible_points.csv', float_precision='round_trip')
keys = ['panel', 'metric', 'method', 'amino_acid']
old = old.sort_values(keys).reset_index(drop=True)
new = new.sort_values(keys).reset_index(drop=True)
assert old[keys].equals(new[keys])
assert np.array_equal(old[['x', 'y']].to_numpy(), new[['x', 'y']].to_numpy())
im = np.asarray(Image.open(ROOT / 'fig2_predictive_capability.png').convert('RGB'))
h, w = im.shape[:2]
rects = [[x, y, .365, hh] for y, hh in [(.703, .199), (.395, .202), (.075, .164)] for x in [.107, .601]]
densities = {}
for k, (x, y, ww, hh) in zip('abcdef', rects):
    crop = im[round((1-y-hh)*h):round((1-y)*h), round(x*w):round((x+ww)*w)]
    fraction = float(np.any(crop < 240, axis=-1).mean())
    densities[k] = fraction
    assert fraction > .015
outputs = {ext: {'path': str(ROOT / f'fig2_predictive_capability.{ext}'),
                 'sha256': sha(ROOT / f'fig2_predictive_capability.{ext}')}
           for ext in ['png', 'pdf', 'svg']}
audit.update(status='passed', caption_characters=len(cap), exact_type_points_equal_v43=True,
             all_selected_endpoints_verified=8, pixel_content_fractions=densities, outputs=outputs)
(ROOT / 'data_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf8')
visual = {'status':'passed', 'figure_png_sha256':outputs['png']['sha256'],
          'review_method':'Final native PNG opened and inspected; all six panels checked.',
          'passes': {'AP0_7':'Verbatim baseline; explicit semantic colors and line styles; no default font or border; legends outside data; vector export.',
                     'CL1_7':'183 x 190 mm; minimum 9 pt source, 8.2446 pt at 167.64 mm; 450 dpi PNG; embedded TrueType PDF and editable SVG.',
                     'VI1_7':'Distinct numerical, ranking and within-type evidence; no data subset or tail clipping; complete PR endpoints; negative type gains visible.',
                     'VV1_5':'No text/data occlusion; metric subtitles separated; legends and colorbar clear; all panels contain visible data.'},
          'repairs':['Separated numerical subtitles from titles.', 'Shifted type comparison legend from titles.',
                     'Removed redundant footer that crowded x-axis labels.', 'Disabled default colorbar rasterization.',
                     'Added all-selected PR endpoints at recall 1 without changing AP or recall grid.',
                     'Increased top color-key margin.'],
          'pixel_content_fractions':densities, 'complete_type_source_inheritance':True,
          'scientific_limits':'Raw Hchi ranking combines type and within-type information; independent within-type results reported separately. Record bootstrap is not a family independence claim.'}
(ROOT / 'visual_QA.json').write_text(json.dumps(visual, ensure_ascii=False, indent=2), encoding='utf8')
lines = ['# Figure 2 statistics and reproducibility', '',
         'This evaluates existing five-fold out-of-fold predictions. No model fitting, new inference, MD simulation, prediction recalibration or threshold optimization is performed. Source identity, fold, residue identities, targets and valid masks are verified by ../analysis/input_alignment_audit.json; original per-record MAE and Spearman values are reproduced for every method.', '',
         'Statistical unit: one input record. All 3,063 records are retained for both targets. RMSF includes 961,321 valid positions; Hchi includes 824,988 positions. Four methods use identical common-valid support. Full is the complete structural model, Intra the independently trained intrachain model, Static the static-feature regression baseline, and AA the training-fold amino-acid-type mean baseline.', '',
         'a/b: Every Full/reference point is retained in prediction_reference_<target>.npz. Point weight is 1 / (3,063 * valid positions in its record). The 90 * 90 histogram stores cell probability mass summing to 1. Both axes cover the same full range and have equal physical scaling. RMSF uses symlog with threshold 1 Angstrom, base 10, linscale 1; Hchi is linear. Colors use a logarithmic probability-mass key with zero cells white. Original units are retained. The ideal diagonal is not fitted. The black curve is a descriptive record-weighted conditional mean of unchanged predictions in 20 MD-reference bins. Its band is a source-stratified record-bootstrap 95% interval from the first 1,000 of 10,000 resamples, seed 20261018. MAE and Spearman annotations are record means, not pooled statistics.', '',
         'c/d: For each record/target, MD >= its 90th percentile (NumPy linear quantile) defines positives, retaining threshold ties. Continuous original predictions rank positions; score ties remain groups. No records have a one-class outcome. Curves show mean precision over records at recall 0,.01,...,1, using the first score group reaching recall. The no-selection endpoint is (0,1); the original all-selected endpoint is appended at recall 1 with mean prevalence. This vertical tail completes the original threshold path without changing AP. Grey lines are mean prevalence (0.1025791 RMSF; 0.1029139 Hchi). AP is the noninterpolated sum of recall increases times threshold precision averaged over records, not area under the mean curve. AUROC is tied-score ROC trapezoid area averaged over records. AP/AUROC intervals use 10,000 source-stratified record bootstrap draws, seed 20261018, using the same draws for methods. PR interval data are in CSV; four overlapping bands are not drawn. No hypothesis tests or multiplicity correction are used.', '',
         '| Target | Method | Mean AP [95% CI] | Mean AUROC [95% CI] |',
         '|---|---|---|---|']
for t in ['rmsf','chi_entropy']:
    for m in ['Full','Intra','Static','AA']:
        v = summary['metrics'][t][m]
        lines.append(f'| {t} | {m} | {v["mean_ap"]:.6f} [{v["ap_ci"][0]:.6f}, {v["ap_ci"][1]:.6f}] | {v["mean_auroc"]:.6f} [{v["auroc_ci"][0]:.6f}, {v["auroc_ci"][1]:.6f}] |')
lines += ['', 'e/f: All 76 type-level observations and their x positions exactly equal old Figure 2c/d. Original 80-row paired-support CSV is copied byte-for-byte; four Hchi A/G rows remain unavailable, not zero gains. There are 40 RMSF points and 36 Hchi points. Two negative Hchi Intra comparisons (P/V) remain. Gain is 100 * (comparator record-mean MAE - Full record-mean MAE) / comparator record-mean MAE, on identical common pair support: ratio of record means, not mean of per-record ratios. Ordering, values, markers, y limits and Hchi symlog threshold (10 percentage points) are inherited. No CIs or tests are drawn for type gains.', '',
          'Traceability: ../analysis contains all numerical grids, trend files and raw threshold points. source_type_visible_points.csv contains complete e/f coordinates; source_archive/source_type_paired_gains_all80_rows.csv contains paired support. data_audit.json records source and output hashes. Recompute by running ../analysis/compute_metrics.py, then plot_predictive_capability.py, then finalize_figure2_qa.py in the bundled runtime.', '',
          'Auxiliary within-amino-acid ranking is documented in the analysis protocol and reported in the manuscript; it does not alter c/d scores or e/f gains.']
(ROOT / 'statistics_and_reproducibility.md').write_text('\n'.join(lines)+'\n', encoding='utf8')
(ROOT / 'figure2_QA.md').write_text('# Figure 2 QA\n\nREADY. Four skill QA passes are complete. Numeric and visual checks are in data_audit.json and visual_QA.json. Every OOF position, complete PR endpoints, unchanged type points, full RMSF tails, source-to-manuscript font floor above 8 pt, and vector PDF/SVG masters are retained.\n', encoding='utf8')
print(json.dumps({'status':'passed', 'caption_characters':len(cap), 'outputs':outputs}, ensure_ascii=False, indent=2))
