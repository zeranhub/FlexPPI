# Academic Figure Skill QA Report

Figure: ITrace frozen-checkpoint full-complex context ablation  
Target: Nature-family standard, double-column (183 mm)  
Backend: Python/matplotlib  
Exports: vector PDF, editable SVG, 300 dpi PNG

## Figure contract

- Core conclusion: the same frozen FlexPPI checkpoint benefits from full-complex
  context overall, but the gain is chain-role dependent and isolated short
  peptides exhibit severe length extrapolation.
- Panel a: paired complex-level Spearman values establish that the aggregate
  gain is observed across individual complexes.
- Panel b: seven dimensionless paired metrics establish cross-metric
  consistency.
- Panel c: chain-role effects establish biological and architectural
  heterogeneity.
- Panel d: absolute RMSF scale establishes the isolated-peptide extrapolation
  failure that qualifies the aggregate absolute-error result.
- Independent statistical unit: PDB complex, `n = 244`.

## Pass 0 - Anti-pattern scan

- [PASS] AP-0: verbatim typography, color and export baselines are present.
- [PASS] AP-1/AP-2: no default, jet, rainbow or HSV palette is used.
- [PASS] AP-3: top and right spines are removed.
- [PASS] AP-4: panel d legend does not occlude the isolated-peptide point;
  panel a uses the caption rather than an internal legend.
- [PASS] AP-5: vector PDF and SVG are exported in addition to PNG.
- [N/A] AP-6: no bar/box summary with fewer than ten observations is used.
- [PASS] AP-7: Arial/Helvetica/Liberation Sans is explicitly configured.

## Pass 1 - Code compliance

- [PASS] CL-1: minimum declared font size is 6.5 pt.
- [PASS] CL-2: canvas width is 183 mm and height is 132 mm.
- [PASS] CL-3: PNG export is 300 dpi (`2137 x 1554` pixels).
- [PASS] CL-4: PDF TrueType fonts are embedded; SVG text remains editable.
- [PASS] CL-5: spine width is 0.6 pt.
- [PASS] CL-6: ticks point outward.
- [PASS] CL-7: PDF, PNG and SVG exports all exist and are non-empty.

## Pass 2 - Visual logic and data integrity

- [PASS] VI-1: panel a is the hero evidence; the diagonal and winner-specific
  shape/color encoding make the paired conclusion visible immediately.
- [PASS] VI-2: full-complex and isolated-chain conditions use red circles and
  blue squares, so color is not the sole cue.
- [PASS] VI-3: no decorative borders, gradients, 3-D effects or redundant
  grids are present.
- [PASS] VI-4: signed effect panels include a zero reference; panel d labels
  its logarithmic axis explicitly.
- [PASS] VI-5: intervals are defined as 2,000 complex-level bootstrap
  percentile 95% CIs; paired tests are two-sided Wilcoxon signed-rank tests;
  Benjamini-Hochberg q-values are stored for the seven metric and five role
  families.
- [PASS] VI-6: panel labels a-d use identical position, size and weight.
- [PASS] VI-7: multi-panel styling is consistent and no known revision-case
  anti-pattern remains.

## Pass 3 - Rendered output verification

- [PASS] VV-1: no data, confidence interval, tick label or legend is occluded.
- [PASS] VV-2: the 2 x 2 grid is aligned and all labels remain within the
  exported canvas.
- [PASS] VV-3: text is legible in the 183 mm manuscript rendering.
- [PASS] VV-4: red, blue, black and grey elements remain visually distinct;
  marker shape provides redundant encoding.
- [PASS] VV-5: all source arrays are finite. Approximate non-white pixel
  densities are 9.44% (a), 4.27% (b), 3.36% (c) and 4.31% (d), within the
  1.5%-50% target. The generic scatter-correlation threshold is not applicable
  to panel a because no x-y regression claim is made; its signal is the paired
  displacement from the equality line. Panel d's between-condition range ratio
  is 0.945.

## Data traceability

- Panel a:
  `paper/data/itrace_context/panel_a_per_complex_spearman.csv`.
- Panel b:
  `paper/data/itrace_context/panel_b_metric_advantages.csv`.
- Panel c:
  `paper/data/itrace_context/panel_c_role_spearman.csv`.
- Panel d:
  `paper/data/itrace_context/panel_d_role_rmsf_scale.csv` and its summary file.
- Complete external-clean metric tables:
  `paper/data/itrace_context/per_complex_metrics_external_clean.csv` and
  `paper/data/itrace_context/per_chain_metrics_external_clean.csv`.
- Statistics and checkpoint provenance:
  `paper/figures/fig_itrace_context_ablation_statistics.json`.
- Input rows retained: 244 complexes, 1,220 chains, 199,668 residues; no
  visualization downsampling was performed.

## Verdict

READY. All required checks passed after two render-review cycles. The final
manuscript also compiled twice with XeLaTeX without warnings, overfull boxes,
underfull boxes or unresolved references, and pages 1, 14-17 and 20-23 were
rendered and visually inspected.
