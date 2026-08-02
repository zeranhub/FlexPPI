# ITrace biological-asymmetry figure QA

## Figure contract

- Core conclusion: Bound-state ITrace MD shows a rigid pMHC contact surface but a relatively flexible TCR contact surface, especially on TCRbeta; full-complex FlexPPI collapses this asymmetry into near-universal rigidification.
- Evidence chain: (a) interface-side asymmetry; (b) contact-local versus non-contact context gain; (c) TCR alpha/beta loop asymmetry.
- Target: Nature-family double-column figure, 183 mm, vector PDF plus 300 dpi PNG.
- Review risk: these are bound-state regional patterns, not apo-to-holo changes; peptide residues are excluded from contact-local ranking because isolated peptides are out-of-distribution short inputs.

## Statistics and traceability

- Source: `E:\Project\FlexPPI\output\itrace_zero_shot_context\raw\common_context_predictions.csv`.
- External-clean cohort: 244 complexes and 199,668 residues after removing 3KXF (833 residues), the sole exact PDB training overlap.
- Panel (a): all 7,370 MHC/TCR interface residues; residues are averaged by complex and chain role before complex-level inference.
- Panel (b): all 11,458 CDR1/CDR2/CDR3 residues; residues are averaged by complex and named loop before complex-level inference.
- Panel (c): all 197,390 non-peptide residues; Spearman correlations are computed separately within each complex and interface stratum.
- Center and interval: mean and 5,000-fold complex-bootstrap 95% confidence interval, n=244 complexes.
- Tests: two-sided paired Wilcoxon signed-rank tests; no multiple-testing claim is attached to the figure.

## Four-pass QA

- Pass 0: PASS. Mandatory typography, palette and export baselines are present; no default/rainbow palette, four-sided border, internal legend or raster-only export.
- Pass 1: PASS. Minimum font size is 6.3 pt, width is 183 mm, PDF fonts use TrueType embedding, and both PDF and 300 dpi PNG are exported.
- Pass 2: PASS. The hero panel carries the asymmetric-interface conclusion; method identity uses both color and marker shape; the zero reference is explicit; panel labels and axes are consistent.
- Pass 3: PASS. The 300 dpi render was inspected at original resolution: error bars and points are unobscured, all three panels align, labels remain legible at the 183 mm target width, method colors and marker shapes are distinct, and every panel contains a visible data signal. The official skill validator also exited successfully.
