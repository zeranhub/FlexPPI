# Figure 2 statistics and reproducibility

This evaluates existing five-fold out-of-fold predictions. No model fitting, new inference, MD simulation, prediction recalibration or threshold optimization is performed. Source identity, fold, residue identities, targets and valid masks are verified by ../analysis/input_alignment_audit.json; original per-record MAE and Spearman values are reproduced for every method.

Statistical unit: one input record. All 3,063 records are retained for both targets. RMSF includes 961,321 valid positions; Hchi includes 824,988 positions. Four methods use identical common-valid support. Full is the complete structural model, Intra the independently trained intrachain model, Static the static-feature regression baseline, and AA the training-fold amino-acid-type mean baseline.

a/b: Every Full/reference point is retained in prediction_reference_<target>.npz. Point weight is 1 / (3,063 * valid positions in its record). The 90 * 90 histogram stores cell probability mass summing to 1. Both axes cover the same full range and have equal physical scaling. RMSF uses symlog with threshold 1 Angstrom, base 10, linscale 1; Hchi is linear. Colors use a logarithmic probability-mass key with zero cells white. Original units are retained. The ideal diagonal is not fitted. The black curve is a descriptive record-weighted conditional mean of unchanged predictions in 20 MD-reference bins. Its band is a source-stratified record-bootstrap 95% interval from the first 1,000 of 10,000 resamples, seed 20261018. MAE and Spearman annotations are record means, not pooled statistics.

c/d: For each record/target, MD >= its 90th percentile (NumPy linear quantile) defines positives, retaining threshold ties. Continuous original predictions rank positions; score ties remain groups. No records have a one-class outcome. Curves show mean precision over records at recall 0,.01,...,1, using the first score group reaching recall. The no-selection endpoint is (0,1); the original all-selected endpoint is appended at recall 1 with mean prevalence. This vertical tail completes the original threshold path without changing AP. Grey lines are mean prevalence (0.1025791 RMSF; 0.1029139 Hchi). AP is the noninterpolated sum of recall increases times threshold precision averaged over records, not area under the mean curve. AUROC is tied-score ROC trapezoid area averaged over records. AP/AUROC intervals use 10,000 source-stratified record bootstrap draws, seed 20261018, using the same draws for methods. PR interval data are in CSV; four overlapping bands are not drawn. No hypothesis tests or multiplicity correction are used.

| Target | Method | Mean AP [95% CI] | Mean AUROC [95% CI] |
|---|---|---|---|
| rmsf | Full | 0.730037 [0.724241, 0.736057] | 0.934391 [0.932318, 0.936495] |
| rmsf | Intra | 0.695543 [0.689315, 0.701854] | 0.921456 [0.918907, 0.924062] |
| rmsf | Static | 0.657131 [0.651249, 0.663247] | 0.910346 [0.908176, 0.912608] |
| rmsf | AA | 0.168053 [0.165473, 0.170710] | 0.625688 [0.622892, 0.628461] |
| chi_entropy | Full | 0.926270 [0.924229, 0.928233] | 0.989238 [0.988920, 0.989536] |
| chi_entropy | Intra | 0.920378 [0.918183, 0.922493] | 0.988566 [0.988225, 0.988894] |
| chi_entropy | Static | 0.906589 [0.904256, 0.908805] | 0.986724 [0.986377, 0.987061] |
| chi_entropy | AA | 0.778923 [0.775299, 0.782532] | 0.974894 [0.974377, 0.975414] |

e/f: All 76 type-level observations and their x positions exactly equal old Figure 2c/d. Original 80-row paired-support CSV is copied byte-for-byte; four Hchi A/G rows remain unavailable, not zero gains. There are 40 RMSF points and 36 Hchi points. Two negative Hchi Intra comparisons (P/V) remain. Gain is 100 * (comparator record-mean MAE - Full record-mean MAE) / comparator record-mean MAE, on identical common pair support: ratio of record means, not mean of per-record ratios. Ordering, values, markers, y limits and Hchi symlog threshold (10 percentage points) are inherited. No CIs or tests are drawn for type gains.

Traceability: ../analysis contains all numerical grids, trend files and raw threshold points. source_type_visible_points.csv contains complete e/f coordinates; source_archive/source_type_paired_gains_all80_rows.csv contains paired support. data_audit.json records source and output hashes. Recompute by running ../analysis/compute_metrics.py, then plot_predictive_capability.py, then finalize_figure2_qa.py in the bundled runtime.

Auxiliary within-amino-acid ranking is documented in the analysis protocol and reported in the manuscript; it does not alter c/d scores or e/f gains.
