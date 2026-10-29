# v45 supplementary computations for the denser Figure 2

`four_method_binned_trend.csv` has 160 rows: two targets x four methods x 20
reference bins. `x_mean` is the mean MD reference and `prediction_mean` the
unaltered mean prediction, each with per-position weight `1 / record_valid_n`.
`record_equal_mass` sums to one per target/method. All four methods share the
same observations and bin edges. Empty reference bins have undefined means;
they are not dropped from the exported 20-row support. Only Full's `ci_low` and
`ci_high` are calculated, using the exact v44 source-stratified record draws.
An interval is shown only where defined; there is no extrapolation across empty
bins. `n_positions` and `n_records` report the observed support per bin.

`within_type_macro_pr.csv` contains 404 rows (101 recall coordinates x four
methods) for task `within_type_hchi`. Precision is first averaged over eligible
amino-acid types within each record, then over records. The 101 recall points
are display samples; do not integrate them to recalculate AP. Use the exact AP
in `v44_summary.json`, which is independently reproduced in `addition_audit.json`.
No confidence band is calculated for this new display.

`within_type_all_selected_endpoint.csv` supplies one genuine full-selection
endpoint for each method. Recall is 1, precision is the mean of each record's
mean eligible-type positive fraction. Append this after the sampled grid's last
point, because the first threshold reaching recall 1 can still leave negative
positions unselected. This preserves the original complete threshold curve.

`raw_within_type_curves_<method>.npz` archives all original threshold points.
Arrays are `type_code`, `threshold`, `precision`, `recall`, `fpr`, `tpr`, in
ascending `type_code` groups and within-group descending threshold order.
`within_type_codes.csv` maps each group to the record, sample, source, fold,
amino-acid identity, valid count, fixed MD q90 threshold, positive count and
prevalence. The leading threshold is positive infinity with precision 1 and
recall 0. Ties in score are selected together. Every group ends at full
selection. `within_type_complete_curve_metrics.csv` stores each original
type AP/AUROC; `within_type_record_macro_metrics.csv` stores the corresponding
within-record type means. `within_type_record_precision_grid_<method>.npz`
stores the 3,063 x 101 matrix before averaging records.

`source_ap_matrix.csv` contains 36 cells: target/task x source x method.
Task names are `rmsf`, `chi_entropy`, `within_type_hchi`; source names are
`atlas`, `dynarepo`, `itrace`. `mean_ap` is the mean original record AP for
that source. `full_minus_method_ap` equals Full's source AP minus that cell's
method AP and is zero in Full rows. `n_records` is the true record count.
Sources are descriptive summaries; no new significance tests are introduced.

Files prefixed `v44_` are byte-for-byte copies of their existing v44 source:
original two-task macro PR curves, full-selection endpoints and summary metrics.
The original v44 density NPZ is consumed directly by the plot without changing
the histogram, edges, original positions or record weights. Every consumed
input is SHA-256 bound in `inputs_sha256.json`. `addition_audit.json` binds code,
protocol, input manifest, generated outputs and all numerical reproduction checks.
