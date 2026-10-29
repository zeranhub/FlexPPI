# Figure 2 design contract

The figure tests numerical agreement, retrieval of high-dynamic positions, and gains beyond amino-acid type. It does not create a new classifier or claim that a ranking curve measures numerical calibration. The user explicitly approved the six-panel revision on 2026-10-18.

Archetype: quantitative grid, three target-matched rows, RMSF left and joint side-chain entropy Hχ right. First row: all common-valid Full OOF predictions versus MD reference. Second row: all four methods' record-equal precision-recall curves for a prespecified reference-defined high-dynamic class. Last row: exact inherited v43 paired within-amino-acid MAE gains, retaining negative gains and unavailable types.

The final source canvas is 183 × 190 mm, with a 9 pt minimum font. At 167.64 mm manuscript width this becomes 8.24 pt. Export is fully vector PDF and editable SVG with a 450 dpi PNG preview. All density cells are vector quadrilaterals; no point subsampling is used. Equal x/y ranges, transforms and physical scales are used in the numerical agreement panels.

## Asset route audit

The skill hub, typography, color, export, journal dimensions, complete checklist, directory map, figure contract, multipanel layout and revision cases were read. Python 3.12.14 and the existing NumPy/Pandas/Matplotlib runtime were verified; no R executable was available. Production companions were visually inspected.

- a/b: `assets/figures/KernelDensity/plot_CorrelationDensity1.R` and PNG were inspected. It produces an unweighted R 2-D KDE, while the approved scientific encoding requires exact record-weighted probability mass and complete tails. Therefore parameter inheritance is used: sequential filled density, clean axes, predicted-versus-observed orientation, external density key. The numerical density is supplied by the separately audited analysis engine; no smoothing or mass truncation is performed.
- c/d: `assets/figures/AUROC/plot_AUROC2.py` and PNG were inspected. ROC's false-positive/true-positive axes are semantically incompatible with precision-recall. This is explicitly a cross-type inheritance, retaining the thin comparison curves, grey reference, clean two-spine axes and method-specific redundant line styles, rather than falsely relabelling a ROC plot. The skill's default rainbow-style colors and large canvas were superseded by its mandatory journal baseline and manuscript font floor.
- e/f: `assets/figures/BarComparison/plot_comparison_GeneRegulatory.py` and PNG were inspected. Its dataset-level bars do not match categorical paired type-level points. The already verified v43 categorical plotting behavior and complete 80-row source CSV are inherited exactly, preserving source ordering, x offsets, y values, point markers, axis limits and Hχ symlog scaling. No new type-level statistics are calculated.

Review risks handled in the encoding: AP is a mean of per-record AP values and never labelled as the area under the macro PR curve; reference-quantile ties remain positive and the plotted no-skill baseline is actual prevalence; the raw Hχ ranking is complemented by the within-type evidence rather than claiming it is type-independent; source-stratified record resampling is described as record uncertainty, not family-independent replication. No hypothesis-test p values or asterisks are plotted.

Before rendering, the anticipated visual risks are density tails hiding the central trend, four-method PR legends crowding curves, and 20 amino-acid labels becoming small. They are addressed with a logarithmic density color key, an equal-axis RMSF symlog transform, external two-column PR legends, and an explicit 9 pt font floor.
