# RMSF baseline comparison figures — QA record

## Figure contract

The core conclusion is that FlexPPI has a broad per-complex ranking advantage on the exact-PDB-clean ITrace external cohort, while the size and even the direction of that advantage depend on the functional region. AlphaFlow is shown separately because only six official ATLAS-overlap systems are available; it is therefore a small descriptive case study rather than another row in the 244-complex ranking. Both figures target a Nature-family double-column layout and are delivered as vector PDF, editable SVG, and 300 dpi PNG.

## Pass 0 — scientific integrity

- **AP0, provenance:** PASS. Every plotted value is loaded from frozen evaluation CSV files; no synthetic or illustrative data are used.
- **AP1, sample definition:** PASS. The main figure uses one equal-weight observation per ITrace complex after exact-PDB exclusion (`n = 244`). The AlphaFlow figure uses one equal-weight protein from the official overlap (`n = 6`, 1,486 mapped residues).
- **AP2, axes and units:** PASS. Correlation panels use Spearman rho; absolute RMSF and MAE use angstroms; standardized residue profiles are explicitly described as independent within-chain, within-method population z-scores.
- **AP3, uncertainty:** PASS. The primary paired mean difference and AlphaFlow convergence curves use 2,000 bootstrap resamples with seed 2026 and show 95% confidence intervals.
- **AP4, statistical test:** PASS. The primary ITrace comparison is the prespecified paired, two-sided Wilcoxon test. Its recomputed value exactly matches the frozen result (`p = 4.16099799218336e-7`). No multiple-comparison correction is applied to this single prespecified contrast. The six-system AlphaFlow Wilcoxon result is presented only as descriptive evidence.
- **AP5, representative examples:** PASS. 4JRY and 8RLV are not extrema: each is the complex closest to the median paired difference within the positive or negative half of the distribution.
- **AP6, train/test scope:** PASS. The contract states a frozen FlexPPI checkpoint, the ITrace exact-PDB-clean external cohort for the main comparison, and the ATLAS development-domain overlap for AlphaFlow.
- **AP7, no pooling across incompatible scopes:** PASS. AlphaFlow is deliberately separated from the full-cohort heatmap and violin plot.

## Pass 1 — technical output

- **CL1, dimensions:** PASS. The main figure is exactly 183 x 151 mm; the AlphaFlow figure is exactly 183 x 126 mm.
- **CL2, raster resolution:** PASS. PNG exports are 2,161 x 1,783 and 2,161 x 1,488 pixels respectively, corresponding to 300 dpi at the declared physical sizes.
- **CL3, vector delivery:** PASS. PDF and SVG versions are present for both figures.
- **CL4, typography:** PASS. Arial/Arial Bold are embedded in the PDFs; mathematical glyphs use the matching DejaVu Sans math face. The minimum explicit font size is at least 5 pt and panel labels are consistent.
- **CL5, line and marker rendering:** PASS. Lines remain distinguishable at print size, ticks point outward, top/right spines are removed where appropriate, and the identity line is visually subordinate to the data.
- **CL6, source data:** PASS. Long-form source-data exports contain 8,679 rows for the main figure and 492 rows for AlphaFlow, with panel identifiers retained.
- **CL7, reproducibility:** PASS. `make_rmsf_baseline_comparison_figures.py` regenerates all graphics, source-data files, and the statistics JSON from the frozen evaluation artifacts.

## Pass 2 — visual and statistical communication

The main conclusion is legible without relying on the caption: panel a shows the full per-complex distribution for seven methods; panel b exposes both wins and losses relative to BackFlip-seq instead of hiding heterogeneity in a mean; panel c shows the region-dependent reversal of the winning method; and panel d anchors the aggregate results in two typical residue-level profiles. Color is redundant with point position, line style, direct labels, and black best-cell boxes, so the interpretation does not depend on red/blue perception alone. The heatmap is centered at zero, uses all 13 by 7 values, and has sufficient signal (`mean-row range = 0.492`, 91 unique values). The paired panel intentionally has low cross-method Pearson correlation (`r = 0.067`); its purpose is not linear calibration but paired distance from `y = x`, supported by a standardized mean difference of 0.365 and a 61.9% win rate.

For AlphaFlow, all six systems remain visible behind the mean convergence curve, direct labels replace a crowded legend, and 6TGK-C is explicitly identified as an absolute-scale outlier (`MAE = 2.471 Å`). This prevents the high rank correlations from being misread as uniformly accurate absolute RMSF calibration. Because `n = 6`, the panel is interpreted descriptively and is not pooled with the ITrace rankings.

## Pass 3 — rendered-paper inspection

PASS. Both source PNGs and the compiled PDF pages were inspected at full-page and figure scale. Axis labels, heatmap cells, method names, chain boundaries, leader lines, panel letters, and captions are readable; no annotation covers a data point needed for interpretation; there is no clipping at the left or bottom edges; and the figures remain internally aligned after LaTeX placement. The final XeLaTeX/BibTeX/XeLaTeX/XeLaTeX build completed without undefined references. The baseline comparison occupies a dedicated float page, which preserves legibility at the journal-width target.

## Key numerical checks

The ITrace mean Spearman values are FlexPPI 0.740, BackFlip-seq 0.704, FastProtFlex 0.685, BackFlip-structure 0.637, RMSF-net-PDB 0.526, ProDy-ANM 0.403, and ProDy-GNM 0.371. FlexPPI minus BackFlip-seq has mean difference +0.0355, median difference +0.0284, bootstrap 95% CI [+0.0238, +0.0475], and 151/244 wins. AlphaFlow minus FlexPPI has mean paired Spearman difference +0.158 across the six overlap systems and wins 6/6; AlphaFlow mean Spearman rises from 0.816 at 10 frames to 0.909 at 250 frames.

## Verdict

**READY.** The two figures pass data-integrity, statistical-annotation, export, source-data, typography, and rendered-paper checks. The only material interpretive caveat is the deliberately separated six-system AlphaFlow scope.
