# Figure 2 relayout QA

READY. Main figure: 183 × 239.44 mm; supplementary S3: 183 × 89 mm. Raster previews are 450 dpi. Main PDF and SVG contain zero raster image objects; original native source paths, fonts and plot geometry are retained.

Pass 0, anti-patterns: inherited manuscript palette, Arial fonts, exterior legends, vector output and retained styles/markers pass. Original density panel has independent upper/right Hchi axes by design. Original grouped bars retain their light major grid and family bootstrap error bars.

Pass 1, code compliance: mandatory typography/color/export blocks present; physical width and height within requested bounds; no text extends outside either newly rendered canvas; smallest newly added type is the 6.8 pt cohort line, above the 5 pt floor. Font embedding is verified in the PDF. Editable text is retained in SVG.

Pass 2, data and logic: original a/b reused at unchanged physical scale. Native bottom group retains all 63 bar means and intervals; only the panel letter changes from d to f. C contains all four original within-type Hchi PR series and their frozen AP values. D/E contain all seven methods and all 14,014 points from the frozen external mean curves. Shared aligned metric table contains independently calculated AP/AUROC for 615 records. S3 retains all 816 vertices for the original RMSF and overall Hchi curves. Main c and S3 retain the original continuous piecewise-linear precision-axis transformations, explicitly documented for captions. No data reduction, training or statistical recalculation occurs.

Pass 3, rendered output: complete main figure and S3 inspected with image viewer. Main c/d/e are aligned in the second row, table is below the axis labels, all seven method line styles and marker keys are visible, and no text covers the curves. Bottom f retains all three metrics and the first axis's 1.0 tick. Original a/b at 160 dpi are pixel-identical to the native source crop; bottom is checked against a reference native PDF placed at the same global coordinates, with only the isolated d/f letter area excluded. Exact comparison values are in figure_QA.json.

Source paths, checksums and statistical definitions are in source_manifest.json and statistics_and_reproducibility.md. To reproduce, run python relayout_fig2.py. The script leaves the regenerated report at visual-review-pending until a human or agent performs a fresh full-image visual inspection.
