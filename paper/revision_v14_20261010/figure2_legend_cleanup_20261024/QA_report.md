# Figure 2 final legend cleanup QA

READY. Final main figure is 183 × 210.44 mm, with 450 dpi PNG, native vector PDF and editable-text SVG. All three curve axes retain their original 42 × 40 mm physical sizes, colors, line styles, markers, axis scales and data. Original a/b and f bar plots retain their physical scales.

The internal four-method legend is placed in the upper right blank region of c. The external seven-method legend is in e's lower right blank region, below the ROC curves and chance diagonal. The old legend below b is deleted using its precisely identified native PDF drawing block (four keys and four encoded text labels), preserving every preceding byte of the original content stream. No geometric redaction is applied to the dense a/b plots.

Removed from the main image at the user's request: c's AP label and four AP values; d/e cohort subtitles; external AP/AUROC metric table; c's Internal controls and 3,063 records lines; and n/family count suffixes from f's three interface legends. External curve titles are now RMSF PR and RMSF ROC. Counts, cohort definitions and metric details belong in the revised caption/body or supplementary Table S10. Supplement S3 is unchanged.

Pass 0–2: mandatory style/export baselines present, source SHA-256 checks unchanged, no point dropped, no recomputation or smoothing, Arial fonts embedded, no text outside the new canvas, all model names/axes/ticks retained. All 408 within-type vertices and 14,014 external mean-curve points are preserved. Original f bar means and error bars are unchanged.

Pass 3: the complete final image was inspected after the final simplification. No clipped titles, missing ticks or curve occlusion were observed. Internal/external legend boxes are wholly inside their assigned axes. At 160 dpi, a/b are pixel-identical to the native source outside the removed legend; bottom bars/axes/method labels are pixel-identical outside the simplified group legend and isolated d/f letter. Native PDF geometry comparisons against v65 verify all four c and all seven d/e series are unchanged except the intended 4 mm upward layout translation. Main PDF and SVG contain zero raster image objects.

Reproduction: python legend_cleanup_fig2.py. Frozen input files are in source_archive, with checksums in source_manifest.json. Automated comparison code is verify_cleanup.py, with detailed results in figure_QA.json. A fresh regeneration leaves the visual-review status pending until the new whole image is inspected.
