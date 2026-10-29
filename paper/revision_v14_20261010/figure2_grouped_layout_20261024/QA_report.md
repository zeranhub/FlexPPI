# Figure 2 grouped layout QA

READY. Final main figure is 183 × 214.44 mm, with 450 dpi PNG, native vector PDF and editable-text SVG. The panel sequence is a, b, c, d, e: external PR/ROC share one d label, and all three bottom bars share e.

Layout: original a/b reused unchanged, including b's independent Full/Intra/Static/AA horizontal legend. C retains its internal legend. The PR axis moves 5 mm right, reducing its separation from ROC to 14 mm compared with the 24 mm c-to-PR separation; both external axes retain 42 × 40 mm size and short titles. ROC retains the shared seven-method legend. The bottom group's original long legend is removed through its precisely identified final native drawing block. A concise three-item legend with the same blue colors and effective font size is placed centrally at local y = 9 mm, just above the bars and separate from e's left-side letter. Counts remain in the caption.

The full final image was inspected. No data occlusion, clipped labels, added small-print statistics, or legend collision was observed. The paired external axes have matching heights. AP/cohort annotation bands and unnecessary counts remain absent. No statistics, axes, line colors, styles, markers or data are changed.

At 160 dpi, a/b including its legend are pixel-identical to the original native crop. All three bottom bar panels, axes, error bars and method labels are pixel-identical outside the relocated group legend and the d/e panel letter. Native PDF path comparisons verify all four c and seven external PR/ROC curves are identical to v65 after only the intended layout translations. PDF and SVG contain zero raster image objects. The archived source-file checksums are unchanged.

Reproduction: python grouped_layout_fig2.py. Source files and manifest are self-contained in this folder. Automated layout/data verification is in verify_grouped.py; detailed results are in figure_QA.json. Supplement S3 is unchanged. External AP/AUROC numbers remain in supplementary Table S10 handled by the document editor.
