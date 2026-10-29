# Artifact contract

Revise v43 Chinese manuscript following the approved six-panel Figure 2 plan. Use all existing five-fold out-of-fold predictions, with no new model training, inference, or MD. Figure 2 presents prediction/reference agreement, high-dynamic-site precision/recall, and within-amino-acid-type error gains. Figure labels remain English and prose remains Chinese.

Preserve all 10 native three-line tables, 36 native math objects, 28 references, comments/anchors, and Figures 1/3/4/5/S1. Do not alter the KRAS case, the recent biological backgrounds, or free-energy validation. Replace only the Figure 2 media part and tightly scoped prose. Use a dedicated new DOCX version; v43 is retained.

Statistical unit is a record. Methods use identical valid support. Ties are handled as a score group. The MD 90th percentile defines positives before evaluating predictions. Report eligibility and undefined cases. Step AP and AUROC are record-macro metrics, separate from a plotted macro PR curve's area. All points contribute to record-equal density; no display-tail trimming or downsampling.

Export using an owned invisible Word instance, render the resulting PDF with the packaged render_docx.py and inspect every final page image. Only publish the delivery pointer after numerical, package, scientific, and visual audits pass.
