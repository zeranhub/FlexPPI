# Figure 2 grouped layout final visual QA

## Final document and renderer

Document: `E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v67_20261024_中文修订_图2组合子图与图例重排.docx`.

SHA-256: `e7a6cec9918daf2250d9273f534766ac29b2b4a12a9c9ed6bd4f40fd65b86b36`.

Final result: **31 pages, all PASS**. The final PDF and canonical `page-N.png` files are in `render_v67_final`. These are internal quality-control outputs.

The workspace wrapper imports the packaged `render_docx.py` without modifying managed files and replaces only DOCX-to-PDF conversion with Word 16.0 read-only PDF export. A separate hidden COM instance opens the specified DOCX read-only, exports PDF, closes that document without saving, and quits only its own application. Bundled Python and bundled Poppler retain the packaged renderer's rasterization and page naming at 150 dpi, approximately 1,241 × 1,754 pixels per A4 page. No installed LibreOffice, process termination, or Word save/rewrite of the manuscript is used.

## Inspection and comparison evidence

All 30 pages of the reference `figure2_relayout_20261024/render_v65_final` had previously received independent full-page PNG visual inspection, recorded in that directory's QA report. Final v67 pages 1–6, 8–10, 15–24 and 26–30 are exactly pixel-identical to the corresponding reviewed v65 pages, including footer page numbers. Those 24 completed independent checks are carried forward by exact pixel identity. The comparison does not rely only on matching text or contact sheets.

The seven changed or new pages, 7, 11, 12, 13, 14, 25 and 31, received new independent full-page PNG inspection. The rendering reviewer viewed pages 7, 11, 13, 25 and 31; the root reviewer independently viewed pages 11, 12, 14 and 31. All passed. `render_v67_final/page_comparison_v65.json` records every page's reference match, text similarity, whole-page PNG identity and optional body comparison excluding only the footer number. Whole-page identity, rather than the optional footer mask, is the evidence used for all carried-forward checks in this final v67 review.

The main figure and caption remain on pages 11 and 12, following the original manuscript's existing image/caption page split. Figure 3 and its complete caption are on page 14; the mostly empty page found in the abandoned v66 pre-final export is resolved. Figure S2 and caption remain together on page 25, Figure S3 and caption remain together on page 26, and Appendix H plus the complete Table S10 caption and all seven method rows fit on page 31. No unexpected blank page was introduced.

Figure 2 uses letters a–e. The b panel has its own four-method legend; c retains the within-type Hχ PR with its four-method internal legend; d groups the external RMSF PR and ROC axes under one letter, with a complete seven-method legend inside the ROC axis; e retains all three interface-stratified bar charts with a shared centered blue-group legend close above them. No panel labels, axis labels, method names, bars or legend items are clipped. The simplified figure contains no old numerical AP/AUROC legend table; the values are fully present in Table S10.

The root reviewer separately verified the a–e references in the document XML and the Table S10 values. The page images show the Section C reference to external PR/ROC as Figure 2d and the Appendix H/Table S10 reference to Figure 2d. The results use Figure 2e for interface-stratified bar charts.

## Individual page log

| Final page | Result | Evidence |
|---|---|---|
| 1 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 1. |
| 2 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 2. |
| 3 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 3. |
| 4 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 4. |
| 5 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 5. |
| 6 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 6. |
| 7 | PASS | New independent PNG view; methods text and headings intact. |
| 8 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 8. |
| 9 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 9. |
| 10 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 10. |
| 11 | PASS | New independent PNG views by both reviewers; complete Figure 2 a–e. |
| 12 | PASS | New independent PNG view by root reviewer; complete Figure 2 caption and results text. |
| 13 | PASS | New independent PNG view; paragraph reflow intact and KIX text remains on this page. |
| 14 | PASS | New independent PNG view by root reviewer; complete Figure 3 and caption, no abnormal empty page. |
| 15 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 15. |
| 16 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 16. |
| 17 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 17. |
| 18 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 18. |
| 19 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 19. |
| 20 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 20. |
| 21 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 21. |
| 22 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 22. |
| 23 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 23. |
| 24 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 24. |
| 25 | PASS | New independent PNG view; Section C, Table S4, Figure S2 and complete caption fit together. |
| 26 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 26; Figure S3 and caption. |
| 27 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 27. |
| 28 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 28. |
| 29 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 29. |
| 30 | PASS | Exact whole-page PNG identity to independently reviewed v65 page 30. |
| 31 | PASS | New independent PNG views by both reviewers; Appendix H and full Table S10 intact. |

No new clipping, abnormal overlap, missing glyphs, broken table rows or footer placement defects were found. Final source hash was checked after export, and all Word instances created for this work were closed through their own COM objects.
