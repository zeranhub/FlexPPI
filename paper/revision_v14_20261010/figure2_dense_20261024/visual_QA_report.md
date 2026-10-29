# Figure 2 dense layout final visual QA

Final document: `E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v70_20261024_中文修订_图2紧凑终版.docx`.

Final SHA-256: `5e8741e4be591e9a3d45aebed460c0bbeab844f4dfc8e14e52312d5505618841`.

Final figure PNG SHA-256 supplied and verified by the root/figure workflow: `3695952ad99dca46f341ab4dc0f80871a4f3f1eaf24b5df9be9b1f182063f9e6`.

**Final result: 31 pages, all PASS.** The final internal rendering is `render_v70_final_r2`. The older `render_v70_final` belongs to the superseded first image and must not be treated as the final manuscript render.

## Renderer and comparison

The workspace wrapper imports the packaged `render_docx.py` without modifying managed files and replaces only DOCX-to-PDF conversion with a separate hidden Word 16.0 COM instance. It opens the manuscript read-only, exports PDF, closes the document without saving, and quits only its own instance. Canonical rasterization and page naming use the bundled Python and bundled Poppler at 150 dpi, approximately 1,241 × 1,754 pixels per A4 page. No installed LibreOffice, process termination or Word save/rewrite is used.

The reference is the fully inspected 31-page v68 export in `figure2_title_spacing_20261024/render_v68_final`. The final r2 export also contains 31 pages. Actual whole-page comparison found changes only on pages 11, 12 and 13. Pages 1–10 and 14–31 are exactly identical to their reviewed v68 PNGs, including body pixels and footer page numbers; completed independent-page visual checks are carried forward for these 28 pages. This evidence comes from full-image pixel comparison, not text matching or contact sheets. `render_v70_final_r2/page_comparison_v68.json` records every page comparison.

## Changed pages

Final pages 11, 12 and 13 were independently viewed at the rendered full-page resolution by the rendering reviewer. The root reviewer also independently viewed final r2 page 11 and confirmed the complete figure and caption. All changed pages passed.

- Page 11 contains the complete compact Figure 2, its complete caption, and the following short results paragraph, all on the same page. The a–e letters, axes, numerical labels, curves, legends and three bar charts are intact. The b legend uses the vacant lower region without covering data; the c and d legends remain legible; the shared interface-group legend and e letter are separated from adjacent row labels. The a title and MD Hχ upper-axis label do not collide. No item is clipped or distorted.
- Page 12 starts with the source-stratified screening results and continues through the interface-stratified and external-comparison results to the beginning of Section 3.2. Paragraphs are fully visible, with no overlap or lost lines.
- Page 13 continues Section 3.2 and includes the complete KIX results text. It ends before the unchanged dedicated Figure 3 page. The lower unused area reflects the existing new-page placement of Figure 3; it does not create an extra page or interrupt a paragraph. Figure 3 remains intact on page 14, which is pixel-identical to the reviewed v68 page 14.

Supplementary figures and tables retain their approved pagination: S2 on page 25, S3 on page 26, and Appendix H/Table S10 on page 31. These pages are pixel-identical to the reviewed v68 counterparts.

## Individual page record

| Final page | Result | Evidence |
|---|---|---|
| 1 | PASS | Exact whole-page PNG identity to reviewed v68 page 1. |
| 2 | PASS | Exact whole-page PNG identity to reviewed v68 page 2. |
| 3 | PASS | Exact whole-page PNG identity to reviewed v68 page 3. |
| 4 | PASS | Exact whole-page PNG identity to reviewed v68 page 4. |
| 5 | PASS | Exact whole-page PNG identity to reviewed v68 page 5. |
| 6 | PASS | Exact whole-page PNG identity to reviewed v68 page 6. |
| 7 | PASS | Exact whole-page PNG identity to reviewed v68 page 7. |
| 8 | PASS | Exact whole-page PNG identity to reviewed v68 page 8. |
| 9 | PASS | Exact whole-page PNG identity to reviewed v68 page 9. |
| 10 | PASS | Exact whole-page PNG identity to reviewed v68 page 10. |
| 11 | PASS | New independent final r2 PNG views by rendering and root reviewers. |
| 12 | PASS | New independent final r2 PNG view; results text reflow intact. |
| 13 | PASS | New independent final r2 PNG view; Section 3.2 and KIX text intact. |
| 14 | PASS | Exact whole-page PNG identity to reviewed v68 page 14. |
| 15 | PASS | Exact whole-page PNG identity to reviewed v68 page 15. |
| 16 | PASS | Exact whole-page PNG identity to reviewed v68 page 16. |
| 17 | PASS | Exact whole-page PNG identity to reviewed v68 page 17. |
| 18 | PASS | Exact whole-page PNG identity to reviewed v68 page 18. |
| 19 | PASS | Exact whole-page PNG identity to reviewed v68 page 19. |
| 20 | PASS | Exact whole-page PNG identity to reviewed v68 page 20. |
| 21 | PASS | Exact whole-page PNG identity to reviewed v68 page 21. |
| 22 | PASS | Exact whole-page PNG identity to reviewed v68 page 22. |
| 23 | PASS | Exact whole-page PNG identity to reviewed v68 page 23. |
| 24 | PASS | Exact whole-page PNG identity to reviewed v68 page 24. |
| 25 | PASS | Exact whole-page PNG identity to reviewed v68 page 25. |
| 26 | PASS | Exact whole-page PNG identity to reviewed v68 page 26. |
| 27 | PASS | Exact whole-page PNG identity to reviewed v68 page 27. |
| 28 | PASS | Exact whole-page PNG identity to reviewed v68 page 28. |
| 29 | PASS | Exact whole-page PNG identity to reviewed v68 page 29. |
| 30 | PASS | Exact whole-page PNG identity to reviewed v68 page 30. |
| 31 | PASS | Exact whole-page PNG identity to reviewed v68 page 31. |

No new clipping, abnormal overlap, missing glyph, broken table, footer defect or extra page was found. Final document hash was checked after the read-only export. All Word instances created for this export were closed through their own COM objects.
