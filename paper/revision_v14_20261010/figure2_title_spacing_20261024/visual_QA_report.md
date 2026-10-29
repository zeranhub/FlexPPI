# Figure 2 title spacing final visual QA

Final document: `E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v68_20261024_中文修订_图2标题间距优化.docx`.

SHA-256: `352e7ac5a351ea5687e90ab07594ed59697bf92acbba2a949336d63af8627269`.

**Final result: 31 pages, all PASS.** Final internal QA files are in `render_v68_final`.

The workspace wrapper imports the packaged `render_docx.py` without editing managed files and replaces only DOCX-to-PDF conversion with a separate hidden Word 16.0 COM instance. The DOCX is opened read-only, exported without saving or rewriting, closed without saving, and only that Word instance is quit. Canonical rasterization and PNG page numbering use the bundled Python and bundled Poppler at 150 dpi. No installed LibreOffice or process termination is used.

The fully inspected 31-page v67 final export in `figure2_grouped_layout_20261024/render_v67_final` is the comparison reference. The new export also has 31 pages. Every page except page 11 is exactly identical to the corresponding v67 PNG, including all body pixels and footer page numbers. The completed independent-page visual checks recorded in the v67 QA report therefore apply unchanged to these 30 pages. This is a full-image pixel comparison, not a text-only or contact-sheet check. Per-page evidence is recorded in `render_v68_final/page_comparison_v67.json`.

Final page 11 was independently inspected as a full-page PNG by both the rendering reviewer and the root reviewer. The a and b titles and their letters sit closer to their plot areas. There is no collision between the a title and the upper MD Hχ label, and the b title remains clear of the violin plot. The whole figure is intact: a–e labels, numerical annotations, axes, legends, PR and ROC curves, and all three bottom bar charts are legible without clipping or overlap. The remaining layout and values match the approved grouped figure.

Figure 2 remains on page 11 and its full caption on page 12. Figure 3 and its caption remain on page 14. Figure S2 remains on page 25, Figure S3 on page 26, and Appendix H with the full Table S10 on page 31. All these unchanged pages are pixel-identical to their reviewed v67 counterparts.

| Final page | Result | Evidence |
|---|---|---|
| 1 | PASS | Exact whole-page PNG identity to reviewed v67 page 1. |
| 2 | PASS | Exact whole-page PNG identity to reviewed v67 page 2. |
| 3 | PASS | Exact whole-page PNG identity to reviewed v67 page 3. |
| 4 | PASS | Exact whole-page PNG identity to reviewed v67 page 4. |
| 5 | PASS | Exact whole-page PNG identity to reviewed v67 page 5. |
| 6 | PASS | Exact whole-page PNG identity to reviewed v67 page 6. |
| 7 | PASS | Exact whole-page PNG identity to reviewed v67 page 7. |
| 8 | PASS | Exact whole-page PNG identity to reviewed v67 page 8. |
| 9 | PASS | Exact whole-page PNG identity to reviewed v67 page 9. |
| 10 | PASS | Exact whole-page PNG identity to reviewed v67 page 10. |
| 11 | PASS | New independent full-page PNG views by both reviewers. |
| 12 | PASS | Exact whole-page PNG identity to reviewed v67 page 12. |
| 13 | PASS | Exact whole-page PNG identity to reviewed v67 page 13. |
| 14 | PASS | Exact whole-page PNG identity to reviewed v67 page 14. |
| 15 | PASS | Exact whole-page PNG identity to reviewed v67 page 15. |
| 16 | PASS | Exact whole-page PNG identity to reviewed v67 page 16. |
| 17 | PASS | Exact whole-page PNG identity to reviewed v67 page 17. |
| 18 | PASS | Exact whole-page PNG identity to reviewed v67 page 18. |
| 19 | PASS | Exact whole-page PNG identity to reviewed v67 page 19. |
| 20 | PASS | Exact whole-page PNG identity to reviewed v67 page 20. |
| 21 | PASS | Exact whole-page PNG identity to reviewed v67 page 21. |
| 22 | PASS | Exact whole-page PNG identity to reviewed v67 page 22. |
| 23 | PASS | Exact whole-page PNG identity to reviewed v67 page 23. |
| 24 | PASS | Exact whole-page PNG identity to reviewed v67 page 24. |
| 25 | PASS | Exact whole-page PNG identity to reviewed v67 page 25. |
| 26 | PASS | Exact whole-page PNG identity to reviewed v67 page 26. |
| 27 | PASS | Exact whole-page PNG identity to reviewed v67 page 27. |
| 28 | PASS | Exact whole-page PNG identity to reviewed v67 page 28. |
| 29 | PASS | Exact whole-page PNG identity to reviewed v67 page 29. |
| 30 | PASS | Exact whole-page PNG identity to reviewed v67 page 30. |
| 31 | PASS | Exact whole-page PNG identity to reviewed v67 page 31. |

No new clipping, overlap, missing glyph, broken table, footer defect or pagination change was found. The source hash was checked after export. The Word instance created for this export was closed through its own COM object.
