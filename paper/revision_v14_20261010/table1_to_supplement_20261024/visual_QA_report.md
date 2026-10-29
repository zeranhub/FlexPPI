# Table 1 relocation final visual QA

Final document: `E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v71_20261024_中文修订_表1对齐图2并移补充.docx`.

Final SHA-256: `4bf956d86b66ba898f62ad48ff4aee0f41c9269d52e7585f4825672601bca132`.

**Final result: 32 pages, all PASS.** The final PDF and canonical PNG pages are in `render_v71_final_r3`. Earlier `render_v71_final` and `render_v71_final_r2` contain superseded pagination and are internal diagnostics only.

## Renderer

The workspace wrapper imports the packaged `render_docx.py` without modifying managed files and replaces only DOCX-to-PDF conversion with a separate hidden Word 16.0 COM instance. The specified manuscript is opened read-only, exported to PDF, closed without saving and only that Word instance is quit. Canonical PNG rasterization uses bundled Python and bundled Poppler at 150 dpi, approximately 1,241 × 1,754 pixels per A4 page. No installed LibreOffice, process termination or Word save/rewrite is used.

## Review evidence

The approved 31-page v70 export in `figure2_dense_20261024/render_v70_final_r2` is the original visual reference. The first v71 export had 33 pages. Its genuinely changed pages 6, 7, 9, 10, 11, 12, 25, 26 and 33 received new independent full-page PNG inspection. The manuscript-audit reviewer checked pages 6, 7, 9, 12, 25 and 26; the rendering and root reviewers checked pages 10, 11 and 33. Unchanged pages were verified by full-image identity or exact body-pixel identity with changed footer numbers. The r1-to-r2 repair removed duplicate OOXML paragraph identifiers; all 33 r2 PNGs were exactly identical to r1, as recorded in `render_v71_final_r2/r1_vs_r2_pixels.json`.

The r2 layout left a large blank area on page 10 after Table 1 was removed. This issue was reported and fixed by moving the intact Figure 2 image/caption group after the associated five results paragraphs. No text, figure data, table values or element formatting was changed by that order-only repair. The final r3 export has 32 pages.

Actual r2-to-r3 full-image comparison found:

- Final pages 1–9 are exactly identical to r2 pages 1–9, including page numbers.
- Final pages 10–12 have new body layouts and were independently viewed at full-page resolution by both the rendering reviewer and the root reviewer.
- Final pages 13–32 have exactly identical body pixels to r2 pages 14–33. Only the centered footer numbers change by minus one. The already completed independent body checks therefore carry forward. All final footer numbers were verified against PDF page indices and visually inspected in a footer-only montage. The printed sequence is continuous 1–32, with no missing, duplicated or misplaced number.

The exact mappings are saved in `render_v71_final_r3/render_run_comparison.json`. Printed-number checks are in `page_number_audit.json`, and the visual footer audit is `footer_numbers_1_to_32.png`. The initial v70-to-v71 mapping is recorded in `render_v71_final/page_comparison_v70.json`.

## Final layout findings

Page 10 is filled by the internal, source-stratified and complex-comparison results; the former removal-induced blank area is resolved. Figure 2 and its complete caption are together on page 11, and the following paragraph is fully visible. Page 12 flows naturally through the remaining results, Section 3.2 and Section 3.3. Figure 3 remains complete on its dedicated page 13.

Supplementary section C now spans pages 24–25. Table S4 continues with a repeated column header and retains its last rows and bottom rule. Figure S2 and its complete caption stay together on page 25. Figure S3 and its complete caption occupy page 26. Appendix H/Table S10 is intact on page 31.

Appendix I/Table S11 fits entirely on page 32. Its caption, two-level RMSF/Hχ column header, all three A/B/C blocks, all 20 table rows, bottom rule and complete explanatory note are visible on that page. The unavailable-value em dashes and Å/nats units are preserved; no cell text or note line is clipped. The root and manuscript-audit reviewers separately verified model order, cohort labels, values and references against the manuscript/data audit.

## Individual final-page record

| Final page | Result | Inspection basis |
|---|---|---|
| 1 | PASS | Exact whole-page identity to r2; approved v70 page 1. |
| 2 | PASS | Exact whole-page identity to r2; approved v70 page 2. |
| 3 | PASS | Exact whole-page identity to r2; approved v70 page 3. |
| 4 | PASS | Exact whole-page identity to r2; approved v70 page 4. |
| 5 | PASS | Exact whole-page identity to r2; approved v70 page 5. |
| 6 | PASS | New independent methods-page view in r1; exact r1/r2/r3 identity. |
| 7 | PASS | New independent methods-page view in r1; exact r1/r2/r3 identity. |
| 8 | PASS | Exact whole-page identity to r2; approved v70 page 8. |
| 9 | PASS | New independent results-page view in r1; exact r1/r2/r3 identity. |
| 10 | PASS | New independent final r3 views; former large blank area resolved. |
| 11 | PASS | New independent final r3 views; complete Figure 2 and caption together. |
| 12 | PASS | New independent final r3 views; results and section flow intact. |
| 13 | PASS | Exact body-pixel identity to reviewed r2 page 14; new footer 13 verified. |
| 14 | PASS | Exact body-pixel identity to reviewed r2 page 15; new footer 14 verified. |
| 15 | PASS | Exact body-pixel identity to reviewed r2 page 16; new footer 15 verified. |
| 16 | PASS | Exact body-pixel identity to reviewed r2 page 17; new footer 16 verified. |
| 17 | PASS | Exact body-pixel identity to reviewed r2 page 18; new footer 17 verified. |
| 18 | PASS | Exact body-pixel identity to reviewed r2 page 19; new footer 18 verified. |
| 19 | PASS | Exact body-pixel identity to reviewed r2 page 20; new footer 19 verified. |
| 20 | PASS | Exact body-pixel identity to reviewed r2 page 21; new footer 20 verified. |
| 21 | PASS | Exact body-pixel identity to reviewed r2 page 22; new footer 21 verified. |
| 22 | PASS | Exact body-pixel identity to reviewed r2 page 23; new footer 22 verified. |
| 23 | PASS | Exact body-pixel identity to reviewed r2 page 24; new footer 23 verified. |
| 24 | PASS | Exact body-pixel identity to independently reviewed expanded-C page 25; new footer 24 verified. |
| 25 | PASS | Exact body-pixel identity to independently reviewed S4/S2 page 26; new footer 25 verified. |
| 26 | PASS | Exact body-pixel identity to reviewed r2 S3 page 27; new footer 26 verified. |
| 27 | PASS | Exact body-pixel identity to reviewed r2 page 28; new footer 27 verified. |
| 28 | PASS | Exact body-pixel identity to reviewed r2 page 29; new footer 28 verified. |
| 29 | PASS | Exact body-pixel identity to reviewed r2 page 30; new footer 29 verified. |
| 30 | PASS | Exact body-pixel identity to reviewed r2 page 31; new footer 30 verified. |
| 31 | PASS | Exact body-pixel identity to reviewed r2 H/S10 page 32; new footer 31 verified. |
| 32 | PASS | Exact body-pixel identity to independently reviewed r2 S11 page 33, plus root's direct final r3 view; new footer 32 verified. |

No new clipping, abnormal overlap, missing glyph, broken table row or footer defect remains. The manuscript's rendering and structural checks are complete. Final document hash was checked after export; all Word instances created for export were closed through their own COM objects.
