# Figure 2 relayout visual QA

## Renderer

The packaged `render_docx.py` was imported without changes to managed files. A workspace wrapper replaces only `convert_to_pdf` with Microsoft Word 16.0 read-only PDF export. The packaged renderer retains its PDF rasterization, page numbering, PNG output and optional PDF persistence behavior. The bundled Python and bundled Poppler binaries are used. No installed LibreOffice is used.

Word export creates a separate hidden COM application, opens only the specified manuscript read-only, exports PDF, closes that document without saving, and quits its own application. No process termination is used. PNG resolution is 150 dpi, approximately 1,241 × 1,754 pixels per A4 page.

## v64 reference

Source: `paper_v64_20261024_中文修订_图2加多链界面子图.docx`.

Reference rendering completed successfully: 29 PDF pages and 29 canonical `page-N.png` files in `render_reference`.

- Figure 2 occupies page 11 and is fully visible, with all original a–d regions intact.
- Figure 2 caption starts at the top of page 12. The image/caption page split already exists in the source manuscript.
- Supplementary section C, Table S4, Figure S2 and its caption share page 25. Figure S2 and its caption are fully visible and together.
- Reference pages 11, 12 and 25 have been individually inspected.

## Final v65

Source: `paper_v65_20261024_中文修订_图2取舍与外部PR_ROC正文版.docx`.

SHA-256: `6831ef35511ee1dac44f55f39d60c9108d23d98ec7590eb3ebfa13ca16223542`.

Final rendering completed successfully: 30 PDF pages and 30 canonical `page-N.png` files in `render_v65_final`. The earlier `render_v65` directory is a pre-final inspection run; use only `render_v65_final` for the final document. The final DOCX was opened read-only and was not rewritten by Word.

All 30 pages passed independent-page visual inspection at the rendered resolution. The primary reviewer checked pages 7–14 and 25–30, including the final page 10 after removal of one repeated sentence. A second reviewer checked each of pages 1–6 and 15–24 as independent full-resolution PNGs. These 16 pages are exactly equal pixel-for-pixel across the reference, initial v65, and final v65 exports. Pages 7–9 of the initial and final exports are also exactly equal; only page 10 changed between the two v65 exports, and final page 10 was inspected again.

| Final pages | Result |
|---|---|
| 1–6 | Pass. Text, figure, equations, and footers intact. PNGs exactly equal to v64. |
| 7–10 | Pass. Methods/results text reflow and Table 1 remain within page bounds; no clipped rows or text overlap. |
| 11 | Pass. Figure 2 entirely visible. Top row a/b preserved; second row c within-type Hχ, d external PR, e external ROC; bottom f retains all three interface-stratified bar charts. Legends, method names, AP/AUROC values and axis labels are legible and not cut off. |
| 12 | Pass. Complete Figure 2 caption at the top of the following page, matching the source manuscript's existing image/caption page split. Following results paragraphs are fully visible. |
| 13–14 | Pass. Reflow adds the final two lines of the KIX paragraph above Figure 3 on page 14; Figure 3 and its complete caption still remain together and within the page. |
| 15–24 | Pass. Tables, figures, figure labels and footers intact. PNGs exactly equal to v64. |
| 25 | Pass. Section C, Table S4, Figure S2 and its caption remain intact on the same page. |
| 26 | Pass. New Figure S3 contains the RMSF and overall Hχ PR panels and full legend. Its complete caption is on the same page. |
| 27–30 | Pass. Original supplementary sections D–G shift by one page; text, Tables S5–S9 and footers are intact. Normalized page text matches the original corresponding pages exactly. |

No new text clipping, object overlap, missing glyphs, broken tables or footer placement defects were found. No unexpected blank page was introduced. Section D retains its original new-page start, so the lower part of the Figure S3 page is intentionally unused.

The report's page references are verified against the final Word-exported PDF. `render_v65_final/page_comparison.json` records reference-to-final text similarity and exact PNG comparisons. `render_v65_final/final_vs_initial_pixels.json` records that final page 10 is the sole changed page since the preliminary v65 export. All Word COM instances created for export have been closed through their own COM objects; no WINWORD process remains from this work.

## Individual page inspection log

| Final page | Result | Inspection basis |
|---|---|---|
| 1 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 2 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 3 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 4 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 5 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 6 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 7 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 8 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 9 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 10 | PASS | Final PNG independently viewed again after sentence removal. |
| 11 | PASS | Final PNG independently viewed; full Figure 2. |
| 12 | PASS | Final PNG independently viewed; full Figure 2 caption. |
| 13 | PASS | Final PNG independently viewed. |
| 14 | PASS | Final PNG independently viewed; full Figure 3 and caption. |
| 15 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 16 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 17 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 18 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 19 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 20 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 21 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 22 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 23 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 24 | PASS | Independent PNG view; pixel identity confirms final equals reviewed preliminary PNG. |
| 25 | PASS | Final PNG independently viewed; S2 remains with caption. |
| 26 | PASS | Final PNG independently viewed; new S3 remains with caption. |
| 27 | PASS | Final PNG independently viewed. |
| 28 | PASS | Final PNG independently viewed. |
| 29 | PASS | Final PNG independently viewed. |
| 30 | PASS | Final PNG independently viewed. |
