# Figure 2 legend cleanup visual QA

## Renderer

The workspace wrapper imports the packaged `render_docx.py` without modifying managed files and replaces only its DOCX-to-PDF conversion function with Word 16.0 read-only PDF export. A separate hidden COM instance opens the specified DOCX read-only, exports PDF, closes that document without saving, and quits only its own application. The bundled Python and bundled Poppler perform canonical page PNG output at 150 dpi. No installed LibreOffice, process termination or Word rewrite of the DOCX is used.

The fully inspected 30-page v65 final render in `figure2_relayout_20261024/render_v65_final` is the reference. Exact PNG identity is used to carry forward completed independent-page visual checks. When page numbering changes, a separate comparison masks only the verified centered footer-number box and retains all body pixels; that comparison is recorded separately from whole-page pixel identity.

## First v66 pre-final render

DOCX SHA-256: `0462fd62ea68d320e68e675f55f39d60c9108d23d98ec7590eb3ebfa13ca16223542`.

The first render completed with 32 pages in `render_v66_final`, but is not approved for delivery because of a new pagination defect.

- Figure 2 and its complete caption are visible on pages 11 and 12 respectively. The new c and e internal legends, simplified titles and interface-group legend are complete and legible.
- The longer Figure 2 caption pushes the four-line final KIX results paragraph to page 14. Figure 3 and its caption then move to page 15, leaving most of page 14 empty. This is a new layout defect relative to v65 and requires correction before final approval.
- Appendix H, the full Table S10 caption and all seven method rows are visible together on page 32 with no truncation or broken table rows.
- Pages 1–6 and 8–10 are exactly pixel-identical to the reviewed v65 pages. The body areas of v66 pages 16–31 are exactly pixel-identical to v65 pages 15–30; only centered footer page numbers differ.

Final review is pending a corrected DOCX and a new read-only export.
