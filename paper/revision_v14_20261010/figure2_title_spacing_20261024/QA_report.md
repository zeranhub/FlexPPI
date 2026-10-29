# Figure 2 title spacing QA

READY. Final dimensions: 183 × 205.44 mm, shortened from 214.44 mm. Only the a/b title and panel-letter placement matrices were edited: each moves down 9 mm, then 9 mm is cropped from the entire figure's top white space. Plot/axis/legend sizes, all curve and bar data, colors and fonts are unchanged. Every other panel simply moves up 9 mm with the crop.

Measured title-bottom to axis-top gap for both a and b: 14.95 mm before, 5.95 mm after. A's main title and its upper MD Hχ axis title retain 12.95 mm horizontal clearance; their text boxes do not intersect. The two main titles retain the same baseline and their letters move with them.

The complete final preview was inspected. No title obscures coordinates or labels, no top whitespace remains unnecessarily, and all a–e grouping, b's independent color key, paired PR/ROC in d, and e's centered blue legend remain intact.

Native integrity: the complete original native XObject streams are retained, with exactly four text-placement matrix substitutions in the top source stream. Reversing those four substitutions reconstructs the original stream byte for byte. No plot is redrawn or rasterized; the resulting PDF/SVG each contain zero raster image objects. The full archived source PDF and SHA-256 manifest are in this folder. Reproduction: python adjust_title_spacing.py. Preview PNG is 450 dpi; the smaller preview is 160 dpi. Detailed measurements are recorded in title_spacing_QA.json.
