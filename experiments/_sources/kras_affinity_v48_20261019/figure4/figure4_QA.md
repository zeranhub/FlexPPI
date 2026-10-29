# Figure 4 QA

Pass 0: Typography, color and export baseline blocks retained; no default palette, rainbow, four-sided border, redundant legend or raster-only export. State colors preserve their meaning across panels. The categorical difference is also encoded by open circle versus filled square markers.

Pass 1: The 183 × 130 mm master is within the target width and height. PDF uses TrueType embedding and SVG retains text. Axes use 0.6 pt left/bottom spines and outward ticks. PNG is 450 dpi. All nonempty visible text uses 9 pt or 10 pt panel letters, with conventional math subscript reduction.

Pass 2: Both state curves retain all positions 1–10 and absolute untransformed model predictions. The y-axis is linear and starts at zero. No curve fit, interpolation beyond straight sequential connecting lines, residue selection, refitting or MD comparison is displayed. State labels are placed at the endpoints in the reserved right margin. The mutation background identifies the chemically altered residue rather than a significance category.

Pass 3: Agent opened and inspected the exported original PNG. All labels and all 20 data markers are visible. Panel a's sequence, binding constants and superposition annotations remain legible without overlaps after increasing fonts. The two profiles remain separated from their right-end labels. Root independent review is pending and recorded in visual_QA.json.
