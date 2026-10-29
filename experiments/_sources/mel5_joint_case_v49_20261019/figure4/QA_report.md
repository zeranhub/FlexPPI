# Figure 4 verification

The final PNG was opened and visually reviewed after the last generation. All visible text is English and at least9pt; the figure is183×185mm, exported as450dpi PNG and editable-text SVG plus vector PDF. The two model panels retain all45 aligned CDR positions and all23 valid, unchanged-amino-acid pMHC interface positions. The predefined30-site interface source is archived in full; its7 residues without valid side-chain chi angles are not plotted in the entropy panel. Data limits were checked independently of text limits.

The structure is the actual existing4JFF model input, with common-numbering betaPhe53 and betaGly54 mapped to inputE54/E55. It is used only to interpret the model's residue localization. Distances and atom labels identify Q72OE1–Phe53N and Q72NE2–Gly54O. No molecular coordinates were edited.

The MD panel has three visible replica points per state per criterion; short vertical marks and right-hand numbers are unweighted means. It uses existing ITrace trajectories and is explicitly labeled ExistingMD. No new production, fitting, model inference, training or calibration was done. No statistical test is applied. See statistics_and_reproducibility.json for n definitions, criteria, provenance and source hashes.

The production library was inspected for LineTrend, PairedBoxScatter and GroupedBarChart, including actual preview images. Existing production scripts use fundamentally different data structures, so their thin open axes, paired connectors and individual-point logic were inherited rather than executing those unrelated examples. The real-coordinate PyMOL script was adapted from the previous manuscript's molecular renderer.
