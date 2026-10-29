# Figure 4 antibody extension first preview

The user requested an actual first figure draft for the 6A77 narrative, to be reviewed before insertion into the manuscript. The existing TCR panels a–d are retained. New panels e–g extend the evidence to the bound B5209B Fab–Robo1 Ig5 complex.

Core claim: residue-level backbone contacts and side-chain state diversity can coexist at the same antibody recognition-loop position. Full-loop predictions and their existing MD reference provide the dynamic context; atom-specific contact analysis explains the distinction at R101. No claim of affinity mechanism, energetic hotspot, mutation optimization, binding-induced rigidification, or a complete membrane BCR.

- e: actual canonical molecular coordinates locate R96/R101 and the distinct local contact environments. Static distances are not dynamic hydrogen-bond occupancies.
- f: all 15 H-CDR3 positions, with RMSF and Hχ on separate aligned axes. Full predictions versus the frozen MD reference; all 12 valid Hχ positions retained, undefined Ala values left blank, valid Pro values retained. R96/R101 are highlighted without omitting other positions.
- g: MD frame fraction with an R101 backbone/light-chain contact versus side-chain/light-chain contact, minimum heavy-atom distance <=4 A, 90.3% versus 2.3%. These are two descriptive proportions of one correlated 1000-frame trajectory, with no error bars, p values, or independent-frame claim.

Archetype: asymmetric mixed-modality. The structure and complete residue profiles carry most area. The small occupancy panel supplies different evidence and does not duplicate the Arg entropy values already visible in f.

Target: 183 mm double-column figure, RGB, Arial, editable vector chart text/lines (PDF/SVG), actual molecular render as embedded raster, 300 dpi PNG. The added row is approximately 80 mm high; old 152 mm TCR figure plus the new row stays below the 247 mm graphic height guide. No new image will be embedded into the DOCX in this revision.

Asset scan: academic-figure-skill LineTrend/plot_trend.py and BarCategorical/plot_correctness_by_category.py and their companion PNGs were inspected. The first takes a hard-coded cumulative monthly publication-count/event dataset; the second a hard-coded multi-method/category correctness matrix. These semantics are incompatible with residue profiles and one-trace contact proportions. Parameter inheritance uses clean left/bottom spines, directly labeled values, and the current Figure 4 scientific style; no template data or statistical decorations are reused. Existing Figure 4 colors and typography are the user-provided visual authority.

Anticipated risks: tiny residue labels, close neighboring R95/R96 tick labels, structure label overlap, and confusion between static Y49 near-contact and dynamic contact with any light-chain residue. Use sparse numeric x ticks while plotting every residue, place R96/R101 labels above profiles, editable structure labels with leaders, and explicit any-light-chain definition in the caption. All masks, atom definitions, thresholds, model hashes and source tables are retained in the source-data/QA reports.

Boundary: R96/R101 entropy differences are also recovered by Intra and are not proof of a Full-exclusive advantage. This is a retrospective local-environment illustration. Existing whole-six-CDR results remain in the supplement; not all six loops improve. The local physical trajectory timescale remains uncalibrated; no three-replicate or convergence claim.
