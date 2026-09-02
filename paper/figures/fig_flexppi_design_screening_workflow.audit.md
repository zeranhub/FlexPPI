# FlexPPI design-screening workflow reconstruction audit

Reference: `codex-clipboard-6fd8e26c-7913-4788-9158-1e2c043371aa.png` (828 × 425 px)

## Style tokens

- Canvas: 828 × 425, solid white background, no border or texture.
- Typeface: Arial/Helvetica-like sans serif; dark charcoal text (`#222222`).
- Design module: three offset mint sheets; front fill `#A8DEB8`.
- Prediction module: three offset blue sheets; front fill `#B1C5E6`.
- Workflow arrows: thin charcoal strokes with compact filled arrowheads.
- Flexibility curves: target mint `#55C995`; FlexPPI prediction purple `#51458B`.

## Visible-element inventory

| id | approximate region / bbox | content / visual | medium | style notes | status |
|---|---|---|---|---|---|
| bg | full canvas, 0,0,828,425 | white page | native | white, no stroke | accepted |
| chart-target | upper-left, 18,28,195,150 | target flexibility profile chart, axes and labels | native | thin gray axes, mint smooth curve | accepted |
| step-1 | 215,64,47,53 | `(1)` and right arrow | native | 21 px label, thin charcoal arrow | accepted |
| design-stack | 273,14,269,155 | offset mint cards | native | two pale rear sheets, saturated front sheet | accepted |
| design-labels | 295,42,247,115 | `Protein Structure Design Models`; general-purpose design subtitle | native | centered title, italic subtitle | accepted |
| step-2 | 553,64,54,53 | `(2)` and right arrow | native | reference-matched geometry | accepted |
| design-candidates | 607,8,164,181 | four green/gray protein cartoons | crop | tight reference crop, white background matched to page | accepted |
| step-3 | 652,197,54,54 | `(3)` and down arrow | native | vertical transition | accepted |
| predictor-stack | 588,263,212,116 | offset blue cards | native | two pale rear sheets, blue front sheet | accepted |
| predictor-labels | 611,282,189,73 | `FlexPPI`; `Flexibility prediction` | native | centered title, italic subtitle | accepted |
| step-4 | 488,288,67,53 | `(4)` and left arrow | native | points from FlexPPI to candidate | accepted |
| candidate | 326,252,96,140 | blue/red candidate protein cartoon | crop | tight reference crop, white background matched to page | accepted |
| colorbar | 434,261,48,121 | red–white–blue flexibility scale with max/flex/min labels | native | stacked editable swatches | accepted |
| candidate-label | 334,393,109,26 | `Candidate` | native | centered 14 px label | accepted |
| step-5 | 230,288,89,53 | `(5)` and left arrow | native | points from candidate to comparison chart | accepted |
| chart-compare | lower-left, 18,252,195,151 | target vs FlexPPI flexibility profiles and legend | native | mint and purple smooth curves | accepted |

## Arrow inventory

| id | source → target | geometry | style | semantic role | status |
|---|---|---|---|---|---|
| a1 | target chart → design models | horizontal right | charcoal, filled classic head | provide target profile | accepted |
| a2 | design models → candidates | horizontal right | charcoal, filled classic head | generate candidate structures | accepted |
| a3 | candidate structures → FlexPPI | vertical down | charcoal, filled classic head | submit candidates for prediction | accepted |
| a4 | FlexPPI → scored candidate | horizontal left | charcoal, filled classic head | map predicted flexibility onto candidate | accepted |
| a5 | scored candidate → profile chart | horizontal left | charcoal, filled classic head | compare predicted and target profiles | accepted |

## Requested semantic substitutions

- Replaced the FlIPS-specific design block with the category `Protein Structure Design Models` and the subtitle `Diverse protein structure generators`.
- Replaced BackFlip with `FlexPPI` and kept `Flexibility prediction` as its role.
- Replaced the lower legend's `Molecular Dynamics` series with `FlexPPI prediction` so the output comparison matches the requested predictor.

## Verification result

- `check_drawio.py`: passed with exit code 0; 55 cells, 39 vertices, 14 edges, 2 embedded image cells, and 19 text cells.
- Draw.io CLI export: passed with exit code 0.
- Full-size visual comparison: all inventory items accepted; no missing element, text clipping, visible crop seam, generic substitute icon, or incorrect arrow direction remains.
