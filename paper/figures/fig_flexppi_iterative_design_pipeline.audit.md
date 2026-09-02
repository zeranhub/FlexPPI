# Visual inventory and audit: FlexPPI iterative design pipeline

## Scope

- Target: `fig_flexppi_iterative_design_pipeline.drawio`
- Preview: `fig_flexppi_iterative_design_pipeline.png`
- Figure type: Solution Overview, circular iterative pipeline
- Canvas: 1400 × 980 drawing units, white background
- Reference style: existing FlexPPI paper figures; flat rounded cards, Okabe–Ito colours, no shadows or gradients

## Visible-element inventory

| ID | Approximate region | Content / visual | Medium | Style notes | Status |
|---|---|---|---|---|---|
| bg-01 | full canvas | White page and faint circular guide ring | native | white background; light-grey dashed ring | accepted |
| hub-01 | centre | “FlexPPI 引导的迭代蛋白设计” and four-part summary | native | white ellipse, purple outline, centred text | accepted |
| step-01 | top | 从头设计起点 | native | blue rounded card; target/epitope icon | accepted |
| step-02 | upper right | 现有工具生成候选 | native | orange rounded card; three candidate symbols | accepted |
| step-03 | lower right | FlexPPI 第一次柔性分析 | native | purple rounded card; residue heat strip | accepted |
| step-04 | bottom | 高柔性候选位点筛选与替换 | native | vermillion rounded card; highlighted sequence substitution | accepted |
| step-05 | lower left | AlphaFold 结构预测与过滤 | native | green rounded card; apo/complex symbol and check | accepted |
| step-06 | upper left | FlexPPI 再次柔性分析 | native | purple rounded card; before/after heat strips | accepted |
| edge-01 | top to upper right | Step 1 → Step 2 | native connector | grey, curved, clockwise, block arrow | accepted |
| edge-02 | right side | Step 2 → Step 3 | native connector | grey, curved, clockwise, block arrow | accepted |
| edge-03 | lower right | Step 3 → Step 4 | native connector | grey, curved, clockwise, block arrow | accepted |
| edge-04 | lower left | Step 4 → Step 5 | native connector | grey, curved, clockwise, block arrow | accepted |
| edge-05 | left side | Step 5 → Step 6 | native connector | grey, curved, clockwise, block arrow | accepted |
| edge-06 | upper left to top | Step 6 → Step 1 | native connector | vermillion dashed feedback, label “下一轮设计” | accepted |
| note-01 | inside step 4 | Protection rule for catalytic, binding, and functional-flexibility sites | native text | explicit caution; no colour-only meaning | accepted |
| outcome-01 | inside step 6 | 达标候选进入备选池；未达标继续迭代 | native pill/text | green acceptance cue plus text | accepted |

## Arrow inventory

| ID | Source | Target | Route | Semantic role | Status |
|---|---|---|---|---|---|
| edge-01 | step-01 | step-02 | outer upper-right arc | initial candidate generation | accepted |
| edge-02 | step-02 | step-03 | outer right arc | candidate structures to first FlexPPI pass | accepted |
| edge-03 | step-03 | step-04 | outer lower-right arc | high-flexibility profile to redesign | accepted |
| edge-04 | step-04 | step-05 | outer lower-left arc | variants to AlphaFold | accepted |
| edge-05 | step-05 | step-06 | outer left arc | predicted structures to second FlexPPI pass | accepted |
| edge-06 | step-06 | step-01 | outer upper-left arc | next design round | accepted |

## Style tokens

- Target / initial design blue: `#0072B2`, tint `#EAF3FB`
- Candidate generation orange: `#E69F00`, tint `#FFF4D6`
- FlexPPI purple: `#CC79A7`, tint `#F4EAF5`
- Redesign vermillion: `#D55E00`, tint `#FBE9E5`
- AlphaFold / accepted green: `#009E73`, tint `#E8F5F0`
- Neutral text: `#1F2937`; neutral connector: `#6B7280`
- Font: Microsoft YaHei; card title 21–22 px; body 17–18 px; centre title 28 px
- No gradients, shadows, 3D effects, decorative icons, or raster images

## Audit result

- Draw.io structural check passed: 63 cells, 54 vertices, 7 edges, and no embedded image/SVG cells.
- PNG preview inspected at full exported size (2802 × 1962 px): all six cards, centre hub, icons, arrows, arrowheads, feedback label, and circular guide ring are visible.
- No text overlap, clipping, cropped arrowhead, background seam, gradient, shadow, or raster element was found after the final two visual iterations.
- The bottom redesign card was widened and its wording reflowed; the AlphaFold title and first FlexPPI body text were adjusted to remove accidental wrapping.
- Minimum text size is 17 px in the 1400-unit source canvas, corresponding to approximately 8.2 pt at 178 mm final width.
- PDF and SVG vector exports completed; PNG is a preview only.

### Figure Designer gate

- Vector format: pass.
- Font size: pass at the specified 178 mm width; recheck if inserted narrower.
- Colour-blind safe: pass; Okabe–Ito colours are paired with numbers, labels, and line styles.
- Self-contained caption: user still needs to place the caption in the manuscript.
- Honest axis range: not applicable.
- No chartjunk: pass.
- Integrity gate 1–5 and 7: pass; gate 6 is not applicable because this is a Solution Overview rather than a Motivated Example.
- Severity: 0 CRITICAL, 0 MAJOR, 1 MINOR. Manual review point: confirm final manuscript insertion width is not below 178 mm without enlarging the source typography.
