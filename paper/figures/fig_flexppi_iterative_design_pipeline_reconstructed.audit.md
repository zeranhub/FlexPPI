# FlexPPI iterative design pipeline reconstruction audit

## Scope and source

- Status: accepted after Draw.io structural validation, 2× PNG inspection, and successful PDF/SVG export.
- Reference: `C:/Users/DELL/AppData/Local/Temp/codex-clipboard-ff83db32-911b-4be2-88e5-a90e7edad978.png`
- Reference dimensions: 1448 × 1086 px.
- Reconstruction canvas: 1448 × 1086 Draw.io units with a white background.
- Fidelity target: reproduce the supplied screenshot, including its Chinese wording, six-card circular composition, color coding, soft shadows, center oval, iconography, mini-plots, and freehand-looking arrows.

## Visible-element inventory

### Global composition

- [x] White 4:3-ish landscape canvas with generous margins.
- [x] Six rounded rectangular cards arranged clockwise around a large central oval.
- [x] Soft grey drop shadows under all six cards and the central oval.
- [x] Card fills are near-white with a very subtle tint matching each card border.
- [x] Six clockwise transition arrows remain visually behind the cards and center oval.

### Center summary

- [x] Large white horizontal oval with thin pale-grey outline and soft shadow.
- [x] Main title: `FlexPPI 引导的迭代蛋白设计` in dark navy and bold.
- [x] Subtitle: `生成 → 诊断 → 定点改造 → 复评` in medium grey.
- [x] Green outlined rounded pill: `降低非期望柔性 · 保留功能性柔性`.

### Card 1 — blue, upper center

- [x] Blue rounded border and pale-blue fill.
- [x] Blue circular badge with white `1`.
- [x] Title: `蛋白 de novo 设计`.
- [x] Bullets: `目标结构与表位`; `功能约束、期望柔性轮廓`.
- [x] Concentric blue target/crosshair icon with radial blue center.

### Card 2 — amber, upper right

- [x] Amber rounded border and pale-amber fill.
- [x] Amber circular badge with white `2`.
- [x] Title: `现有工具生成候选`.
- [x] Bullets: `RFdiffusion / BindCraft`; `ProteinMPNN 等序列设计工具`; `得到多样化初始蛋白`.
- [x] Three short blue, orange, and green horizontal rounded bars aligned on the right.

### Card 3 — purple, lower right

- [x] Purple rounded border and pale-purple fill.
- [x] Purple circular badge with white `3`.
- [x] Title: `FlexPPI 柔性分析`.
- [x] Bullets: `预测 RMSF、侧链熵`; `定位高柔性候选位点`; `同时保留预测不确定性`.
- [x] Right-side five-bin vertical red–white–blue scale with the red label `高`.

### Card 4 — orange, bottom center

- [x] Orange rounded border and pale-orange fill.
- [x] Orange circular badge with white `4`.
- [x] Title: `高柔性候选位点筛选与替换`.
- [x] Bullets: `Top-k 高柔性残基 → 约束性替换`; `排除催化位点、关键界面残基`; `以及功能性铰链`; `得到一组备选变体蛋白`.
- [x] Upper grey sequence box: `A  G  L  V  Y`.
- [x] Small orange downward replacement arrow.
- [x] Lower orange outlined sequence box: `A  G  I  V  Y`, with `I` highlighted in orange.

### Card 5 — green, lower left

- [x] Green rounded border and pale-green fill.
- [x] Green circular badge with white `5`.
- [x] Title: `AlphaFold 结构预测与过滤`.
- [x] Bullets: `预测备选蛋白及复合物结构`; `检查置信度、clash 与界面保持`; `保留结构可信的候选`.
- [x] Protein ribbon rendering and overlapping green check badge at lower right.

### Card 6 — magenta, upper left

- [x] Magenta rounded border and pale-magenta fill.
- [x] Magenta circular badge with white `6`.
- [x] Title: `FlexPPI 再次柔性分析`.
- [x] Bullets: `比较改造前后的柔性轮廓`; `同时检查不确定性与结构一致性`; `达标候选进入备选池`.
- [x] `前` and `后` labels, magenta-to-blue before/after flexibility mini-plots, and a small right-pointing magenta arrow.
- [x] Bottom magenta outlined pill: `未达标：继续迭代`.

## Arrow inventory and implementation

- [x] 1 → 2: dark-blue clockwise arc across the upper right.
- [x] 2 → 3: amber downward arc along the right edge.
- [x] 3 → 4: purple long arc sweeping left along the lower right.
- [x] 4 → 5: green long arc sweeping left along the bottom.
- [x] 5 → 6: green upward arc along the left edge.
- [x] 6 → 1: magenta dashed arc across the upper left, labelled `下一轮设计`.
- [x] Implementation constraint: all six arrows are an embedded transparent SVG overlay made from explicit Bézier paths plus triangular SVG markers. No Draw.io connector, `endArrow`, or built-in connection arrow is used.

## Reconstruction medium decisions

- Native Draw.io: all editable text, card panels, numbered badges, central oval, pills, target icon, colored tool bars, heat scale, and sequence boxes.
- SVG: the six transition arrows and the two mini flexibility curves, because these require hand-drawn curve geometry rather than connector routing.
- Raster crop: the detailed AlphaFold-style protein ribbon and check badge, because this is complex artwork and exact visual fidelity is more important than redrawing it as a generic protein icon. The crop is kept local in the reconstruction asset directory.

## Visual QA

- [x] `.drawio` structural checker passes.
- [x] PNG export succeeds at 2× scale.
- [x] PDF and SVG exports succeed.
- [x] Full-size visual inspection confirms all text is legible and contained.
- [x] Full-size visual inspection confirms arrows do not cross through card content.
- [x] Full-size visual inspection confirms arrowheads and dashed iteration arc closely match the reference.
- [x] Full-size visual inspection confirms the protein crop has no conspicuous neighboring text or border contamination.

## Remaining risks

- No blocking visual issue remains. The protein ribbon/check artwork is intentionally retained as one embedded raster crop, so that small artwork region is not decomposed into independently editable vector residues; the cards, labels, badges, scales, sequence boxes, and all other diagram content remain editable.
