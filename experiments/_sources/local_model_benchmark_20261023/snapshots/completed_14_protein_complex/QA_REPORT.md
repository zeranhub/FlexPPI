# 当前横评图 QA

状态：`PASS_CURRENT_SNAPSHOT_WITH_DOCUMENTED_VALIDATOR_LIMITATIONS`。本报告只对应下列已保存文件的哈希，不代表尚未运行完的模型或未来重绘图已经通过检查。

原 skill validator 已通过 import 显式执行，共 17 条原始 findings：{'PASS': 13, 'FAIL': 1, 'WARN': 3}。原文保留在 JSON 和 `figures/qa_raw_skill_findings.json`。补充文件/数据检查 37/37 通过。

AP-0 的配色 FAIL 来自精确空白匹配：代码写成 `DIVERGING   = [`，检查器只接受 `DIVERGING = [`；AST 已验证变量及颜色值完整存在。CL-2 无法解析 `183/25.4`，AST 给出两个画布的名义宽度均为 183 mm；同时独立测量实际导出的 PDF/PNG，不用名义值代替实际尺寸。AP-4 的 30 字符 legend 搜索与 AP-6 对注释中 bars 的词法搜索限制也保留并解释。

PDF 字号检查已修正旋转文字的解析误判：pdfplumber 的 size 对竖排坐标标签可能是字形宽度，7 pt 的 i 会被读成 1.554；本次直接读取 PDF 的 Tf 字号，并考虑文字与图形矩阵缩放，原始 bbox size 也留档。

绘图来源绑定：plot_audit 的 generated_utc 以及 plot_script_sha256 / summary_sha256 / scoring_audit_sha256 / training_exposure_audit_sha256 全部逐项校验；训练暴露审计SHA与 references/training_exposure_protocol_sources_sha256.json 实际文件逐字节哈希核对。当前绑定状态为 True。任何旧字段缺失、哈希不符或UTC时间无效均要求重绘，不作为通过。

数据检查：主图 14 个完整方法/配置，两组为 {'protein': 2203, 'protein_complex': 860}，总计 3,063 条；source.csv 共 78 行。逐条分类全部明确，核对全部方法×组×适用指标、有效残基覆盖和summary数值。相对分数未写成绝对 Å MAE。

| 文件 | PNG像素 / DPI | PDF实际宽×高 mm | PDF最小字号 | 字体嵌入 |
|---|---|---|---|---|
| benchmark_current | [2161, 1665] / (299.9994, 299.9994) | 183.000 × 141.000 | 6.30 pt | True |
| benchmark_distributions_current | [2161, 921] / (299.9994, 299.9994) | 183.000 × 78.000 | 6.30 pt | True |

实际 PDF 均检查了字体文件流、文本字号和字符边界；画布固定183 mm，并独立测量PDF与PNG实际尺寸。

视觉检查与哈希：

- `benchmark_current.png`：`inspected_current_hash_with_notes`，SHA-256 `01a4936a15611e3af38e603cc95cc1cf6ac46143638a663c03d71667b4edd39f`。
  - Viewed the actual 2161 x 1665 PNG at original resolution after the final render. Two biological rows are Proteins (2203 records, 1500 families) and Protein complexes (860 records, 261 families), with panels a-f and 14 fixed method/configuration rows in each group.
  - Blue, ochre and teal are visibly retained. The plot is not rendered in grayscale; direct labels, circle/diamond/triangle/square markers and fixed row positions supplement color.
  - The final MAE axes show 0.0, 0.5 and 1.0 ticks. The previous awkward endpoint label is removed. No clipped labels or overlap among method names, axis text, group headings or the footer/legend was observed.
  - Points and family-cluster intervals fit within each axis. Protein-group intervals are often smaller than their markers; this is not a claim that intervals are absent. Larger complex-group intervals remain visible.
  - All six nonphysical-MAE cells for FastProtFlex, ANM and GNM visibly say relative score. No physical error values were fabricated for them.
  - The footer visibly states 14 completed entries, same 3063 records and mixed training exposure. It is on separate lines above the bottom color legend. This does not claim equal training or completion of deferred heavy models.

- `benchmark_distributions_current.png`：`inspected_current_hash_with_notes`，SHA-256 `36709cab40a2eb34b13c1be7cb973937e94599167c0869542f55920458d2000f`。
  - Viewed the actual 2161 x 921 PNG at original resolution after the final render. Two panels use the same two biological labels and n=2203/860.
  - The six displayed configurations are FlexPPI, BackFlip-2.1 assembly, BackFlip-2.1 chains, DynaProt-M, FastProtFlex and ANM. All six names fit in the external two-row legend.
  - Blue and teal curves and four ochre curves retain named solid/dashed/dotted/dash-dot patterns. ANM has its separate dash-dot-dot pattern; this supplements color even where some protein-group curves overlap.
  - Titles, axis labels, tick labels, legend and training-exposure footer are visible without observed clipping or overlap. Negative-rank observations remain within the data-driven full x range; no local zoom or nonlinear axis is used.
  - The distribution shape is descriptive across all records; the visual does not imply that the 244 related TCR-pMHC records are independent sequence families.

限制：未做色觉模拟或纸张打印检查；为控制内存没有读取大体积逐记录/逐残基评分文件。未来图像哈希变化时，本脚本自动把视觉复核改为 pending，不会沿用旧图的通过记录。
