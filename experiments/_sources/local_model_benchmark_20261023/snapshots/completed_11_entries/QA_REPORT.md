# 当前横评图 QA

状态：`PASS_CURRENT_SNAPSHOT_WITH_DOCUMENTED_VALIDATOR_LIMITATIONS`。本报告只对应下列已保存文件的哈希，不代表尚未运行完的模型或未来重绘图已经通过检查。

原 skill validator 已通过 import 显式执行，共 17 条原始 findings：{'PASS': 14, 'FAIL': 1, 'WARN': 2}。原文保留在 JSON 和 `figures/qa_raw_skill_findings.json`。补充文件/数据检查 31/31 通过。

AP-0 的配色 FAIL 来自精确空白匹配：代码写成 `DIVERGING   = [`，检查器只接受 `DIVERGING = [`；AST 已验证变量及颜色值完整存在。CL-2 无法解析 `183/25.4`，AST 给出两个画布的名义宽度均为 183 mm；同时独立测量实际导出的 PDF/PNG，不用名义值代替实际尺寸。AP-4 的 30 字符 legend 搜索与 AP-6 对注释中 bars 的词法搜索限制也保留并解释。

PDF 字号检查已修正旋转文字的解析误判：pdfplumber 的 size 对竖排坐标标签可能是字形宽度，7 pt 的 i 会被读成 1.554；本次直接读取 PDF 的 Tf 字号，并考虑文字与图形矩阵缩放，原始 bbox size 也留档。

绘图来源绑定：plot_audit 的 generated_utc 以及 plot_script_sha256 / summary_sha256 / scoring_audit_sha256 全部逐项校验；当前绑定状态为 True。任何旧字段缺失、哈希不符或UTC时间无效均要求重绘，不作为通过。

数据检查：主图 11 个完整方法/配置，三个互斥组为 2,193 / 626 / 244 条，总计 3,063 条；source.csv 共 93 行。核对了全部方法×组×适用指标、有效残基覆盖、summary 数值和完成状态；相对分数方法没有被写成绝对 Å MAE。未将 incomplete 方法放入主图。

| 文件 | PNG像素 / DPI | PDF实际宽×高 mm | PDF最小字号 | 字体嵌入 |
|---|---|---|---|---|
| benchmark_current | [2182, 1579] / (299.9994, 299.9994) | 184.824 × 134.028 | 5.40 pt | True |
| benchmark_distributions_current | [2154, 743] / (299.9994, 299.9994) | 182.356 × 62.875 | 6.00 pt | True |

实际 PDF 均检查了字体文件流、文本字号和字符边界；名义 183 mm 经 tight crop 后的实际差异已记录，未隐藏。

视觉检查与哈希：

- `benchmark_current.png`：`inspected_current_hash_with_notes`，SHA-256 `501fde5cf189547e565acbaa98bde08382d065f95fa9170888f56437df814572`。
  - Viewed the actual current 3-column by 3-row dot/interval figure. Eleven method/configuration rows are present in consistent order, including assembly and chains versions of both BackFlip-2.1 and BackFlip-1.
  - Blue/ochre/teal semantic colors are visible. Direct method labels and marker shapes supplement color.
  - No title, row-label, axis-label, or external-legend overlap/clipping was observed at the inspected image resolution.
  - ANM/GNM absolute-MAE cells visibly say relative score; TCR-pMHC points are shown without independent-family confidence intervals.
  - The longer BackFlip row labels fit without observed collision or clipping. Some intervals are smaller than their markers; this is not evidence that intervals are absent.
  - This is an eleven-completed-entry progress snapshot. Partial FastProtFlex and other unfinished models are absent from the primary figure, rather than being compared on smaller surviving subsets.

- `benchmark_distributions_current.png`：`inspected_current_hash_with_notes`，SHA-256 `c0e4de60093aee48d75e2155b0c068fdb119ae75d5877f559c364e3f4fbf26b3`。
  - Viewed the actual current three-panel empirical-CDF figure. Titles, axes, and the external bottom legend are visible without observed clipping.
  - The five displayed methods are FlexPPI, BackFlip-2.1 assembly, BackFlip-2.1 chains, DynaProt-M, and ANM. The five-entry, two-row external legend fits without observed clipping.
  - The learned-baseline curves share an ochre family color but have solid/dashed/dotted styles visible in the actual figure. Direct legend labels identify each; ANM is teal with a distinct dash pattern.
  - Single-chain curves partly overlap; no local zoom or nonlinear axis is used. No color-vision simulation or paper printout was inspected.

限制：未做色觉模拟或纸张打印检查；为控制内存没有读取大体积逐记录/逐残基评分文件。未来图像哈希变化时，本脚本自动把视觉复核改为 pending，不会沿用旧图的通过记录。
