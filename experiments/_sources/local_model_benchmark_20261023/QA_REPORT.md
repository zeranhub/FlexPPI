# 当前横评图 QA

状态：`AUTOMATED_CHECKS_PASS_VISUAL_REVIEW_PENDING`。本报告只对应下列已保存文件的哈希，不代表尚未运行完的模型或未来重绘图已经通过检查。

原 skill validator 已通过 import 显式执行，共 17 条原始 findings：{'PASS': 14, 'FAIL': 1, 'WARN': 2}。原文保留在 JSON 和 `figures/qa_raw_skill_findings.json`。补充文件/数据检查 33/33 通过。

AP-0 的配色 FAIL 来自精确空白匹配：代码写成 `DIVERGING   = [`，检查器只接受 `DIVERGING = [`；AST 已验证变量及颜色值完整存在。CL-2 无法解析 `183/25.4`，AST 给出两个画布的名义宽度均为 183 mm；同时独立测量实际导出的 PDF/PNG，不用名义值代替实际尺寸。AP-4 的 30 字符 legend 搜索与 AP-6 对注释中 bars 的词法搜索限制也保留并解释。

PDF 字号检查已修正旋转文字的解析误判：pdfplumber 的 size 对竖排坐标标签可能是字形宽度，7 pt 的 i 会被读成 1.554；本次直接读取 PDF 的 Tf 字号，并考虑文字与图形矩阵缩放，原始 bbox size 也留档。

绘图来源绑定：plot_audit 的 generated_utc 以及 plot_script_sha256 / summary_sha256 / scoring_audit_sha256 / training_exposure_audit_sha256 全部逐项校验；训练暴露审计SHA与 references/training_exposure_protocol_sources_sha256.json 实际文件逐字节哈希核对。当前绑定状态为 True。任何旧字段缺失、哈希不符或UTC时间无效均要求重绘，不作为通过。

数据检查：主图 14 个完整方法/配置，三个互斥组为 2,193 / 626 / 244 条，总计 3,063 条；source.csv 共 117 行。核对了全部方法×组×适用指标、有效残基覆盖、summary 数值和完成状态；相对分数方法没有被写成绝对 Å MAE。未将 incomplete 方法放入主图。

| 文件 | PNG像素 / DPI | PDF实际宽×高 mm | PDF最小字号 | 字体嵌入 |
|---|---|---|---|---|
| benchmark_current | [2182, 1867] / (299.9994, 299.9994) | 184.824 × 158.382 | 5.40 pt | True |
| benchmark_distributions_current | [2154, 743] / (299.9994, 299.9994) | 182.356 × 62.875 | 6.00 pt | True |

实际 PDF 均检查了字体文件流、文本字号和字符边界；名义 183 mm 经 tight crop 后的实际差异已记录，未隐藏。

视觉检查与哈希：

- `benchmark_current.png`：`pending_actual_visual_review`，SHA-256 `bb280e7981284ff0d03f1bcb728ffdc4766a439b506cf5c6fbd78be04b8a074e`。

- `benchmark_distributions_current.png`：`pending_actual_visual_review`，SHA-256 `4186b1c838a915c12516c6d0f287c69b4ddddbe01e534d8c943e90535d9810e0`。

限制：未做色觉模拟或纸张打印检查；为控制内存没有读取大体积逐记录/逐残基评分文件。未来图像哈希变化时，本脚本自动把视觉复核改为 pending，不会沿用旧图的通过记录。
