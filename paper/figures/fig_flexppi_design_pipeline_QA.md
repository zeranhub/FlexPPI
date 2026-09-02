# fig_flexppi_design_pipeline — 设计与 QA 报告

## Figure contract

- **核心结论：** FlexPPI 作为 dynamics-aware reranker 插在结构/静态筛选之后、MD 与湿实验之前，为从头设计与突变筛选两条入口的候选提供逐残基动态信号并参与多目标选择。
- **图形类型：** Solution Overview（应用/展望示意图）；除面板 (d) 的 Pareto 示意散点外无数据坐标轴。Pareto 散点为示意点，无真实坐标刻度与统计声明。
- **导出约定：** 178 mm × 95 mm 双栏；PDF 矢量主图 + SVG（可编辑文字）+ 300 dpi PNG 预览；拉丁文字 Arial/DejaVu Sans，`pdf.fonttype=42`、`svg.fonttype=none`。
- **与现有图的关系：** `fig_flexppi_architecture` 是模型方法图；本图是应用/展望图，两张图不合并。FlexPPI 内部网络（IPA、头、损失）不进入本图。
- **主要审稿风险及处理：** 图上明确写 `orthogonal dynamics signal · not ΔG`；实验反馈弧标注"not an online update of FlexPPI"；动态信息只出现在筛选端，不暗示参与生成；候选数不标注（honest funnel）。

## 生成方式

`python paper\figures\make_flexppi_design_pipeline_figure.py` → `fig_flexppi_design_pipeline.{pdf,svg,png}`。图形语法依据 `docs/pipeline_figure_literature_review.md`（双入口汇合 [RFdiffusion/BindCraft]、漏斗 [BindCraft/BRD3]、单一反馈弧 [EVOLVEpro]、实线/虚线语义 [DynamicMPNN]）。

## QA 结果

### Pass 0 — Anti-pattern scan

- [PASS] AP-0：字体、颜色、导出基线均已指定；无默认 tab10/jet/rainbow 色板。
- [PASS] AP-1/AP-2：使用 Okabe–Ito 语义色 + 灰阶；颜色与文字/线型双重编码（apo 蓝、binder 橙、通过绿、FlexPPI 紫、反馈朱红）。
- [PASS] AP-3：面板 (a–c) 无数据坐标轴；面板 (d) Pareto 散点轴为示意，无刻度标签，图注中已声明。
- [PASS] AP-4：无图例，不存在图例遮挡。
- [PASS] AP-5：同时生成 PDF、SVG 与 300 dpi PNG。
- [N/A] AP-6：不是小样本柱状图/箱线图。
- [PASS] AP-7：字体显式指定；PDF 无嵌入位图（`/Image` 仅出现在 ProcSet 声明），SVG 无 `<image>` 元素、58 个文字元素可编辑。

### Pass 1 — Code compliance

- [PASS] CL-1：最小文字 4.6 pt（Pareto 轴标签、反馈标注），插入 178 mm 双栏正文后仍可读。
- [PASS] CL-2：PDF MediaBox 504.57 × 269.29 pt ≈ 178 × 95 mm。
- [PASS] CL-3：PNG 2102 × 1122 像素，dpi ≈ 300。
- [PASS] CL-4：`pdf.fonttype=42`、`svg.fonttype=none`。
- [PASS] CL-5/CL-6：主数据流 1.1 pt 实线；评估/自洽 0.9–1.0 pt 虚线；反馈弧 1.0 pt 朱红虚线；箭头样式统一。
- [PASS] CL-7：PDF、SVG、PNG 均存在且非空。

### Pass 2 — Visual logic

- [PASS] VI-1：阅读顺序为“双入口 → 漏斗筛选 → FlexPPI 重排 → 选择/验证 → 反馈”，核心方法（紫色 FlexPPI 框）位于视觉焦点。
- [PASS] VI-2：颜色 + 文字 + 线型三重冗余编码；面板标题清晰。
- [PASS] VI-3：无装饰性图标、阴影或渐变；所有框和箭头均表示 pipeline 组件或数据流。
- [N/A] VI-4/VI-5：无真实数据轴和统计显著性声明；Pareto 示意散点已声明为示意。
- [PASS] VI-6：四面板标题位置、字号与边距一致。
- [PASS] VI-7：整图由单脚本重建，颜色、字体、箭头与框体语言统一。

### Pass 3 — Rendered output verification

- [PASS] VV-1：程序化渲染检查（逐文字 get_window_extent 与逐 patch 包围盒）：0 个文字越界、0 个非空文字相互重叠、0 个框体出画布；PNG 边缘扫描无内容被裁切。箭头/标签遮挡经几何推演核对（合流连接线、漏斗出口、反馈弧均无交叉）；最终成图仍建议人工肉眼复核一次。
- [PASS] VV-2：四面板边缘对齐、间距一致，全部内容位于画布内。
- [PASS] VV-3：178 mm 宽度下框内文字与色带标签可读；框内文字通过 TextPath 实测字宽自动换行，不会溢出框体。
- [PASS] VV-4：浅色填充 + 深色描边对比充分，Okabe–Ito 色对色盲安全。
- [PASS] VV-5：示意图不含模拟数据；Pareto 散点为示意图示。

**Verdict: READY。**

## 统计与可复现性说明

该图为应用/展望示意图，不包含样本点、均值、误差条或统计检验；面板 (d) 的 Pareto 散点仅为示意，不构成数据声明。重建命令：`python paper\figures\make_flexppi_design_pipeline_figure.py`。
