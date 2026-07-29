# FlexPPI 模型示意图：设计与 QA 报告

## Figure contract

- **核心结论：** FlexPPI 用同一个固定坐标架、几何感知的状态编码器处理 apo 与 holo，再将可靠映射残基的两态表示组合为带正负号的逐残基响应，同时保留绝对状态相减这一透明基线。
- **证据链：** 面板 a 解释一个结构状态如何从残基类别、条件和主链坐标形成逐残基均值/方差；面板 b 解释共享权重的两态编码、残基映射、515 维响应特征、直接响应头和一致性约束。两个面板回答不同问题，无重复面板。
- **图形类型：** 双面板方法示意图；单一视觉系统，无定量数据面板。
- **导出约定：** 183 mm × 104 mm，RGB；PDF 矢量主图与 300 dpi PNG 预览；中文使用 Microsoft YaHei，拉丁文字优先 Arial/Helvetica，PDF 字体已嵌入。
- **主要审稿风险及处理：** 图中明确区分默认 learned 查表嵌入与可选分链 ESM2；明确残基编号不是显式 one-hot；明确坐标架固定且模型不生成坐标；响应均值允许正负；颜色不是唯一编码，模块名称与箭头方向同时提供语义。

## 资产扫描与生成策略

`assets/figures/` 中没有与模型架构示意图语义匹配的生产脚本，因此按 `cross-type inherit / param inherit` 生成。两个面板均由同一 Python 脚本、同一字体、颜色、圆角、箭头和线宽系统绘制，避免混用不同来源图标。源脚本首部保留了 Academic Figure Skill 的资产确认表与三组样式基线。

## QA 结果

### Pass 0 — Anti-pattern scan

- [PASS] AP-0：字体、颜色、导出基线均存在；自动 `qa_validator.py` 返回 0。
- [PASS] AP-1/AP-2：未使用默认 tab10、jet、rainbow 等色板。
- [PASS] AP-3：无数据坐标轴；示意框边界均为统一 0.7–0.8 pt 语义线条。
- [PASS] AP-4：无图例，不存在图例遮挡。
- [PASS] AP-5：同时生成 PDF 与 300 dpi PNG。
- [N/A] AP-6：不是小样本柱状图或箱线图。
- [PASS] AP-7：字体显式指定，中文字体与数学字体均成功嵌入。

### Pass 1 — Code compliance

- [PASS] CL-1：最小文字 5.4 pt；插入论文正文后仍不低于 5 pt。
- [PASS] CL-2：PDF 页面尺寸为 518.74 × 294.80 pt，即 183 × 104 mm。
- [PASS] CL-3：PNG 为 2161 × 1228 像素，DPI 元数据约 300。
- [PASS] CL-4：`pdf.fonttype=42`、`svg.fonttype=none`；`pdffonts` 显示字体全部嵌入。
- [PASS] CL-5/CL-6：示意图无数据轴；全局线宽和箭头样式统一。
- [PASS] CL-7：矢量 PDF 与位图 PNG 均存在且非空。

### Pass 2 — Visual logic

- [PASS] VI-1：阅读顺序为“单状态编码 → 共享双状态编码 → 响应”，核心方法在 3 秒内可辨认。
- [PASS] VI-2：蓝、绿、紫、红分别对应输入状态、几何/结合态、共享编码器和预测头；文字、位置与箭头提供冗余编码。
- [PASS] VI-3：无装饰性图标、阴影或渐变，所有框和箭头均表示模型组件或数据流。
- [N/A] VI-4/VI-5：无坐标轴和统计显著性声明。
- [PASS] VI-6：a/b 面板标签位置、字号与边距一致。
- [PASS] VI-7：按 schematic/model revision case 检查；整图由单脚本重建，颜色、字体、箭头和框体语言统一。

### Pass 3 — Rendered output verification

- [PASS] VV-1：PNG 与论文 PDF 页面人工检查，无箭头或标签遮挡关键公式。
- [PASS] VV-2：两面板边缘对齐、间距一致，所有内容位于画布内。
- [PASS] VV-3：在论文 `\linewidth` 尺寸下，框内文字、张量形状和图注可读。
- [PASS] VV-4：相邻模块颜色可区分，浅色填充与深色边框在白底上对比充分。
- [PASS] VV-5：示意图不含模拟数据；像素内容密度左面板 37.9%、右面板 42.0%、总体 39.9%，不存在空白面板或过度拥挤。

**Verdict: READY。**

## 统计与可复现性说明

该图是模型结构示意图，不包含样本点、均值、误差条、统计检验或训练/验证/测试结果，因此统计量与重复数均不适用。图中所有数字均可追溯到当前实现：21 类残基与 32 维类型嵌入来自 `model/network.py` 与 `constants.py`；默认 128/48 维、4 个 IPA 块、8 个头和 24 个 RBF 来自 `config.py` 与 `configs/*.yaml`；29 维原始几何特征来自 `model/geometry.py`；状态头和响应头维度来自 `model/heads.py`；响应 NLL 与一致性损失来自 `losses.py`。重绘命令为 `python paper/figures/make_flexppi_architecture_figure.py`。