# FlexPPI 主图设计规范

## 1. 图类型

- 类型：Solution Overview（方法总览图）。
- 推荐范式：Multi-layer 多阶段范式，内部采用从左到右的数据流。FlexPPI 同时包含绝对状态预训练、apo/holo 配对微调和静态结构推理，三个阶段使用同一套状态编码器，但输入组合、监督信号和输出不同；纯线性流水线无法清楚表达这种共享关系。
- 不采用 System Architecture，因为模型不存在并发服务、控制器或反馈代理；不采用纯 Pipeline，因为那会让读者误以为 apo 与 holo 经过不同网络，或误以为 MD 轨迹是推理输入。

## 2. 画布与总体布局

- 最终尺寸：双栏通栏，178 mm × 94 mm；白色背景；导出 PDF 和 SVG 两份矢量文件。
- 四个面板从左到右排列，其中面板 (c) 和 (d) 在最右侧上下堆叠：
  - (a) Unified Dynamic States：x=0–42 mm，宽 42 mm。
  - (b) Shared Invariant State Encoder：x=46–112 mm，宽 66 mm。
  - (c) Two-stage Learning：x=116–178 mm，y=0–62 mm。
  - (d) Static-structure Inference：x=116–178 mm，y=66–94 mm。
- 面板之间保留 4 mm 空隙。主要数据流使用 1.1 pt 实线箭头，训练监督使用 0.9 pt 深灰虚线箭头，残基对应使用 0.9 pt 点划线。所有箭头均从左向右或从上向下，不允许交叉。

```text
┌──────────────────┐   ┌───────────────────────────────┐   ┌──────────────────────────┐
│ (a) Unified      │   │ (b) Shared state encoder θ   │   │ (c) Two-stage learning   │
│ dynamic states   │──▶│ Sequence ─┐                   │──▶│ Stage I: absolute state  │
│                  │   │           ├▶ Fixed-frame IPA  │   │ Stage II: paired response│
│ apo receptor     │   │ Geometry ─┘       │            │   │ + consistency             │
│ holo + partner   │   │               residue hᵢ      │   └──────────────────────────┘
│ MD labels ⋯⋯⋯⋯⋯⋯│   │               State Head      │   ┌──────────────────────────┐
│ residue map π    │   │               mean + variance │──▶│ (d) Inference            │
└──────────────────┘   └───────────────────────────────┘   │ apo, holo, induced Δ     │
                                                           └──────────────────────────┘
```

## 3. 面板 (a)：Unified Dynamic States

面板标题使用 “(a) Unified dynamic-state representation”。上半部分画一个蓝色目标受体，标签为 “Apo state: target receptor only”；下半部分画同一个蓝色受体与一个橙色伙伴接触，标签为 “Holo state: receptor + partner context”。不要用两个不同形状的受体，以免读者以为它们是不同蛋白。每个状态旁放四个小输入标签：“Sequence”，“Static N–CA–C/O coordinates”，“Target-chain mask”和“Condition (T, window, source)”。其中 N–CA–C 三个原子用三个小圆点和两条短线表示，明确局部坐标架来自观测骨架而不是模型生成坐标。

在两个状态下方画三条逐残基色带，分别标为 “RMSF (Å)”、“Side-chain entropy Sχ”和“SASA (Å²)”。左下角放一个灰色 MD trajectory 图标，向三条标签色带发出虚线箭头，并写 “training labels only”；MD 图标不得向编码器发箭头。apo 与 holo 之间画一组点划线残基对应，标为 “sequence-based residue map π: iapo ↔ jholo”。这同时解释了为什么变化量是在对应受体残基上计算，而不是对两个不同长度张量直接相减。

## 4. 面板 (b)：Shared Invariant State Encoder

面板标题使用 “(b) Shared invariant state encoder”。最外层使用浅灰边框，右上角标注 “same parameters θ for apo, holo, monomers, and complexes”。内部严格使用代码中的真实模块名，从左到右分成三段。

第一段上方是 “Chain-wise sequence encoder”，副标题为 “ESM or learned amino-acid embedding”；下方是 “Condition MLP”，输入为温度、时间窗口和数据来源。第二段是 “Invariant pair geometry”，内部只列四项：“Cα distance RBF”，“local direction”，“same-chain flag”和“sequence separation”。第三段是主模块 “Fixed-frame IPA blocks × L”，下方用一行小字写 “N–CA–C frames remain fixed; representations update, coordinates do not”。IPA 输出用一排蓝色小方块表示逐残基隐藏表示，并标为 “context-conditioned residue states hᵢˢ”。

从隐藏表示向上分出 “Shared absolute-state head”，输出两个并排的小框：“mean ŷᵢˢ = (RMSF, Sχ, SASA)”和“uncertainty σ̂²ᵢˢ”。不要画 AlphaFold recycling、torsion-angle head、structure generation 或坐标更新，因为这些模块不在 FlexPPI 实现中。

## 5. 面板 (c)：Two-stage Learning

面板标题使用 “(c) Absolute-state pretraining and paired-response fine-tuning”。用上下两个浅色区域区分训练阶段。

上半区标题为 “Stage I — absolute-state pretraining”。左侧画 “All labelled states: monomers, complexes, and unpaired systems”，右侧画绝对状态头的均值和方差；MD 标签通过虚线只进入损失框。损失框写 “Noise-aware absolute NLL  L_abs”，下方写 “model variance + replica/window label variance”。这里不要声称三个数据集都具有三个标签；可在标签旁以掩码符号写 “available labels only”，表示 ATLAS 当前主要提供 RMSF，而缺失的 Sχ/SASA 不作为零值处理。

下半区标题为 “Stage II — paired apo/holo fine-tuning”。画两个来自共享编码器的表示 hᵢᵃᵖᵒ 与 hπ(i)ʰᵒˡᵒ，一起进入 “Paired response head”。在响应头上方用五个小标签列出真实输入：“apo state”，“holo state”，“representation difference”，“element-wise product”和“absolute-state subtraction”。输出写为 “direct signed response Δŷᵢ and uncertainty”。右侧放三个损失框：“L_abs: retain absolute states”，“L_resp: supervise ΔRMSF, ΔSχ, ΔSASA”和“L_cons: direct response ≈ ŷholo − ŷapo”。底部只保留一条总目标：

\[
\mathcal L = \lambda_{\mathrm{abs}}\mathcal L_{\mathrm{abs}}
+ \lambda_{\mathrm{resp}}\mathcal L_{\mathrm{resp}}
+ \lambda_{\mathrm{cons}}\mathcal L_{\mathrm{cons}}.
\]

## 6. 面板 (d)：Static-structure Inference

标题使用 “(d) Static-structure inference”。左侧只放两个静态 PDB 图标：“apo receptor PDB”和“holo complex PDB”，二者均进入同一个共享编码器。右侧用三条水平逐残基色带输出：“pre-binding state ŷapo”，“post-binding state ŷholo”和“partner-induced response Δŷ”。在第三条色带旁同时标出两个很小的比较标签：“direct response head”和“state subtraction baseline”，说明模型保留了可检验的简单基线。整个面板不出现 MD trajectory，从视觉上明确推理不需要轨迹。

## 7. 颜色、字体与图形语法

- 使用 Okabe–Ito 色盲安全色板：目标受体蓝 `#0072B2`；伙伴橙 `#E69F00`；绝对状态/状态头绿 `#009E73`；直接响应朱红 `#D55E00`；IPA 紫 `#CC79A7`；数据与不确定性深灰 `#6B7280`；面板浅底色 `#F7F8FA`。
- 颜色必须与文字或线型双重编码：apo/holo 必须直接写标签；MD 监督必须使用虚线；残基映射必须使用点划线；共享权重必须写 “shared θ”，不能只依靠相同颜色暗示。
- 最终缩放后的面板编号 10 pt 粗体，面板标题 9 pt，模块名 8.5 pt，公式与辅助注释 8 pt；字体与论文正文一致。模块框使用 0.8 pt 深灰描边、2 mm 圆角，不使用阴影、渐变、立体效果或拟物数据库图标。
- 三条动态标签色带使用同一条以零为中心的蓝—灰—朱红发散色带来表示变化量；绝对值色带使用 Viridis 顺序色板。除变化色带外，不使用渐变填充。

## 8. 建议图注

英文图注：

> **FlexPPI uses one shared, rigid-motion-invariant state encoder to predict absolute residue dynamics in any structural context and learns binding-induced responses only from aligned apo–holo pairs.** Each state contains a target receptor, optional partner context, static backbone coordinates, and experimental conditions; MD trajectories provide masked RMSF, side-chain entropy, and SASA supervision during training but are not inference inputs. Absolute-state pretraining learns context-conditioned means and uncertainties, whereas paired fine-tuning adds a signed response head and constrains its prediction to agree with the difference between the two absolute-state predictions. At inference, two static structures yield the pre-binding state, post-binding state, and residue-wise partner-induced response.

中文图注：

> **FlexPPI 使用同一个刚体运动不变的状态编码器预测任意结构环境中的逐残基绝对动态，并仅通过残基对齐的 apo–holo 样本学习伙伴诱导响应。** 每个状态包含目标受体、可选伙伴环境、静态骨架坐标和实验条件；MD 轨迹只在训练期间提供带掩码的 RMSF、侧链熵和 SASA 监督，而不是推理输入。绝对状态预训练学习环境条件化的均值与不确定性，配对微调进一步增加带正负号的响应头，并约束其预测与两种绝对状态预测之差一致。推理时，两份静态结构即可输出结合前状态、结合后状态和逐残基伙伴诱导变化。

## 9. 工具与交付格式

- 首选 draw.io 制作第一版，因为该图需要精确对齐的模块框、共享箭头和多阶段分区；版式稳定后可用 Figma 做颜色、间距和蛋白图标的最后润色。
- 论文提交使用 PDF；同时保存 SVG 和可编辑 `.drawio` 源文件。不要把截图、PNG 或 AI 生成的蛋白示意图嵌入最终图中。
- 建议输出两个版本：正文完整版本 178 mm × 94 mm；答辩或幻灯片版本删除损失细节、保留四个面板主路径。

## 10. 质量控制与完整性门

- Vector format：设计通过，成图时必须确认 PDF/SVG 内没有嵌入位图。
- Font size：设计通过，成图后以论文最终 178 mm 宽度检查所有文字不低于 8 pt。
- Colour-blind safe：通过，使用 Okabe–Ito 并配合文字与线型双编码。
- Self-contained caption：通过，首句说明共享编码器与配对响应这一核心方法结论。
- Honest axis range：不适用，本图没有数值坐标轴。
- No chartjunk：通过，明确禁止阴影、3D、装饰渐变和无语义图标。
- Integrity gate：Gate 1–5 与 Gate 7 通过；Gate 6 不适用，因为本图不是 Motivated Example。成图后的字体、矢量属性和实际颜色仍需视觉复核。
- 严重度：0 CRITICAL，0 MAJOR，2 MINOR。两个成图后检查项是：确认最终缩放字体，以及确认导出的 PDF 没有把蛋白示意图栅格化。

