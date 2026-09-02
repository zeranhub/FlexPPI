# FlexPPI：RMSF 与 Hχ 的介绍、文献和可复用图片指南

## 核心结论

FlexPPI 中的 Hx 应规范写为 \(H_\chi\)，读作 “H chi”，表示侧链 \(\chi\) 二面角联合旋转异构状态的 Shannon 熵；它不是 hydrogen exchange（氢交换，HX），也不等于完整热力学构象熵。RMSF 和 \(H_\chi\) 分别概括“主链位置晃动的幅度”和“侧链在不同取向之间分散的程度”，是互补的逐残基动态统计量。

## 两个指标应怎样解释

对每一帧先消除蛋白整体平移和旋转。设残基 \(i\) 的 Cα 坐标在第 \(t\) 帧为 \(\mathbf r_{it}\)，其时间平均位置为 \(\bar{\mathbf r}_i\)，则

\[
\mathrm{RMSF}_i=\sqrt{\frac{1}{T}\sum_{t=1}^{T}\|\mathbf r_{it}-\bar{\mathbf r}_i\|_2^2}.
\]

RMSF 的单位为 Å。高 RMSF 表示该残基的位置云更宽、围绕平均位置的摆动更大；低 RMSF 表示位置更集中。它不包含运动方向、残基间协同运动、跃迁速率或完整自由能景观。若未先做刚体对齐，蛋白整体漂移和转动会被误算成内部柔性。

对侧链，将每个可定义的 \(\chi_1\) 至 \(\chi_4\) 二面角划分为旋转异构区间，并把同一帧中的多个 \(\chi\) 角组合成联合状态 \(c\)。若状态 \(c\) 在轨迹中的频率为 \(p_{ic}\)，则

\[
H_{\chi,i}=-\sum_{c:p_{ic}>0}p_{ic}\log p_{ic}.
\]

FlexPPI 使用自然对数，所以单位是 nats。若侧链几乎总停留在一个状态，\(H_\chi\) 接近 0；若在多个状态间较均匀地切换，\(H_\chi\) 更高。最大可能值依赖残基具有多少个 \(\chi\) 角，因此不能把所有氨基酸的绝对值不加条件地视为同一量尺。Gly 和 Ala 没有可定义的侧链 \(\chi\) 角，应标为缺失而不是 0。

## 最值得使用的图片

| 优先级 | 文献与图 | 图能解释什么 | 许可与建议 |
|---|---|---|---|
| 1 | Vander Meersche et al., ATLAS, NAR 2024, Figure 3 | 同一页面同时呈现 MD、逐残基 RMSF profile、重复轨迹和按 RMSF 着色的三维结构；最适合讲“轨迹如何压缩成逐残基柔性图” | CC BY 4.0。建议裁取 3B 的 RMSF 曲线和 3E 的结构着色，图注写 Adapted from，并注明裁切/重排 |
| 1 | Park & Saitou, ROTAS, BMC Bioinformatics 2014, Figure 3 | 球棍模型标出侧链 \(\chi\) 角，并用 Newman 图显示三种典型 \(\chi_1\) rotamer；最适合解释 \(H_\chi\) 的离散状态来自哪里 | CC BY 4.0。可直接改编；再自行增加低熵/高熵概率柱图和 Shannon 公式 |
| 2 | Song et al., RMSF-net, Nature Communications 2024, Figure 1 或 Figure 4 | Figure 1 给出 RMSF 公式及 MD 标签到预测的流程；Figure 4 把逐残基曲线和结构着色并列 | CC BY 4.0。定义 RMSF 时优先 Fig. 1；展示结构上的柔性定位时优先 Fig. 4 |
| 2 | Clark et al., PLOS Computational Biology 2019, Figure 3 | 展示 apo、holo 以及合并状态的 \(\chi_1\) 角范围分布，说明结合可把侧链推入新的旋转异构势阱 | CC BY。它测量的是 \(\chi_1\) range，不是 \(H_\chi\)，只能作为生物学动机 |
| 3 | Cole & Warwicker, Protein Science 2002, Figure 1 和 Figure 6 | Figure 1 将蛋白表面按侧链构象熵着色并显示 PPI 界面；Figure 6 展示单残基对界面熵的贡献 | 与 FlexPPI 的 PPI 场景最贴近，但原文显示传统版权声明；论文公开发表前应申请许可，或仅借叙事后重绘 |
| 3 | Scouras & Daggett, Protein Science 2011, rotamer distributions | 从 807 个蛋白 MD 中展示侧链 rotamer occupancy 和 transitions，为“侧链不是一个固定取向”提供强证据 | 适合作为方法学引用；历史 Protein Science 版权页未明确给出 CC BY，不建议未经核实直接转载 |

## 推荐的最终定义图

不建议把两篇论文的图生硬拼在一起。最稳妥的办法是自行重画一张两栏图：左栏画多帧对齐后的半透明主链，在某一残基处画坐标点云、平均位置和距离箭头，下方放 RMSF profile；右栏画一个具有 \(\chi_1/\chi_2\) 的侧链、三个 rotamer 势阱以及两组概率柱图，其中 \(p=(1,0,0)\) 对应低 \(H_\chi\)，\(p=(1/3,1/3,1/3)\) 对应高 \(H_\chi\)。底部用一句话收束：RMSF measures positional fluctuation amplitude; \(H_\chi\) measures side-chain rotamer-state diversity。视觉语言可分别参考 ATLAS Figure 3 和 ROTAS Figure 3，但所有元素重新绘制并引用来源。

若需要第三栏，可画四种组合：低 RMSF/低 \(H_\chi\)、低 RMSF/高 \(H_\chi\)、高 RMSF/低 \(H_\chi\)、高 RMSF/高 \(H_\chi\)，直观说明“主链晃动”和“侧链换向”不是同一个量。Clark et al. 报告的主链位移与 \(\chi_1\) 范围相关性很弱，可作为这一区分的文献支撑。

## 建议的中文介绍

蛋白质并不是一座静止的雕塑，而是在热涨落下不断访问一组相近构象。完整分子动力学轨迹包含大量坐标，因此需要把它压缩成可以逐残基比较的统计量。RMSF 衡量一个残基在对齐后的轨迹中围绕平均位置移动了多远，单位为 Å；数值越高，说明该位置的空间波动幅度越大。\(H_\chi\) 则关注侧链而不是主链位置：它把侧链 \(\chi\) 二面角划分为若干旋转异构状态，再用这些状态在轨迹中的占据概率计算 Shannon 熵；数值越高，说明侧链在更多取向之间分散。两者分别回答“这个残基的位置晃动多大”和“它的侧链取向有多丰富”，因此应共同使用，而不能相互替代。

## 核心文献阅读路径

1. McCammon, Gelin & Karplus, Dynamics of folded proteins, Nature, 1977：蛋白质全原子 MD 的奠基工作。
2. Karplus & McCammon, Molecular dynamics simulations of biomolecules, Nature Structural Biology, 2002：MD 与生物分子动态的经典综述。
3. Vander Meersche et al., ATLAS, Nucleic Acids Research, 2024，DOI 10.1093/nar/gkad1084：标准化重复 MD 和逐残基 RMSF 数据。
4. Song et al., RMSF-net, Nature Communications, 2024，DOI 10.1038/s41467-024-49858-x：直接从结构/密度预测 RMSF，并在正文明确给出 RMSF 公式。
5. Scouras & Daggett, The Dynameomics rotamer library, Protein Science, 2011，DOI 10.1002/pro.565：MD 中的侧链 rotamer 分布、替代构象和切换。
6. Park & Saitou, ROTAS, BMC Bioinformatics, 2014，DOI 10.1186/1471-2105-15-307：用图清楚定义侧链 \(\chi\) 角和 rotamer 状态。
7. Cole & Warwicker, Side-chain conformational entropy at protein–protein interfaces, Protein Science, 2002，DOI 10.1110/ps.0222702：侧链熵与 PPI 界面的直接先例。
8. Clark et al., Inherent versus induced protein flexibility, PLOS Computational Biology, 2019，DOI 10.1371/journal.pcbi.1006705：说明主链位移与侧链 rotamer 变化需要分开描述。
9. Frederick et al., Conformational entropy in molecular recognition by proteins, Nature, 2007；Tzeng & Kalodimos, Protein activity regulation by conformational entropy, Nature, 2012：用于说明构象熵与识别/调控的重要性，但不能把其 NMR 构象熵与 FlexPPI 的 MD \(H_\chi\) 直接画等号。

## 署名模板

原样转载可写：Reproduced from [authors], [journal, year], under CC BY 4.0。

裁切、改色、加标注或重排可写：Adapted from [authors], [journal, year], under CC BY 4.0; cropped, recolored, and annotations added。

CC BY 图片仍必须给出作者、来源、许可链接并说明改动。对没有明确开放许可的图，不因“能在 PMC 免费看到”就自动拥有转载权；此时应申请许可或重绘概念。

## 本地图片包

上述候选图已下载到 `docs/literature_figures/rmsf_hchi_intro/`，文件格式、像素尺寸和图片内容均已检查。文件清单、官方原图链接及逐图许可说明见 [本地图片索引](literature_figures/rmsf_hchi_intro/README.md)。
