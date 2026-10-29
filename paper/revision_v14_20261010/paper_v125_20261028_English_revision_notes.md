# FlexPPI 全文英文化说明

英文稿文件为 paper_v125_20261028_English_NC_GenomeBiology.docx。本次依据指定仓库 [nature-polishing](https://github.com/Yuan1z0825/nature-skills/blob/main/skills/nature-polishing/SKILL.md) 及其共享规范，完成主文、图注、补充材料、表格中文标签和原有批注的英文化。采用一致的英式拼写，将中文并列短句改写为主语和逻辑关系清晰的英文句子，并分别表述模型预测、MD参照、原子接触解释和文献中的SPR实验证据。原始中文文件未修改。

MEL5叙述保留整段CDR2β/CDR3β的连续响应、未发生氨基酸替换的残基上的主链变化，以及界面对侧HLA的侧链响应。“未替换”统一按 unmutated residues/positions 理解，避免 unchanged backbone 被误读为坐标不变。CDR3β的0.786相关明确指全部13个位点，0.850指其中9个未替换位点；这些数值来自原稿。补充方法进一步明确“成熟态减MEL5”的差值方向、两组氢键共同使用的几何阈值，以及按记录–氨基酸类型组判断有效位点数的计分门槛。

| 中文概念 | 统一英文表达 | 使用边界 |
|---|---|---|
| 主链均方根波动 | backbone root mean square fluctuation (RMSF) | 逐链对齐后的Cα位置波动，不包含链整体刚体运动 |
| 侧链联合状态熵 | joint side-chain state entropy (Hχ) | 固定χ角区间及观测窗口内的多角联合占据熵 |
| 完整／链内模型 | Full / Intra | 分别为全部有效链共同编码和逐链独立编码 |
| 静态特征／类型均值对照 | Static / AA | static-feature baseline / amino-acid-type mean baseline |
| 体系留出 | system-held-out evaluation | 原五折协议，不能等同于全部组成链家族隔离 |
| 全组成链家族留出 | family-held-out evaluation based on all constituent chains | 与原体系划分分开描述 |
| 同类型比较 | within-amino-acid-type comparison | 控制残基类型及侧链自由度差异 |
| 低相对自由能角区总宽度 | total width of low relative free-energy angular regions | MD角度占据的经验投影，不等同于结合自由能 |

主文摘要压缩为144词，集中呈现双读数、总体误差改善及MEL5的区域定位能力。原摘要中的CaM、泛素和83条记录的详细数值仍完整保留在结果、图注和补充材料中。本次保留原有章节顺序、所有结果和方法细节，没有把结果移入补充材料，也没有为压缩篇幅删除不利对照或限定条件。主要调整为句内和段内重写，核心论证顺序未变。第3.1节建立预测表现，第3.2节分解类型先验与几何环境的贡献，第3.3节解释读数的构象含义，第3.4节定位MEL5的识别环与跨界面响应。

按空格分词估算，引言660词、方法2,776词、结果正文2,823词、讨论491词、结论130词，四个主图图注合计1,001词。引言、结果正文、讨论和结论合计4,104词；该口径不包含摘要、方法、图注、参考文献和补充材料。中文原稿与英文词数不可直接比较，因此不把跨语言长度变化报告为压缩比例。本次沿用原稿的统计数值及其报告位置，没有重新计算实验成绩。

已核查213个翻译条目的覆盖性；正文及原有批注中无未处理的中文文本。6张表、8幅实际插入文档的图片和36个原生Word公式均保留；图片二进制内容与公式结构逐项检查一致。除摘要有意压缩外，所有段落中的小数数值经程序对比无变化；保留原参考文献内容和编号。使用本机Word导出36页并逐页检查页面图，未发现文字裁切、重叠、公式缺字或表格溢出。为使大幅组合图与完整图注保持在一起，部分页面保留自然留白。

## 期刊适配与原稿待补信息

本文件是完整英文语言稿，沿用原稿的章节顺序、参考文献样式和主文与补充材料合并形式。最终选定期刊后，应按相应栏目调整文件拆分、章节顺序和声明。Nature Communications 的[投稿说明](https://www.nature.com/ncomms/submit/how-to-submit)强调面向跨领域读者解释术语，并将补充材料单独提供；其最终Word格式采用双倍行距和左对齐。本次未把保留原排版的语言稿称作已完成全部投稿格式要求的文件，也未机械采用技能库内可能需要更新的期刊字数口径。

Genome Biology 的栏目要求不同：[Methodology](https://link.springer.com/journal/13059/submission-guidelines/methodology)当前要求不超过100词的非结构式摘要，而[Research](https://link.springer.com/journal/13059/submission-guidelines/research)要求不超过250词，并分为Background、Results、Conclusions三段。以下附一版适用于Methodology摘要长度的替换文本。两栏目的方法位置和Declarations要求也应在选定栏目后落实，不能仅根据英文语言风格视为已经适配。上述网页于2026年10月28日核查。

原稿代码与数据段落列出了复现材料，但没有实际的公开仓库地址、版本、DOI或许可证，因此英文稿保留事实描述，没有虚构发布链接。原稿也没有独立的局限性段落；已有的关键边界包括有限窗口与固定χ角区间、来源间观测设计差异、逐链对齐移除刚体运动，以及32条外部记录只属于一个关联家族。这些事实均已保留，但其是否需要在Discussion中再集中论述，应由作者在投稿定稿时判断。原稿6条批注已译为英文并保留；部分“待补”批注反映较早版本，可能与当前新增的单／联合目标或家族隔离结果存在时间差，需作者按当前研究状态确认。

## Genome Biology Methodology 短摘要备选

FlexPPI jointly predicts backbone fluctuations and joint side-chain state entropy from amino-acid identities and static backbone geometry. Across 3,063 records in five-fold system-held-out evaluation, mean absolute errors decreased by 14.8% and 11.3% relative to static-feature baselines. Within-type analyses supported geometry-dependent positional information. Predictions tracked side-chain angular occupancy diversity and, in family-held-out MEL5 complexes, recovered continuous recognition-loop constraints and responses at unmutated residues and opposing human leukocyte antigen side chains. Structures and replicate molecular dynamics simulations provided contact-level explanations. These complementary residue-level readouts support local dynamics analysis in single chains and complexes.
