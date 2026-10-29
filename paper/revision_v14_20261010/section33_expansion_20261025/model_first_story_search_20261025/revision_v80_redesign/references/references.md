# 图3重绘的正式期刊参考与看图记录

检索日期：2026-10-25。已读取 agent-reach 的 SKILL.md、references/search.md 和 references/web.md；前序 Exa 免费额度已用尽，因此本次采用 web 搜索/官方期刊网页与官方媒体原图。以下四张图均由工具实际打开图片检查，不是根据搜索摘要推测。图片仅供内部版式研究，不会复制进本论文，不将它们当作本研究的实验依据。参考图文件与正文数据隔离保存。

## 首选：共用残基坐标的多轨道信息与结构对应

Arbesfeld JA et al. **Mapping MAVE data for use in human genomics applications.** Genome Biology. 2025;26:179. DOI: 10.1186/s13059-025-03647-x。

官方文章：https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03647-x

实际查看：**Fig. 2**；官方原图：https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-025-03647-x/MediaObjects/13059_2025_3647_Fig2_HTML.png 。已打开本地原图 GB_MAVE_Mapping_Fig2.png。图A将蛋白同一坐标下的变异、二级结构、可及面积、置信度、接触和功能注释纵向对齐；图B用完整连续的结构将上方信号落到三维位置。图中上下联线清楚地区分数据与结构含义。

可借鉴：CaM不必把每条曲线、每种接触类别都拆成一个独立稀疏panel；可将 ΔRMSF、ΔHχ、初始接触类别和EF-loop区域在成熟编号4–148共用横轴上组合，结构中的M72/M109/M144与曲线中的标记对应。借鉴“信息对应”的关系，不复制网页UI、彩虹氨基酸条、巨大的黑色折线箭头或这篇文章的科学结果。

## 配套：全体样本证据与实验/模型结构并排比较

Tekpinar M, David L, Henry T, Carbone A. **PRESCOTT: a population aware, epistatic, and structural model accurately predicts missense effects.** Genome Biology. 2025;26:113. DOI: 10.1186/s13059-025-03581-y。

官方图页：https://link.springer.com/article/10.1186/s13059-025-03581-y/figures/6

实际查看：**Fig. 6**；官方原图：https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-025-03581-y/MediaObjects/13059_2025_3581_Fig6_HTML.png 。已打开 GB_PRESCOTT_Fig6.png。左侧全数据性能点图与右侧同一视角的实验/模型着色结构构成两大panel，右侧三结构共用一条色标；图注明确PyMOL。

可借鉴：CaM中完整122个侧链的预测对应要先于局部例子，并明确全体点的数量与误差；若用两结构作视觉比较，应保留一致语义、近似视角、相同结构表现方式。不能把不同状态的真实变化绘成单纯好坏打分，也不应套用其致病性蓝黄颜色含义。图6左侧较长标签和大图例不适合直接照搬。

## 配套：密集定量面板旁安排足够大的完整结构

Ruperti F et al. **Cross-phyla protein annotation by structural prediction and alignment.** Genome Biology. 2023;24:113. DOI: 10.1186/s13059-023-02942-9。

官方文章：https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02942-9

实际查看：**Fig. 2**；官方原图：https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-023-02942-9/MediaObjects/13059_2023_2942_Fig2_HTML.png 。已打开 GB_MorF_Fig2.png。左侧是一张密集分组点图，右侧两个大尺寸结构叠合图；结构轮廓完整，没有为了局部强调而裁成碎片。蓝黄语义与左侧注释关联。

可借鉴：结构面板需要足够面积，不能被细小独立面板挤成装饰；CaM的整体状态形状与伙伴走向应完整可读，关键残基用sticks和简单指引突出。采用本论文的蓝、暖赭与青绿，而非机械照搬蓝黄。

## 补充：统一类别顺序、语义颜色贯穿不同统计读数

Kwon S et al. **Genomics 2 Proteins portal: a resource and discovery tool for linking genetic screening outputs to protein sequences and structures.** Nature Methods. 2024. DOI: 10.1038/s41592-024-02409-0。

官方文章：https://www.nature.com/articles/s41592-024-02409-0

实际查看：**Fig. 4**；官方原图：https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41592-024-02409-0/MediaObjects/41592_2024_2409_Fig4_HTML.png 。已打开 NM_G2P_Fig4.png。上方分组热图，下方分组柱图与小提琴图保持相同类别顺序和语义配色，数据类型改变但阅读顺序保持稳定。

可借鉴：CaM接触类别、全部9个Met在所有位置保持固定顺序；MD、Full和Intrachain分别用固定线型/符号与现稿语义色，避免同一颜色在结构中表示伙伴、在曲线中又表示模型而无解释。其小字密集布局不应照搬。

## 图3具体执行建议

新图建议由4个信息较完整的大panel组织：A为两个连续完整的CaM–伙伴结构与相应紧凑程度/关键位点标记；B为共用成熟编号横轴的主链变化、侧链变化、接触与EF-loop轨道；C为全部122个侧链的预测–MD变化对应，并将四种接触分组统计直接嵌在同一区域；D将关键位点的真实结构局部、MD/模型双状态数值、全部9个Met变化放在相邻区，避免单个局部结构只有一根棒图。具体空间由最终字号可读性决定，不能为压缩panel牺牲数据可读性。

主链RMSF和Hχ单位、范围不同，必须保持独立纵轴或独立热图色标；不能共享一个数值色标。图中保留主链linker方向预测错误和全部122/145个位点，不为视觉简化排除不吻合值。现稿Full蓝、Intrachain暖赭、辅助青绿优先于参考图原色。

## 新增原始结构论文书目信息

2YGG：Köster S, Pavkov-Keller T, Kühlbrandt W, Yildiz Ö. **Structure of human Na+/H+ exchanger NHE1 regulatory region in complex with calmodulin and Ca2+.** Journal of Biological Chemistry. 2011;286(47):40954–40961. DOI: 10.1074/jbc.M111.286906。官方JBC PDF：https://www.jbc.org/article/S0021925820504456/pdf 。PDB官方关联：https://pdbj.org/mine/summary/2ygg 。正式PDF首页核实作者、卷期、页码与DOI。

4LZX：Chen LT, Liang WX, Chen S, Li RK, Tan JL, Xu PF, Luo LF, Wang L, Yu SH, Meng G, Li KK, Liu TX, Chen Z, Chen SJ. **Functional and molecular features of the calmodulin-interacting protein IQCG required for haematopoiesis in zebrafish.** Nature Communications. 2014;5:3811. DOI: 10.1038/ncomms4811。官方RCSB：https://www.rcsb.org/structure/4LZX 。该页明确列出全部作者、论文题目、5:3811与DOI，并将4LZX命名为Ca2+-free CaM complex；同论文另一个结构4M1L为Ca2+-bound，二者不可混淆。

结构论文只是给输入结构的出处与背景。本研究轨迹和模型得到的定量结果仍以本地数据为准，不能用原始晶体结构的钙状态替代轨迹拓扑的离子审计。

## 检索失败记录与图像证据

PRESCOTT Fig.4、Fig.5和MaveDB2019 Fig.2的官方图页/图注可读取，但某些原图请求返回3038字节的HTML拦截页，未把这些当成“已实际看图”的证据，已从本文件夹清理。PDF截图对部分页面超时，未据此声称检查成功。以上四张被计入参考的PNG均已用view_image实际展示检查，且有真实PNG文件头、像素尺寸和SHA256记录见reference_image_audit.json。
