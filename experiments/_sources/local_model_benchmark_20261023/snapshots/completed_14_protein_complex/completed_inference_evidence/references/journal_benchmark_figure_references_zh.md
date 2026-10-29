# 真实横评图板与可复用构图（2026-10-23）

以下均已阅读原论文图注，并把原始PDF对应页渲染后逐张查看；没有仅凭论文标题推测图形。Nature/Science三篇满足主刊范围，RMSF-net作为同领域补充。图像为研究参考，不能把原图数值用于本轮模型成绩。

## 1. AlphaMissense — Science 2023，Fig. 2A–C

[原文](https://www.science.org/doi/10.1126/science.adg7492)；本地原始版面见 [Fig. 2 页](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/references/alphamissense_fig2-05.png)。PDF镜像是Science排印原文，原图文字与出版方图注交叉核对；PDF第5页对应印刷文章第4页。

真实图形是三个并排的**横向排序柱图**，不是点图。A为ClinVar总体分类，B为逐基因平均auROC，C为另一临床benchmark；各行保留方法名与95% bootstrap区间，少数模型跨分面使用固定强调色，灰色训练组有明确语义。它在一排中展示多个独立benchmark，读者无需辨认十余条几乎重叠曲线。

用于FlexPPI时，可把这种“方法为行、来源/任务为列”的结构改成横向**点加区间**以进一步节省墨色；原图并未使用这种改法，因此这是本次适配建议。建议各列固定相同模型行顺序，分别给MAE、Spearman ρ、AP的真实尺度，标n与计分方向；公开权重的训练暴露或单链输入限制应有文字标识。不要照搬原图截到0.5的柱基线，也不要把全部对照都变灰。

## 2. AlphaFold 3 — Nature 2024，Fig. 1c

[原文](https://www.nature.com/articles/s41586-024-07487-w)；[独立图页](https://www.nature.com/articles/s41586-024-07487-w/figures/1)；[本地原图页](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/references/af3_fig1.png)，PDF第2页、印刷页494。

真实图形是同一个c子图内横排四个任务块：配体、核酸、共价修饰、蛋白。每个任务块用少量分组柱、误差线、明确n和适用指标，AF3保持蓝色；任务之间由空间组织而非额外装饰区分。图注明确均值、置信区间、样本单位与采样协议。

适合把RMSF幅值、RMSF排序、高RMSF筛选、同类型Hχ筛选组织成清楚的任务分面。外部RMSF方法不支持Hχ时应在相应面板缺席或写不适用，不把别的输出强塞进去。模型多于约6个时，不宜复制密集分组柱，改用对齐的行式点区间图仍保留其任务分面的优点。各任务指标不同，不能把它们混成单一“总性能百分数”。

## 3. BioEmu — Science 2025，Fig. 2左列与Fig. 3A右部

[原文](https://www.science.org/doi/10.1126/science.adv9817)；[作者公开接收稿PDF](https://www.microsoft.com/en-us/research/wp-content/uploads/2024/12/accepted_paper.pdf)；[Fig. 2原图页](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/references/bioemu_fig2.png)与[Fig. 3原图页](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/references/bioemu_fig3.png)，接收稿PDF第16、17页。

Fig. 2左列按三种运动任务分行，用pretrained/finetuned两种模型的覆盖曲线对照；右侧才放代表性结构/自由能实例。它是**模型阶段消融**与任务分层，不能称为大量外部模型排行榜。Fig. 3A右上把计算成本与各蛋白的误差并排展现，成本采用对数纵轴，右侧“Errors & Scores”保留逐蛋白变化与总体均值。

适合借鉴“左侧紧凑总体统计、右侧少数有明确解释目的的实例”以及“性能与耗时分开、共享对象顺序”的组织方式。本文本机横评可加CPU秒/链或秒/有效残基的速度面板，与MAE/AP同列比较；速度必须来自本轮同机器实测，不能抄A100或文献GPU小时。无需复制BioEmu的大量蛋白结构图、彩虹热图或完整自由能版式。

## 4. RMSF-net — Nature Communications 2024，Fig. 2c–g

[原文](https://www.nature.com/articles/s41467-024-49858-x)；[独立图页](https://www.nature.com/articles/s41467-024-49858-x/figures/2)；[本地原图页](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/references/rmsfnet_fig2.png)，PDF第5页。

Fig. 2d、f按相关系数区间统计每个方法覆盖多少蛋白，展示平均成绩背后的分布；2c、e用**同一蛋白、两个模型的分数做配对散点**，以对角线明确显示哪一个更好。2g是按某模型分数先过滤样本后的结果，不能照搬该选择策略来证明本文优势。

适合为总体点区间主图增加一块“逐记录配对改善”图，例如选预先定义的最强对照绘制AP差或MAE差的分布，并保留零差线。模型太多时逐模型配对散点会重复占空间；应使用有共同记录集合的模型×来源分面，或只对预先指定的关键比较做配对分析。不可挑选Full胜出的样本才作图。

## 本次横评的构图落地建议

建议采用行式总体对照作为主板：行是模型，列是All、ATLAS、DynaRepo、ITrace及外部集；MAE、ρ和AP分成可直接比较的指标块，各自保留相同模型顺序和样本数。颜色以本文蓝色、暖色和青色基调为主，外部模型可按输入类型用同色系深浅并辅以形状，不能堆满不同彩虹色，也不能全部灰化。参考图中的区间形式要遵从本文真实估计方法，不能把记录间IQR改叫置信区间。

单个主图还可保留一块逐记录分布和一块关键PR曲线，但不必将每个新模型都叠到PR图里；完整曲线与来源分层表放补充，主图突出共同集上的幅值、排序、筛选三方面。先锁定成功推理的队列、每模型coverage、单位与训练暴露，再定最终模型行，避免画好布局后为凑行数混入不同协议。

## 访问与核验边界

使用 agent-reach 的 Exa 搜索和 Jina 网页路由、已安装GitHub connector检查官方接口。Jina在本次网络中超时，Exa后期出现429；随后用出版方/作者原始PDF读取与本地Poppler渲染完成视觉核验。PDF与PNG仅作为参考来源保存在本目录。Agent Reach版本检查返回GitHub API限流，未报告或声称已获得更新结果。
