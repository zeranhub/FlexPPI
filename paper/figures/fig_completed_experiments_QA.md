# 已完成实验汇总图：统计与质量检查

## 图形契约

核心结论是：绝对状态预训练已经学到稳定的逐残基柔性排序信号，并在共同静态
RMSF 任务上稳定优于 FastProtFlex；现有配对响应仍然较弱，TCR3D 只提供中等
强度的外部 B-factor 代理证据。图采用 2×2 quantitative grid；（a）仅显示
RMSF 与侧链构象熵，按确认要求不显示 SASA。

主要审稿风险是：（a）—（c）使用内部验证或模型选择结果，不能称为完全独立
测试；（c）只有 7 个配对体系，区间较宽；（d）的 B-factor 受晶格、分辨率和
精修影响，不能等同于 MD-RMSF 真值。

## 数据与统计口径

- （a）来自 `state_validation_pretrain.json` 和
  `state_validation_finetune.json`。中心值为蛋白等权 Spearman，区间为体系
  bootstrap 95% CI；RMSF 的 n=276，侧链构象熵的 n=85。
- （b）来自 ATLAS 和 DynaRepo 的 FastProtFlex 比较报告。中心值为蛋白或
  复合物等权 Spearman，区间为体系 bootstrap 95% CI；n 分别为 191、31 和
  54。这里只比较两种方法都能回答的相对排序任务。
- （c）来自 `pair_validation.json`。中心值为 7 个体系等权的 MAE，区间为
  体系 bootstrap 95% CI。图中不报告 p 值，也不声称统计显著。
- （d）从 `tcr3d_residue_predictions.csv` 的全部 129,268 个残基重新计算
  每链 Spearman。没有抽样或降采样；没有残基因非有限值被删除。1,189 条已
  评测链中有 1,187 条可计算相关性，另外 2 条因序列过短或变量为常数而排除。
  重新计算的链等权均值 0.3886278959 与归档报告一致。橙色点为链等权均值，
  灰色线为归档的链 bootstrap 95% CI [0.3777164215, 0.3968168522]。

本图没有假设检验或多重比较校正。误差线只表示 95% 置信区间，不表示标准差或
标准误。用于绘图的 1,203 行数据已写入
`fig_completed_experiments_source_data.csv`，完整来源和信号检查写入
`fig_completed_experiments_statistics.json`。

## 可复现性与资产继承

生成脚本为 `make_completed_experiments_figure.py`。当前环境没有 R，而生产
资产库中的配对散点模板为 R 脚本；同时本图使用的是归档点估计与 bootstrap
区间，不是该模板需要的原始重复观测，因此（a）—（d）均采用参数继承：继承
生产资产的颜色、点形、误差线、坐标轴和 violin 视觉语言，用 Python 重建与
本数据结构相符的图层。PDF 是矢量主文件，PNG 为约 300 dpi 预览，绘图过程
未对数据做降采样。

## 四遍质量检查

1. 数据与图形：所有 JSON/CSV 均成功读取；脚本校验了 TCR3D 的链数和均值；
   所有面板均具有非退化的数值范围，误差线没有负长度。
2. 布局：四个面板使用同一网格、字体、配色与线宽；标题、图例、面板标签和
   数据标记无重叠；图例不遮挡数据；长数据集名称采用换行。
3. 可读性：使用蓝—橙—紫—灰且辅以圆、方、菱形三种点形；不依赖红绿差异；
   坐标轴标明 Spearman、MAE 和单位；置信区间含义在图注中说明。
4. 导出：`qa_validator.py` 退出码为 0；矢量 PDF 为单页，PNG 尺寸
   2374×1369 像素且元数据分辨率约为 300 dpi；PDF 与 PNG 已逐面板视觉核验。

最终判定：**Ready**。该判定表示图形和统计说明可进入当前稿件，不改变上述
证据等级限制。
