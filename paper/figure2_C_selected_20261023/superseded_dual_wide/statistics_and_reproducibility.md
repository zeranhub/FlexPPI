# 选定C双坐标叠加统计与可复现说明

沿用原五折体系留出预测，不训练、不重新推理、不校准；统计单位为记录。3,063条记录来自ATLAS 1,929、DynaRepo 890、ITrace 244。各记录内预测与参考的支持和阳性定义均由冻结协议确定，统计细节保留在`source_archive/PRESPECIFIED_PROTOCOL.md`及`PRESPECIFIED_ADDITIONS.md`，原统计说明也逐字节保留。

密度源为`prediction_reference_rmsf.npz`与`prediction_reference_chi_entropy.npz`，直接使用所有有效reference/prediction及原histogram_mass。各条记录总权重相同，残基在本记录内等权，二维格为质量而非原始计数。Full原趋势和区间来自`four_method_binned_trend.csv`；区间用原来源分层记录bootstrap的前1,000次抽样。图中MAE和ρ直接读`summary.json`，并与`per_record_metrics.csv`记录均值核对。

PR源为`macro_pr_curve.csv`、`within_type_macro_pr.csv`和对应`all_selected_endpoint.csv`/`within_type_all_selected_endpoint.csv`。每条曲线保留101个显示网格值及真实全选终点，不平滑、不抽样、不对显示网格积分估计AP。原10,000次来源分层bootstrap的种子为20261018；本次主图不额外增加AP区间、检验、星号或多重比较。方法差异的大小沿用原真实值。

b–d的分布采用原始每任务每方法的3,063个AP，无记录筛选。`all_record_AP_distribution_points.csv`保存全部36,756值，`AP_distribution_summary.csv`保存均值、中位数、25/75分位及极值。核密度使用同一0.20带宽，128个求值位置，轮廓不延伸到实际值以外；横向确定性错开不修改AP。分布四分位表示记录间差异，不是均值估计的不确定性。

`source_ap_matrix.csv`包含全部36个来源/任务/方法结果，用独立伴图保留。没有将来源分层评估称为来源留出。原v61全部方法分箱趋势与来源热图另在`original_v61_validation.pdf`完整保留，本图件a为同一绘图区的双横轴/双纵轴叠加，b–d为分布，e–g为PR；正式Word尚未替换。


a仅改变显示方式：两个物理坐标变换完全重合，蓝色下/左轴对应RMSF，暖赭上/右轴对应Hχ，原始范围、网格、分箱、记录权重和轴变换均未改变。原带权质量采用相同LogNorm范围，以恒定色相和0.015–0.35透明度映射叠加，零质量格完全透明；没有混合两个指标的统计量，也没有做跨指标归一化。四方法颜色键仅作用于b–g。
