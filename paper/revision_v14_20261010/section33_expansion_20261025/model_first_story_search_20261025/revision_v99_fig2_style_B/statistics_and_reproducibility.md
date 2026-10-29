# 数值与复现说明

本次为B版空间布局与图2样式统一，不新增模型训练、推理、MD、筛样或显著性检验。已冻结模型结果和参考的协议、fold、seed、checkpoint等字段仍见各原始源表及原文；本次未改动这些字段。统计结果与上一版逐字段完全相同。

b是同一CaMKI肽段的25个残基，不是25个独立实验；MD、Full、Intra各25值，共75值。两模型MAE分别为0.1148592012Å、2.383203344Å，按25位点绝对误差均值计算。源列为cam_all25.csv的CaMKI_native_position、target_rmsf、pred_rmsf、intra_rmsf。

c是泛素72个序列位置中的66个有效位点，MD和Full各66值，共132值。原始6个A/G无定义位点仍断线；三个Pro的合法零参考保留。Full MAE为0.07953050649754886nats。源表ubiquitin_72_position_profile_with_original_mask.csv及其valid_hchi掩码完整未动。

d全部6个Gln各有1个MD和1个Full Hχ，共12纵值，另使用原始6个sidechain_SASA_A2横坐标。无回归、无因果推断，不新增显著性检验。源表ubiquitin_all6_Gln_Hchi_SASA_plot.csv。空心点为MD、蓝色实心点为Full；只有同一残基的竖向配对线。

e按chi1–chi4分组的位点数为66/52/23/10，总151个点、302个坐标。双方先减去各自的氨基酸类型均值。横轴是同一联合Hχ预测，纵轴为相应MD单角fine-bin entropy；不是四个独立预测头。Spearman rho沿用per_chi_correlations.csv：0.4975052713、0.4471864059、0.5247035573、0.1272727273（显示两位小数）。原始全部点来自defined_angles_long.csv。未画回归线或新增P值。

f的样本单位为record，每组83值，共4组332值。record内部先按氨基酸类型校正，再计算Full预测与MD占据描述量的Spearman rho。记录级点的横向抖动只用于显示、不参与统计。均值及bootstrap95%区间沿用retained_panel_statistics.csv；未重新采样或改变区间。细箱熵Hχ均值0.5098929841，CI0.4907300073–0.5274475608；RMSF均值0.3205688592，CI0.3052365755–0.3359174296。低F*角宽Hχ均值0.3995555843，CI0.3818306871–0.4160917128；RMSF均值0.2385388632，CI0.2206801593–0.2554976387。原始记录值来自retained_absolute_correlations_all498.csv，其中Full两预测量的四组用于主图，Intra源值仍随原表保留。

主图冻结长表all_plotted_values.csv有853行，哈希233dcbf5e9222b248d46db9f5b377bde605c0905877209bfb8b375543b1a6b21，与v98字节相同。SASA横坐标和均值区间也完全保留。两PyMOL原图、PSE、PDB及相机清单字节不变。

样式取自当前论文图2figure2_dense_20261024：Arial，轴标/刻度/图例9pt，普通标题10pt，面板字母11pt，轴线0.6pt，Full#2166AC实线1.25pt，Intra#B87936虚线1.1pt。按图2相同Word宽度167.6mm显示时，普通文字约8.24pt。图3原生画布183×202mm，PNG450dpi，PDF嵌入Arial字体，SVG保留可编辑文本；分子为原始真实栅格渲染。Word未替换。
