# 图2 外部模型高RMSF位点筛选曲线

建议正文使用 fig2e_external_PR.pdf（独立子图e），补充材料使用 figS_external_ROC.pdf。external_RMSF_PR_ROC.pdf 是两图并排查看版。三者均另有600 dpi PNG及保留可编辑文字的SVG，PDF嵌入TrueType字体。单图89 × 110 mm，并排图183 × 110 mm。原始Word稿件未在本次绘图中改写。

本图要回答的问题是：对同一批蛋白复合物，模型能否优先找出真实主链波动最高的残基。PR横轴Recall表示真正高波动位点中已找回的比例，纵轴Precision表示筛选结果中真正高波动位点的比例；ROC横轴为其余残基被误选的比例，纵轴与Recall相同。它们评价RMSF排序所支持的筛选用途，不是蛋白界面或是否发生相互作用的分类。

## 数据与统计口径

样本为图2d完全相同的615条DynaRepo蛋白复合物记录，覆盖261个序列家族和260,340个有效残基。1个、2至3个、至少4个界面的记录分别为414、144和57条。同一家族可以出现在多个界面层，因此不能将各层的家族数直接相加当作总数。七种方法使用完全相同的残基位置，未因预测缺失删除任何有效位置。

每条记录以MD参考RMSF的第90百分位及以上定义阳性，真实值并列全部保留；模型分数并列时作为一整个阈值组。PR在共同召回率上取该记录首次达到该召回率的阈值组精确率，再对615条记录等权平均，不对精确率作单调包络。Recall=0使用首个阈值组的真实精确率，Recall=1保留首次达到1的精确率，并额外画出连接到全部残基入选的末端。横向灰虚线为逐记录等权阳性率10.139574%。ROC保留每条记录完整的原始阈值点和垂直段顺序，在共同假阳性率上分段线性插值，再等权平均；另保留原点。两类曲线共同网格均为0至1的1001点，只用于显示，不改变标量指标。

图例中的AP来自每条记录原始阈值下的“召回率增量乘精确率”之和，AUROC也先逐记录计算，两者均再对记录等权平均。不能将显示曲线的视觉面积当作图例AP；标量AUROC也以独立的逐记录计算为准。重算AP与冻结图2d逐记录数据的最大绝对差为4.44e-16。FlexPPI沿用已有五折折外预测，外部方法沿用冻结的作者权重推理输出；本次没有训练或调整模型，没有分数校准，也没有跨记录直接合并残基。

summary_metrics.csv给出所有方法均值及95%区间。区间按261个家族整簇、有放回重抽2,000次，种子20261024；重抽样中的均值仍按记录数加权分母，同一组重抽权重用于七种方法。各指标的n=615条记录；残基数表示覆盖范围，不作为独立重复数。默认图不画置信带，避免七条曲线相互遮挡。mean_curves.csv额外保存FlexPPI逐点区间，可通过绘图参数开启；该区间不是整条曲线的同时置信带。本图未进行显著性检验，不报告p值或多重比较结论。

该结果只对应上述615条DynaRepo复合物，与论文全3,063条记录的总体AP不是同一口径；观察到的均值排名不等于所有蛋白体系、所有训练协议下的统一优势。

## 文件与复现

- fig2e_external_PR.pdf / .svg / .png：正文候选子图e。
- figS_external_ROC.pdf / .svg / .png：补充材料候选ROC图。
- external_RMSF_PR_ROC.pdf / .svg / .png：两图并排查看版。
- summary_metrics.csv：七种方法的平均AP、AUROC及家族重抽样95%区间。
- per_record_metrics.csv：615 × 7条记录级指标及阳性率。
- mean_curves.csv：七种方法PR和ROC的显示网格均值；共14,014行。
- cohort_records.csv：同图2d的615条成员清单。
- raw_threshold_curves.npz：每条记录、每种方法的全部原始阈值曲线，使用offsets恢复各段。
- compute_external_curves.py：从本地冻结预测与真值重新计算全部数据。
- plot_external_curves.py：只需本目录的结果CSV和审计文件即可重画图件。
- calculation_audit.json：输入文件的绝对路径、SHA256、阈值与统计规则。
- QA_report.md / plot_layout_audit.json：数值与版式检查记录。

本机运行环境为 E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe，已含numpy、pandas、scipy、scikit-learn、matplotlib。源数据无需复制到本图目录；计算脚本依赖calculation_audit.json记录的本地原始文件。绘图脚本及结果CSV可以整体搬到其他已安装numpy、pandas、matplotlib的环境中运行。

```powershell
& 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe' '.\compute_external_curves.py'
& 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe' '.\plot_external_curves.py'
```

## 建议图注

图2e：在所比较的方法中，FlexPPI在DynaRepo蛋白复合物的高RMSF位点筛选中取得最高的平均AP。七种方法在与图2d相同的615条记录、260,340个有效残基上评价；每条记录以MD参考RMSF第90百分位及以上定义阳性，并保留并列值。PR曲线在共同召回率处取各记录首次达到该召回率的阈值组精确率，再对记录等权平均；曲线末端保留全部入选状态，灰色横虚线为实际平均阳性率10.14%。图例为原始阈值曲线逐记录计算后平均的AP，不是显示曲线的积分；相同预测分数作为同一阈值组。BackFlip-2.1和BackFlip-1.0均采用assembly输入模式。

Figure 2e. FlexPPI achieves the highest mean average precision among the evaluated methods for high-RMSF residue retrieval in DynaRepo protein complexes. All seven methods are evaluated on the same 615 records and 260,340 valid residues used in Figure 2d. Positive residues have MD-reference RMSF at or above the record-wise 90th percentile, with ties retained. Precision at each common recall is taken from the first threshold group reaching that recall in each record and then averaged equally across records. The all-selected endpoint is retained, and the gray horizontal line indicates the mean positive fraction (10.14%). Legend values are record-mean AP calculated from the original threshold groups, not integrals of the displayed mean curves. Equal prediction scores enter together; both BackFlip versions use assembly inputs.

补充ROC图：七种方法在同一615条DynaRepo复合物上的高RMSF位点ROC曲线。阳性定义和有效位置同图2e；各记录保留完整阈值顺序，在共同假阳性率上插值后等权平均。图例AUROC先逐记录计算再等权平均；灰色对角线表示随机参考。全部方法的AP与AUROC家族整簇重抽样95%区间见配套源表。
