# 用模型定义的高低组验证 MD 自由能占据曲线

先用冻结 Full Hχ 预测在每条记录、同一种氨基酸的位置内排序，复用 v26 完全相同的高、低半组及残基索引；之后才读取这些位置的真实 MD 角度。全部 83 条记录、全部 31 个自由能阈值均保留，没有根据 MD 结果另选位点或阈值。每组总计 16250 个位置，642 个奇数组的中间位置仅在高低分组中略去。相关分析仍使用 v26 的全部 33142 个同类型可比较位置。

横坐标 E 是相对自由能的截断值。纵坐标 W(E) 是满足 F*=−ln(p/pmax)≤E 的角度箱的总长度。p 是 MD 帧落在一个 5° 角度箱中的比例，pmax 是该角概率最大的箱；因此 E=1 表示保留概率至少达到最高箱概率 1/e 的角区。箱可以分散在几个位置，所以宽度可以跨越多个取向区域。每角先平均 0、1.25、2.5、3.75° 四个分箱原点，再平均该残基所有已定义侧链角，包含 Pro。没有平滑、插值或加虚拟计数。使用每条记录的 replica1，在 0–100 ns 范围取 501 帧，间隔 200 ps。

先在每条记录内按同类型氨基酸等权平均高低组宽度，再将 83 条记录等权平均，得到两条群体曲线。它们是模型预测分组对应的 MD 实测曲线，不是模型输出的自由能函数。模型验证来自预测高组在真实 MD 中覆盖更宽的低自由能角区；逐角曲线的形状由 MD 给出。

正文图使用 `population_curve_summary_frozen_primary_indices.csv`。它复用 v26 已冻结的主要高低组统计所保存的同一组 10000 次记录抽样索引，所有 31 个曲线点以及高、低两组共用这些索引。每次分别包含 ATLAS 30、DynaRepo 23、ITrace 30 条记录，阴影是每个阈值点处的 95% 百分位区间，不是整条曲线同时覆盖的区间。E=1 时低组 35.412503°、高组 43.170245°、配对差 7.757742°；差值区间 7.261125–8.262658° 与 v26 主结果和原有 10000 次抽样结果完全相同。

最初另行冻结的 seed20261014 方案及结果保留在 `analysis_plan.json`、`population_curve_summary.csv`、`primary_paired_bootstrap_curves.npz`，没有删除或覆盖。展示曲线复用原主证据链的索引是后续明确决定，记录在 `frozen_primary_indices_addendum.json` 中；它不是根据区间宽窄选择随机种子。E=1 差值沿用 v26“各同类型组高减低后再平均”的浮点数值，其与“两条平均曲线相减”只有约 1e-15 的累加顺序差；全部位置、分组和记录的实际量已先行核验一致。

`all_per_residue_widths.npz` 保存 33493 个原有共同位置在 31 个阈值下的宽度矩阵和身份；其中 `energy_cutoffs` 是 31 个 E 值，`width_degrees` 每一行对应一个位置、每列对应一个阈值，`eligible_sameAA` 标明是否属于原有至少 4 个同类型位置的可比较组。`prediction_defined_groups.csv` 保存原预测排序定义的高低组索引，`per_AA_group_curves.csv` 保存每个同类型组的真实高低曲线，`per_record_curves.csv` 保存 83×31 行真实记录曲线。`frozen_primary_paired_bootstrap_curves.npz` 保存展示用原索引、高低均值曲线及差值曲线的全部 10000 次结果。

另一份计算实现以独立的直方图累加和显式负对数公式复算全部 83 个缓存、33493×31 宽度矩阵、所有预测分组、氨基酸组与记录曲线，以及全部保存的曲线抽样，10452 项检查通过。核验报告为 `second_implementation_verification.json`，复算脚本为 `verify_population_width_curves.py`。展示用原主索引的 930000 个曲线抽样值另由 `verify_frozen_primary_indices.py` 全部重放，报告为 `frozen_primary_indices_verification.json`。这是同一实验代理使用另一份实现的复核，不宣称为另一个代理的审计。

复算原计划及展示索引的命令：

```powershell
& 'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -X utf8 'E:\Project\FlexPPI\paper\revision_v14_20261010\experiment_completion\free_energy_curves_v27_20261014\population_curves\run_population_width_curves.py' --phase compute
& 'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -X utf8 'E:\Project\FlexPPI\paper\revision_v14_20261010\experiment_completion\free_energy_curves_v27_20261014\population_curves\reuse_frozen_primary_bootstrap_indices.py'
```

全部原来源文件在运算前后核验 SHA256。没有新增 MD、训练、模型推理或 GPU 工作。
