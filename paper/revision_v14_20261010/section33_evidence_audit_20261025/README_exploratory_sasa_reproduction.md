# 泛素全 66 位点起始结构 SASA 探索核验

输入结构：`E:/Project/FlexPPI/paper/figures_v59/03_KIX与泛素局部动态验证/source_archive/ubiquitin/inputs/A003R_input_structure.pdb`。Hχ 参考与留出预测：`E:/Project/FlexPPI/experiments/_case_validation_sources/data/ubiquitin_fold0_all_residues.csv`。必须按 `chain == A` 筛泛素：该源表实际同时包含 A/B 两条链，共 119 行；泛素有 72 个位置，其中 66 个 Hχ 有效。结构包含 1,921 个原子，去氢后完整双链复合物 935 个重原子均参加遮蔽计算。

复算命令（PowerShell）：

```powershell
& 'E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe' 'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_evidence_audit_20261025\reproduce_exploratory_sasa.py'
```

脚本在自身所在目录下的 `reproduced/` 写入复算 CSV 与 summary，不改输入文件或当前保存的原探索结果。实际初次计算 2.6 秒。环境：MDTraj 1.11.0、NumPy 1.26.4、SciPy 1.14.1。采用 Shrake–Rupley、1.4 Å 水探针、每原子 960 球面点；侧链排除 N/CA/C/O/OXT。MDTraj 的长度单位 nm，面积结果由 nm² 乘 100 转为 Å²。完整复合物参加 SASA 计算，但分析仅取泛素 A 链。

| 分析子集 | n | 侧链 SASA 与 MD Hχ 的 Spearman ρ | 模型预测与 MD Hχ 的 Spearman ρ |
|---|---:|---:|---:|
| 全部有效位置，原值 | 66 | 0.519082 | 0.979375 |
| 全部有效位置，各氨基酸内部中心化 | 66 | 0.401374 | 0.838417 |
| 类型至少 2 位点，各氨基酸内部中心化 | 63 | 0.404374 | 0.848159 |
| 类型至少 4 位点，各氨基酸内部中心化 | 50 | 0.488222 | 0.880432 |
| 全部 6 个谷氨酰胺 Q | 6 | 0.142857 | 0.942857 |

类型至少 4 位点的子集包括 D/E/I/K/L/Q/T/V，数量分别为 5/6/7/7/8/6/7/4。中心化是各列分别减去当前蛋白同氨基酸位点的均值；这里不是减训练集类型均值。原始源表中减训练集类型均值后的模型/MD 相关另行复现为 0.868192，与正文一致。

全部 6 个 Q 的起始侧链 SASA / MD Hχ / 预测 Hχ 分别为：Q2，70.106 Å² / 3.008879 / 2.989756；Q31，80.261 / 2.438781 / 2.461743；Q40，71.308 / 2.186479 / 2.554097；Q41，0 / 1.890488 / 2.000081；Q49，69.534 / 3.003505 / 2.945453；Q62，100.712 / 2.631066 / 2.683872。Q2/Q49 与 Q40 的起始暴露度相近而 Hχ 差约 0.82 nats，说明单一暴露度不是一一对应的解释；模型方向正确但对 Q40 偏高，压缩了这组差异。

这些是一个蛋白、一个起始坐标的描述性探索，没有统计显著性推断，也没有训练及独立检验 SASA 基线。不能把两列相关直接写成正式的模型对基线性能提升。SASA 是起始静态结构的暴露度，Hχ 来自 MD 侧链联合占据；二者不是同一个测量。正式论证应补跨蛋白、类型控制、按蛋白划分的 SASA/接触数基线，并使用与模型相同留出样本评价。
