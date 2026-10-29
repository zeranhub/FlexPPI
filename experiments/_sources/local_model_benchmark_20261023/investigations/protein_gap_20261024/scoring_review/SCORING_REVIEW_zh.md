# 蛋白组差距：评分、原始输出、单位和目标定义的独立复核

结论是：在本次已冻结的参考与模型输出上，FlexPPI Full 与 BackFlip-2.1 的差距确实存在，没有找到由 RMSF/熵列拿反、预测均值/方差拿反、残基顺序错位、双方有效残基不同、缓存取错、单位乘错或者把全局相关当逐蛋白相关造成的证据。该结论是输出文件和评分的核验，不是模型能力差距的因果解释。本审计没有重新推理、重新训练、修改标签、修改冻结文件或按误差挑选样本；读取了全部 3,063 个 Full 原始 NPZ、BackFlip-2.1 两种输入方式的全部 6,126 个原始 NPZ，并只读检查了五个 Full checkpoint 和一个 BackFlip checkpoint 的元数据。

所有 Full 原文件的 SHA 与冻结预测清单一致；原文件 `mean[:,0]` 与最终用于评分的 Full 数组逐位相同，`target[:,0]`、`valid[:,0]`、氨基酸、链名和原始残基编号均与共同参考逐位相同。BackFlip 的每个输出均有成功 status 和匹配的输出 SHA；其 `prediction`、`target`、`valid`、链名及残基编号也与共同参考和最终数组逐位相同。双方两个 BackFlip 输入方式均覆盖 3,063 条、961,321 个有效位置。本次共同参考 CSV 中 961,321 行的 `valid_rmsf` 全部为 1，原始输出的 RMSF valid 列也完全一致；没有为某个方法另删难预测位置。独立重算 Full、BackFlip-2.1 assembly 和 chains 的全部 9,189 条记录的 Spearman、Pearson、MAE，与发布表的最大绝对差为严格的 0。最差样本例子 1l2w_I、4ke2_C、1ng6_A、5dbl_A 也包含在这些全量原始文件校验中。

Full 的实际含义是五折 OOF（每条记录只取对应未参与训练折的模型）而不是目前某个部署 checkpoint，也不是五个模型预测平均。五个文件均为 `trained_checkpoints/full/fold_0...fold_4/epoch_040.pt`，实际 checkpoint 内均记录 epoch 40、seed 7、模型类型 full、输出指标顺序 `['rmsf','chi_entropy']`；五个 SHA 与原始 `prediction_manifest.jsonl` 的对应折完全相同，记录数依次为 612、612、611、614、614。导出脚本直接保存预测头的 mean，RMSF 为第 0 列，没有漏掉反归一化；训练保留源码的 StateHead 对均值使用 softplus，而 dataset 原样读入物理标签，没有把标准化分数当 Å 的步骤。这排除了拿错模型输出的解释，却同时限定了横评结论：这里比较的是固定 epoch 40 的原五折 Full 与发布模型，而不是穷尽优化后的 FlexPPI 最佳配置。

BackFlip-2.1 的实际权重是官方 ATLAS tag，SHA `1af6e13a0e0161b2b92fa6fb9cddd2dee8f3e71c284f5059c77d7a30b9c527b9`，与本机执行记录一致；不是 joint 或 mdCATH 版本。只读 checkpoint 记录 epoch **221**、global_step **81,918**。官方 config 允许 max_epochs 500、按 valid/loss 保存 best，但 500 是上限，不能写成这份权重实际训练了 500 个 epoch。Full 与发布权重的训练数据、更新次数和 checkpoint 选点规则不同，这是核实的实验条件差异；不能只凭 epoch 数多寡断言差距由训练不足造成。BackFlip 协方差输出按 `10*sqrt(trace(C))` 转为 Å，与官方 `batched_rmsf_from_covar` 完全相同；Full 的 mean 已是 Å。统一的正比例单位变换本来就不会改变同一记录内的 Spearman、Pearson 或 AP，因此即使怀疑幅度单位，也解释不了这里的排序差距。

本图采用的是每条记录先计算相关性、再把记录等权平均，不能和把所有蛋白残基串在一起计算一次相关性混用。独立复算结果如下，BackFlip 采用 chains 配置：

| 当前蛋白组来源 | 记录数 | Full 逐记录平均 ρ | BackFlip 逐记录平均 ρ | Full 串接残基后 ρ | BackFlip 串接残基后 ρ |
|---|---:|---:|---:|---:|---:|
| 全部蛋白 | 2,203 | 0.791201 | 0.897902 | 0.610959 | 0.905086 |
| ATLAS | 1,929 | 0.788367 | 0.905829 | 0.619623 | 0.913469 |
| DynaRepo 蛋白 | 274 | 0.811153 | 0.842089 | 0.569007 | 0.857999 |

串接残基相关性还混入了不同蛋白间柔性高低是否校准一致的信息，且长蛋白含更多残基；它不是本图的替代评分。这次改成串接残基也没有消除差距，反而使 Full 的差距更大，因此本次问题不来自“某一方按全局算、另一方按逐蛋白算”。蛋白组实际由 ATLAS 的 1,929 条单输入链、DynaRepo 的 264 条单输入链及 10 条多输入链构成；ATLAS 占蛋白记录的 87.6%。BackFlip 对 ATLAS 的逐记录优势较大，对 DynaRepo 蛋白的优势较小。这一来源分布本身不能证明训练泄漏；官方 BackFlip-2.1 的具体训练 CSV 未取得，已冻结训练暴露审计仍将其划分记为未知，不能把全组称为各模型相同且严格未见的测试集。

目标定义存在真实但尚不能定量归因的差别。FlexPPI 保留标签对每条链独立做刚体拟合，去掉链整体相互移动；BackFlip 官方预处理把输入的 CA 轨迹整体拟合到参考，对完整多链体系时两者定义不同。不过蛋白组只有 10 条多输入链，链级/整体拟合差别不能单独解释占绝大多数的 ATLAS 单输入链上的差距。另一个容易忽略的区别甚至可以影响单链：FlexPPI ATLAS 标签从三份预计算 `RMSF_R1/R2/R3` 先分别得到 RMSF，再逐残基平均；公开 BackFlip 数据生成脚本先合并多个 `R*.xtc`，再从协方差开方导出 RMSF。“先分别开方后平均”与“合并后求协方差再开方”通常不完全相同，合并也可能保留副本之间平均结构的差异。官方还提供相对首帧计算协方差的可选分支；公开脚本主入口默认不启用该选项，但发布权重对应的原始训练 NPZ/生成参数没有完整绑定，不能断言训练标签就是其中某一分支。以上说明是在共同参考上的迁移比较，不代表已找到对 Full 不利的评分错误，更不能据此否定 BackFlip 在相同冻结参考上的实测优势。

下一步值得验证的假设包括固定训练预算、来源采样、双任务目标以及 NLL/标签方差对大幅度残基的影响；这些属于优化机制假设，需要单独预先固定实验或进一步分析训练日志才能成立。本报告没有对它们作因果断言，也没有修改主横评分数。root 正另行完成预测幅度、端点/核心及最差蛋白剖析，避免将相同的输出诊断重复当作独立验证。

## 代码与数据定位

共同参考哈希和输入读取见 [evaluate_benchmark.py](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/evaluate_benchmark.py:48)，Full 源列和标签身份检查见同文件 [第 76 行](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/evaluate_benchmark.py:76)，双方共用 valid mask 和逐记录统计见 [第 226 行](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/evaluate_benchmark.py:226)，单记录统计定义见 [第 31 行](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/evaluate_benchmark.py:31)。冻结缓存来源见 [prepare_frozen_cache.py](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/prepare_frozen_cache.py:4)；本次直接复核了缓存和原始 Full 输出，未仅依赖缓存文字声明。

Full 固定 epoch/seed 及导出列见 [训练保留导出脚本](E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/evaluate_corrected_models.py:45)、[保存 mean/valid 的位置](E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/evaluate_corrected_models.py:62)；直接标签加载见 [dataset.py](E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/src/flexppi/data/dataset.py:49)，正值预测均值见 [heads.py](E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/src/flexppi/model/heads.py:20)。全部 checkpoint 元数据与哈希另见本目录 `full_checkpoint_identity.json`。

BackFlip 权重和 adapter 绑定见 [run_benchmark.py](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/backflip/run_benchmark.py:93)，官方 frames 接口和单位转换见 [第 140 行](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/backflip/run_benchmark.py:140)，输入逐残基坐标/顺序映射检查见 [prepare_inputs.py](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/backflip/prepare_inputs.py:56)，参考 valid 定义见 [第 88 行](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/backflip/prepare_inputs.py:88)。BackFlip 官方转换公式见 `backflip/upstream_archive/backflip-887c6a73b216e3fe2ff9f551ae6e2c983c21b292/backflip/data/flexibility_utils.py:564–576`，已实际加载的 epoch 元数据见本目录 `backflip_checkpoint_identity.json`。

独立链 Kabsch 及围绕平均坐标计算 RMSF 见 [prepare_corrected_labels.py](E:/Project/FlexPPI/experiments/core_corrected_pipeline/prepare_corrected_labels.py:106)；ATLAS 三副本先算后平均见 [保留 atlas.py](E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/src/flexppi/data/atlas.py:81)。BackFlip 官方公开数据生成分支见 `backflip/upstream_archive/backflip-887c6a73b216e3fe2ff9f551ae6e2c983c21b292/scripts/dataset/generate_dataset.py:96–104`，合并与整体拟合见同一官方源码下 `backflip/data/flexibility_utils.py:79–113`，平均协方差与首帧分支见 `126–143`。

本目录 `scoring_identity_audit.json` 绑定共同 truth、最终数组、原始清单、adapter、scorer 和本次验证脚本的 SHA；两个 raw identity CSV 逐条保存所验证原始 NPZ 的 SHA；`independently_recomputed_record_scores.csv` 是所有独立重算分数；`aggregation_definition_comparison.csv` 区分等权逐记录、按长度加权的逐记录和串接残基三种定义；`group_input_chain_composition.csv` 保留来源组成。所有这些都是派生诊断，主分组与主评分仍保持冻结版本。
