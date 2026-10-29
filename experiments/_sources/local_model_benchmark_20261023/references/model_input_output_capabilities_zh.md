# 横评模型的输入与输出能力

状态快照：2026-10-24 00:46（北京时间）。全量指同一3,063条记录；“完整结构输入”只说明本轮保留伴链几何，不能据此称为原生多链学习。显式链身份指进入预测计算的链关系特征，不包括仅用于选链、保存编号的标签。

| 模型 | 输入 | 本轮结构上下文 | 显式链身份 | RMSF单位 | 同定义联合侧链熵Hχ | 本轮状态 |
|---|---|---|---|---|---|---|
| Full / Intra [1] | 主链结构＋残基序列 | Full完整结构；Intra各链独立 | Full有同链关系；Intra按链拆分 | 绝对Å | **是** | 均全量；复用冻结五折OOF |
| BackFlip-2.1 [2] | 主链几何＋残基身份，无pLM | 完整结构、逐链各一套 | 本轮无；跨链连续位置编号 | 协方差nm²→RMSF Å | 否 | 两套均全量完成 |
| BackFlip-1 / 1-seq [3] | 1：主链几何；seq再加残基身份 | 完整结构、逐链各一套 | 本轮无；跨链连续位置编号 | 直接Å | 否 | 1两套全量；seq各1,318/3,063，运行中 |
| DynaProt-M [4] | 主链几何＋残基身份 | 4,973链独立预测后拼回记录 | 选链标签不作为跨链特征 | 直接Å | 否 | 全量完成，4,973/4,973链 |
| FastProtFlex [5] | 完整蛋白重原子坐标图 | 保留结构内全部蛋白链 | 无，图由原子距离构建 | 标准化相对分数 | 否 | 2,697/3,063，运行中 |
| ANM / GNM [6] | Cα坐标弹性网络 | 完整结构，跨链接触纳入 | 无，网络由坐标距离构建 | 相对幅度，非校准Å | 否 | 两者全量本机重现；GNM有简并模式数值差异 |
| PLANET-MD full / seq / mini [7] | full：序列ESM3＋结构；seq/mini：序列ESM3 | 已备逐链协议；本轮尚无前向实测 | 选单链；无跨链身份输入 | 原始nm，拟×10为Å | 否；SHP≠Hχ | 静态检查完成，推理待部署/延期 |
| RMSF-net（PDB-only）[8] | PDB原子生成模拟密度→3D CNN | 三条完整结构烟测，含五链1AO7 | 无，密度没有链标签通道 | 相对分数 | 否 | 3/3,063烟测成功；全量延期 |
| Flexpert-3D / Seq [9] | 3D：序列ProtT5＋ANM；Seq：序列ProtT5 | 官方入口指定单链；本轮未实测 | 选链标签不作为跨链特征 | 原始nm，需×10为Å | 否 | 权重/环境准备，未成功前向 |

证据（链接定位本地源码或本轮审计）：[1] [同链特征](E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/src/flexppi/model/geometry.py:88)、[逐链编码](E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/src/flexppi/model/information_ablation.py:178)、[双输出配置](E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/configs/full_fold0.yaml:31)。[2] [实际输入与单位](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/backflip/run_benchmark.py:145)、[运行审计](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/backflip/results/backflip-2.1/run_metadata.json:1)。[3] [官方结构入口](E:/Project/FlexPPI_baselines/BackFlip/backflip/deployment/utils.py:267)、[输入/单位/状态说明](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/backflip/README.md:1)。[4] [官方单链接口](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/vendor/dynaprot-main/dynaprot/inference.py:109)、[Å输出](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/vendor/dynaprot-main/dynaprot/inference.py:273)。[5] [原图模型和归一化](E:/Project/protein_flexibility_models/FastProtFlex/predict.r:24)。[6] [全量及数值差异审计](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/rmsfnet_nma/README_zh.md:5)。[7] [逐链协议/三版本检查](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/planet/README.md:1)、[SHP定义](E:/Project/protein_flexibility_models/PLANET-MD/README.md:153)。[8] [实际密度适配与输出](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/rmsfnet_nma/run_rmsfnet_local.py:55)、[三条烟测审计](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/rmsfnet_nma/rmsfnet/smoke/summary.json:1)。[9] [官方单链选择](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/vendor/Flexpert-main/predict.py:68)、[单位更正](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/vendor/Flexpert-main/README.md:78)。

仅完成全部共同记录的方法进入主榜；部分运行、烟测和未推理候选均单列。复合物标签按链独立去除刚体运动；不同模型的训练数据暴露和目标定义并未统一。
