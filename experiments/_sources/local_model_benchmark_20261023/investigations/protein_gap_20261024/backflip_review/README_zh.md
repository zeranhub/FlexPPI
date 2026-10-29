# BackFlip-2.1 静态输入与直接信息泄漏审查

本次全量只读审计为 **PASS：3,063/3,063 条记录通过，异常 0**。BackFlip 的源结构哈希逐记录同时匹配横评 cohort、BackFlip 输入清单、Full 的冻结预测导出清单以及训练冻结清单；从相同源结构生成的完整结构和逐链输入，在保留残基的序列、链和原始残基编号映射、N/CA/C/O 坐标及缺失原子掩码方面均一致。共检查 8,036 个规范化 PDB、7,680,628 个原子，B 因子非零原子为 0，occupancy 不等于 1 的原子为 0；共同的完整几何位点为 961,321，完整结构与逐链上下文合计核对 1,922,642 个坐标位置。源 PDB 中实际读到的 3,840,854 个主链原子，其 B 因子也全部为 0。因此，本轮未发现 BackFlip 因输入含有 MD 波动标签、B 因子或换用另一套静态结构而获得额外信息的证据。

这里需要区分“PDB 读取器读取过 B 因子”和“模型使用 B 因子”。官方读取器确实请求 `b_factor` 附加字段，但随后建立模型特征时仅从 N、CA、C 坐标计算每个残基的位置和朝向，并提取残基身份。本轮调用传入的位置、朝向及残基身份，与官方 `predict_from_frames` 内部建立的残基掩码和默认残基序号一起构成实际输入；B 因子、occupancy、参考 NPZ、RMSF/侧链熵标签、来源数据集编号均未传入网络。实际配置启用氨基酸身份、关闭 ESM 特征，故 BackFlip-2.1 应描述为“主链几何与序列残基身份输入”，不能写成仅几何输入。961,321 个保留位点全部为标准 20 种残基，官方未知残基转 ALA 的回退触发数为 0，与 Full 的残基身份分歧数为 0。

本地推理脚本在网络已输出协方差并转换为 RMSF 之后，才读取单独保存的 `reference.npz` 以写出对照值和评估掩码；参考文件没有进入前向计算。RMSF 换算采用协方差三个对角元素之和的平方根，再乘 10，将 nm 转为 Å，与官方实现一致。权重为发布的 `backflip-2.1.ckpt`，SHA256 为 `1af6e13a0e0161b2b92fa6fb9cddd2dee8f3e71c284f5059c77d7a30b9c527b9`；配置 SHA256 为 `acb6eb8683370fd6f2da3c2899a7025e524309e79b646b3bb94d79c1781e9ba8`。本轮实际使用的 6 个核心官方源码文件逐字节匹配已下载的官方源码归档，提交号为 `887c6a73b216e3fe2ff9f551ae6e2c983c21b292`，没有数值网络层替换的证据。

Full 的冻结结构解析器只读取坐标，不读取 B 因子；其实际 `forward_state` 读取残基身份、主链坐标和掩码、链编号、残基序号等特征，不读取标签。数据整理时标签可以与特征存在同一个 batch 中，但这不表示网络消费标签；实际冻结 Full 配置也关闭了 condition 特征。Full 导出代码和训练冻结清单的哈希均匹配原评估审计，5 个所检查的模型/解析源码均匹配冻结执行记录。本次证明的是源结构相同、BackFlip 保留位点与源结构坐标一致；两个模型仍有不同的特征构造、坐标变换及缺失几何处理，不能据此声称内部数值输入完全相同。

在多链处理上，完整结构上下文保留跨链空间几何，但本地调用没有显式传递链身份，官方接口默认使用拼接后的残基序号。逐链上下文则将每条链分别推理后按冻结映射拼回原记录。2,193 条单链的两种上下文共享同一次实际前向结果，状态日志的 6,126 个最终记录均为 complete。因此，“可以把复合物坐标交给网络”不等于“原生学习了链身份和复合物动力学”。同时，新 protein 组为 2,203 条，其中 2,193 条为单链，其余 10 条为已核定同源多聚体或同一蛋白多片段体系，不能仅用多链支持差异解释该组的整体性能差距。

以上审查排除了本轮可检验的直接输入通道问题，但没有证明发布权重不存在训练暴露。BackFlip-2.1 的精确训练成员清单仍不能与当前发布 checkpoint 完整绑定，应继续报告为未知；未知不能作为存在泄漏的证据。其发布配置体现 ATLAS/单体数据设置、协方差预测目标，与 Full 的多来源、RMSF 与侧链熵联合任务确有训练设置差异，但这些只能作为待验证解释，不能由本次源码审查直接推断为性能差距的原因。模型优劣及差距来源应由同记录配对结果、来源分层和已冻结训练暴露敏感性分析共同支持。

证据文件为本目录的 `static_input_audit.json`（全量结论与源文件哈希）、`static_input_per_record_audit.csv`（逐记录结果）、`source_feature_audit.json`（实际输入字段、冻结执行源码绑定和残基身份审查）、`audit_static_inputs.py`（本次全量检查脚本）。`evidence_sha256.json` 冻结这些新增证据及本说明的哈希。审查不执行模型、不读取本次性能数值、不访问网络，不修改原结构、标签、权重、冻结源码或现有结果。

主要源码定位如下，路径均为本机已存版本；`ROOT` 指 `E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023`，`CORE` 指 `E:/Project/FlexPPI/experiments/core_corrected_pipeline`，`BF_SOURCE` 指 `ROOT/backflip/upstream_archive/backflip-887c6a73b216e3fe2ff9f551ae6e2c983c21b292`：

- 规范化输入与 B 因子置零：`ROOT/backflip/prepare_inputs.py:39`、`:44`、`:59`、`:76`、`:80`；参考数组单独保存：`:86`。
- 实际前向与参考读取顺序：`ROOT/backflip/run_benchmark.py:140`、`:145`、`:154`、`:162`。
- B 因子被文件读取器加载：`BF_SOURCE/backflip/deployment/utils.py:200`；模型框架只由坐标和残基身份构建：`:260`、`:266`、`:273`、`:314`、`:320`。
- 官方推理特征白名单：`BF_SOURCE/backflip/deployment/inference_class.py:169`、`:189`、`:194`、`:198`；嵌入层实际消费字段：`BF_SOURCE/backflip/models/protein_embedder.py:75`、`:91`。
- 训练数据标签配置与推理输入开关：`ROOT/backflip/weights/backflip-2.1/config.yaml:12`、`:22`、`:122`、`:126`、`:140`。配置中的训练标签列不能解释为本轮推理通道。
- 官方 RMSF 单位换算：`BF_SOURCE/backflip/data/flexibility_utils.py:564`。
- Full 冻结解析及前向输入：`CORE/training_payload/src/flexppi/data/structures.py:61`、`:97`；`CORE/training_payload/src/flexppi/model/network.py:68`；`CORE/training_payload/configs/full_fold0.yaml:29`。
- Full 导出前向及随后取标签：`CORE/training_payload/evaluate_corrected_models.py:61`、`:64`、`:72`。
