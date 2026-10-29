# FlexPPI v2 五折横评：本地基线盘点

审计日期：2026-09-06。此次只读检查了本地代码、说明、权重文件、旧预测表，以及轻量运行环境导入；没有下载、安装、训练、运行基线全量推理或启动 Docker 服务。本文件是运行准备清单，不代表以下所有模型都已完成新版横评。

## 路径纠正与当前计算环境

用户给出的 `E:\Project\FlexPPI\_baselines` 与 `E:\Project\protein\_flexibility\_models` 在当前工作站不存在。实际相邻目录是 **`E:\Project\FlexPPI\baselines`** 和 **`E:\Project\protein_flexibility_models`**。前者主要保存 RMSF-net、AlphaFlow 的部署说明和 Dockerfile，并不是所有模型代码、权重的集合；后者包含 FastProtFlex、RMSF-net、PLANET-MD、DynaMine、OPUS-Bfactor、DEFMap 等资源。项目自己的第三方源码另在 `E:\Project\FlexPPI\tmp\third_party`。

`E:\Project\FlexPPI\.venv-md\Scripts\python.exe` 实测 PyTorch `2.8.0+cpu`、`torch.cuda.is_available() == False`，可导入 NumPy、SciPy、pandas、MDAnalysis；未发现 `prody`、`backflip`、TensorFlow。系统可见 AMD Radeon HD 5450、Intel UHD Graphics 770 与虚拟显示适配器，没有可见 NVIDIA GPU。Docker 客户端 `29.4.3` 存在，但连接 `dockerDesktopLinuxEngine` 管道失败，服务当前未运行，因此本次无法重新核实 Docker 镜像、容器及容器内权重是否仍存在。不能把旧部署说明中的镜像名称当成当前已验证可调用状态。

## 主榜 RMSF 基线

| 模型 | 真实本地资源与运行入口 | 当前可运行性及证据 | 横评限制 |
|---|---|---|---|
| FastProtFlex | 源码 `E:\Project\protein_flexibility_models\FastProtFlex`；现有批处理 `E:\Project\FlexPPI\scripts\fastprotflex_batch.R`；协调脚本 `E:\Project\FlexPPI\scripts\compare_fastprotflex.py` | `C:\Program Files\R\R-4.6.0\bin\Rscript.exe` 实测可运行，`bio3d`、`orca` 的 `requireNamespace()` 均为 TRUE。现有 CPU 入口环境可用，但此次未执行新预测。其线性系数直接保存在批处理脚本中，不需要神经网络权重文件。 | 输出相对柔性分数，不列绝对 Å MAE。适配器以 ORCA 的等价轨道计数替代上游 netdist 入口；单次整体计算对超过 10,000 个蛋白重原子的结构明确报错，不会自动换成分块近似。新全量队列需先审计这一限制。 |
| BackFlip-seq / BackFlip-structure | 适配器 `E:\Project\FlexPPI\scripts\run_backflip_baseline_batch.py`；旧记录标签 `backflip-1.0-seq`、`backflip-1.0`；旧预测缓存见下节 | 原流程调用容器中的 `backflip.deployment.inference_class.BackFlip.from_tag()`。当前主环境无 backflip，限定盘点目录未找到独立主机源码/权重；Docker 服务未运行，容器内实际权重路径未获复核。旧逐残基缓存可重计指标；新增样本需恢复并检查原容器或另行配置环境。 | 适配器读取 `global_rmsf`；上游训练标签的对齐定义不能仅凭名称认定与当前标签相同。多链输入必须按既定方案逐链预测，并注明不接收伙伴上下文。旧静态复合物结果必须采用 `*_static_complex_chainwise.csv` 替换原整体 PDB 预测，不能直接拼接导致重复。 |
| ProDy-ANM / GNM | 适配器 `E:\Project\FlexPPI\scripts\run_prody_baseline_batch.py` | 旧协议记录 ProDy `2.6.1`，Cα 节点、15 Å 截断、gamma 1、20 个最低非平凡模态；当前主环境没有 ProDy，旧缓存可用，新增样本需恢复环境。ProDy 是物理模型，不需要学习权重。 | 旧缓存为 `calcSqFlucts()` 返回的平方涨落。沿用旧榜须保持变换口径；若改成平方根涨落，需对全部方法统一说明并重新统计，因为 Pearson、标准化误差会改变。只比较相对指标，不给未经校准的绝对 Å MAE。 |
| RMSF-net-PDB | 源码 `E:\Project\FlexPPI\tmp\third_party\RMSF-net` 及 `E:\Project\protein_flexibility_models\RMSF-net`；PDB-only 权重各自 `model_cpu\sim_model.pth`，项目副本大小 26,857,801 字节；适配器 `E:\Project\FlexPPI\scripts\run_rmsf_net_baseline_batch.py`；部署 `E:\Project\FlexPPI\baselines\rmsf_net\Dockerfile` | 旧固定版本 `eb353b0d9d909a6c46d2808467a9d0f16e2bbda7`，记录镜像 `flexppi-baseline/rmsf-net:pdb-eb353b0`。源码、权重、CPU 适配器和旧缓存均在，但 Docker 当前不可连接，尚未验证原容器或直接 CPU 加载整条路径。 | 使用 PDB-only CPU 分支，不引入实验 cryo-EM 密度。适配器用 NumPy 重现文档 molmap 参数，属于适配器实现。输出按相对分数评价；不把 PDB B-factor 字段或网络原始分数误作已经校准的 RMSF Å。 |
| PLANET-MD / RocketSHP（新增候选） | 源码 `E:\Project\protein_flexibility_models\PLANET-MD`；环境 `E:\Project\protein_flexibility_models\envs\planet_md_py312\Scripts\python.exe`；入口 `E:\Project\protein_flexibility_models\PLANET-MD\planet_md\interface\predict.py` | 实测 `import planet_md` 成功，定位到上述源码；PyTorch `2.12.0+cpu`，CUDA False。有结构+序列版与序列版检查点，以及 ESM3 缓存，但本次**没有加载权重或完成单样本推理**，因此仍是候选而非已完成基线。 | 官方 CLI 默认 `cuda:0`，序列提取只取第一条链。必须显式 CPU、离线本地缓存、逐链/残基编号映射，并分别调用正确的序列与结构特征，再验证输出单位与训练集重合，才能进主榜。CPU ESM3 成本需通过小样本测量，不能保证全量耗时。 |

PLANET-MD 的两套已发现检查点位于 `E:\Project\protein_flexibility_models\hf_cache\hub\models--samsl--rocketshp\snapshots\7831d318182d4f4fc1acbdafdbf23862d8bc1657\checkpoints`，分别为 `rocketshp_pretrained_20250505_v1.ckpt`（356,157,598 字节）和 `rocketshp_pretrained_20250507_v1_seq.ckpt`（351,949,258 字节）；另一快照 `0c8357ea44f8f6c511dc967d0a41a9f1bd275029` 中也发现同名文件。ESM3 缓存位于 `E:\Project\protein_flexibility_models\hf_cache\hub\models--EvolutionaryScale--esm3-sm-open-v1\snapshots\47f0545b2b6daf26a93439a3cd610f4f7f3d5478`，包含主模型与结构编码器/解码器权重。本次没有做各快照内容一致性校验，也没有发现 mini 版已下载检查点。

## 旧缓存覆盖：520 个系统，不是新版 3,068 个系统全覆盖

可复用原始基线缓存统一位于 `E:\Project\FlexPPI\artifacts\evaluation\additional_baselines\raw`。冻结旧队列为 ATLAS 单链 191、DynaRepo 单链 31、多链复合物 54，即静态 276 个系统；另有 ITrace 244 个复合物，总计 **520 个系统**。新模型五折覆盖比旧队列更大，因此现有旧缓存只支持“旧队列复测子榜”，不能声称所有基线已覆盖全部新版折外样本。

| 原始缓存 | 实际读取行数 | 唯一系统数 |
|---|---:|---:|
| `backflip_seq_static.csv` / `backflip_structure_static.csv`（各自） | 67,077 | 276 |
| `backflip_seq_static_complex_chainwise.csv` / `backflip_structure_static_complex_chainwise.csv`（各自；用于替换对应复合物） | 19,601 | 54 |
| `prody_static.csv` | 67,064 | 276 |
| `rmsf_net_static.csv` | 67,077 | 276 |
| `backflip_seq_itrace.csv` / `backflip_structure_itrace.csv` / `prody_itrace.csv` / `rmsf_net_itrace.csv`（各自） | 199,668 | 244 |
| `alphaflow_atlas_official_overlap.csv` | 1,486 | 6 |

ProDy 静态原始行数比另几个基线少 13 行，不能以原始 CSV 行数相同为匹配条件，也不能静默取交集。旧统计脚本 `E:\Project\FlexPPI\scripts\evaluate_additional_rmsf_baselines.py` 按 `sample_id + chain + author_residue_id` 唯一键先合并基线，再对冻结真值表左连接，并在缺失冻结残基时硬报错。新版应保留并增强这种覆盖审计，使用同一有效残基集合；“有某个系统的部分预测”不等于该系统已完成。旧指标表位于 `E:\Project\FlexPPI\artifacts\evaluation\additional_baselines\results`，新统计应写入独立目录，避免覆盖原图表证据。

FastProtFlex 旧原始预测与共同残基数据另在 `E:\Project\FlexPPI\artifacts\evaluation\fastprotflex_comparison` 及 `E:\Project\FlexPPI\output\itrace_benchmark\raw`。评估器默认仍读取旧的 `flexppi_pretrain_rmsf_angstrom` 列，并且输出路径固定到旧目录，所以**不能不修改入口就直接重跑这个旧脚本生成所谓 v2 结果**。

五折模型横评应给每个系统只选其未参与训练的那一折模型预测，并记录 fold ID；不能把五个模型的预测平均后在它们训练过的样本上报测试成绩。已训练外部基线的同输入预测可按该系统的 fold ID 分组重计，无需为了“五折”把确定性预训练推理重复五遍；如果新拟合了任何校准器或简单统计基线，则必须每折只用其训练部分拟合。跨链、跨系统加权及置信区间单位也要与旧榜清楚区分。

## 补充榜或待配置资源

| 资源 | 已核实本地情况 | 定位与阻碍 |
|---|---|---|
| AlphaFlow MD+Templates | 源码 `E:\Project\FlexPPI\tmp\third_party\alphaflow`；权重 `E:\Project\FlexPPI\tmp\third_party\alphaflow\params\alphaflow_12l_md_templates_distilled_202406.pt`；作者公开轨迹 `E:\Project\FlexPPI\tmp\third_party\alphaflow\samples\alphaflow_12l_md_templates_distilled_202406`；评估器 `E:\Project\FlexPPI\scripts\evaluate_alphaflow_ensembles.py` | 官方入口依赖 CUDA/OpenFold，当前无 NVIDIA、Docker 未启动，不能本机重新生成官方集合。可继续用作者 250 帧集合中的旧 ATLAS-191 重合 6 个系统做独立案例，不与全队列或 ITrace 混榜。若要拓展到新 ATLAS 队列，要重新审计公开集合重合及训练重合，不能沿用“6”推断全部覆盖。 |
| DynaMine / b2bTools | 源码 `E:\Project\protein_flexibility_models\DynaMine_b2bTools_source\b2btools-3.0.8`；环境 `E:\Project\protein_flexibility_models\envs\b2btools_py312\Scripts\python.exe` | 实测标准 `import b2bTools` 失败，堆栈为 DisoMine 的 `runpsipred_single.py` 抛 `NotImplementedError: Platform not supported yet: win32`。需要 Linux 环境或经过验证的独立 DynaMine 入口；未在此次尝试绕过。预测的是序列内禀动力学性质，可做代理指标补充，不与 RMSF Å 或 χ 转子熵直接等同。 |
| OPUS-Bfactor | 源码及入口 `E:\Project\protein_flexibility_models\opus_bfactor_standalone\opus_bfactor\run_opus_bf.py`；权重 `models\opus_bf_seq.h5`、`models\opus_bf_struct.h5`；环境 `E:\Project\protein_flexibility_models\envs\opus_bfactor_py39\Scripts\python.exe` | 实测 torch `2.2.2+cpu` 与 TensorFlow `2.10.1` 导入成功；但完整流水线未验证，原入口含上游硬编码 Linux 路径、DSSP 入口、固定测试集清单和 56 进程。需要适配后做小样本验证，不应原样运行。输出 B 因子，只做代理补充榜。 |
| B-factor codes repository | `E:\Project\protein_flexibility_models\B_factor_codes_repository`；权重 `checkpoint\b_factor_model.pth`、`checkpoint\b_factor_all_features.pth` | README 给的是训练/测试 notebook 与预处理数组流程，未找到已完成的本项目批量残基映射。作为 B 因子代理候选，不认定即刻可调用。 |
| DEFMap | `E:\Project\protein_flexibility_models\DEFMap`；`model\model_res5A.h5`、`model_res6A.h5`、`model_res7A.h5` 存在 | README 指定 cryo-EM 密度输入及旧 Keras/TensorFlow GPU、EMAN2/HTMD 工具链；其输入信息与 PDB-only 主榜不同。仅可在有相应实验密度及匹配数据时建额外输入条件榜。 |
| Flexpert、BioEmu、CABS-flex 等 | 当前限定路径盘点没有核实出可用部署；旧 `additional_baselines\README.md` 也未将其列作完成全队列结果 | 不将文献清单或名称出现计为本地可运行，更不填充推测成绩。后续需要单独确认代码、权重、输入定义、环境和训练集重合。 |

现有主要横评资源可以支持先完成新版折外预测在旧 520 系统上的统一复测与图表更新；全量 3,068 系统横评、新增 PLANET-MD，以及需要恢复容器的基线，应在总体评估协议确认后分批执行，并明确标注失败、缺失、部分覆盖与缓存复用状态。
