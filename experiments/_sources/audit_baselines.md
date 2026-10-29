# 新增基线与全链家族划分审计（2026-10-10）

本审计以 `E:\Project\FlexPPI` 的实际基线实现，以及旧工作树 `C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\experiments\corrected_pipeline` 的校正后冻结实验为依据。当前工作树7471不作为实验数值来源。所有新增横评共用3,063条记录、4,973条组成链、961,321个有效RMSF残基和824,988个有效Hχ残基。原520条横评及早期3068条清单不可直接拼入新主表；五条长度排除保持原主实验定义。

## 实际可运行模型与边界

|方法|现有实现/权重|新增资产与真实状态|可比较目标与限制|
|---|---|---|---|
|ProDy ANM/GNM|ProDy2.6.1，无训练权重；原脚本`E:\Project\FlexPPI\scripts\run_prody_baseline_batch.py`|`baselines/run_prody_3063.py`；3,063条全部真实完成，零失败，CPU计时合计676.294秒|完整蛋白组装体CA输入；20个低频非零模式、15Å截断、γ=1；取均方波动平方根。尺度为相对值，只比较排序和独立标准化曲线误差，不赋予Å误差。|
|RMSF-net PDB|上游commit `eb353b0d9d909a6c46d2808467a9d0f16e2bbda7`；`model_cpu/sim_model.pth` 26,857,801字节，SHA256 `cd12dce9489ed57c5b0b8d5400161e1b4e165963622dc78323fdeb07651d8c9b`|`baselines/run_rmsf_net_3063.py`及既有PDB密度adapter；5条烟测全部通过。全队列由先前授权启动，尚无全量完成证明。|保留完整组装体重原子；4Å Gaussian密度、1.5Å网格、12Å边缘填充，是原Chimera molmap操作的独立复现。网络架构/权重不变，但不能称上游完整原生前处理。保留归一化前相对分数，无Å MAE。|
|PLANET-MD full/seq（RocketSHP）|上游commit `917999482343c6400bef51180c817a8af01ee871`；full权重356,157,598字节，SHA256 `b31a5eb353e773c346a413fbca7e63735a6ab53249a628073a7c01ceee393bf9`；seq权重351,949,258字节，SHA256 `54f8ee6df8ed0c45e33995431b720cb7676b38d7521a5101602fae58f574e243`|`baselines/run_planet_md_3063.py`及冻结离线payload；二模型各2链CPU烟测通过，已由方法代理启动完整CPU队列，尚无全量完成证明。|上游CLI只取第一链，因此使用4,973条独立链并逐残基映射回3,063条完整记录。ESM3特征与结构特征层级从checkpoint读取。原MDTraj未缩放RMSF标签为nm，输出固定乘10换为Å；不做误差拟合。SHP不同于本文χ状态熵，不能用于Hχ主比较。|
|FastProtFlex|`E:\Project\protein_flexibility_models\FastProtFlex`及`scripts/fastprotflex_batch.R`，既有R实现|未并入本次3063完成结果|相对刚性/灵活性分数，R/bio3d/ORCA依赖，重原子超过10,000会失败；须先同队列重跑及覆盖审计。|
|BackFlip2.1|校正冻结推理/映射/评分已有真实全量产物|主稿已有其结果，无需重复占GPU|已完整对齐3,063条、961,321残基；整体RMSF优于FlexPPI的负结果必须保留。Released权重训练重叠和模拟协议与本文新训练划分不同。|
|DynaProt、Flexpert、PEGASUS、PIPENN|本次本地资源检查未找到可直接部署的前3者完整实现/权重；DynaProt/Flexpert有文献文件。PIPENN在稿件中为序列界面分类器。|不能记作已完成新增基线|这里只报告本地部署证据，不推断全球代码不可得。PIPENN目标不同，不能直接放进RMSF/Hχ回归横评。|

PLANET-MD原包要求Python3.11与numpy<2；本机esm3.2.3要求Python3.12，不能直接混用。已固定作者指定ESM fork commit `f578182d51a7c78c68cc32ee96a41037938baa8e`（支持Python>=3.11,<3.13）。根代理设置独立`planet_env`，numpy1.26.4优先，沿用已冻结torch2.6.0+cu124。离线payload总文件3,666,680,986字节，tar文件3,666,780,160字节，SHA256 `4ac97db36445ec96788028546408e205549f9e87dc63f3177fd54b7ecb0ce9e9`；包括ESM3序列权重、结构encoder、两预测权重及小型元数据，省略推理未调用的function/structure decoder。默认每条输入300K，930条原记录为310K，这属于跨条件协议差异，不可据此选择性排除样本。

## 同队列、同残基与完成证明

`baselines/prepare_baseline_cohort.py`逐个核对原训练冻结输入/标签SHA256，生成完整结构和独立链PDB。只序列化第一模型已接受的标准蛋白重原子、删除重复altloc/氢；每条完整结构的序列、链、作者残基号、主链缺失掩码和全部主链坐标与原结构逐项一致。独立链只将链标识和编号规范化，原记录索引/作者号保留在冻结truth。输出主manifest不读取新基线误差。Linux实际cohort冻结SHA256为 `b58795275fa2843a502fdbb68f417f9a2b0b80f76e28d226ad5d58c9bf430a91`；structure_manifest为 `c2b7c7b28fb0bff98b3609380eedb5fab8464432e42660a2186cb21bcfdbac30`。本地生成cohort与Linux序列化文件换行不同，路径/hash各自冻结，不混用其manifest哈希。

独立的`baselines/verify_additional_baselines.py`不导入PyTorch或模型、不调用GPU。它验证冻结四文件哈希及全部成员数量、原truth残基身份、每条status JSON完整状态、实际PDB输入哈希、NPZ输出哈希、残基形状/有限性、RMSF-net原作者号/链顺序，以及PLANET逐链序列和nm×10转换。只有全部3,063记录或4,973链成功才写`run_audit.json`并exit0；失败/缺失列入`verification_failures.csv`及`verification_audit.json`并exit1。先前完成证明若存在，重验前改名保留，以免失败检查后仍留下过期的complete。检查器SHA256 `7353a155f4be8a79656cf81309cb2a60fe05329c1bd77abeda17d82fcdd49d18`。真实负向测试目录只有5条RMSF-net输出时，5条通过而3,058条缺失，exit1且无完成证明；在独立副本注入NaN并更新文件哈希后仍被有限性检查拒绝。

`baselines/score_common_residues.py`已有`--rmsf-net`、`--planet-full`、`--planet-seq`入口。必须先独立核验，再使用`--require-complete`；评分脚本自身不替代status/文件hash核验。所有方法按同一冻结有效RMSF残基比较，每条记录一票，组成链在记录内合并。对于相对分数，报告Spearman/Pearson、分别减去各自均值并除以自身标准差后的曲线误差，以及前10%残基检索精度；这种标准化只比较曲线形状，不从真实值拟合预测Å尺度。对PLANET及FlexPPI保留原单位MAE/RMSE。缺失结果必须先披露，不能静默取交集后称全量。

ProDy同一完整队列复算结果：ANM、GNM、Full、intrachain平均Spearman分别0.715、0.663、0.796、0.777；标准化MAE分别0.450、0.487、0.371、0.399。所有方法覆盖全部3,063条、961,321个有效残基。ATLAS中ANM rho0.785接近Full0.788，而DynaRepo0.635对0.774、ITrace0.453对0.938，应该保留这种来源差异，不能只宣称普遍大幅超过弹性网络。

## 全组成链家族划分与新外部去同源

原MMseqs搜索的3,770条长唯一序列中，低复杂度A009Z链C（21aa，`GKPSRPPRPSRPPPPTPRRPA`）没有任何查询行。不能把它等同于无同源。已在hit-gpu CPU独立重跑全库关闭mask/comp-bias搜索，原TSV完全保留，新结果另存`family_folds/all_chain_unmasked_similarity.tsv`。参数为`--min-seq-id 0.3 -c 0.7 --cov-mode 0 -s 7.5 --mask 0 --comp-bias-corr 0 --max-seqs 3770 --threads 2`，14列格式记录在日志。新TSV SHA256为 `d68b6d479cdce81071c1a2dd92dff1737e3e0ee6b6c8f111fe7532c946354d8f`；原TSV `d7427d96738d2e80ec153568eba29ce1388ceea4c3eb9199eb225900fd7278de`。输入FASTA3,770条与全部chain_inventory序列一一匹配且每条长序列有自比覆盖；275条训练短链（174种唯一短序列，含外部）用exact相同序列连接，未按missing忽略。

`baselines/prepare_family_folds.py`合并原/补搜长链边（identity>=0.30、整数对齐跨度/原链长度>=0.70，双向均满足），全部长/短链exact序列、同来源biological_system，将完整复合体记录按传递连通组整体分折。外部32记录参与依赖图，但不进入主3063折。合格非self独立边61,191条，比原增加544条；训练有1,662家族组，最大训练组633条/含外部665条。五折记录数633/608/608/607/607，多链记录数462/112/105/92/99。每一条exact/长链边及生物组都核验未跨折。分配只看成员数和来源作为平衡，未读取标签或预测；Full+intrachain共10fit，从零训练，seed7、40epoch，与原预算一致。

全部244条ITrace记录被共同组成链连入同一大组，均为fold0测试；该折训练无ITrace，另外4折测试无ITrace。因此ITrace结果是一次严格家族/来源外推，不可填写虚构的五折均值±SD。全队列应使用记录macro汇总，同时保留各折/来源结果。序列搜索的相似性阈值不是检测全部进化同源的数学保证；exact短肽/共有伙伴可保守合并不同结构域家族，组大小和来源集中程度必须公开。

补搜将A009Z纳入外部相连组，新的严格外部删除633条、保留2,430条（原协议为632/2,431）。新family fold0训练集合恰等于新外部去同源集合，因此其Full/intrachain最终40epoch权重可直接复用预测固定32条ThermoPCD，无需另加2fit。`family_folds/evaluate_external_2430.py`已准备并上传，旧2431结果不覆盖。对应Static/AA CPU已真实独立重拟合，训练625,577个RMSF有效残基、533,074个Hχ残基，总用时22.999秒；输出`baselines/static_external_2430`。Static RMSF MAE0.471/rho0.645，Hχ MAE0.211/rho0.886；AA RMSF MAE0.669/rho0.279，Hχ MAE0.243/rho0.855。外部仍为一个连接家族，不提供独立家族CI或p值。

## 执行与资源

GPU由根代理串行调度神经训练和MD；新增基线已全部改为CPU。`baselines/gpu_additional_steps.json`因根调度器兼容而保留原文件名，实际只有5个CPU步骤：等待RMSF原进程退出、加锁检查原wrapper缓存并独立核验；等待PLANET-full CPU队列并只读复验；等待PLANET-seq CPU队列并只读复验；原划分共同评分；新家族划分共同评分。没有GPU推理命令，也不再次启动PLANET全队列。RMSF通过`/proc`核对实际wrapper及目标输出，`flock`锁会继承到子进程防止调度器重复启动；既有失败记录不加`--retry-failed`自动抹去，独立核验仍将拒绝失败。两PLANET状态文件为`planet_full_cpu_queue_state.json`和`planet_seq_cpu_queue_state.json`，只有原队列已完成独立verifier、原float32/精确10倍转换全扫以及最终冻结输入检查后才可complete；等待脚本每20秒检查、最多24小时，failed状态或已观测失败记录立即exit1。它逐条只读检查当前status/NPZ/metadata哈希，保留原`run_audit.json`证书，以免重写含时间信息的证书破坏队列保存的哈希。

每个调度步骤经`execute_frozen_baseline_step.py`检查独立`additional_cpu_execution_contract.json`的SHA256及16个已冻结文件，再由`check_additional_baseline_freeze.py`核对所有3,063个完整结构和4,973个链输入、cohort四文件、RMSF源码/adapter/权重、PLANET源码与53个payload文件（总3,666,680,986字节）、CPU `SMOKE_PASSED`及对应审计/真实输出、`PLANET_ENV_READY`。实际强检已在CPU通过，用时约6.3秒，无模型推理。负向测试在独立fixture证明：failed状态、失败记录和缺失状态超时均exit1，已有伪旧complete证书不能绕过；故意改变预期cohort哈希也被拒绝。负向报告为`baselines/guard_negative_tests_20261010/test_audit.json`。执行contract SHA256为`419e132ac514e519d10105fd28194d13f688381b46140eb0a52ee6a862e81fa7`，调度JSON SHA256为`d547800de05c1b9230d1334d2551faa7ccd698aac7f93c212f196841862b9fa4`。新脚本只新增独立检查/调度文件，未改变仍在运行的模型wrapper。

Full/seq completion文件分别为`/root/flexppi_completion_20261010/baselines/planet_full/run_audit.json`和`planet_seq/run_audit.json`；两套共同评价分别在`additional_original_oof`、`additional_family_oof`。评分前必须再次检查两CPU state及当前完成证书，不能因为旧证书存在就忽略failed状态。JSON本身只是冻结执行协议，不能作为已完成实验结果。

RMSF-net为CPU网络，烟测5条耗时22.546、10.833、3.484、16.158和13.234秒，不据此对复杂结构作未经测量的精确总时长承诺；网格大小由完整组装体空间跨度决定。PLANET的ESM3序列权重约2.804GB，预测权重约0.35GB/模型，并有结构encoder和输入长度相关开销。离线CPU烟测415/130残基两链，full耗时10.093/2.740秒、seq耗时7.485/1.694秒，原始float32权重不变、无训练或校准；冻结最大链长1,024。方法代理2026-10-10 07:36UTC启动Full/Seq各2线程队列，父/子PID73982/73988和73983/73989；07:38UTC已有22/24链成功、零失败，状态仍running。启动验收确认无GPU上下文、真实预测与烟测输出哈希一致。这些局部计时不能充当全队列耗时或覆盖结果，长链及异常结构失败将如实保留、由verifier拒绝。旧严格训练实测Full约1.018小时、intra约1.263小时，仅可作原机器、原2431集合的历史成本证据，新家族10fit成本随来源/长度混合变化。服务器计费单价未提供，本审计不编造人民币/美元费用。

精确CPU核验模板：`/root/flexppi_completion_20261010/cpu_env/bin/python /root/flexppi_completion_20261010/baselines/verify_additional_baselines.py --cohort /root/flexppi_completion_20261010/baselines/cohort --output /root/flexppi_completion_20261010/baselines/rmsf_net --model rmsf-net --expected-cohort-freeze-sha256 b58795275fa2843a502fdbb68f417f9a2b0b80f76e28d226ad5d58c9bf430a91 --weights /root/flexppi_completion_20261010/baselines/RMSF-net/model_cpu/sim_model.pth`。PLANET替换output/model/weights为JSON中对应项。完整的原划分/家族划分评分命令亦保留在JSON中，必须在所有前置完成证明有效后执行。

环境与权重证明见`baselines/cpu_env_freeze.txt`、`weight_inventory.json`、`planet_payload_freeze.json`；全量ProDy输出和同残基评分见`baselines/prody_common_evaluation`；新Static拟合证据见`baselines/static_external_2430/run_audit.json`；家族分折输入/配置哈希与规模见`family_folds/training_freeze.json`。未把未完成的RMSF-net/PLANET或新神经训练数值写成实验结果。
