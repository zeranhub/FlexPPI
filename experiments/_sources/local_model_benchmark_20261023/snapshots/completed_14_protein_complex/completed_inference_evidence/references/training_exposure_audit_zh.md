# 本轮横评的训练暴露与OOF划分审计

本审计只读取本地已冻结队列、官方小配置、公开ID清单及预测导出清单，没有加载模型权重张量、运行推理或训练、访问网络、重做序列聚类，也没有读取模型性能表。当前3,063条主队列及预测、排名完全保留，仍应称为**已发布权重与共同参考标签的评价**，不能称为所有模型采用相同训练集的独立留出竞赛。逐条敏感性分析membership先于性能读取冻结；后续须注明为本次审计衍生的敏感性分析，完整报告train / val / test / unlisted / unresolved / unknown，不按结果挑选子集。

“精确ID”必须同时满足同一PDB编号和原始作者链编号：仅统一PDB大小写及点号/下划线分隔符，保留链编号大小写。同PDB不同链不计为精确同样本；ID相同也不证明坐标、轨迹、标签或历史优化批次完全相同。没有ID命中不等于排除同源序列、其他训练数据或预训练暴露。

## 外部权重与官方清单的绑定强度

| 本轮权重 | 能确认的绑定 | 本队列ATLAS 1,929条的精确PDB＋链命中 | 可采用的表述 |
|---|---|---|---|
| DynaProt-M | 随权重README与data_config指向`atlas_chains.npy`；官方Dataset追加`_train/_val/_test`并解析包内同名清单 | train **1,262**；val **39**；test **82**；unlisted **546** | 发布checkpoint配置所记载的训练清单暴露（documented training-list exposure）；不是独立核实权重历史 |
| BackFlip-1.0 | 随权重config与README确认Flexpert split名称；原`flexpert_*.csv`本地未找到 | 官方候选Flexpert ID清单：train **1,108**；val **137**；test **138**；unlisted **546** | checkpoint绑定split名称，精确ID名单来自官方交叉来源；尚未绑定原checkpoint CSV的逐字节清单 |
| BackFlip-1.0-seq | 与1.0相同的Flexpert路径声明；有独立权重哈希 | 同上 | 同上；不能把候选清单重叠升级为已核实每次训练读取历史 |
| BackFlip-2.1 | 配置只有`alphaflow_train/val/test.csv`路径名，本地没有可绑定的ID文件 | **未知，不给出重叠数** | 不能因“AlphaFlow”名称相近而套用DynaProt的atlas清单 |

DynaProt发布atlas三组清单原规模为1,266 / 39 / 82；其总表有1,390个ID，三组并集1,387个，未分配的`6j33_B / 6qfk_A / 6sup_A`均不在本队列。Flexpert三组原规模为1,112 / 139 / 139：DynaProt仓库的三份`flexpert_chains_*.npy`与Flexpert官方仓库`data/atlas_splits.json`在规范化分隔符后逐组集合完全相同。这支持“官方Flexpert候选划分”核对，但缺失的BackFlip发布CSV仍单独保留限制。val、test命中分别报告，不能都称为训练集命中。

ATLAS的原PDB＋链可直接由`sample_id`解析；其冻结结构链名`_`不能拿来替代原作者链。DynaRepo的`A00xx`编码由冻结`biological_system`恢复PDB；ITrace也由该字段和`itrace_`来源核对原PDB。但现有`chains.csv`保存的是MD输入链名，先前冻结注释已明确原始PDB链可能重命名，未找到经认证的逆向作者链映射。因此这1,134条记录的3,044条组成链均保留原链解析未解决状态，绝不把当前A/B/C等链名直接提升为原始PDB链。

其中七条DynaRepo记录与公开清单有同PDB命中，需要继续保留为**unresolved**，不计精确同样本：

| 记录 | PDB | DynaProt atlas清单的同PDB角色 | Flexpert候选清单的同PDB角色 |
|---|---|---|---|
| A000E | 1B2S | train | val |
| A001B | 1L2W | train | train |
| A0028 | 1WQJ | train | test |
| A002X | 2GH0 | train | train |
| A005S | 3OWT | train | train |
| A00BI | 5E3E | train | train |
| A00AX | 6FUB | val | train |

其余883条DynaRepo和全部244条ITrace的原PDB均不出现在这两套公开清单中，故在这两份特定ID清单内归为unlisted，并继续注明原链映射未解决。这个PDB层面的排除不依赖猜测链名，但不能表述为“无任何训练暴露”。每个外部方法的3,063条membership中，DynaProt为train 1,262 / val 39 / test 82 / unlisted 1,673 / unresolved 7；两版BackFlip-1的候选清单为train 1,108 / val 137 / test 138 / unlisted 1,673 / unresolved 7；BackFlip-2.1全部为unknown。

## 本文Full / Intra也不满足序列家族留出

对`neural_exports/evaluation/prediction_manifest.jsonl`进行了流式核查，实际为**6,126行**，Full、Intra各3,063行，身份不重复。每条记录的顶层fold、嵌套metadata.fold均与`records.csv.original_fold`一致；每条记录都在对应冻结test清单中，**没有一条出现在该模型对应的train清单中**。五对冻结train/test文件的并集均为完整3,063条、交集为空。原始2,294个生物学分组没有跨原训练折，支持记录/生物学分组OOF的表述。

然而，现有1,662个family中，**63个跨original_fold**；Full与Intra各有**904条测试记录**的同family其他成员出现在对应训练清单中。新`family_fold`确实把这些family完整放在同一折，但本轮复用的预测不来自这套新划分，不能把family_fold或family bootstrap的存在解释为已经按家族隔离训练。244条TCR记录在实际original_fold 0–4中的数量为49 / 48 / 49 / 49 / 49，而family_fold全部为0；故它们属于记录OOF评价，不能称为TCR序列家族外推。预测折号与family_fold数字碰巧相同的590条也不改变该结论。全部导出行声明`source_split_version=state_rmsf_hchi_v2_performance_first`。

因此，**family成组bootstrap只处理本组分数之间的依赖，不能补足训练时的家族隔离**。这项限制与外部权重的已知/候选训练清单重叠应同时披露，不能只讨论外部模型。`strict_external_evaluation`是另一项32条外部记录、仅以source/family-purged 2,431条训练的独立实验，已有清单明确其范围；它没有用于替换当前3,063条Full/Intra预测，也不能作为当前主榜已满足家族留出的证据。

## 冻结交付与复现

`training_exposure_membership_frozen.csv`共**18,378行**（DynaProt、BackFlip-1、BackFlip-1-seq、BackFlip-2各3,063行，加Full/Intra各3,063行），主键`method + sample_id`。核心列为`method, checkpoint_sha256, sample_id, dataset, group, documented_split, binding_strength, exact_identifier_match, identifier_resolution`；另含原PDB、可解析原链、公开split家族、同PDB候选、original/family fold、本文同family训练成员标志及快照时间。Full/Intra的`documented_split=record_biological_group_OOF`不冒充外部公开train/test角色。

`training_exposure_chain_details.csv`有19,892行，保留4个外部权重×4,973条链的可解析性及未提升为精确匹配的当前MD标签候选；`training_exposure_record_details.csv`有12,252行外部逐记录汇总；`training_exposure_official_split_ids.csv`保留全部公开ID及源哈希；`training_exposure_own_oof_details.csv`保留6,126条导出逐条核对；`training_exposure_own_family_details.csv`保留1,662组跨折记录；`training_exposure_summary.csv`提供按模型/来源汇总。协议、40份本地来源SHA及全部CSV的SHA在`training_exposure_protocol_sources_sha256.json`，membership及协议的最终SHA另存`training_exposure_delivery_sha256.json`。复现脚本为`audit_training_exposure.py`，不导入torch，不读取性能表。

来源定位：[DynaProt权重随附映射](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/weights/dynaprot/dynaprot-M/README.md:3)、[数据配置](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/weights/dynaprot/dynaprot-M/data_config.yaml:7)、[Dataset分组解析](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/vendor/dynaprot-main/dynaprot/data/datasets.py:32)、[BackFlip-1官方训练划分说明](E:/Project/FlexPPI_baselines/BackFlip/README.md:97)、[随权重Flexpert路径](E:/Project/FlexPPI_baselines/model_cache/backflip/backflip-1.0/config.yaml:18)、[BackFlip-2.1路径](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/backflip/weights/backflip-2.1/config.yaml:18)、[主队列导出协议](E:/Project/FlexPPI/experiments/core_corrected_pipeline/neural_exports/evaluation/evaluation_audit.json:1)、[独立strict实验协议](E:/Project/FlexPPI/experiments/core_corrected_pipeline/neural_exports/strict_external_evaluation/evaluation_audit.json:1)。
