# Protein组差距：来源、折外划分与训练暴露的只读核查

核查日期：2026-10-24。本报告以主图的 **BackFlip-2.1逐链预测后合并（chains）** 为主要对照，assembly结果保留在源表作为补充。结论是：protein组的差距主要集中在ATLAS，而且是排序、最高10%位点AP和标准化形状误差同时落后，不能只用绝对尺度不一致解释。2,203条protein中，1,929条来自ATLAS，占87.56%；ATLAS贡献了protein组Spearman均值差的96.39%、AP均值差的96.59%和绝对MAE差的95.72%。这些比例是来源组成对观察差距的算术分解，**不是已经证明导致差距的训练因果机制**。Full的记录折外预测是真实的，但这批预测并非全组成链家族隔离训练；BackFlip-2.1的训练ID名单仍未知，现有本地证据不能把优势归咎于泄漏。

## 1. 测试蛋白究竟来自哪里

按冻结的实际生物组成成员表，protein包括ATLAS 1,929条以及DynaRepo 274条。274条DynaRepo中，212条来自antigen_single、52条来自test53_single，另10条是本轮逐链身份审计确认的同源/同蛋白拆域多链输入（dynabench_complex 1条、pdbbind_complex 7条、test53_complex 2条）。protein中没有ITrace记录；全部244条ITrace都在protein_complex。完整来源及有效位置数见`source_composition.csv`，此处没有按性能重新选样本。

| protein来源 | 记录数 | Full Spearman | BF2 chains Spearman | Full AP | BF2 AP | Full MAE Å | BF2 MAE Å |
|---|---:|---:|---:|---:|---:|---:|---:|
| ATLAS | 1,929 | 0.788367 | 0.905829 | 0.746420 | 0.871646 | 0.692618 | 0.390281 |
| DynaRepo | 274 | 0.811153 | 0.842089 | 0.675533 | 0.706681 | 0.634676 | 0.539402 |
| 全部protein | 2,203 | 0.791201 | 0.897902 | 0.737603 | 0.851129 | 0.685412 | 0.408828 |

protein整体Spearman差为0.106701，其中ATLAS的加权贡献为0.102853，DynaRepo为0.003848。ATLAS内部的Spearman差为0.117462，明显大于DynaRepo protein的0.030936。BF2在1,898/2,203条protein的Spearman上优于Full，其中ATLAS为1,727/1,929；因此该现象不是只由少量特别差的个案拉动。标准化MAE从Full的0.366319降至BF2的0.236723，ATLAS对应为0.359697与0.218231，说明把每条记录分别去均值并标准化后，差距仍存在。精确分解、所有原始配对差及各source_group结果分别见`protein_gap_source_contribution.csv`、`paired_record_differences.csv`、`source_paired_comparison.csv`。

作为上下文，protein_complex 860条的Full/BF2 chains Spearman为0.808959/0.719627，AP为0.710653/0.509732，MAE为0.299659/0.864586 Å。因此把两组混在一起会掩盖ATLAS与复合物上的不同表现；这并不能反向证明任一组结果错误。

## 2. Full确实是OOF，但不是本轮家族折上的OOF

本次重新只读核对了十份实际训练/测试JSONL与6,126行预测导出manifest。五个fold中，train/test交集均为空、并集均为全部3,063条；test集合与`original_fold`完全一致。Full/Intra各3,063条预测均来自对应测试折模型，顶层fold与嵌套metadata.fold相同，**没有一条被评分的记录出现在生成它的对应训练清单中**。导出审计明确写明“corresponding heldout fold only, no five-model OOF ensemble”，固定epoch40、seed7。这排除了“现在误拿Full的训练集预测，或误把五模型OOF集成当单模型”这两种解释。

但是，当前复用的是`original_fold`，不是后来冻结的`family_fold`。protein组有274条记录的同一连通家族其他记录进入相应训练折，其中ATLAS 139条、DynaRepo 135条；它们的原生物体系没有跨折。全部3,063条记录中同家族训练成员暴露为904条，原63个家族跨original_fold。family bootstrap只改变置信区间的抽样单位，不会把已有模型重新变成家族隔离训练。逐来源核对见`own_OOF_and_family_overlap.csv`，实际train/test组成见`actual_fold_source_counts.csv`。

| protein来源 | Full同家族成员是否在该折训练 | n | Full ρ | BF2 chains ρ | BF2−Full |
|---|---|---:|---:|---:|---:|
| ATLAS | 否 | 1,790 | 0.785195 | 0.906448 | 0.121254 |
| ATLAS | 是 | 139 | 0.829219 | 0.897855 | 0.068636 |
| DynaRepo | 否 | 139 | 0.794324 | 0.844030 | 0.049706 |
| DynaRepo | 是 | 135 | 0.828481 | 0.840090 | 0.011610 |

两类均存在BF2优势，ATLAS无同家族训练成员的子集差距更大。这个后验分层只能说明Full原划分的已知家族成员暴露没有造成观察到的落后，也不能用它断言严格新家族OOF重训练后的成绩。此处没有改成员、重做聚类或重跑模型。

## 3. 外部训练暴露：哪些能查明，哪些不能

`training_exposure_membership_frozen.csv`对BackFlip-2.1的全部3,063条均标为unknown；本机发布权重配置只给出`./alphaflow_train.csv`、`./alphaflow_val.csv`、`./alphaflow_test.csv`路径，没有绑定实际ID文件。本次在`E:/Project/FlexPPI_baselines`及本轮benchmark目录用这些文件名只读查找，仍未找到。故**不能将ATLAS记录分成BF2 seen/unseen，也不能声称BF2在这些测试蛋白上有训练泄漏**。名单缺失也不证明没有暴露。DynaProt的ATLAS名单和BF1的Flexpert候选名单不属于BF2，不能因为AlphaFlow/ATLAS名称相近就移用。

其他发布方法的证据可说明整个横评并非“所有方法同训练集”的比赛：ATLAS中DynaProt-M绑定公开名单的train/val/test/unlisted为1,262/39/82/546；BF1两版官方候选Flexpert名单为1,108/137/138/546，候选名单仍未绑定BackFlip发布CSV的逐字节版本。这些类别及所有数据来源的计数在`exposure_status_counts.csv`；`other_model_published_exposure_strata.csv`完整保留三种已知/候选名单模型的各类别配对成绩。DynaProt在其公开test及unlisted上仍高于Full，但这些类别既不证明独立于所有同源/预训练数据，也不能用于跨模型的test榜，更不能被改称BF2的暴露分层。

## 4. 训练目标与资源配置存在已记录差异，但目前只能作为解释候选

Full原训练payload配置采用从头训练、两个等权目标RMSF与Hχ、固定40轮、不设验证集、固定使用第40轮权重；采样为`source_system_balanced`。实际代码先让ATLAS、DynaRepo、ITrace各占总采样概率的1/3，再在来源内部对体系等权，有放回抽样，每轮抽样数等于训练记录数。每折训练集中ATLAS约占63%的原始记录，却只获得1/3的预期抽样量；按各折实际清单计算，每条ATLAS记录平均每轮约被抽到0.529次，40轮约21.2次。这是配置的期望值，不是已重建的逐步访问日志。ITrace虽只占约8%的训练记录，也获得1/3采样量。由此可提出可检验的解释：当前训练设计更重视多来源与双目标平衡，而protein评价却由ATLAS主导。

BF2发布config的`dataset_type=ATLAS`、`allowed_oligomers=[monomeric]`，训练文件名是alphaflow三份名单，提供的训练器上限为`max_epochs=500`，checkpoint按`valid/loss`保留top5。这里的500是**允许的最大轮数，不能写成已证实实际训练500轮**；发布checkpoint的完整历史与ID名单仍未知。ATLAS单体上的任务/数据侧重点、Full每记录预期训练量、双目标权衡、优化与模型能力差异，都可能共同影响本次结果。当前没有受控重训练对照，不能给这些机制分配因果百分比，不能直接说“只要加epoch或把ATLAS权重调高就能补平差距”。各折预期采样表见`configured_training_source_sampling.csv`。

## 5. 为什么原论文总表看起来没有protein组差这么大

v63论文表1A汇总的是**全部3,063条记录的记录等权均值**，不是protein子集、不是每个来源等权，也不是每个残基等权：Full RMSF MAE 0.577、Spearman 0.796；原BackFlip-2采用逐链2.1，MAE 0.537、Spearman 0.848。本轮现有逐记录表重新聚合全部记录得到Full 0.577104/0.796187、BF2 chains 0.536791/0.847847，与论文舍入值一致；因此在这里没有发现“论文总分来自另一个更有利子集”的证据。新的protein组剔除了Full较占优势的860条复合物，只剩主要为ATLAS的2,203条，均值差自然更大。

论文表1B的32条ThermoPCD是另一外部集合，不能混作protein组或当前3,063主队列。论文补充材料D的全组成链家族隔离评价也是另一次训练，不是本轮Full缓存。正文2.5及表1注记对strict训练样本写2,430，而较早审计写2,431，属于独立外部实验的版本差异，当前protein差距分析没有读取或使用这32条预测，也没有据此改分。论文方法与聚合相关段落原文摘录保存在`manuscript_scope_evidence.csv`；更广泛的论文文本一致性不在本报告范围。

另需避免输入上下文混比：当前BF2 assembly全3,063记录的Spearman为0.828062，chains为0.847847，protein分别为0.897110和0.897902。ATLAS实际单链，两者成绩一致；主要上下文变化发生在多链体系。论文原表1用chains，本报告主结论也用chains，assembly不作为重新挑选的更有利默认口径。

## 复现与文件完整性

运行`C:/Users/DELL/AppData/Local/Programs/Python/Python39/python.exe -X utf8 reproduce_split_review.py`即可由现有表、元数据与日志重新生成本目录CSV和`review_audit.json`。脚本只写本调查目录，无torch导入、不加载权重张量、不联网、不推理；读取结束后重新校验全部输入SHA256，确认其在本次调查期间未变。逐记录支持完全配对，未新增显著性检验、未按模型结果删记录。本报告可支持“ATLAS来源主导了当前protein差距，训练口径与外部暴露不完全可比”的结论；它没有证明BF2泄漏，也没有单凭来源分解证明FlexPPI落后的内部因果机制。
