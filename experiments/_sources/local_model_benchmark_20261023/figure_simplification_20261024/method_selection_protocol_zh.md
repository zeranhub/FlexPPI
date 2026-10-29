# 主图方法精简协议（准备版，等待样本范围澄清）

主图按用户要求移除 `Intra`，保留 `Full`；BackFlip 固定保留 `backflip-2.1 (chains)` 和 `backflip-1.0 (chains)`。这两项是两个主版本在相同“逐链预测、按原记录拼回”的输入上下文下的结果，选择依据是版本代表性及输入上下文一致，未查看成绩后挑选。不同版本的训练数据、架构及残基身份等输入特征仍可能不同，不能把其差异解释成纯粹的版本号或单一模块效应。所有 assembly 与 1.0-seq 结果仍保留在完整 14 项原始结果和补充材料中；本协议只规定展示子集。

主图固定顺序为 `Full`、`backflip-2.1 (chains)`、`backflip-1.0 (chains)`、`DynaProt-M`、`FastProtFlex`、`ANM`、`GNM`、`Static`、`AA`，共 9 行。两行生物分组、三列指标不变，即 a–c 为 protein 组的 ρ、high-RMSF AP、Å MAE，d–f 为 protein_complex 组的相同指标。每组有 24 个真实数值点，合计 48 点；FastProtFlex、ANM、GNM 的 MAE 仍以明确的 relative score 空格表示，两组共 6 格，不做单位拟合或填补。初始排版可将 183×141 mm 缩为约 183×116 mm，保留 6.3–7 pt 字号、3 种语义颜色、直接方法标签及外置图例，再以实际 PNG 复核间距。旧 14 行图、其图注和 QA 均不覆盖。

ECDF 同步固定为 `Full`、`backflip-2.1 (chains)`、`backflip-1.0 (chains)`、`DynaProt-M`、`FastProtFlex`、`ANM`，两组各一轴。移除原 ECDF 的 BF2 assembly，补入 BF1 chains；原 ECDF 本来没有 Intra。每条曲线使用该最终分组的全部记录，不抽样、不平滑。按名称固定线型：Full 实线，BF2 chains 虚线，BF1 chains 实线，DynaProt 点线，FastProtFlex 点划线，ANM dash-dot-dot；三种语义颜色保持蓝 `#31688E`、赭 `#B78038`、青 `#348980`。BF1 与 BF2 用同色不同线型表示同类发布模型，图例写完整版本和 chains 上下文。

实现时必须分开 `RAW_COMPLETE_METHODS`（冻结的 14 项完整结果）与 `DISPLAY_METHODS`（以上 9 项），先验证原始全量覆盖与源 SHA，再按展示集读取既有统计。不能直接把现脚本 `ORDER` 改成 9 项却继续要求 `set(complete_methods)==set(ORDER)`；该脚本目前第 144 行把展示顺序与完整性门控耦合。展示视图的 coverage 为 18 个 method×group 单元，主图数值单元为 48；完整原始数据的 28 个覆盖单元和 78 个主指标单元仍原样保留。`Full` 和 `Intra` 必须使用精确 key 区分，禁止按字符串 FlexPPI 模糊过滤。图源 CSV、图注、footer、plot audit 和 QA 的行数/方法名/ECDF 清单均需同步，不能只删画面标签。图例改为单个 FlexPPI 名称，避免继续称作多个 FlexPPI configurations。

protein 样本范围仍等待用户澄清，当前不能绘图。若只改名称而保留 2,203 条，沿用当前冻结二分组及其既有均值和 family bootstrap 区间，仅应用展示方法过滤；全队列仍为 3,063 条、961,321 有效位点。若排除单链，protein 组只剩预先审计的 10 条多链同蛋白体系，complex 组仍为 860 条，合计 870 条；必须另存明确的记录 membership/范围协议，从既有逐记录分数重新做等记录权重聚合与 family bootstrap，不能将 2,203 条的统计直接改名为 n=10。此时全量源数据仍须通过 3,063 条完整性门控，展示子队列改用新的真实记录数、有效位点数与 family 数，不能继续硬编码 3,063/961,321。ECDF 同样使用各组实际全部记录，10 条会自然呈阶梯形。

两种范围都不重推理、不改预测、真值或掩码，不按分数改变纳入记录。最终另存新脚本、图及 source CSV，绑定当前 plot/scoring/summary/membership/training-exposure SHA；caption 保留发布权重与自身 OOF 的混合训练条件、family bootstrap 不等于训练家族隔离的说明。“9 项”是展示的模型/配置条目数，不能写成 9 个独立模型。本准备阶段未读取性能均值、未运行模型、未生成图，也未修改现有源码、分类、原始统计或冻结文件。
