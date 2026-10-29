# 旧分类字段语义审计：蛋白与蛋白复合物

用户已明确新定义：同一种生物学蛋白的单链、片段或同源多聚体归“蛋白”，不同蛋白伙伴结合归“蛋白复合物”。本笔记只审查既有字段及来源，不创建最终分组，不读取性能表，也不修改冻结 cohort、当前图表或结果。后续应采用独立的实际 MD 链到官方蛋白身份映射，不能把本笔记的诊断数当成最终组数。

当前 `records.csv` 有 3,063 条，SHA256 为 `5c1153ea9717e4c54a790e39dfda82e6d5e48f482ca0cb6a30e1e35cb903f08e`。旧字段的实际统计如下，详见 [counts.json](counts.json)。

| 字段 | 实际取值与数量 | 对新定义的证据力 |
|---|---|---|
| `ppi_classification_original` | high_confidence_ppi 368；not_ppi_single_protein_chain 13；needs_manual_review 3；not_audited 2,679 | 几何接触筛选，未区分同源与异源 |
| `explicit_high_confidence_ppi_or_tcr` | True 612；False 2,451 | 只是前述 368 加 ITrace 244，不能全部直接归复合物 |
| `source_designated_complex` | True 870；False 2,193 | 只是来源标签后缀，不能代替实际输入的生物学身份 |
| `observed_distinct_chain_sequences` | 1：2,194；2：447；3：135；4：40；5：247 | 不同序列可能只是缺失片段、突变或编号解析后的同一蛋白 |
| 实际输入链数 | 单链 2,193；多链 870 | 用于输入范围核对，不能单独判定异源复合物 |

`source_group` 实数为 ATLAS 1,929、pdbbind_complex 384、test53_complex 38、test53_single 52、antigen_single 212、dynabench_complex 204、itrace_complex 244。来源标签的 870/2,193 与多链/单链的数目恰好一样，但成员不同：`A00KQ/3D5O` 标为 pdbbind_complex，冻结输入只有一条 171 残基链；`A00K8/7CQC` 标为 antigen_single，冻结输入有 65、14、107 残基三条链。

旧 PPI 分类真正来自 [audit_pdbbind_ppi.py](C:/Users/DELL/.codex/worktrees/2788/FlexPPI/scripts/audit_pdbbind_ppi.py:80)：至少两条各不少于 20 残基的链，且某链对在 Cα 距离不超过 8 Å 时有至少 5 对接触、每侧至少 3 个界面残基，便记 high_confidence_ppi。它不读取 entity、UniProt 或序列同源性，且只审 PDBbind 来源。`not_ppi_single_protein_chain` 的实际意思是“不足两条达到 20 残基的链”，13 条中 12 条仍有两条输入链。当前全部 384 条 PDBbind 记录与原 [audit.tsv](E:/Project/FlexPPI/data/splits/state_rmsf_hchi_v2_performance_first/pdbbind_ppi_audit.tsv) 的逐 ID 标签一致；原 summary 的 390/406 是更大原审计，不能替代本轮的 368/384。标签经 [prepare_performance_first_v2.py](C:/Users/DELL/.codex/worktrees/2788/FlexPPI/scripts/prepare_performance_first_v2.py:141) 读入 metadata，再由冻结脚本照抄。

来源分组则由 [dynarepo.py](C:/Users/DELL/.codex/worktrees/2788/FlexPPI/src/flexppi/data/dynarepo.py:14) 按 description 的 test53、PDBbind、antigen 等字串生成。冻结脚本第 149 行只检查 `_complex` 后缀，第 150 行只数实际链序列种类。已有 [grouping_protocol.json](E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/cohort/grouping_protocol.json:14) 明确提示：MD 链可能重新命名，旧作者链标签与 entity 的匹配仅作诊断；序列不同不能证明异源。旧诊断的 102 个 one_entity、3 个 one_source_entity、591 个 multi_entity、174 个 unresolved，均不能直接纳入论文二分组。

反例已整理为 [key_cases.csv](key_cases.csv)，每例绑定本地官方缓存路径。`A000A/1AK4` 被旧诊断列为 one_entity，但缓存描述是 Cyclophilin A 与 HIV-1 capsid，UniProt 为 P62937/P12497，说明碰巧相同的作者链字母会误导。`A001Y/1T63` 被旧规则判 not_ppi，仅因第二链只有 14 残基；来源却指向 androgen receptor 与 NCOA2，不能因短肽而自动当作单蛋白。`A003H/2OIN` 的 NS4A peptide 与 polyprotein 共用 P27958，是“相同 accession 未必相同成熟蛋白”的边界；GFP–nanobody 或 DARPin 伙伴缺少第二个 UniProt 也不能自动当成同源多聚体。

仅凭本轮已核对旧资料，可建立一个保守证据下限：2,193 个单链输入加 `A00C7/3H6G` 的两个完全相同 383 残基链、同一来源 entity/P42260，可支持至少 2,194 个“蛋白”输入；ITrace 的 244 个 TCR–pMHC 输入有明确复合物来源支持。其余 625 条是“仍需实际链—蛋白身份映射”，不是宣称它们都无法判定，也不是最终可用于作图的第三组。不同蛋白复合物的实际数量显然还包括该 625 条中的许多记录，应由新的映射审计补齐，不能从这里删去它们。

最小可执行规则是：先对实际输入链做可追溯的官方来源序列/实体映射并处理链改名；映射到同一种成熟生物学蛋白的拷贝、突变或片段归蛋白；映射到至少两种不同蛋白伙伴归复合物；无法唯一映射、仅有截短小肽或成熟产物含义不清者显式待审。不要新增任意序列相似度阈值来强行填满两组。唯一需要在协议中进一步写清的语义是“蛋白身份以成熟功能蛋白为准”：病毒多聚蛋白产生的不同成熟蛋白不能只因共享 UniProt 而合为同源多聚体；同一成熟蛋白的断裂或缺失片段则不能只因序列不同而拆成异源复合物。该说明与用户已给的定义一致，无需重复询问已明确的同源多聚体归属。

所有读取来源的 SHA 与前后未变化检查见 [source_sha256.json](source_sha256.json)。这些新文件是分组解释的独立审计记录，不替代新的最终分组 membership，也不对任何模型重新推理或评分。
