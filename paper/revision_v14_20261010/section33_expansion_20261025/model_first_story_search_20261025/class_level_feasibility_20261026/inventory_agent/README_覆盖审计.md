# 类别分析的数据覆盖审计（2026-10-26）

只读既有预测、原始注释、链序列、接触标记及家族分组；未修改论文、既有数据，未运行任何模型训练、推理或分子模拟。所有类别均在检查模型表现前按结构和注释纳入；本结果只回答能否做分析，不据此宣布模型已经找到类别规律。

原始完整非 TCR 队列为 `E:\Project\FlexPPI\experiments\non_tcr_case_screen_20261025\all_616_record_screen.csv`，实际616记录、570个不同PDB、261个 all-chain connected-component family。其 `all_chain_metrics.csv` 实际1535链，提供 sample_id/pdb_id/chain/description/identity/n_residues，以及各链 Full/Intra 的 RMSF/Hχ MAE、相关系数。`E:\Project\FlexPPI\experiments\_sources\baselines\cohort\chain_manifest.csv` 共4973链，提供序列；`...\family_folds\record_family_folds.csv` 共3063记录，提供 component_id 与 family holdout fold；`feature_baseline\features_all.parquet` 共961321残基，提供跨链 CA 8Å 接触及两指标有效掩码。这里家族是全链相似性连通分量的划分单位，不等于蛋白功能类别：不少抗体相关体系会经连通关系被并入同一大分量。

`family_folds\evaluation\per_record_metrics.csv` 共12252行 = 3063记录 × 2模型 × 2指标。其同目录 evaluation_audit.json 明确训练完成10个 family OOF fit，对每个记录仅用对应家族留出的折预测。以上指标来自原始数据标签，不能冒充统一MD观察时窗的结果。需要注意，`case_selection_audit_20261026\all_616_record_screen.csv` 这个副本现有只有3条CaM记录；它的名字不能作为616样本的证据。

按5–40 aa短链、同记录受体至少70 aa、身份不同、双方均有跨链 CA 8Å 接触且两指标有效位点至少1个，得到105条短链、92记录、90个PDB、67个家族；其中严格双链体系56记录、46家族，不存在短链其实只接触第三条链的歧义。具体成员见 `short_peptide_5to40_candidates.csv`。例包括CaM–CaMKI、CaM–IQCG、雄激素受体–NCOA2、vinculin–talin、TroponinC–TroponinI、PP1α–PPP1R10、KAT6A–HistoneH3、KEAP1–NRF2等。此筛选得到的是短肽介导识别的结构候选，不能仅因PDB短链就叫天然无序蛋白，仍须核查构建体、天然片段边界和生物作用。它与现有83记录统一0–100ns队列的record交集为0。

同一受体完全相同序列且至少70 aa，要求至少两个不同PDB、至少两种不同伙伴身份签名（排除仅同PDB不同record），得到23个序列组、7个家族；严格双链子集15个序列组、5个家族。成员表 `exact_receptor_context_members.csv`，分组表 `exact_receptor_different_partner_candidates.csv` 和 `two_chain_exact_receptor_different_partner_candidates.csv`。例包括RNaseA、两种溶菌酶、trypsin、H-Ras、泛素、ISG15、xylanase inhibitor及UDG inhibitor。后两者是同一抑制剂面对不同酶，不宜叙述成受体面对不同抑制剂。不同命名可能来自同源伙伴而非不同功能伙伴，必须核查UniProt/构建体；多链体系不同伙伴签名也可能只是增加附属组分，故双链子集更直接。15个完全相同序列组不等于15个独立家族，尤其不同长度的泛素在同一受体家族。

单靠链名称中的酶和抑制剂词汇得到31记录、10个家族，见 `enzyme_protein_inhibitor_description_candidates.csv`。该清单是人工核查起点而非完整确证队列；不能用此粗计给类别最终下结论。

现有统一观察窗口83记录队列为 `matched_complementarity\cohort83_record_family_inventory.csv`，83记录、42家族。其中Atlas30记录全单链；Dynarepo23记录中17多链、6单链；ITrace30记录均多链，但同一个家族。整体多链仅47记录、13家族。30个TCR体系可以做该家族内部分析，但不能称30独立家族；83总体也不应全部称PPI复合物。

所有数字和字段已保存于 `inventory_summary.json`，重现脚本为 `inventory.py`。本审计没有因模型表现好坏删样本或分组。
