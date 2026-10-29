# 预先固定的抗原候选特征

`extract_features.py`只读取结构、推断NPZ及固定协议，**不读取SPR结果、训练标签或MD值**。主模板2BNQ/1AO7已在协议中固定；每个候选在全部同体系模板上独立计算，输出76行。输出的`candidate_id`是生物学候选ID（与标签一致），`input_id`才是模板前缀的76个推断输入ID。

主文件为`candidate_features.csv`，列为`candidate_id,panel_id,template_pdb,input_id,peptide_sequence,parent_candidate_id,parent_input_id,template_is_primary,is_panel_parent,sequence_n_mutations,BLOSUM62_deficit,delta_TCR_pMHC_BSA_A2,delta_peptide_TCR_contact_residue_pairs_4p5A,R,H,static_eligible,preparation_qc_pass,dynamic_qc_pass,feature_region_id,n_U_RMSF,n_U_Hchi,n_U_all,checkpoint_sha256,...sourcehash`。主分析取`template_is_primary=True`共30行；不能按closest选择或跨模板择优。76行不是76个实验候选。

第一步`freeze-regions`从每模板自身的完整母肽复合物确定U：HLA重链及TCRα/β任一重原子与任一肽重原子距离≤8Å，排除β2m。R使用原始固定主链有N/CA/C的U；H进一步要求至少一个化学χ、排除Gly/Ala、完整canonical重原子身份。U及距离、作者编号、角色、AA和两个有效掩码永久保存为`regions/<PDB>_U8A.csv/json`与`fixed_regions.json`；每个候选完全复用同一清单。

第二步`extract`读取推断`run_audit.json`及每输入`prediction.npz`，强制核模型SHA、完整76输入、原子顺序、AA/主链对应和物理尺度metric。R/H均先作每个U残基的候选减同模板母肽，再对各自固定范围取算术平均；既不取绝对值，也不假定负号代表更高亲和力，不对softplus后mean再次expm1，不使用方差权重。整个肽及U逐残基曲线、HLA/TCR分区均值另存archive文件，不新增为可调主特征。

静态计算仅用root准备的完整重原子结构，不用删侧链的模型输入冒充静态突变结构。Biopython ShrakeRupley的探针1.4Å、每原子960球面点，固定元素半径C1.70/N1.55/O1.52/S1.80Å；BSA=[SASA(TCRα+β)+SASA(pMHC含β2m)-SASA(完整复合物)]/2；接触数为肽/TCR残基唯一配对中至少一个非氢原子距离≤4.5Å。每个候选都减相同模板下母肽值。所有常规重原子完整性和原始NCAO坐标一致性均核查，失败保留行并记QC。缓存以完整PDB SHA和固定设置为键，不按实验结果或图形效果换构象。

要求manifest：fixed_backbone_manifest.json中candidates包含candidate_id（推断ID）、biological_candidate_id、template_pdb、roles（五种角色到链ID）、prepared_input_pdb、peptide_sequence、panel_reference_peptide、is_panel_parent、is_primary_template。static_complete_manifest.json可用candidates/records列表，以相同推断ID索引，每条提供`static_complete_pdb`（也接受complete_pdb/output_pdb/prepared_pdb）。如果字段不同，应只修schema适配，不修改特征定义。

执行分两步，依赖使用现有本地包目录，脚本不安装软件、不联系远程：

```text
python extract_features.py freeze-regions --runtime-path <Biopython/SciPy的现有包目录>
python extract_features.py extract --prediction-audit <正式推断输出目录>/run_audit.json --runtime-path <同上>
```

需要时重复`--runtime-path`。2026-10-21已在root授权后完成五模板U冻结、全76份完整结构静态计算及正式特征提取，`feature_manifest.json`为passed，76份静态/动态/准备QC均通过，主模板候选30个，五个母肽的六项差值均严格为0。特征CSV的SHA-256为`f28d14f1b549cea6bd2c21f883b8854336cd84d95cf36527636a2f259c36fe3e`。提取程序未读取SPR，未做候选结果评价或改变预先定义的范围。实际使用`structures/static_complete_inputs.json`作为完整结构的不可变路径索引，正式预测audit为`predictions/fixed_backbone/run_audit.json`；调用时分别以`--static-manifest`和`--prediction-audit`指定。

首轮冻结的AA字典漏LEU错误已在任何新候选预测前修复，记录见`BUGFIX_AA_MAPPING_20261021.md`。独立测试`input_contract_test.json`核了20AA映射/χ数及全76份60956残基的官方parser身份、编号、原始主链和静态canonical完整性；`static_calculation_selftest.json`核了孤立原子解析面积、无接触BSA=0及唯一残基对计数。所有写入仅限本feature_extraction目录；`feature_manifest.json`记录协议、输入、U、NPZ、checkpoint和输出文件SHA。
