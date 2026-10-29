# 图 3 图注草稿（CaM–CaMKI 版本，v4，游离/结合两态）

**图 3 | CaM–CaMKI 肽段的伙伴依赖折叠。** a，两态示意：左侧为游离 CaMKI 肽的无序系综（示意图），中间箭头表示 CaM 结合，右侧为 1MXE 复合物的 PyMOL cartoon；游离态是模型预测的示意，结合态有 MD 参考。b，肽段逐位点 Cα RMSF：游离态预测（灰虚线，同一完整模型移除 CaM 后得到）、结合态 MD（黑）、结合态完整模型（蓝）、结合态链内模型（橙虚线）。c，结合态界面（n=43）与界面外（n=126）的 RMSF MAE；蓝为完整模型，橙为链内模型。d，肽段平均 RMSF：游离态预测 4.22 Å，结合态 MD 0.79 Å、完整模型 0.86 Å、链内模型 2.00 Å；箭头标出结合带来的 −3.36 Å 变化。e，伙伴移除对照：折叠伙伴（1XWD、3QT2、2ZXX）变化 0.04–0.23 Å，短链/肽段（1MXE、3WWQ、2HWN、3P8M、6A7V）变化 1.19–4.73 Å；这说明模型学到的是“短链有序性依赖伙伴”的上下文规则，而不是某条肽独有的识别。

**证据层级要写清楚**：结合态的蓝色/橙色曲线和 c 面板是模型对 MD 标签的留出预测；游离态的灰色曲线和 d 面板的 4.22 Å 是同一模型在“移除 CaM”输入下的反事实预测，**没有游离态 MD 参考**。1MXE 原论文只提供定性的“游离肽无固定结构、结合后形成螺旋”，所以我们把游离态标成“model”而不是“MD”，不能把它写成“模型验证了游离态”。如果以后补了游离肽的 MD 或 NMR，才可以把游离态也升级为定量验证。

## 面板与数据对应

- a：游离态为示意图；结合态为 `cam_structure.png`（PyMOL）；渲染脚本 `render_structure_pymol.py`。
- b：`cam_camki_per_residue.csv` 的 `full_no_partner_rmsf`（游离态预测）、`target_rmsf`（结合态 MD）、`full_rmsf`、`intra_rmsf`。
- c：`case_metrics.csv` 中 `model=main_oof, metric=rmsf, subset=interface_ca8/outside_interface`。
- d：`peptide_ablation_summary.csv` 的 `mean_rmsf`。
- e：`partner_removal_controls.csv` 的 `delta`。

## 需要保留的边界

- 游离态没有 MD 标签；灰色曲线是“移除 CaM”的输入遮挡预测，不是游离肽的 MD 轨迹。
- 1MXE 的肽只有 25 个残基；其他短链移除伙伴后也会柔性化，因此 e 面板必须保留短链对照。
- χ 熵只在界面附近优于链内模型；本图以 RMSF 为主，χ 数据放补充表。
- 伙伴移除是输入遮挡（单链或 100 Å 平移），不是重新训练一个删除伙伴的模型。
- MD 是 DynaRepo 的 100 ns 窗口模拟参考（源记录声明 3 个重复，本地只保留 1 条轨迹）；游离态无序来自 1MXE 原始论文，不是本数据直接测量。
