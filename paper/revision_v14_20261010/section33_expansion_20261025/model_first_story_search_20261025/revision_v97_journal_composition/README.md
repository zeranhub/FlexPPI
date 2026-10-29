# v97 期刊参考组图重排

179×200 mm，来源 v96。两个结构集中顶部，对应 b/c 两列序列图；d 全六 Gln，e 完整2×2，f 全83记录。853个绘图数值逐位与前版完全一致，MAE由同一冻结表计算。模型、MD、掩码和结构坐标未变。

figure3/compose_fig3_v97.py 是自包含重绘脚本，使用既有科学 Python。compose_body_v97.py 是布局源；prepare_v97.py 用于从同目录 v96 初始化，普通重绘无需运行初始化。pymol/ 提供原生PNG、可编辑PSE、精确PDB和相机元数据，不含已删除局部图。build_v97_document.py 使用 bundled Python 和旁侧 v96 Word，仅更新主图3及图注。Word已检查29页，11页为唯一变化页，已实际打开；另28页与已验证v96完全相同。

参考与设计取舍.md 和 reference_notes/ 保存9篇11张实际观图记录及统计口径。发表图仅作设计参考，未复制为本文素材。源包不含原始MD完整轨迹或训练权重。
