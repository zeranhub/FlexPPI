# 图3：KIX与泛素的局部动态验证

`fig3_kix_ubiquitin.svg`是最终可编辑排版；结构面板内嵌，图中曲线、坐标轴和文字保留矢量。直接运行`python compose_figure.py`会使用`source_archive/data/`中的冻结结果与`source_archive/kix/`、`source_archive/ubiquitin/renders/`中的真实结构渲染，生成`reproduced/fig3_kix_ubiquitin.{png,pdf,svg}`，无需旧实验目录或原始轨迹。也可从项目根目录运行`python paper/figures_v58/reproduce.py 3`。

KIX完整结构及已保存的局部视角位于`source_archive/kix/KIX_case_validation.pse`，用PyMOL打开即可编辑。`reproduce_kix.pml`按文件自身位置重建场景；其`reproduce_kix.py`和`whole_protein_crop.py`使用随包真实PDB及数值源表重绘。泛素整体结构位于`source_archive/ubiquitin/ubiquitin_2QHO.pse`，Q41与Q62另有独立PSE和PML入口；`render_from_bundle.py`默认从本目录`inputs/`和`representative_frames/`读取已归档的真实参考结构与MD代表帧，重新渲染时需要PyMOL。`renders/`保留全图与局部裁切，`crop_manifest.json`保存固定裁切框与相机。日常编辑、打开PSE和重绘图板不需要重新抽取MD帧。
