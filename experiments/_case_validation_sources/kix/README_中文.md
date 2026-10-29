# KIX 实例：可复现的 PyMOL 图像包

双击 `KIX_case_validation.pse` 即可打开全部真实结构及五个已经保存的场景。PSE 内已包含原子坐标、颜色和相机，不需要原来的 E 盘文件。在 PyMOL 命令栏输入 `scene KIX_context, recall` 看整体三元复合物；输入 `scene KIX_L12_MD, recall` 看 MD 柔性着色局部；输入 `scene KIX_L12_FlexPPI_fold1, recall` 看当前 fold1 预测着色局部。对应的 `KIX_L12_MD_whole` 和 `KIX_L12_FlexPPI_fold1_whole` 是局部截图相同朝向的完整复合物全景。所有场景都显示完整 A25、B87、C19 共131个残基的 cartoon；局部会话缩小视野也能看到全部蛋白。

从源文件重建：用 PyMOL 的 File → Run Script 打开 `reproduce_kix.pml`，或输入 `@完整路径/reproduce_kix.pml`。脚本会根据自身位置读取同目录文件，因此整个 kix 文件夹可以移动。重建及裁切需要 PyMOL、Python 标准库和 Pillow；如需重新渲染，先运行脚本，再在 PyMOL 输入 `python kix_namespace["render_kix"](cmd)`。随包的 `whole_protein_crop.py` 先完整显示蛋白并 ray 渲染，再按固定像素框直接裁切。

整体图是 1800×1600、300 dpi。两个局部图分别来自 PyMOL 原生 ray tracing 的 6000×4800、600 dpi 完整复合物截图 `kix_l12_md_whole_source.png` 和 `kix_l12_fold1_whole_source.png`，按同一像素框直接裁切，没有重新缩放。`*_crop_box_preview.png` 用红框标出裁切区域；`crop_manifest.json` 记录原图、裁切框、结果尺寸、全景/局部相机和 SHA256，并验证局部图与原图对应像素完全一致。整体结构链 B 为 KIX（灰），链 A 为 cMyb（绿），链 C 为 MLL（紫），完整 L12 残基614–621为橙色。局部使用完整三元复合物，KIX全部87位按MD/模型RMSF百分位着色，结合肽保留绿色/紫色；L12显示完整残基的实际主链与侧链键，主链细棒半径0.07 Å、侧链0.12 Å，614与621标注边界。全部56个重原子及57条真实键形成一个连续片段，包含七条肽键和两处脯氨酸闭环。cartoon关闭环区平滑及片层压平，源原子坐标保持不变。`stick_connectivity_QA.json`记录连接检查。白—浅蓝—深蓝共用0–100百分位色标。

`prepared_2AGH_model1_ternary.pdb` 是真实 prepared 2AGH model 1 三元复合物，坐标未修改；氢原子只是不显示。`kix_both_fold1_source.csv` 包含完整 87 位的 MD 和模型百分位，`kix_l12_fold1_audit.csv` 用于核对完整 8 位 L12 的区域均值。当前使用 `model_fold=1` 的单个模型检查点，没有使用旧 ensemble_mean 上色。MD 是一个匹配轨迹，残基不是独立实验重复。`scene_config.json` 保存各图完整 18 元相机参数，`render_provenance.json` 保存数值来源、SHA256、场景名、像素范围和整体图的链/环标签定位坐标。各 panel 的相机和来源由该文件逐一记录。

本图预测统一使用原稿v13的 chemical_identity_corrected_v1 完整模型、种子7、fold1、第40轮权重。`kix_v13_refresh_QA.json`记录检查点SHA及四状态源文件，`kix_l12_fold1_audit.csv`只汇总该模型的L12。四状态87位排序相关为0.896/0.881/0.729/0.846，top18命中14/16/11/14；MD与模型的局部图共享源结构、相机及0--100百分位色标。
