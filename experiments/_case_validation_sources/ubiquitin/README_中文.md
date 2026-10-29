# 2QHO 泛素–EDD1 UBA：完整蛋白截图与局部裁切的 PyMOL 图包

f 子图保留完整泛素–UBA 复合物：参考输入的 A 链 72 个泛素残基和 B 链 47 个 UBA 残基共 119 个残基均以 cartoon 表示。局部图的获取顺序为在 PyMOL 中表示完整蛋白并保留真实 MD 侧链与各帧骨架、按关注残基的方向渲染整蛋白 6000×5000 像素原图，再只对原图执行固定像素框裁切。邻近残基以浅灰细棒显示完整重原子键，包括 N/CA/C/O 及侧链的连接；关注位点每个实际 MD 谷氨酰胺显示完整 9 个重原子和 8 条真实键，主链用灰色细棒，侧链用橙色或紫色粗棒。各代表帧自己的整条泛素 72 个 Cα 灰色细迹线显示侧链归属，彩色根部小球标出该帧的真实 Cα。不同帧采用不同透明度便于分辨。全部原子保持随包真实帧坐标，详见 `stick_connectivity_QA.json`。

`q41_whole_source.png` 与 `q62_whole_source.png` 是完整复合物的原始 PyMOL ray 图；`q41_crop_box_preview.png` 与 `q62_crop_box_preview.png` 在整体图中用红框标出裁切位置；`q41_zoom.png` 与 `q62_zoom.png` 是供论文 f 子图使用的原图像素裁切结果。`crop_manifest.json` 逐项记录完整图与局部图的尺寸、裁切框左/上/右/下像素坐标、全景与焦点相机、完整原子投影边界、源图/局部图 SHA256、119 个 CA 的 cartoon 覆盖数，以及裁切像素与源图对应区域完全相等的核查。

双击 `f_Q41_local.pse` 或 `f_Q62_local.pse` 可直接打开对应焦点视图；整个参考复合物始终显示，缩小视角即可继续观察完整蛋白。每份 PSE 保存五个场景：`01_complex` 为带位点与链标签的整体定位，`03_Q41_whole_for_crop`、`03_Q62_whole_for_crop` 为各自局部图的完整截图视角，`02_Q41_MD_observed`、`02_Q62_MD_observed` 为相应裁切焦点视角。通过右侧 scene 按钮，或在命令栏输入 `scene 03_Q41_whole_for_crop, recall`、`scene 03_Q62_whole_for_crop, recall`，可以从局部直接切回同方向的完整蛋白截图。完整与焦点场景使用同一相机距离和完整蛋白的近远裁切面，仅调整图像平面位置和视野角；焦点场景没有再次按局部厚度截取蛋白。

`f_Q41_local.pml` 与 `f_Q62_local.pml` 可通过 PyMOL 的 File → Run Script 打开；脚本按自身位置载入 `ubiquitin_2QHO.pse` 并恢复焦点场景，能够从任意工作目录执行。结构坐标、实际 MD 代表帧、完整 protein 表示和五个视角均嵌入 PSE，不需要原来的 D 盘轨迹。根目录的 `打开f_Q41局部_PyMOL.cmd`、`打开f_Q62局部_PyMOL.cmd` 是本机的双击入口；换电脑时直接打开 PSE 即可。

从源文件完整重绘，在本目录执行 `D:\Anaconda\envs\pvb_cpu\python.exe render_from_bundle.py`；也可以运行 `rebuild_relative.pml` 后在 PyMOL 输入 `python _render_ns["render"](_bundle)`。完整重绘会重新 ray 渲染整蛋白截图、生成精确像素裁切、保存五个场景与两份独立 f PSE。日常在 PyMOL 用 File → Run Script 运行 `rebuild_relative.pml` 时，默认只重建完整蛋白表示和五个场景，并校验已有截图与裁切的相机、边框和文件哈希，不重复耗时渲染。脚本通过自身目录读取完整参考 PDB、六个 MD 代表帧、`whole_protein_crop.py` 和源表，不依赖当前工作目录；只需要具有 NumPy、Matplotlib、Pillow 的 PyMOL Python 环境。`views.json` 保存全部五个完整的 18 项相机参数，`render_metadata.json` 记录版本与渲染尺寸，`pse_scene_audit.json` 对每个场景检查所有 119 个 CA 的 cartoon 表示与未改变的源坐标，`f_local_entry_QA.json` 记录直接打开入口。

输入是 DynaRepo A003R 的 `structure.pdb`，结构身份为 2QHO 泛素–EDD1 UBA。实际输入链 A 有 72 个泛素残基，编号 1–72；链 B 有 47 个 UBA 残基，DynaRepo 采用串接后的全局编号 73–119。Q41 和 Q62 确认为 A 链 GLN 41/62，对应缓存零起始索引 40/61。`residue_mapping.csv` 记录全部 119 个残基的输入编号、链内序号、缓存全局索引。包内 PDB 始终保留实际输入编号，没有把 UBA 的 DynaRepo 编号冒充原生作者编号。

`representative_frames.csv` 是每根彩色棒的来源表：包括轨迹零起始帧号、时间、联合角度类别、三个 χ 角、全部标签帧占据率、叠合残差和叠合残基。每个位点按全部 5000 个标签帧的类别占据率排序，取最高三个类别，并从该类别内选择三个 χ 角最接近圆周平均角的**实际帧**，未计算平均构象、未人工摆动侧链。原轨迹共 5001 帧；标签按前 5000 帧分成 5 个 1000 帧窗口，最后一帧不参与标签。Q41 展示类别 16/9/15，占据率 26.50%/21.70%/20.38%，实际帧 2397/401/1171；Q62 展示类别 3/4/13，占据率 19.12%/12.56%/12.56%，实际帧 924/1143/1467。完整 27 类频率和全部 5001 帧角度保留在 `inputs`，不存在删除未显示类别的分析。

每个真实帧先按该位点参考结构 9 Å 周围的泛素残基的 N、CA、C、O 主链原子叠合到输入结构。此处 9 Å 指任意残基原子与所选 GLN 任意原子的参考结构空间距离；Q41 使用 151 个主链原子，Q62 使用 92 个。所有代表帧的所有重原子一同做同一个刚性旋转和平移，侧链内部几何没有改变。Q41 三个局部主链 RMSD 为 0.700/0.764/0.960 Å，Q62 为 0.511/0.534/0.630 Å。叠合选择是两个位点各自固定的参考邻域，未根据代表帧结果调整。

这些彩色结构表示 **MD 观测**，模型预测的是 RMSF 和侧链熵，并未生成这些侧链构象。角度类别采用原标签的三角联合分箱，部分不同类别处在 ±180° 或其他分箱边界附近，实际空间朝向可能接近；因此不能把三个类别直接叫作三个物理上完全分离的 rotamer。结构图用于观察 Q41 周围较多的蛋白邻域与 Q62 向外的局部环境；全部 27 类频率及覆盖曲线才是侧链分布宽度的定量依据。`reference_environment_contacts.csv` 记录实际输入中侧链到邻近残基重原子的最小距离，属于结构描述，不直接等同于结合能或对熵差异的因果证明。`reference_solvent_accessibility.csv` 用完整复合物 935 个重原子、1.4 Å 探针和 PyMOL 3.0.0 的 `get_area`（`dot_solvent=1`, `dot_density=4`）计算输入结构的溶剂可接触面积：Q41 侧链为 0.0 Å²，Q62 为 100.67 Å²。这为“埋藏侧链与面向溶剂侧链”提供同一输入结构上的可复算依据；不是 MD 平均值。

`frame_geometry_audit.csv` 核查每个已打包真实帧的三个 χ 角与原始轨迹归档角度，6 个代表结构全部通过；最大差异 0.0849°，来自 PDB 坐标保留三位小数。可用 MDAnalysis 运行 `validate_representative_frames.py` 重做这项检查。`pse_scene_audit.json` 核查了保存 PSE 的五个场景均显示完整的 119 残基主链、启用正确对象并恢复对应相机，不依赖外部结构路径。

包中原始输入、每帧结构和占据率均可审计；`provenance.json` 包含原始 PDB 与 XTC 的 SHA256 校验、来源路径和处理规则。重新抽帧才需要原轨迹和 MDAnalysis，命令为工作区脚本 `scripts/case_validation_20261009/render_ubiquitin.py --extract`；日常打开、调整角度或重画只需要本文件夹。

当前数值使用原稿v13的 chemical_identity_corrected_v1 完整模型、种子7、fold0、第40轮留出预测，泛素A链为72个RMSF及66个侧链熵有效位置；全部A/B两链119位的标签、掩码与预测均保存于data源表。`v13_ubiquitin_refresh_report.json`记录来源与训练折类型均值。Q41/Q62的真实MD联合朝向、窗口熵及所选代表帧保持一致；RMSF采用v13处理标签，源表保留独立重算值与实际差异。
