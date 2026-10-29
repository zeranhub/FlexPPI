# v13原稿新增结构实例与PyMOL复现

本图包对应恢复后的原稿 paper_v13_20261008.tex。原稿的19张图及5张表完整保留，新增a--h结构实例组合图在结果部分补充主链活动区域定位与侧链状态多样性的用途。数值统一采用原稿v13的chemical_identity_corrected_v1完整模型、种子7、第40轮：KIX使用固定fold1，泛素A003R使用fold0留出预测。`case_model_protocol.json`记录模型与缓存SHA。

`fig_local_dynamic_cases.pdf`为183×225mm组合图，SVG保存曲线和文字，PNG为300dpi预览。KIX四状态含全部87个位置；泛素曲线含A链全部72个RMSF和66个侧链熵有效位置。曲线按真实作者编号连接有效观测，MD空心圆标出数值位置，轴下短线标出无效侧链位置。`data/ubiquitin_fold0_all_residues.csv`保留A/B两链全部119位，作图时筛选A链。`compose_figure.py`从本包的完整数值与PyMOL图重新拼图，运行`python compose_figure.py`即可。

双击`kix/KIX_case_validation.pse`、`ubiquitin/ubiquitin_2QHO.pse`可以打开完整蛋白、配色和相机。f可直接打开`ubiquitin/f_Q41_local.pse`或`ubiquitin/f_Q62_local.pse`，同名PML恢复对应局部视图。根目录四个中文CMD调用本机PyMOL；复制到其他电脑时直接用自己的PyMOL打开PSE即可。两个主会话各含五个整体与局部场景，切换scene按钮或缩小视角可看完整蛋白。

局部制作顺序为完整蛋白表示、完整图ray渲染、按固定像素框裁切。KIX三条链保留131个Cα cartoon，泛素--UBA两条链保留119个。L12完整主链与侧链真实键连续显示；Q41/Q62每个实际MD代表帧保留该残基的9个重原子和8条真实键，彩色根部小球标出该帧的真实Cα，并叠加自己整条泛素72个Cα的灰色细迹线。所有源坐标保持不变。`*_whole_source.png`保存全景，`*_crop_box_preview.png`标出裁切范围，`crop_manifest.json`保存相机、像素框及校验值，局部PNG逐像素取自全景。

在PyMOL选择File→Run Script打开`kix/reproduce_kix.pml`或`ubiquitin/rebuild_relative.pml`可以从随包源文件重建，不依赖当前工作目录。重新ray渲染的命令见各案例中文README。Q41/Q62的联合朝向来自真实MD，分别取占据率最高三类的代表帧；图h包含全部27类，占据率覆盖90%需要七类和十二类。模型输出RMSF与Hχ数值，侧链帧号、角度和原PDB随包保存。

`参考图与画法说明.md`记录在本地文献库查看的原图及表示层次。`delivered_session_checks.json`核对全蛋白表示、相机和像素裁切；`connected_sticks_delivery_QA.json`核对真实键、每帧自身骨架和不变坐标；`portable_connected_sticks_QA.json`记录移至另一目录后实际PML重建的结果。请重新打开当前PSE以载入更新后的场景。
