# V91 定量结构映射资产

最终推荐拼版使用 `1MXE_Full_RMSF_focus.png` 与 `1MXE_Full_Hchi_focus.png`。二者是同一刚体相机、同一像素比例和同一裁剪矩形的原生透明 PNG，尺寸 1638×841 px（比例 1.9477）；`*_focus_white_preview.png` 是按标准 alpha 在白底合成的审查图，适合查看实际白底效果。原始完整画幅 2600×1600 与 PSE 均保留。透明 PNG 在部分查看器中会忽略部分透明度而显得 CaM 很重，正式 matplotlib/PIL 合成应保留 alpha。

肽链全长 294–318 的 25 个位点全部显示；RMSF 采用 `pred_rmsf`，范围 0–2.3 Å，映射至主链 cartoon；只有 W303/M316 两侧链用细灰 sticks 辅助定位。Hχ 采用 `pred_chi_entropy`，范围 0–4.3 nats，全部 21 个有效侧链以 sticks 显示，A299/A306/A309/A311 的 Hχ 未定义，保持中性灰而非伪造零值。两套标度均为原生 matplotlib cividis，色标由 root 合图器绘制。Hχ 连接主链为 #9DA4AA，透明度 0.35；CaM 为淡灰 cartoon，透明度 0.84，不参与数值着色。这里 Full/Intra 是模型身份，本资产只显示 Full，绝无游离/结合状态转换含义。

选择沿肽长轴滚转 90° 的相机，是因为肽链水平、N 端在左/C 端在右，W303 芳环和周围侧链分散度较好。相机旋转没有改变坐标；没有把结构拉直或人为展开。局部裁剪由所有 211 个真实靶肽重原子的二维投影定义并加至少 70 px 边距，扩充高度至约 1.95 横纵比；没有裁掉任何靶肽。上下的远端灰色 CaM 被有意裁去，因此本资产是完整靶肽及邻近 CaM 支架的局部展示，不是完整复合物全景。

`render_manifest.json` 包含最终两个资产、所有 25 位点的 CA 和侧链重原子几何中心在原图/裁剪图中的像素坐标、原始三维坐标、相机、crop、色标范围。外部标签可直接读 `delivery_assets["Hchi"]["focus_named_landmarks"]["K296"]["sidechain_center"]` 等；RMSF 标签通常连接 `CA`。W303、M316、K296、N297、N308 以及所有其余位点均可查。

独立 `geometry_data_color_audit.json` 检查了源 CSV 25/21 个值与 manifest 一致、PSE 全部原子名/原生编号与源 PDB 一致、保存坐标仅有从三位小数 PDB 读入浮点数的舍入误差（<5×10⁻⁶ Å）且渲染前后加载坐标完全相同。Hχ 逐侧链原子核对颜色；RMSF 逐残基 N 原子核对标量基础颜色，renderer 同时把相同颜色显式赋给该残基 cartoon。未进行新模型预测或 MD。

复绘顺序：`render_quantitative_structure.py`（D:/Anaconda/envs/pvb_cpu/python.exe；环境 V91_PREVIEW=0, V91_ANGLES=90），然后 `verify_and_finalize.py`，`make_shared_focus_crop.py`，`add_named_landmarks.py`。默认不设置 V91_ANGLES 会输出四个试验相机；试验资产只是视觉选择记录，最终使用 90°。只需使用两个 focus PNG 和对应 PSE 即可继续拼版，勿把试验图与最终图混用。

视觉审查：已实际查看两张最终白底 crop。Hχ 的中灰主链现在能够把 21 个彩色侧链连接成连续调控肽；RMSF 核心整体呈深色来自其真实连续低值，保留固定色标，不对微小差异人为拉伸。图上不加机制箭头、活性图标、假构象，也没有用 Hχ 球体大小伪装原子重要性。