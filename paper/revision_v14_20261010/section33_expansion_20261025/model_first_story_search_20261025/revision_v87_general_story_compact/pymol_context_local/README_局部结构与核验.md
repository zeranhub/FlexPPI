# v87：1A06 自抑制结构的 M316 局部

主图直接使用 `M316_contacts_local_trim.png`（889 × 953，透明背景）。图内没有文字；M316、L45、S99 的裁边图像像素坐标保存在 `render_manifest.json` 的 `assets.local.trim.landmarks` 中。建议配对 v86 的 `autoinhibited_trim.png` 整体图。两者相机旋转的九个元素完全相同，仅改变相机位置、放大倍率和裁边，没有重新旋转局部，也没有改变任何原子坐标。

调控段 312–316 为紫色 cartoon；M316 为青色棒；L45 和 S99 为灰色棒。为使 24–28 mm 的局部保持可读，激酶中与这三个侧链无关的长 cartoon 片段不显示。两条虚线仅表示真实原子对之间的距离：M316 SD–L45 CD1 为 4.143 Å，M316 CE–S99 OG 为 3.161 Å；不将它们指定为氢键。此图是 1A06 静态结构的接触依据，没有在 1A06 上运行新的模型预测或 MD，也没有将它描述为实测的转换途径。

输入为 v86 已保存的 `1A06_rigidly_aligned_to_1MXE.pdb`，原先只经过整分子的刚体叠合。该 PDB 保留 2221 个已观测原子；本批渲染前后坐标完全一致。其坐标按 PDB 格式保留到 0.001 Å，因此重新计算的距离与 v86 内存坐标的末尾小数略有不同，显示到 0.01 Å 时相同。1A06 中未解析的 317–318 没有补建。

复现：`D:/Anaconda/envs/pvb_cpu/python.exe render_M316_local.py`。脚本只依赖同目录的 inputs/ 数据，自动设置 PyMOL 的 DLL 目录，输出透明 PNG、裁边 PNG、PSE 及 manifest，不依赖临时补丁或另一个文件夹中的脚本。