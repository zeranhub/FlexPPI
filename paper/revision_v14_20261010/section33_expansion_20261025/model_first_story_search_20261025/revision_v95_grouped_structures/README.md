# v95 结构集中与局部精简

两个完整结构并列于顶行 a/b；删除 Q41/Q62 局部放大。c 是全部25位点 RMSF，d 是泛素完整侧链谱，e 是全部六个 Gln，f 是逐 χ 的2×2关联，g 是83记录分析。878个绘图数值与 v94 完全相同，仅重排和修改子图引用。S4未改。

figure3/compose_fig3_v95.py 为可运行绘图脚本；source_data 为全部面板源表。pymol/ 含两个原始整体结构渲染、可编辑 PSE、精确结构输入及相机记录。此包不含完整原始 MD 轨迹或训练权重。prepare_v95.py 是从既有 v94 初始化的脚本，复现本包成图可直接运行 compose_fig3_v95.py，无需再次初始化。build_v95_document.py 需同目录体系中的 v94 Word 作为保真输入，使用 Codex bundled Python。

科学绘图运行环境：E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe。audit_v95.py 使用 bundled Python（包含 pandas 和 lxml）。

Word 已渲染29页，6张改变页已逐页打开检查，其余23页与已验证 v94 的渲染逐字节相同。主图与图注位于第11页，无裁切或重叠。
