# 7XQ8：hmu-gpu 单条 100 ns 膜 MD

用户于 2026-10-29 明确选择先做 1 次 × 100 ns。运行目录为 `hmu-gpu:/ranze/runs/flexppi_7xq8_md_20261029`，SSH 端口为 30352。登录凭据仅在交互式会话中使用，不保存到本目录。服务器日志使用 UTC，换算北京时间需加 8 小时。

当前已完成原子修复，已启动脱离 SSH 会话的 `launch_pipeline.sh`，依次建膜、最小化、分级平衡、生产采样。初始建膜任务 PID 2300 在 OpenMM 默认缩放插膜阶段出现 NaN 并停止；改用完整尺寸删重叠脂质后重启的 PID 为 2542。初始与重试记录见 `audits/pipeline_launch_20261029.json`、`audits/pipeline_fullsize_launch_20261029.json`；实时状态以服务器 `pipeline.log`、`replica_1/progress.json`、`production_entry_gate.json` 为准。启动不表示已产出 100 ns 数据；只有 `replica_1/complete.json` 存在且轨迹核验通过后才能报告生产完成。

体系是 7XQ8 已解析的六条蛋白链 C、L、v、R、A、B，共 1,889 个残基，保留 author 编号与插入码。没有补建未解析末端、表达信号肽、标签或胞内长尾；去除全部 NAG，使用无糖 POPC 简化膜。六链截断端采用带电 NH3+/COO−，这一假设在 Igβ 跨膜段附近尤其需要注意。本体系用于探索性柔性参考，不能支持天然膜、糖调节或未解析胞内信号尾部的机制结论。

原子修复使用 PDBFixer，只补缺失原子、末端 OXT 和氢；原始重原子在初次补全时保持不动，之后做有骨架约束的局部几何最小化。最终蛋白含 29,126 个原子，全部 23 对二硫键长为 2.019–2.052 Å，最大 Cα 位移为 1.856 Å，逐残基位移保存在 `audits/protein_preparation_audit.json`。原 PDB 的 11 条 SSBOND 显式进入拓扑；缺失 Fab 二硫键参考 3NGB 抗体部分，额外可变区二硫键通过序列对应 C42–C108/v42–v108；IgM 重轻链连接依据 UniProt P01871 的恒定区 C14 与轻链连接注释，精确映射到 C137/v137，并连接 L216/R216。UniProt 恒定区 C28–C88 对应 7XQ8 C151–C211，其他恒定区 Cys 也逐一映射，不能用 IgG 的重链 C216 编号替代 IgM 的 C137。重建连接属于有参考证据的模型处理，不是新增的原实验 SSBOND 记录。

参数使用 OpenMM 8.4、CUDA mixed、CHARMM36m（`charmm36_2024.xml`）及其配套 `charmm36_2024/water.xml`，POPC 膜、310 K、0.15 M NaCl，并额外中和净电荷；四跨膜段共同几何轴对齐 z，胞外域朝负 z。`Modeller.addMembrane` 同时建膜和水，最小边缘 padding 为 1.5 nm，实际盒长、原子数和离子数在建膜后记录。最终体系用 PME、1.2 nm 截断、1.0 nm LJ switch、HBonds 约束；使用 LangevinMiddle（摩擦 1/ps）和膜专用 MonteCarloMembraneBarostat（1 bar、表面张力 0、XY 同比例变化、Z 独立变化）。这些是本项目方案，不是相关 2025 年论文的逐项复现。

默认的半尺寸 InflateGRO 类插膜在 OpenMM 内部膜放松时产生 NaN。本项目仅在运行进程内替换该方法为完整尺寸蛋白删重叠脂质，阈值 0.30 nm，保持上下层删除数相等；初始最小化最多 2000 步，建膜内部低温放松步长改为 0.5 fs。保留原算法的膜贴片复制、加水及加离子逻辑；生成的完整方法源代码和 SHA256 保存在服务器 `addMembrane_fullsize_implementation.py` 与 `membrane_insertion_method.json`，不修改已安装 OpenMM 库。这个内部建膜放松不是正式平衡，更不计入生产时间。

流程首先全体系有约束最小化，然后 3 × 125 ps NVT（1 fs），3 × 500 ps NPT（2 fs），再 10 ns 无约束预平衡，最后单条 100 ns 生产期（2 fs）。骨架/侧链约束依次为 4000/2000、2000/1000、1000/500、500/200、200/50、50/0，单位 kJ/mol/nm²，预平衡与生产为零；NVT/NPT 初期对脂质 P 头基也逐步释放位置约束。生产期之前检查最后 5 ns 的平均温度、体积与膜面积波动/漂移；失败则停止，不自动将异常平衡计为正式数据。这些检查不能证明热力学收敛。

每 10 ps 保存蛋白 DCD、能量/温度/体积，生产期另每 100 ps 保存全体系 DCD；每 50 ps 保存交替二进制 checkpoint 和原子提交的 progress，阶段末保存 XML state。重启时新建 DCD 分段，物理起始时间写入 DCD 和进度 manifest；不能忽略 manifest 直接拼接重启产生的潜在重叠尾帧。生产前至少有 11.875 ns 的独立平衡时间，不计入 100 ns。后续分析应使用生产段并按链、author 编号、插入码、氨基酸匹配，不能凭行号连接旧预测。

服务器的系统 `libcuda.so.1` 虽能让 nvidia-smi 显示显卡，却无法初始化计算。实测使用进程独立的 `LD_LIBRARY_PATH=/usr/local/cuda/compat/lib` 后，OpenMM CUDA 小体系能量计算通过。启动脚本仅为本任务进程设置该路径，没有修改系统驱动或其他任务环境。

运行或断点续跑使用服务器目录中的 `bash launch_pipeline.sh`；脚本用 flock 防止重复启动。已有后台进程活跃时不要重复启动。作业失败需先读 `pipeline.log` 和最后 checkpoint，再决定修复；不要删除现有轨迹。只有完整轨迹产生后，才能做 RMSD、RMSF、Hχ、膜面积与厚度、前后半程采样检查。单条轨迹无法估计独立重复间变异。

参考：[7XQ8](https://www.rcsb.org/structure/7XQ8)、[3NGB](https://www.rcsb.org/structure/3NGB)、[UniProt P01871](https://www.uniprot.org/uniprotkb/P01871/entry)、[OpenMM 力场与运行文档](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)。相关 BCR 文献的方法原文已保存在 `reference/published_methods.txt`；其使用另一 Fab 构建体和复杂膜，不能把本次简化模拟称为复现。
