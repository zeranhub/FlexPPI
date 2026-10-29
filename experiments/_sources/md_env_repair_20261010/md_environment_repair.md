# KIX 新批次 MD 隔离环境修复审计（2026-10-10）

环境已经就绪。实际解释器为 `/root/flexppi_completion_20261010/md_conda/bin/python`，就绪标记为 `/root/flexppi_completion_20261010/MD_ENV_READY`。`md_environment_audit.json` 记录了实际导入版本、主机 glibc、驱动文本、每份力场文件哈希、官方包哈希、插件载入结果及明确的版本差异。这里仅执行了 CPU 安装、Python 导入、平台名称枚举、合成数值自检和 CLI 帮助检查，未创建 GPU Context，未运行烟测或 MD。根代理统一负责后续 GPU 烟测及生产调度。

原 pip 安装失败来自 ABI 不兼容：主机 glibc 为 2.31，而 [OpenMM 8.4.0.post2 官方 PyPI CPython 3.11 Linux 包](https://pypi.org/project/OpenMM/8.4.0.post2/#files) 为 `manylinux_2_34_x86_64`，要求 glibc ≥2.34。官方包 SHA256 为 `12ffcd82d596bded1382e30af55907754ee481aafc0fc4a921a97de2ea7a8c55`。更换 pip 镜像不能修复这一 ABI 限制，因此保留根代理已经增加镜像、180 秒超时与 `PIP_CONFIG_FILE=/dev/null` 的原安装脚本，另外创建兼容环境。

采用 [conda-forge 官方 8.4.0 发布元数据](https://api.anaconda.org/release/conda-forge/openmm/8.4.0) 中的 `linux-64/openmm-8.4.0-py311h1c77f12_0.conda`。该构建声明 glibc ≥2.17，实际下载包 SHA256 与官方元数据的 `8542b87f2b86a567322db41a2f28616f9a57227b0e7a8f98326a8f6cb9417b7f` 完全一致。实际 OpenMM distribution 版本为 8.4.0，runtime full_version 为 `8.4.0.dev-4768436`，与历史 PyPI 8.4.0.post2 的发行方式不同；这里没有宣称二进制逐字节一致。所有九份 CHARMM36 XML 与历史 `prepared_READY.json` 的 SHA256 逐一完全一致，故力场文件没有换成 CHARMM36m，也没有改动温度、溶剂、盐浓度、截断、时间步长、约束、系综和 10+100 ns 物理参数。

实际导入版本为 Python 3.11.17、NumPy 2.2.6、SciPy 1.15.3、MDAnalysis 2.9.0、PDBFixer 1.12.0、OpenMM 8.4.0。PDBFixer 从 [官方 v1.12 tag](https://github.com/openmm/pdbfixer/releases/tag/v1.12) 对应固定 commit `94cfa4c0ca551cdc5f13320f9a658efd59f2b881` 的源码包安装，源码包 SHA256 为 `eade9cc4abf6ebb44bde64b25e9dca3fdde71836876b6c337a0bfc905177bb6c`。核心版本已经实际导入验证，`pip check` 无破损依赖。conda 求解初次安装的 SciPy 1.17.1 随后由精确的 PyPI wheel 替换为 1.15.3，故 `md_environment_conda_explicit.txt` 仅代表 conda 基层；重建时必须继续执行 `finalize_md_conda.sh` 的 pip 覆盖层，不能仅凭 conda explicit 文件宣称实际 SciPy 为 1.17.1。`md_runtime_versions.txt` 记录实际 distribution 版本，`md_wheel_manifest.json` 保留下载轮子的哈希。

CUDA runtime/NVRTC 固定为 12.6/12.6.85，避免求解器默认采用 12.9。已枚举的平台为 Reference、CPU、CUDA。四个 AMD HIP 插件因缺少 `libhiprtc.so.6` 未载入，错误完整保留在审计 JSON；本批次使用 NVIDIA CUDA，不安装 AMD 驱动来掩盖记录。主机 NVIDIA 驱动为 550.78；[NVIDIA 官方小版本兼容文档](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html) 给出 CUDA 12.x 的驱动兼容范围，同时明确 PTX 动态编译有附加限制。因此 CUDA 插件注册成功不能替代实际四状态 CUDA 烟测，烟测若报 PTX 版本错误，应由根代理选择任务内兼容方案并再次记录，不能修改宿主驱动或省略门禁。

此次安装只写入新任务 `md_conda`、`md_conda_pkgs`、`md_conda_envs`、`md_install_sources`、`md_env_repair_20261010` 及对应审计/日志文件，没有向旧 base 安装包，也没有修改旧项目资产或根代理原安装器。初次未固定 CUDA 的新环境下载在事务提交前中止，改用下列固定命令完成环境创建，已完成的安装证据为 `md_conda_install_cuda126.log` 和 `md_environment_finalize.log`。

```bash
CONDA_PKGS_DIRS=/root/flexppi_completion_20261010/md_conda_pkgs \
CONDA_ENVS_PATH=/root/flexppi_completion_20261010/md_conda_envs \
/root/flexppi_completion_20261010/python311/bin/conda create \
  -p /root/flexppi_completion_20261010/md_conda \
  --override-channels \
  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge \
  --strict-channel-priority --solver libmamba \
  python=3.11 openmm=8.4.0=py311h1c77f12_0 numpy=2.2.6 \
  cuda-version=12.6 cuda-nvrtc=12.6.85 pip -y
bash /root/flexppi_completion_20261010/md_env_repair_20261010/finalize_md_conda.sh
```

# CLI 与阶段门禁复核

新批次先执行 `run_md.py prepare --root "$BATCH"`，然后执行 `audit_rebuilt_prepared.py --root "$BATCH"`，再由根代理执行 `run_md.py smoke --root "$BATCH"`。重建审计要求四个原始输入哈希与历史相同、保留重原子只发生统一平移和 PDB 舍入（残差 <0.002 Å）、体系电荷中和、五个 KIX 组氨酸为指定 HID、周期重原子间隙不少于 2 nm、力场哈希和 system/protein 文件哈希匹配，以及 1.2/1.0 nm 截断/切换。烟测覆盖四个状态各 40 ps，检查温度、DCD 40 帧/1 ps/原子数、检查点回读和几何，写入 `smoke/PASSED.json`，但不能据此宣称平衡或收敛。

生产必须使用 `queue_triplicate.py`，因为它会同时核验 `protocol.json`、prepared 审计及哈希、Full/Static 冻结预测及输入哈希、全部 KIX 87 个位点、所有有效侧链位点的有限值，以及 smoke 标记。`run_md.py queue` 是旧调度入口，只检查 smoke，不能作为新批次生产门禁。`queue_triplicate.py --verify-only` 可以在生产前检查全部上述门禁，不启动 GPU。以下变量 BATCH 表示批次目录、PY 表示已验证解释器。三个质量命令需在对应四状态生产及 per-repeat 分析完成后分别执行，保留每次独立输出；不要使用跳过 DCD 哈希或跳过原生接触的选项。默认质量输出目录必须为空，失败后的重跑请指定新的、可审计的 output-dir，不能删除已有结果后覆盖。

```bash
BATCH=/root/flexppi_completion_20261010/kix_md_triplicate_rebuild_20261010
PY=/root/flexppi_completion_20261010/md_conda/bin/python
"$PY" "$BATCH/code/queue_triplicate.py" --root "$BATCH" --verify-only
# 生产由根代理统一指定 CUDA_VISIBLE_DEVICES 后启动：
"$PY" "$BATCH/code/queue_triplicate.py" --root "$BATCH" --states apo myb mll both --repeats 1 2 3
"$PY" "$BATCH/code/quality_one_repeat.py" --root "$BATCH" --repeat rep01
"$PY" "$BATCH/code/quality_one_repeat.py" --root "$BATCH" --repeat rep02
"$PY" "$BATCH/code/quality_one_repeat.py" --root "$BATCH" --repeat rep03
"$PY" "$BATCH/code/analyze_triplicate_delta.py" --root "$BATCH"
```

生产每个 1 ns 段用 JSON 原子写入作为提交标记，未提交文件在续跑时保留到 interrupted 目录，检查点验证后继续；完整单次运行要求 10 个平衡段、100 个生产段、55,000,000 步，per-repeat 指标要求 100,000 个生产帧。新队列只有这些条件与指标 CSV 哈希通过才记录 job_completed，全部请求的作业通过才记录 worker_completed。worker_completed 明确仍要求质量与重复性分析。最终 delta 脚本要求全部 12 个新批次作业及各段 DCD/检查点/拓扑/指标哈希一致，历史 rep01 不合并。delta 入口本身不检查三份 quality/repXX/summary.json，因此根调度完成条件应另行要求三份质量报告完成并审阅；质量输出中的 passed 表示数据与全帧复算通过，并不是统计收敛证书。整个分析不使用 CI/p 值，不把累积窗口、全部位点或九个跨状态重复组合视为独立重复。

`analyze_job.py --self-test` 与 `quality_one_repeat.py --self-test` 已在实际新环境通过，包含角度边界、联合离散熵、在线 RMSF、合成 DCD 回读/分块、一致的 chi 定义哈希及质量窗口定义；它们仅是合成数值夹具。run_md、queue_triplicate、audit_rebuilt_prepared 和 analyze_triplicate_delta 的 CLI help 已全部通过实际解释器，完整输出见 finalize 日志。未生成或冒充新的 MD 科学结果。

# CPU prepared 重建结果

根代理随后授权启动纯 CPU prepare。实际执行时该 SSH 进程的 `CUDA_VISIBLE_DEVICES` 为空、`OPENMM_CPU_THREADS=4`，脚本的 PDBFixer 与 addHydrogens 都显式选 CPU。四状态 prepare 和 audit_rebuilt_prepared 均退出 0，`prepared/READY.json` 与 `prepared/input_audit.json` 已写入；审计状态为 `all_four_rebuilt_inputs_passed`。保留原始重原子在扣除装盒的统一平移之后，最大坐标差异分别为 apo 6.19×10⁻¹³ Å、myb 6.82×10⁻¹³ Å、mll 6.26×10⁻¹³ Å、both 4.53×10⁻¹³ Å，在 PDB 写出精度下完全一致。apo/myb 没有新增重原子；mll/both 仅增加共同修复的 C841 CG/OD1/OD2 和 C857 OXT。原始 ternary 共同重原子修复记录的既有原子位移为 0 nm。

四个溶剂化体系总原子数分别为 45,816、45,803、45,057、45,931；全部净电荷误差小于 10⁻¹¹ e，HID、周期重原子间隔、1.2/1.0 nm 截断/切换和历史力场哈希审计通过。prepared 整个目录（system.xml、solvated.pdb、protein.pdb、metadata、READY 与审计）及两份 CPU 日志已同步到本地新批次目录。初始 PAYLOAD_MANIFEST 是派发前快照，新的阶段状态另记于 PREPARED_HANDOFF.json，不用初始 prepared=false 状态描述目前阶段。此次仍未创建 GPU Context；GPU smoke 与生产完成状态由根代理随后记录。
