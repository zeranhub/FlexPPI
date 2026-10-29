# FastProtFlex 生命周期降峰候选：尚未运行 R

此目录只准备独立的候选文件和验证工具，当前状态为 `prepared_static_only_R_not_executed`。官方源码、正在运行的 Python/R 包装器、旧预测和原有 RAM 门槛均未修改；没有启动第二个 R 进程。Python 静态检查已通过，不能把这当作实际数值等价或峰值降低已经通过。

`original/FUNCTION_GDV.r` 是官方文件的原字节副本，SHA256 为 `c3efdd8b31b65d33ea13ff70b6f777ca394f5566a19aa947d05c865a259a9b60`；`candidate/FUNCTION_GDV.r` 的 SHA256 为 `5f7327911f1bcbf2b3c3c737c703a7970faeecda072c5b9d591b4642965f8c34`。完整逐行差异见 `FUNCTION_GDV_lifecycle_only.diff`，原始和候选哈希、源文件行号、原 RDS 哈希及验收条件见 `lifecycle_provenance.json`。静态工具确认全部 11 条增删行都只涉及 `rm(...)` / `gc()`；去掉这些生命周期语句后，原程序和候选程序的其他文本及先后顺序完全相同，阈值、图构造、轨道计数、平滑、标准化和 15 个系数均未改动。

候选只调整 `GDVallatonce`：原第 14 行生成矩阵后，提前释放第一次的 `DD`；原第 16 行建图返回后释放 `adj`；原第 19 行轨道计数后释放 `g`；原第 22 行完成 2.1 Å 阈值操作后释放 `DDsave`；原第 23 行得到平滑邻接矩阵后释放第二次的 `DD`。每次释放后调用 `gc()`，并移除原第 18 行已被拆开的重复清理。原分块函数完全不变，冻结队列本身也没有超过 10,000 重原子的输入。最大 8,058 原子时，一个压缩距离对象约 247.66 MiB、完整邻接矩阵约 495.39 MiB；不同阶段的释放量不能相加当作整体峰值下降。首次矩阵或图构造仍可能决定全局峰值，R 分配器也可能保留已回收内存，因此实际能降多少必须测量。

验证工具默认只检查副本、差异和三条原始完整精度 RDS 的哈希，**没有 `--run` 不会启动 R**：

```powershell
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/memory_lifetime_v2/validate_candidate.py'
```

等 root 明确通知可以运行后，使用以下命令串行验证三个来源；此命令尚未执行。每次测量需要新的 `--run-tag`，不会覆盖之前的结果：

```powershell
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' -u 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/memory_lifetime_v2/validate_candidate.py' --run --run-tag smoke_01 --memory-wait-seconds 600 --case-timeout-seconds 180
```

顺序固定为 `2o70_F`（1,385 原子）、`A007V`（3,072 原子）、`itrace_1qrn`（6,620 原子）。只要原 FastProtFlex Python 包装器还在运行，或有任意 R 进程存在，工具就阻止启动；也完整保留原烟测推导的可用 RAM 门槛，不因预期降峰而下调。等待超过 600 秒会写 blocked 并退出 3，任何模型/数值/版本不一致退出 2，全部通过才退出 0。每条最多 180 秒，只有本工具自己启动的进程树会在超时后终止。不会终止现有其他任务或绕过它们的内存保护。

每条输入在独立工作目录运行。`predict.r` 和 `capture_full_precision.R` 都是原字节副本，保留原来的相对 `source(...)` 调用；由于当前目录是隔离的验证目录，它们只会读取该目录中的候选 `FUNCTION_GDV.r`。原输入 PDB 字节也与原成功记录的输入哈希核对。`capture_and_validate.R` 在模型运行后读取新产生的完整精度 `atom_predictions_exact.rds`，与原成功运行的 RDS 进行全部原子比较：原子表列及顺序一致，全部身份列一致，所有 x/y/z 坐标一致，完整原子表一致，全部原子预测向量 `identical(...,num.eq=FALSE)`，数值全部有限，最大绝对差恰好为 0。不是只比较 Cα，也不使用 PDB B-factor 的两位小数。任何条件失败都不允许采用候选。

Python 每约 0.1 秒采样一次 R 进程树 RSS、CPU 时间和线程数，保持数值库单线程、同一逻辑 CPU，并记录启动前可用内存及原门槛。R 版本和包版本还必须与参考运行一致。输出位于 `runs/<run-tag>/<sample_id>/`，每条包含 `comparison.json`、完整精度 RDS/CSV、日志和 `validation.json`；总结果为 `runs/<run-tag>/validation_audit.json`。测得的墙钟时间包括输出捕获与 RDS 比较，历史 baseline 和新候选不在同一时刻运行，报告明确这些限制。即便三个来源数值全部一致、峰值确实降低，也不会自动改动现有推理任务或内存门槛；应另行审核极值结构的资源验证并创建独立的新版本。

## 独立续跑入口：已准备，尚未启用

`../run_local_lifetime_v2.py` 是新的独立入口，原 `run_local.py`、原捕获包装器和官方文件均未改动。默认命令只检查全部 3,063 条调度的身份、分组和来源哈希，不加载 NumPy/Pandas，也不启动 R、归档或写入协议。此命令已实际运行通过；静态检查记录为 `resume_entrypoint_static_audit.json`，其状态明确不包含数值或内存验证：

```powershell
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/run_local_lifetime_v2.py'
```

以下续跑命令仅供 root 在三来源实际验证完成、审阅资源结果且原 Fast/R 已停止后执行，目前没有执行。`--validation-audit` 必须指向真实通过的三来源总审计及其仍在磁盘上的原子 RDS、逐条验证和比较文件；没有真实审计的 `--run` 会立即失败，也不会伪造校准、成功状态或启动 R：

```powershell
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' -u 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/run_local_lifetime_v2.py' --run --validation-audit 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/memory_lifetime_v2/runs/smoke_01/validation_audit.json' --memory-wait-seconds 600
```

续跑沿用原 `full_schedule.csv` 的完整 3,063 条、原顺序、原 `records/<sample_id>/prediction.npz` 路径和全部数组键。每条已完成的缓存都核对输入、来源脚本、全部已记录输出的 SHA256，以及原子身份验证标志和 NPZ 中残基身份、顺序、有效位。完整缓存只读复用，不重算、不改原预测或原 `validation.json`；缓存损坏会报错停止，不能以重算掩盖。仅缺失或未完成记录可以使用候选；新记录仍从整个结构输入并核对全部原子坐标和身份，最终只取 Cα 原子的相对分数。标签副本只按已记录的 CSV 末位解析容差核对，实际评分仍使用冻结原真值，不把相对分数当 Å，也不生成 Hχ。

内存门槛只来自实际测得的峰值：以“原子数的平方”为横轴拟合非负斜率，再抬高截距以覆盖每条观测，最后保留原 15% 安全量和 512 MiB 余量。脚本没有预填峰值。**对于大于实际测量最大原子数的任何缺失记录，一律写 `blocked_unvalidated_atom_count_requires_root_resource_review` 并退出 3；没有手填最大原子数或降低门槛的覆盖选项。** 没有下述真实极值资源扩展审计时，三来源最大仍只有 6,620 原子。范围内若 RAM 或既有 R 不满足，最多等 600 秒，写阻塞状态并退出 3；失败退出 2，全 3,063 条完整才退出 0。

只有真正启用续跑时，才逐字归档原 `run_protocol.json` 为 `run_protocol.official_v1.archived.<原SHA前16>.json`，然后写新协议。新协议用 `execution_entrypoint_path` 和 `execution_entrypoint_sha256` 绑定真实 `run_local_lifetime_v2.py`，`wrapper_path` / `wrapper_sha256` 为同值兼容字段；绝不把新 SHA 绑定到原 `run_local.py`。另明确记录 `official_algorithm_source_path/sha256`、`candidate_function_path/sha256`、`capture_wrapper_path/sha256`、`lifecycle_validation_audit_path/sha256`、`resource_gate_version`、完整 `resource_calibration` 以及原协议归档路径/SHA；缓存策略明确允许已验证 `official_v1` 结果与后续 `lifetime_v2` 结果共存。每条新成功验证也写入这些来源、真实验证审计 SHA 和资源门槛版本。汇总兼容原 `full_results.csv`、`full_audit.json`，另保留 `full_audit.lifetime_v2.json`；尚未采用时不会修改这些文件。控制器若将来采用此入口，应同时识别 `run_local_lifetime_v2.py` 和 `memory_lifetime_v2/validate_candidate.py` 的进程签名，避免原 runner 退出就提前开始下一阶段。

续跑读取真值时仅投影确实消费的七列 `sample_id, residue_index, amino_acid, chain, author_residue_id, target_rmsf, valid_rmsf`，仍保留全部 961,321 行和原默认浮点解析方式，不加载 Hχ 或训练元数据。小样本核查记录见 `truth_usecols_small_sample_audit.json`；这是表格读取兼容性检查，不是候选 R 模型的数值验证。

## 独立 8,058 原子资源探针：尚未运行

`probe_extreme_resources.py` 只处理冻结的 `7zh9_A`（8,058 重原子、1,024 残基），不接受替换体系或手填峰值。默认仍只读检查，不启动 R：

```powershell
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/memory_lifetime_v2/probe_extreme_resources.py'
```

仅在 root 明确通知、三来源真实逐值验证通过并停止原 Fast/R 后，才可运行下述尚未执行的命令。资源启动门槛使用三来源实测二次包络外推至 8,058 原子，保留 15%＋512 MiB；当下可用 RAM 必须满足。这里外推的是资源需要量，不是模型预测分数。内存不满足则等待，600 秒后退出 3；不降低门槛、不重试模型、不改分块规则。单次 R 最长 180 秒，超时仅处理本探针启动的进程树。

```powershell
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' -u 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/memory_lifetime_v2/probe_extreme_resources.py' --run --validation-audit 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/memory_lifetime_v2/runs/smoke_01/validation_audit.json' --run-tag extreme_01 --memory-wait-seconds 600 --case-timeout-seconds 180
```

探针只在 `resource_probes/<run-tag>/7zh9_A/` 保存输入/源副本、日志、原精度 CSV/RDS/PDB、独立 NPZ 和 `integrity.json`，不写原 `records`。全部原子身份、顺序、坐标和预测有限性，以及 Cα 与冻结真值身份及顺序均必须通过。这里没有该极值的原版完整精度参考，审计明确写 `original_full_precision_reference_available=false` 和 `original_equivalence_compared=false`，不能写成原版逐值一致。实际三来源逐值相同审计仍是单独硬要求。

只有实际成功才生成状态 `extreme_resource_and_integrity_check_complete_no_reference_equivalence` 的 `resource_probes/<run-tag>/resource_extension_audit.json`。它记录真实峰值、单 CPU、输入/源码/包装器/三源审计 SHA，以及全部新输出文件的大小和 SHA；若等待则也绑定等待文件。`base_three_source_resource_calibration` 完整保留原三来源校准及其自哈希，`gate` 绑定这一基准；`combined_resource_envelope` 合并三个原峰值和这一个实测极值，重新抬高包络，最大已测量原子数才成为 8,058。失败或阻塞审计没有采用资格，不会产生成功校准。

实际探针完成并经 root 审阅后，续跑入口可额外提供显式扩展参数。入口重新核验两份独立审计、全部探针输出哈希、源哈希和拟合结果；只提供扩展审计不能代替三来源精确等价审计。实际采用时，协议和每条新记录额外写 `resource_extension_audit_path/sha256` 与 `resource_extension_probe_path/sha256`，资源门槛版本变为 `lifetime_v2_measured_three_source_plus_extreme8058_quadratic_upper_v1`。没有完成探针，以下命令会拒绝，不会放宽资源门槛：

```powershell
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' -u 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/run_local_lifetime_v2.py' --run --validation-audit 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/memory_lifetime_v2/runs/smoke_01/validation_audit.json' --resource-extension-audit 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/fastprotflex/memory_lifetime_v2/resource_probes/extreme_01/resource_extension_audit.json' --memory-wait-seconds 600
```

准备阶段仅有静态检查记录 `extreme_probe_static_audit.json`，不存在真正的资源扩展成功审计、三来源候选成功审计或测得峰值。控制器将来还应等待 `probe_extreme_resources.py` 的进程签名；当前队列文件仍不修改。
