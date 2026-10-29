PLANET-MD 的公开 full 与 seq 权重已经在隔离 `planet_env` 中各完成两条冻结链的真实 CPU 预测，输入为 16pk_A（415 个残基）和 1a62_A（130 个残基），全部四个输出均为 complete。审计核验输出序列和长度、全部数值有限、预测 dtype 为 float32、NPZ/输入哈希，以及 `rmsf_angstrom = 10 × rmsf_native_nm` 逐元素完全一致。full 的预测时间为 10.093 和 2.740 秒，seq 为 7.485 和 1.694 秒；这些是模型已加载后的单链预测时间，不包含权重加载，也不是完整 3,063 条记录的运行时间或性能评测结果。

派发前冻结清单中全部源码、离线缓存和公开权重共 3,666,680,986 bytes 的逐文件 SHA256、payload.tar SHA256、cohort 四份主要文件以及全部 4,973 条 chain PDB 的 SHA256 均匹配；完成 smoke 后再次核验源码、缓存及两组权重没有变化。公开 full checkpoint SHA256 为 `b31a5eb353e773c346a413fbca7e63735a6ab53249a628073a7c01ceee393bf9`，seq 为 `54f8ee6df8ed0c45e33995431b720cb7676b38d7521a5101602fae58f574e243`，ESM3 权重为 `5ead5a135c658068db6a4f1b933e72d6110992c4668822e1c0e2dcc53e38acd9`，结构 encoder 权重为 `467acbaee703ba3ccde6e75241a912a316952e5ff071355f85c1d33c68704f40`。源 commit 为 `917999482343c6400bef51180c817a8af01ee871`，ESM3 snapshot 为 `47f0545b2b6daf26a93439a3cd610f4f7f3d5478`，rocketshp snapshot 为 `7831d318182d4f4fc1acbdafdbf23862d8bc1657`。

实际 checkpoint 中 full 的 `struct_features=True`、`struct_stage='encoded'`；seq 的 `struct_features=False`，结构阶段元数据也为 encoded，但 seq 不加载结构 encoder。适配器读取权重中的结构阶段，正确避开默认 quantized 设置。两组 checkpoint 的所有浮点张量为 torch.float32，默认温度状态值均为 300 K。适配器先在 CPU 加载 ESM3，再转至所选设备；此次选择 cpu，未改变 dtype、权重或输入温度，不训练、不校准。单位来自冻结 `trajectory.py` 对 MDTraj `md.rmsf` 的直接调用，以及预处理/数据读取没有另行乘倍率；故 native nm 转为 Å 的十倍换算成立。这个推理环境不安装 MDTraj，因为适配器的推理链不调用轨迹预处理；smoke 实际成功验证了这一点。

此次始终设置 `CUDA_VISIBLE_DEVICES` 为空、`HF_HOME=/root/flexppi_completion_20261010/baselines/hf_cache`、`HF_HUB_OFFLINE=1`，OMP/MKL/OPENBLAS 及 Torch 线程数均为 2。顺序执行 full，再执行 seq，没有 GPU Context、训练或并发 GPU 操作。报告为同目录 `payload_audit.json` 和 `smoke_audit.json`，标记为 `SMOKE_PASSED`；两组输出与完整日志已同步到本地相邻的 `planet_full_cpu_smoke_20261010`、`planet_seq_cpu_smoke_20261010` 文件夹。这只能证明离线原权重和输入适配器能够实际预测，不能用两条链代替完整评测。该基线按独立链输入，未编码相互作用伙伴；公开权重可能涉及 ATLAS/mdCATH，不能当作家族去重后从头训练的模型；SHP 不是本文联合离散 χ 状态熵，因此不拿 SHP 当 Hχ 预测。

可重复命令如下，两个命令需依次执行，输出目录使用新目录以免混入既有结果。这里 cpu 是显式指定的设备；wrapper 仅设 Torch 两线程，然后执行未修改的 root 适配器。

```bash
export CUDA_VISIBLE_DEVICES=
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export HF_HOME=/root/flexppi_completion_20261010/baselines/hf_cache
export HF_HUB_OFFLINE=1
PY=/root/flexppi_completion_20261010/planet_env/bin/python
BASE=/root/flexppi_completion_20261010/baselines
"$PY" "$BASE/run_planet_cpu_smoke.py" --cohort "$BASE/cohort" --source-root "$BASE/PLANET-MD" --weights "$BASE/weights/rocketshp_pretrained_20250505_v1.ckpt" --variant full --output "$BASE/planet_full_cpu_smoke_20261010" --device cpu --limit 2
"$PY" "$BASE/run_planet_cpu_smoke.py" --cohort "$BASE/cohort" --source-root "$BASE/PLANET-MD" --weights "$BASE/weights/rocketshp_pretrained_20250507_v1_seq.ckpt" --variant seq --output "$BASE/planet_seq_cpu_smoke_20261010" --device cpu --limit 2
```
