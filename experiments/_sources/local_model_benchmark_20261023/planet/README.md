# PLANET-MD 本机三个版本横评

`run_local.py` 已完成语法检查、三份 checkpoint 的无张量元数据检查，以及三条烟测输入的重原子身份和坐标静态检查。当前尚未加载 ESM3，也没有把准备好的脚本当成已完成的模型结果。

full 和 seq checkpoint 均为 `struct_stage=encoded`、`struct_dim=1024`；full 使用结构特征，seq 不使用。官方 mini checkpoint 也是**序列模型**，使用更小的下游预测网络。三份 checkpoint 的默认温度均实际读出为 300 K。完整检查记录在 `checkpoint_metadata_static.json`、`inspect_protocol.json`、`input_static_validation.json`。

## 执行命令及成功条件

```powershell
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' -u 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/planet/run_local.py' --mode inspect
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' -u 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/planet/run_local.py' --mode smoke
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe' -u 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/planet/run_local.py' --mode full
```

`inspect` 不导入 torch、不分配模型张量，也不运行推理。`smoke` 对 ATLAS、其他多链、TCR–pMHC 各一条冻结链运行三个 checkpoint，总共验证九次。每个版本的 RMSF 路径都必须与**同一官方模型的完整 forward 输出逐元素完全相等**。ATLAS 烟测还将 seq 结果与先前官方 CPU 输出比较，以独立检查 ESM3 的节省内存加载方式；序列完全相同、权重哈希相同后才比较，容许不同 CPU 数值库的微小浮点误差。

烟测成功要求进程退出码为 0，并且 `smoke_audit.json` 中 `status=complete`、`chains_requested=3`、`official_forward_equivalence_checks=9`、`sources_verified=[atlas,dynarepo,itrace]`、`heads_verified_against_official_forward=[full,seq,mini]`。full 模式会强制检查这些条件，以及权重、源码、脚本、队列和运行环境与烟测完全对应；条件不满足时不能启动全量。任何链或预测头失败时整轮退出码为 2；其他致命异常退出码为 1。不会用失败子集冒充全队列成功。

## 保持原算法的资源处理

CPU 推理始终使用原生 float32、两线程、300 K，原始 nm 输出乘以 10 转为 Å，不裁剪、不拟合、不训练。4,973 链含 3,856 条不同序列，共 770,439 个不同序列残基，因此按序列哈希缓存一次完整精度 ESM3 特征，再交给三个官方预测头。每条链仍独立计算自己的结构特征和残基映射，序列相同不会让结构结果互相替代。

模型使用官方构造函数和原始权重，先在 meta 设备建立形状，再以严格匹配的 `load_state_dict(assign=True)` 关联 mmap 只读权重，避免初始化权重与加载权重同时占用两份内存。只调用官方公式重建不需要训练的 rotary 位置常数。所有模型参数检查为 CPU float32。

官方完整 forward 在 RMSF 之外还构建两两残基输出；最长链的额外中间数组可能超过 2 GiB。此脚本原样执行 `_transform`、附加默认温度和 `rmsf_head`，省去与 RMSF 无依赖关系的其他输出，并用烟测与完整官方 forward 的逐元素比较作为全量前置条件。SHP 并非本论文的侧链联合状态熵 Hχ，不能用它代替 Hχ。

权重稳态约 3 GiB，预计整体进程峰值约 4–5 GiB，尚待烟测实际测量。脚本默认在启动模型前等待至少 6 GiB 可用物理内存，不足时每 30 秒输出一次等待状态。等待上限默认为 600 秒，可用 `--memory-wait-seconds 600` 明确设置；超时会将 `memory_status.json` 和对应 `{mode}_audit.json` 写为 `blocked_insufficient_RAM_before_model_load`，退出码为 **3**，不会无限挂起队列。未到门槛时不会加载模型。所有其他推理退出后，可以仅将启动门槛改为 `--minimum-free-gib 5` 做有界烟测，这不改变任何模型计算；4.5 GiB 门槛的余量偏紧，静态估计不足以确认安全。

真正执行时，`smoke_resource_usage.json` / `full_resource_usage.json` 以约 0.1 秒间隔观察进程峰值 RSS 和系统最低可用内存，并每 5 秒落盘，可据此评估 full 的启动门槛。致命异常保留在 `fatal_error.json`。当前仅有早期 CPU 历史测量，全量粗略预计 4–10 小时；需本次两线程烟测后修正，不能把估计作为最终运行时间。

## 评估接口

`planet/{full,seq,mini}/{input_id}.npz` 保存 `rmsf_native_nm`、`rmsf_angstrom`、`sequence`、`input_id`、`input_index`、原始 `chain`、`author_residue_id`、`amino_acid`、`record_residue_index`。推理脚本不读取标签，评估端通过冻结共同真值与有效位置掩码评分。每条链都校验规范化链 PDB 哈希，并与完整结构中的原始链核对所有重原子、坐标和残基顺序。

同名 JSON 的 `status=complete` 和预测文件 SHA 是权威成功标记。`full_progress.json` 提供当前进度，`full_audit.json` 只在全部条目尝试完毕后生成。ESM 共享开销、结构特征时间、单独 RMSF 预测头时间及烟测额外验证时间分开保存，横评运行时间时应说明共享特征缓存与线程数。
