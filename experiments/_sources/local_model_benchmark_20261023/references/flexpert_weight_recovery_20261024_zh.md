# Flexpert 权重恢复审计（2026-10-24）

本次仅检查已有文件、请求官方 HEAD、执行最多 8 MiB / 60 秒的小型 Range 探针，未启动完整权重下载、未加载蛋白语言模型、未运行推理，也未改动正在进行的其他模型任务。完整机器可读结果为 `flexpert_weight_recovery_audit_20261024.json`，每次探针的请求、响应头、时间、接收字节数与内容 SHA 保留在 `weights/flexpert/probe_*`。

官方服务器已恢复连接。3D 和 Seq 的 HEAD 均返回 HTTP 200，并声明支持 `Accept-Ranges: bytes`。来自官方源码 `download_flexpert_weights.sh` 的地址分别是 `https://data.ciirc.cvut.cz/public/projects/2025Flexpert/flexpert-weights/flexpert_3d_weights.bin` 与 `https://data.ciirc.cvut.cz/public/projects/2025Flexpert/flexpert-weights/flexpert_seq_weights.bin`。本地官方源码已绑定此前核实的 commit `66c5adb348b53183a1e75e7f1ad5fdb395b6ba00`，下载脚本的 SHA 也在 JSON 中保存。权重服务器本身是 CIIRC 的静态端点，提供 ETag、Last-Modified 和长度，没有提供 Hugging Face 式不可变 revision、x-repo-commit 或 LFS SHA256；不能把 ETag 写成官方整文件 SHA256。

| 项目 | Flexpert-3D | Flexpert-Seq |
|---|---|---|
| 官方完整长度（字节） | 4,843,266,070 | 4,842,603,885 |
| 官方 ETag | `"120ae6016-62f0c141fb800"` | `"120a4456d-62f0c13872180"` |
| 官方 Last-Modified | 2025-02-26 14:06:24 GMT | 2025-02-26 14:06:14 GMT |
| 已有文件逻辑长度（字节） | 4,242,281,755 | 12,013,568 |
| `.aria2` 控制文件 | 存在，756 字节 | 不存在 |
| ZIP 中央目录完整 | 否 | 否 |
| 能否当作完整模型加载 | 否 | 否 |

3D 的 4.24 GB 逻辑文件长度不能解释成已经下载了约 88%。aria2 可以先写入相距很远的分段，中间存在空洞或零填充；旧失败日志显示当时只传输了约 28 MiB。因此必须把 `.bin` 与 `.bin.aria2` 作为一对保留，不能按 `.bin` 逻辑长度使用 `curl -C -` 追加，这会跨过尚未收到的中间片段。Seq 的文件开头 64 KiB 和现有前缀末尾 64 KiB 均与官方对应 Range 逐字节相同，但这是采样核对，未证明整个 12 MB 前缀。保守恢复方案将 Seq 重新下载到独立的 `flexpert_seq_weights.verified.part.bin`，保留原先的 12 MB 文件，不覆盖、不误认完整前缀。

已有部分文件 SHA256（仅标识这些部分文件，不是官方完整权重摘要）：

| 文件 | SHA256 |
|---|---|
| `flexpert_3d_weights.bin` | `051d846eaf0dde0d8e2f597318aea6594938d93c3c05ea732d6e4ecc77bc2ff7` |
| `flexpert_3d_weights.bin.aria2` | `caaaa328793b4b6d8227a2cdbb9dab124bbeada02b9e5e429a788c9cac8774dc` |
| `flexpert_seq_weights.bin` | `c35490e0e247964cf565cbbc1fead70faed28cac515ed4c23b669f1cf2254119` |

吞吐仍不足以称为稳定恢复。每个模型的 8 MiB 探针都使用 `If-Range` 锁定上述 ETag，服务器均返回匹配实体的 HTTP 206 和正确 Content-Range，但在 60 秒上限前均未接收完整范围。3D 实收 2,543,849 字节，平均 42,389 B/s；Seq 实收 1,495,849 字节，平均 24,924 B/s；两者均以 curl 28 超时退出，JSON 中的 `complete_exact_range_and_ETag_verified` 明确为 false。TLS/首字节时间也已保留：3D 首字节 1.828 秒，Seq 6.976 秒，因此不是用 64 KiB 请求的握手时间冒充稳态速度。若机械地延续这次平均速度，单个完整权重分别约需 31.7 小时和 54.0 小时；这只是当前慢速探针的数量级外推，不能当作可靠 ETA。探针后 E 盘仍有 102,569,369,600 字节可用（95.53 GiB），磁盘容量充足。

已经准备 `weights/flexpert/resume_verified_flexpert.ps1`，只做了 PowerShell 语法检查，没有执行。它默认仅执行 HEAD 指纹校验和输出命令预览；必须显式加 `-StartDownload` 才会开始传输。每次恢复重新核对原 URL 的 ETag 和长度，并用 `If-Match` 防止服务器实体改变；使用一个连接、有限重试、最多一小时，保留 aria2 控制文件。当前速度很慢，建议先保留为可恢复后补项，待网络稳定后再决定是否运行有界下载。

```powershell
# 仅验证与预览，不开始完整下载
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/weights/flexpert/resume_verified_flexpert.ps1' -Variant 3d

# 以后明确决定恢复时，最多运行一小时；本次没有执行此命令
& 'E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/weights/flexpert/resume_verified_flexpert.ps1' -Variant 3d -MaximumSeconds 3600 -StartDownload
```

完成的判断必须同时满足传输成功退出、完整文件长度吻合、ETag 未改变、aria2 无剩余控制标记、整文件 SHA 留档，以及 ZIP CRC 和内容目录校验。官方端点没有提供整文件密码学摘要，所以 HTTP ETag 只能用作版本实体一致性条件，不能替代官方 SHA256。两份 `.partial.json` 标记已经写在原文件旁，明确禁止把目前的文件当作完整模型。
