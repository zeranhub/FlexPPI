hmu-gpu 的 SSH 服务可达，但当前非交互认证未通过。两次有限探测都在约 0.3 秒收到 `root@hpc.hrbmu.edu.cn: Permission denied (publickey,password)`，并非网络超时。`ssh -G hmu-gpu` 的生效配置是 root、hpc.hrbmu.edu.cn、30352 端口、指定现有 `C:\Users\DELL\.ssh\id_ed25519`、IdentitiesOnly=yes。一次 verbose 诊断确认连接到 172.16.1.18，服务器的 ED25519 host key 与既有 known_hosts 匹配，公钥被实际提供但服务器未接受；服务器仅列出 publickey/password，当前没有 ssh-agent socket，BatchMode 不进行交互密码认证。

由于认证在执行远程命令之前失败，服务器目录没有被检索，未获取任何 KIX/kix/2LXT/allostery 文件、prepared 溶剂体系、四状态生产轨迹、运行参数、分析报告或哈希。因此原始 system.xml/溶剂坐标是否仍在 hmu-gpu、是否存在另一套历史模拟、以及与本地或 hit-gpu 的文件是否同源，目前都只能记为未知；不能把认证拒绝写成“没有找到历史 MD”，更不能据此断言没有独立重复。此次没有修改服务器文件、环境、进程、GPU 或生产队列，也未改稿件。

本地九月四状态各 100 ns 与全部帧质量报告已由根代理找到，可作为历史批次复用；这一事实与本次 hmu-gpu 的访问限制分别记录，不因为访问失败否认本地数据。对 hmu-gpu 的后续复用必须先取得可用的已授权认证方式，再只读读取 metadata、DONE、分析 summary、system/protein/solvated 哈希及轨迹目录；若生产 DCD 哈希与本地九月批次相同，应归为同一批副本，不能作为新独立重复。要确认另一套独立重复，需对应独立种子和完整运行记录，且与已有轨迹的字节哈希不同，不能仅根据目录名中的 rep02/rep03 判断。

本地 README、各状态分析 summary 和队列快照一致记录历史目录 `/ranze/runs/flexppi_kix_md_20260914`、环境 `/ranze/envs/flexppi-kix-md-20260914`。这是可供下一次成功登录后优先验证的实际路径线索，但现有文件没有提供足够证据把它归属为 hmu-gpu；此次无法进入服务器，故不将该路径写成“hmu-gpu 上已找到”。

完整客观状态见同目录 `connection_audit.json`。本次仅执行两次 bounded SSH 探测，以及读取生效 SSH 配置和公钥身份的诊断，未读取或输出私钥内容，也未进行反复尝试或扩大无关目录搜索。
