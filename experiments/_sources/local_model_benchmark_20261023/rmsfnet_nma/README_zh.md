# RMSF-net 与弹性网络模型：本机运行与数值审计

本目录只处理正式冻结队列的3,063条结构、961,321个残基位置。输入为完整原始组装体，保留原始链、作者残基编号及插入码；在每个记录上核对PDB字节哈希、CA数量和原顺序。它与把4,973条组成链分别送入单链接口的运行条件不同。用户确认的主图分组为单链2,193、其他多链626、TCR–pMHC 244；本目录不重新分组，不改变参考标签。

`run_nma_scipy.py` 在本机重新计算ANM和GNM。由于当前Windows环境没有可编译ProDy扩展的MSVC工具链，采用保存于 `upstream_prody/` 的ProDy原始方程及求解流程进行NumPy/SciPy重现，并在运行全队列之前完成三种来源的逐残基数值验证。参数固定为15 Å截断、弹簧系数1、最多20个最低非零模式，零模式阈值严格沿用ProDy的1e-6。每个记录的模式数与原冻结脚本一致，取 `min(20,N-1,3N-6)`。三来源共1,433个残基的ANM最大绝对误差8.54e-13、GNM最大绝对误差1.33e-15；原始对照值仅用于验证实现，不用于拟合、缩放或替换新输出。完整结果位于 `nma_scipy/predictions/`，每条记录为 `record_<sample_id的SHA256前24位>.npz` 与同名JSON。NPZ中的 `anm`、`gnm` 是相对柔性幅度，不能报告成经过校准的Å单位绝对误差；其余字段保存平方幅度、链/残基身份、CA坐标及顺序。

全量本机计算已完成3,063/3,063条记录、961,321个位置，失败0，总计1,892.30秒。全位复核发现ANM最大误差1.4071e-7，全部低于1e-6；GNM初版有10条记录、151个位点超过该阈值，不能声称与原ProDy输出逐位完全一致。继续审查原始C源码发现，ProDy默认GNM的CKDTree把坐标及距离累计转换为float32，而ANM默认dense路径保持float64。对全部3,063条记录扫描两种邻接精度后，只有4f55_A的一条临界接触不同：作者编号19与92之间的距离平方在float64为225.00000000000003，在原CKDTree float32路径为225，边数由2,823变为2,824。已据官方实现的精度规则进行独立v2重算，和原ProDy该记录的最大差降到9.71e-16；这次修正由源代码和邻接矩阵决定，不按预测表现选值。初版全部原文件保留。

最终GNM读取规则为 `gnm_ckdtree_v2/resolved_prediction_index.csv`：文件完整列出3,063条记录的数组、状态、来源与输入/输出哈希，4f55_A读 `gnm_ckdtree_v2/predictions/record_87a406a20cf12f92855d78ca.npz` 的 `gnm`，其余3,062条继续读本轮新计算的v1；所有ANM仍读v1。等价的程序入口是读取 `gnm_ckdtree_v2/overlay_manifest.json` 的 `overlays`，仅替代列明的记录，未列明时回落v1。复现顺序是 `run_nma_scipy.py`、`correct_gnm_ckdtree_precision_v2.py`、`finalize_gnm_v2_audit.py`。v2源码哈希与官方永久链接分别保存在overlay manifest和 `source_refs.json`。

修正后GNM仍有9条记录的40个位点超过1e-6，最大差0.06469；它们在第20与第21个模式之间出现2–4重简并，谱隙仅约3.6e-15–7.5e-14，两种邻接精度的边完全相同，零模式数也均为1，总平方波动差仅约1e-14。这意味着把同一组完全相同频率的运动截取前20个时，不同数值库可在该组中选出不同方向，局部柔性值随之变化。保留这些新本机结果并明确标注为算法重现，不以旧缓存回填，也不宣称GNM全位等价。十条诊断记录的rho与AP影响、边数、谱隙和零模式数见 `gnm_ckdtree_v2/rho_ap_metric_impact.csv`，其中最大rho差约0.01149，最大AP差约0.01246；这只是所诊断记录的差异，正式横评仍重新统计全部3,063条记录。完整边扫描、临界pair、逐记录数值比较和最终解释保存在同目录审计文件中。

`run_rmsfnet_local.py --smoke` 已完成ATLAS、DynaRepo、ITrace各一条完整组装输入：16pk_A、A000A、itrace_1ao7，共3条记录、1,433个残基、211个密度块。官方网络代码与PDB-only权重均保持原样；使用既有冻结PDB密度适配器重现其文档中的4 Å Gaussian molmap处理，1.5 Å格距、12 Å外边界、40³输入块和10³输出核心。这里没有运行Chimera本体，因此不声称与Chimera二进制逐浮点一致。保存网络输出完成10^x−1反变换后、结构内归一化前的分数，按相对柔性输出评价，不将其解释为Å。三条原始smoke的总计时为157.99秒，包含输入解析、密度构建及前向，排除权重下载和首次导入。

RMSF-net尚未启动3,063条记录全量。精确几何统计得到全队列153,284个需预测密度块，按三来源原始batch2的实测吞吐估计约31.77小时；这只是耗时估计，绝不是完成的全量结果。`profile_rmsfnet.py` 仅检查float32的batch1/2/4和channels-last内存布局，只有小幅加速，不能据此当作完整队列性能。后续三条布局等价验证在内存不足时依root要求暂停，见 `rmsfnet/channels_last_validation/deferred.json`，不能把该验证标记为完成。现有3条smoke必须在覆盖率表中显示3/3,063，不能纳入满覆盖主榜。

运行与复核均使用隔离解释器 `../runtime_py312/Scripts/python.exe`。命令分别为 `run_nma_scipy.py --smoke --verify-reference`、`run_nma_scipy.py`、`audit_nma_full_reference.py`；最后一条在全量完成后对全部961,321个位置与旧官方ProDy缓存逐位比对，保存每记录误差、全局最大误差和超过1e-6的位置数。对称简并模式在20个模式截断边界可能使不同数值库返回不同基向量，审计将如实报告差异，不通过替换缓存或改模式数掩盖。RMSF-net的原始续跑入口为 `run_rmsfnet_local.py`（不带 `--smoke`），目前没有自动启动。

`summarize_execution_inventory.py` 只读取各模型的运行状态和已完成摘要，不加载模型，输出 `execution_inventory/model_execution_inventory.csv` 与 `timing_table.csv`。它区分冻结OOF复用、新本机全量、仅部分smoke和未运行候选。不同方法的线程、输入条件及计时范围并不一致，例如FastProtFlex的计时包括每条记录R进程启动且限制为单线程，而神经网络通常复用已加载模型；这些秒数不能直接当作统一速度排名。科学指标及最终覆盖率仍以root的逐位身份验证与统一计分流程为准。
