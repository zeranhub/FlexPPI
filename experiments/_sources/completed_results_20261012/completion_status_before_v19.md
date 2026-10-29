当前可审阅稿件为 `E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v18_20261010_中文修订_复用KIX单轨迹_三线表.docx`。稿件共40页、10张三线表，已完成Word实际渲染核查；原批注、批注锚点与既有图片保留。当前文件SHA256为 `1ee8fed276ee8f86e854f79e242790af86def9d1352c9bc2881b2052e66485b1`，排版检查记录为 `v18_layout_audit.json`。新增补充材料F、表S7/S8，使用2026年9月四状态各一条100ns的完整KIX生产轨迹；400段DCD当前哈希及全部质量分析核验通过。独立内容审查通过。稿件如实区分活动区定位与响应预测：六组L12平均RMSF状态差值均与单轨迹MD参考方向相反，且前后半程存在时段依赖，不宣称收敛、独立重复不确定性或变构机制验证。此为已完成分析的阶段性中文稿，单目标与家族隔离模型训练和新增基线全队列评价尚未完成。

已经可以用于稿件的新增结果包括：冻结3,063条记录上的方差覆盖诊断与分组事后校准；官方RCSB条目元数据补充（3,060条记录解析成功，3条保留未知）；ANM/GNM完整3,063条记录及961,321个共同有效RMSF位置的比较；预先抽取90条记录、依据原始时间戳排除7条短轨迹后，83条记录的100ns共同窗口标签重算与全部冻结对照的配对评价。旧新标签的残基对应与有效掩码一致。ANM/GNM输出只用于相关和标准化剖面误差，没有用测试标签拟合Å幅值。标签子集分析不等同于全队列重训。

2026-10-12北京时间早晨已核对 `hit-gpu:/root/flexppi_completion_20261010`：10次单目标训练和10次全链家族隔离训练均已完成，单目标折外预测、配对分析及家族折外评价均有完成审计。10次单目标训练结束于服务器UTC2026-10-10 23:23，家族训练及评价结束于UTC2026-10-11 15:35。队列随后在 `evaluate_external_2430` 步骤退出，日志显示两目标累计数组与三通道标签掩码的形状不一致。错误已修复并实际通过全部32条外部记录的Full/Intra评估；两目标通道依据标签名称显式对应，两模型与静态对照逐残基标签和掩码一致，有效位置18,748/16,265，实测15.8秒。新父队列PID19962已恢复剩余CPU独立核验与评分；已完成训练不需重跑，冻结权重和输入保持，12条新MD取消控制不变。具体恢复结果以新的服务器状态、修复审计和最终完成门控为准；不能把已完成20个检查点等同于论文全部补充实验已经完成。

2026-10-12早晨快照显示CPU预测输出已齐：RMSF-net 3,063/3,063条记录，PLANET full与seq各4,973/4,973条链。PLANET两项队列均complete且有独立run_audit；RMSF-net全部逐项预测complete但尚无最终独立核验run_audit，剩余主队列应完成它的核验和共同残基评分。所有这些数量是文件/完成状态核查，正式共同评分仍须使用通过完整性门控的输出。PLANET保留公开float32权重、默认300K和原训练范围；RMSF-net使用PDB模拟密度适配，并非真实冷冻电镜输入的直接复现。

训练之后的基线协议已冻结为 `baselines/gpu_additional_steps.json`，实际步骤都是CPU等待、独立核验及共同残基评分。PLANET只有独立核验全部4,973条链且float32、哈希、序列和单位检查通过才会标记complete；主流程只读复核同一完成证书，避免再次生成证书导致哈希指针失配。RMSF-net先等待原进程退出、持独占锁检查缓存，再作全部3,063条记录的独立核验。共同评分分别写 `baselines/additional_original_oof` 和 `baselines/additional_family_oof`，不可混用两个划分。失败或缺失会阻止后续步骤，不能通过删除失败样本变成完整成绩。

用户已将KIX策略改为优先完成论文、复用本地及hmu-gpu已有MD，每种条件的一次轨迹即可。`kix_md_triplicate_rebuild_20261010`中的四种体系保留已完成的输入修复、力场哈希、重原子坐标核查和40ps技术验证；12条尚未启动的新生产轨迹及其后三重复响应差值、逐重复质量分析已取消。取消已在运行时生效：只替换监控父进程71135为76963，新父读取`pipeline_control.json`、重新持有原`gpu_queue.lock`并等待原训练队列71188；训练进程75948及所有CPU基线进程的PID、启动时间和命令保持一致。证据为`pipeline_md_cancellation_audit_20261010.json`。主队列不会再启动这批新生产模拟，已有单轨迹按实际来源、长度及质量审计后用于结构实例。hit-gpu目前只查到四个40ps技术DCD，每个40帧、1ps保存间隔；它们不能替代生产MD，也不能作为三次独立重复或用于宣称收敛。

查看进度可以执行 `ssh hit-gpu '/root/flexppi_completion_20261010/python311/bin/python /root/flexppi_completion_20261010/collect_remote_progress.py'`，它只读取日志、状态和GPU管理器，不加载权重或启动CUDA。完整快照写入服务器 `completion_progress.latest.json`；本地同名文件只是最后一次同步快照。主流程日志为 `pipeline_resumed.log` 与 `pipeline_logs/`，单目标训练日志为 `single_target/logs/`，CPU基线日志为 `baselines/rmsf_net_full.log`、`baselines/planet_full_cpu_queue.log` 和 `baselines/planet_seq_cpu_queue.log`。先确认现有进程和锁，再考虑恢复，不能重复启动正在运行的队列。

服务器计算完成时，主状态将成为 `computations_complete_pending_manuscript_update`。这个状态表示计算已完成，仍需取回实际源表、解释结果、更新正文/三线表/科学图、核查旧632/2,431方案与新633/2,430方案的表述，并重新渲染整个Word文档。未取得完整训练、预测和全队列核验输出时，不填写相应成绩；KIX单轨迹分析已经完成，不再等待新的重复轨迹。目前没有创建每小时自动唤醒本对话的任务；服务器上的持久队列会继续执行，稿件后续更新需继续本对话或另行授权自动检查。


本轮历史资产搜索已完成可访问范围：本地权威轨迹在 `E:/Project/FlexPPI/artifacts/external/kix_md_20260914/{apo,myb,mll,both}/production/`，各100段、100,000帧、100ns，共约8.78GB。原始指标及全帧质量分析在 `artifacts/kix_md_20260914/pilot_results_20260916/`。400段当前SHA、帧数、蛋白原子数、文件大小及五份质量分析表哈希均已核查。模型只使用 `paper/figures/case_validation_v13_20261009/data/` 下修正的第1折结果；`paper/data/case_validation_20261009/` 同名旧表不用于v18。hmu-gpu的SSH服务可达但现有公钥认证被拒绝，未能搜索其目录，是否保存另一套轨迹或原溶剂体系仍为未知；访问限制不会影响本地轨迹复用。详细证明见 `kix_existing_reuse_20261010/reuse_audit.json`、`hmu_audit/connection_audit.json` 和 `v18_content_review.json`。
