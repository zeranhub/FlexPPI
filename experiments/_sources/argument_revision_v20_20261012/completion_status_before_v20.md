当前交付为 `E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v19_20261012_中文修订_补充实验完成_三线表.docx`。中文稿件和补充实验整合已经完成，共45页、14张原生可编辑三线表。文件SHA256为 `b13c4f1e1de8458609fe6ab8d0eb17ceea81b0f0992d1a5d3d0728eb056bdea6`。原始6条批注、批注锚点、原图片与原生数学公式均保留；新增外部验证图4与正文、表2采用相同的2,430条训练方案。独立数值及引用检查通过，全篇实际Word渲染及新增/变化页面视觉检查通过。原v13文件未覆盖。

实际实验已全部完成：单目标十次及全组成链家族十次40轮训练、两组完整折外预测、单目标/联合配对分析、新严格外部32条记录评价，以及RMSF-net完整3,063条结构、PLANET full/seq各4,973条链的独立核验和共同残基评分均通过。服务器完成日志为2026-10-12T00:32:46Z；服务器钟落后本地约11分39秒，校正到北京时间约08:44（校正估计）。20个最终权重的哈希再次核验一致。完成证明见 `computation_completion_status_readonly_20261012_1158.json` 和 `external_adapter_repair_20261012/final_completion_verification.json`，不因再生成完成证书改变既有指针。

稿件已如实更新结果：联合训练相对于单目标的平均MAE收益约0.8%和1.8%，RMSF排序未改善；家族隔离后完整模型整体MAE为0.616 Å/0.173 nats，ITrace仅第0折一次家族/来源外推，其RMSF排序约0.680。新严格外部删除633条、保留2,430条，完整模型RMSF MAE0.474 Å略高于链内及静态对照，Hχ误差相对静态小幅改善约1.5%。PLANET full总体RMSF MAE0.568 Å/ρ0.828优于FlexPPI，但RMSE1.186 Å较高；RMSF-net相对分数不赋予Å。发布权重与同语料训练对照分别解释。旧632/2,431条协议的图S1、S10、S11和表S12明确作为历史敏感性保留。

KIX按用户最新要求复用2026年9月四种状态各一条100ns生产轨迹；400段DCD当前文件哈希、帧数、原子数、大小及全部质量分析源表已再次核验。12条新模拟及其后三重复分析保持取消，没有启动。正文及表S7/S8保留活动区定位与结合状态差值方向不一致、前后半程采样依赖及无独立重复统计的限制。hmu-gpu现有公钥认证失败，是否有另一个服务器副本未知，不影响本地轨迹复用。

本轮复现源表、全部审计和新图的矢量PDF/SVG及PNG位于 `completed_results_20261012/`。最终20个检查点、4份逐残基压缩表、训练配置和原完成证书归档为 `completed_scientific_sources_20261012.tar.gz`，成员清单及SHA256另存 `completed_scientific_sources_manifest.json`，本地备份核验记录另存 `local_scientific_backup_audit.json`。原始服务器队列状态仍标为 `computations_complete_pending_manuscript_update`，这是服务器计算阶段的最终状态；本地 `delivery_status.json` 是论文实际整合和核验完成的权威记录。没有创建自动唤醒或新增任务。
