# ERGO-II 公开预训练模型对照

采用作者官方仓库 https://github.com/IdoSpringer/ERGO-II 的两个公开AE checkpoint，McPAS版本1meajht与VDJdb版本1veajht，同时报告而未按本例结果挑选。仓库提交85d320ab03ade33460be9a81ea3a51b8d37cd998，论文Springer、Tickotsky、Louzoun，Frontiers in Immunology (2021)，doi:10.3389/fimmu.2021.664514。官方预测说明 https://tcr2.math.biu.ac.il/help 规定输出为0–1的结合评分、CDR3以C开始F结尾，允许可选信息缺失。

输入由root预先冻结的binding_panel_inputs.csv映射到官方CSV列；12行全部保留，包含四个晶体核心状态、原论文成熟受体的其他丙氨酸扫描以及明确SPR阴性对照。CDR3由ANARCII 2.0.8编号取得，保留C104/F118；MHC为HLA-A*02:01，细胞类型CD8。TRAV12-2与TRBV30有独立注释，J基因未知保留空白，不捏造等位基因。McPAS词表识别这两个V名称；VDJdb词表仅有带等位基因的形式，按官方函数落入UNK。两模型都识别该MHC。ERGO-II未直接读取工程化CDR1/2序列，因此不能完整感知所有成熟受体替换。

主结果在panel_results/all_ergo_scores.csv，逐模型输入、预测、权重与源码SHA及训练重叠检查在同目录。McPAS四核心态3HG1/4JFF/4JFD/4JFE评分分别8.526699e−7、1.178444e−6、5.201456e−6、1.460777e−6；VDJdb分别0.0746754、0.0860479、0.0889360、0.103508。G4A与I7A均未表现为相对成熟ELA更低的结合评分。这只是本panel在所述配置和缺失特征条件下的结果，不代表ERGO-II总体性能，也不把score比例解释为KD比例。

官方预测接口没有给可移植到SPR测定的固定二元阈值。本次额外保存未调参的0.5参考分类，不能称为经过作者验证或针对该实验校准的阈值。按此参考，McPAS全部12行低于0.5；VDJdb仅E1A、阴性对照SLLMWITQC和YLEPGPVTA高于0.5。原四核心态均为实验可测结合，G4A/I7A不能标为实验negative。所有连续原始分数和原实验标签均保留，未以阈值搜索制造失败案例。

权重由官方Predict.py按默认严格方式加载，所有层与forward/step/validation_step直接使用原始Trainer.py与Models.py。旧pytorch-lightning 0.7.1仅在推理中用作nn.Module容器，故运行适配器提供nn.Module基类及不执行训练的导入占位，未更改网络、参数或张量运算；torch.load只将权重放到CPU。官方AE初始化路径TCR_Autoencoder与发布Models/AE目录不一致，因此将发布的两个AE权重逐字节复制到预期路径，并核验SHA相同。两模型12行批推理与逐行推理最大绝对差分别9.31×10⁻¹⁰和1.19×10⁻⁷。官方example另外跑通，但不是本例预测或性能验证。

训练重叠检查显示：McPAS没有本panel任一αβ–肽或β–肽精确重叠；VDJdb含野生型β–ELA的2条记录，αβ–ELA组合没有精确重叠。ELA肽本身在McPAS和VDJdb训练样本分别出现1161与10274行，不能称此例为unseen peptide，也未完成同源家族隔离证明。原训练文件随官方下载保留，逐行计数见审计。

复现命令：使用requirements_runtime_frozen.txt列明的Python环境，执行run_ergo.py --input panel_official_schema.csv --output-dir panel_results；脚本直接调用download/ERGO-II-master下的官方推理函数。环境复用了本地前次亲和力对照的torch CPU环境，并添加pandas 2.2.3及其依赖，未改变torch/numpy版本。原始官方代码、模型文件和93MB下载存档都保留在本地；论文轻量结果包可只纳入源代码、输入输出与审计，不打包权重和训练样本。
