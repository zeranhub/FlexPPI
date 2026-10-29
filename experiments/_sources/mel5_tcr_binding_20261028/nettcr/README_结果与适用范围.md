# NetTCR-2.2：MEL5结合面板的实际运行结果

本次对统一的12行文献面板运行NetTCR-2.2公开pan模型，全部肽采用同一模型模式和同一官方BLOSUM50编码，不混用已知ELA的肽特异模型与突变肽的pan模型，也不加TCRbase缩放。六CDR由root以ANARCII2.0.8在原晶体SEQRES上统一编号；CDR1、2、3分别取IMGT27–38、56–65、105–117，不含C104/F118。原生MEL5六环与作者训练表中的已知MEL5记录完全一致，构成输入定义的直接核验。十肽长度在官方12位上限内，未截短或改写肽。

模型来源为 https://github.com/mnielLab/NetTCR-2.2 ，冻结提交7cead3fe6dcb539ff8e2d9121586dafca1e059c2。论文为Jensen与Nielsen（2024），Enhancing TCR specificity predictions by combined pan- and peptide-specific training, loss-scaling, and sequence similarity integration，https://doi.org/10.7554/eLife.93934.3 。所有官方文件来源及SHA见download_audit.json，20个TFLite权重均保持原样。

训练重叠在预测前已检查：最终limited训练表含原生MEL5六环及ELA完全相同阳性记录（index26046）；也含成熟受体两条CDR3完全相同但其CDR1/2仍为原生的ELA阳性记录（index14497）。这两组CDR3配对连同交换生成的阴性行共11行，全部在partition0。主对照预先固定为t.0.v.1、t.0.v.2、t.0.v.3、t.0.v.4四个官方模型的均值：它们将分区0作为test而非training/validation，排除了已确认的精确CDR3配对训练和验证重叠。不能由此声称全部近缘家族严格隔离。另完整保留官方20模型均值作敏感性对照；它包含已知重叠，不当作独立验证。20个模型不是20个生物重复。

## 连续结合分数

| 条件 | 实验结合 | 主4模型均值 | 全20模型均值 |
|---|---|---:|---:|
| 3hg1 / ELAGIGILTV | 可测结合 | 0.962290 | 0.965600 |
| 4jff / ELAGIGILTV | 可测结合 | 0.218565 | 0.203049 |
| 4jfd / ELAAIGILTV | 可测结合 | 0.178530 | 0.143416 |
| 4jfe / ELAGIGALTV | 可测结合 | 0.078191 | 0.057257 |
| mature_E1A / ALAGIGILTV | 可测结合 | 0.185265 | 0.138510 |
| mature_I5A / ELAGAGILTV | 原论文NB/阴性对照 | 0.075607 | 0.057914 |
| mature_G6A / ELAGIAILTV | 可测结合 | 0.160013 | 0.111786 |
| mature_L8A / ELAGIGIATV | 可测结合 | 0.088224 | 0.057312 |
| mature_T9A / ELAGIGILAV | 可测结合 | 0.175630 | 0.153643 |
| mature_control_ILAKFLHWL / ILAKFLHWL | 原论文NB/阴性对照 | 0.124873 | 0.135606 |
| mature_control_SLLMWITQC / SLLMWITQC | 原论文NB/阴性对照 | 0.077204 | 0.080086 |
| mature_control_YLEPGPVTA / YLEPGPVTA | 原论文NB/阴性对照 | 0.061115 | 0.052312 |

在成熟受体背景中，G4A与I7A相对ELA的分数均下降，故本模型恢复了这两项肽替换的大方向，不能写成模型普遍连下降方向也无法预测。原生MEL5的分数约0.962，成熟ELA约0.219，其大小次序与SPR亲和成熟方向相反；但这属于特异性分数与亲和力排序的比较，二者仍均为实验可测结合，不是一个阳性/阴性误判的直接定义。I7A可测结合与I5A的NB分数约0.0782与0.0756，较为接近；阴性ILAKFLHWL的分数约0.1249还高于I7A/L8A的阳性分数，因此此面板的连续分数不能被任何一个固定高分为阳性的阈值完全分开。

分数是特异性分类模型的sigmoid输出，不是KD，也不解释或定量恢复亲和力倍数。本次没有找到官方适用于全部pan肽的通用二分类阈值，因此主结果保留连续分数，不调节阈值迎合本面板。CSV另列score≥0.5的探索性结果，只是预先固定的常用数值切分，不宣称其为作者推荐临床或实验阈值。ELA为训练26肽之一，Ala替换及若干负对照为未见肽；作者明确提醒新肽外推较弱。人工亲和成熟受体还可能偏离以自然TCR为主的训练输入分布，相关局限不能改写成柔性机制得到证明。

## 实施与核验

run_nettcr_pan.py直接导入原始keras_utils的BLOSUM50编码函数，以与官方predict.py相同的动态张量形状、输入名映射和TFLite解释器执行；未更改网络或权重，未重训。共享的独立Python3.9环境为相邻imrex/venv，使用TensorFlow/Keras2.15.0和CPU，官方环境文件原用TensorFlow2.11。运行版本差异已实测核验：官方CV表中上述11个分区0记录，以四个t.0模型复算的最大绝对误差为1.4823×10⁻⁸；再对本12行直接调用未改动的官方predict.py与包装程序比较同一t.0.v.1模型，最大误差为2.4144×10⁻⁸。均在浮点转换与输出精度范围内。

主要文件：nettcr22_pan_predictions.csv为12行汇总；predictions_per_model.csv为全部240个逐模型原始分数；standardized_predictions.json可按case_id合并其他模型；known_CDR3_pair_training_overlap.csv保留11条训练重叠；CV_reproduction_validation.csv与official_cli_validation/run.log为两项实际核验；inference_audit.json记录模型分区、输入、适用范围和运行版本。

复现命令：用imrex/venv/Scripts/python.exe运行run_nettcr_pan.py --input ../binding_panel_inputs.csv。需要保留upstream/src、20个官方TFLite权重及用于审计的limited训练表/CV参考表。依赖完整冻结见requirements_runtime_frozen.txt。本文仅为同一家族文献条件的小面板，不作总体预测性能推断，不把模拟帧或模型成员计作生物重复。
