# KRAS实例链内对照的内部有效性审计

本文件为内部分析记录，不作为论文新增讨论或限制段落。最终编辑决定为删除原拟3.5第四段链内比较及方法2末尾相应句，保留完整模型的全10肽与去G6/D6后全部9个未替换位点的直接MD验证。正文和图表不需为这一撤回增加防御性说明。

## 已核对的来源

`tcr_case_v40_20261017/family_predictions/prediction_manifest.jsonl`中7OW5与7OW6的两个Intra条目均对应`variant=intrachain`、`fold=0`，检查点SHA-256为`6b069d9d3521e95f5f2d63a3f19d64fb4071862efa4682473387b0ca45ff73c3`。路径为`/root/flexppi_completion_20261010/family_folds/evaluation/intrachain/itrace_7ow5.npz`及对应7ow6。该来源是全组成链家族第0折模型，与主Full使用相同训练划分。

`family_folds/configs/intrachain_fold0.yaml`明确`information_ablation: intrachain`、`metrics: [rmsf, chi_entropy]`、seed7、40轮、对应fold0训练与测试清单；条件特征关闭。`family_folds/evaluate_family_models.py`按检查点config构建`build_information_ablation(..., variant)`，断言variant、epoch与seed，`strict=True`加载权重，直接从`prediction['mean'][0]`导出NPZ的mean。未见nm/Å单位换算、列顺序交换或额外倍率。源模型`IntrachainFlexPPI`逐链调用同一IPA编码器并将输出回填到原残基顺序，没有额外输出尺度转换。

## 数值问题及解释边界

Intra对7OW5肽V1/V2/V3的绝对RMSF为约7.048/6.519/5.822 Å，而MD与Full处于约0.4–0.7 Å范围。即使去本位9位的变化MAE 0.503741 Å及rho −0.083333算术无误，这项对照也同时带有严重绝对幅值失准。不能仅凭这一差分成绩将Full的优势解释为已经单独确认伙伴链信息的作用。其原因可能涉及短肽上的外推与校准，但本轮未进一步验证，不将推测写成诊断结论。

因此撤回该实例中的环境贡献强推论。全队列完整/链内对照的既有评价不受本次局部编辑决定影响。本文实例直接论证的能力为：同一Full模型从相近复合物静态结构辨认单点抗原改变伴随的主链变化，包括未替换的全部9个位点。没有新增训练、模型推断或MD。
