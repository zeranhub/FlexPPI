from pathlib import Path
import json
r=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\imrex')
(r/'结果说明_ImRex.md').write_text('''# ImRex：官方预训练模型的12例核查

使用作者仓库 https://github.com/pmoris/ImRex 的固定提交d3975b744c030c4bc67c942c82a32b1c83857c4e和两个随仓库发布的预训练H5模型，原始src/scripts/predict/predict.py完全未修改。默认主分析采用README明确推荐的2020-07-24 downsample模型；2020-07-30未下采样模型同时保留作敏感性核查，未按本例分数选择模型。论文：Moris et al., Current challenges for unseen-epitope TCR interaction prediction and a new perspective derived from image classification, Briefings in Bioinformatics, doi:10.1093/bib/bbaa318。

统一输入为root在预测前冻结的binding_panel_inputs.csv，12行均由官方CLI成功输出；CDR3β含保守C104/F118，来自原晶体SEQRES的ANARCII统一编号。ImRex只接收CDR3β与肽序列，未利用TCRα、CDR1/2、HLA或三维结构；输入CDR3β长度10–20、肽8–11，全部12行在适用长度内。主模型与敏感性模型使用官方四种理化特征、绝对差值及居中零填充到20×11×4。原生MEL5与成熟受体之间，模型只能看到CDR3β变化，不能把此模型当作完整受体环境模型。

默认downsample模型对3HG1、4JFF、4JFD、4JFE的score分别为0.755036、0.904950、0.819486、0.727221。四个状态均高于固定0.5参考阈值，符合原文均可测SPR结合；G4A和I7A的分数也均低于成熟ELA。因此这组结果不支持“ImRex在核心两种肽突变的大方向均失败”，必须保留其正确结果。I5A原文NB但得分0.631041；无关肽YLEPGPVTA原文未检测到结合但得分0.610823，属于该阈值下的两个假阳性。其余两个无关肽ILAKFLHWL、SLLMWITQC分别0.277511、0.117412。12行在0.5参考阈值下TP=8、FN=0、TN=2、FP=2。

未下采样full模型对四核心态分别0.657078、0.751317、0.419080、0.764587：G4A低于0.5，I7A略高于成熟ELA；I5A=0.231631低于0.5，YLEPGPVTA=0.632741仍高于0.5。12行TP=7、FN=1、TN=3、FP=1。该checkpoint差异说明不能选择性展示某一版本以制造失败或成功。0.5仅是与作者trainer中BinaryAccuracy()/Precision()/Recall()默认行为相符的固定sigmoid参考口径，作者CLI本身只输出连续score，没有为本SPR测定校准的生物学结合阈值；本次未按本例数据调参。score不是KD，也不能从其差值推算五万倍亲和力效应。

实际训练CSV已随模型审计。两版均含野生MEL5 CDR3β CAWSETGLGTGELFF与ELAGIGILTV的完全相同positive pair，因此3HG1绝不能称独立测试；成熟CDR3β CAWSETGLGMGGWQF在两版训练CSV中均未出现。ELA肽在默认模型训练中有400个positive（含negative共835行），full中960个positive（含negative共2059行）；G4A和I7A肽在两个训练CSV中均未见。精确逐行重叠见inference_audit.json与audited_scores_*.csv。12例均是同一案例家族，分类计数仅为描述性核查，不表示模型总体性能。

执行环境为Windows CPU、Python3.9、TensorFlow/Keras2.15、numpy1.26.4、Biopython1.79、pyteomics4.7.5；冻结依赖见requirements_frozen.txt。旧pyteomics4.2因pandas版本元信息接口报错，发生于任何有效预测之前，升级后使用原CLI即成功；该部署错误日志单独保留，未计为预测失败。未改架构、模型权重、推理代码或再训练。用同一官方特征与TF数据构造再作API/CLI核对，所有24个分数与CLI保存值之差小于3×10⁻⁸，差异只来自CSV浮点显示精度；全精度结果见audited_scores_downsample.csv、audited_scores_full.csv及all_imrex_scores.csv。

复现顺序：prepare_imrex_inputs.py → run_official_cli.py → audit_imrex.py。upstream/为选定提交的官方源代码与原预训练目录，weights/是同一H5权重的短路径副本，training/是同一原始训练CSV的短路径副本，用于绕过Windows路径长度限制。原archive.zip与Python环境无需随小结果包分发；source_file_hashes.json、inference_audit.json保留代码、输入、训练集和权重SHA256。
''',encoding='utf-8')
(r/'deployment_issues.json').write_text(json.dumps({'non_prediction_failures':[{'stage':'git clone','error':'github.com:443 connection timeout; recovered with fixed-commit codeload official archive'},{'stage':'archive extraction','error':'Windows maximum path length in long original checkpoint folders; recovered with extended filesystem path extraction and short identical-byte weight/train copies'},{'stage':'initial import','error':'pyteomics4.2 referenced removed pandas.version metadata; recovered using pyteomics4.7.5, original inference source unchanged'}],'valid_predictions':24,'models_changed':False,'inference_source_changed':False,'training_performed':False},ensure_ascii=False,indent=2),encoding='utf-8')
print('README and deployment audit saved')
