"""Combine independently checked local evidence without changing the benchmark."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
CORE = ROOT.parents[1] / "core_corrected_pipeline"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


profile = pd.read_csv(OUT / "prediction_profile_by_dataset.csv", float_precision="round_trip")
attributes = pd.read_csv(OUT / "gap_by_record_attributes.csv", float_precision="round_trip")
regions = pd.read_csv(OUT / "region_summary.csv", float_precision="round_trip")
names = pd.read_csv(OUT / "all_2203_test_proteins_with_names.csv", float_precision="round_trip")
influence = pd.read_csv(OUT / "outlier_influence_diagnostic_only.csv", float_precision="round_trip")
audit = json.loads((OUT / "profile_audit.json").read_text(encoding="utf-8"))
contribution = pd.read_csv(OUT / "split_review/protein_gap_source_contribution.csv", float_precision="round_trip")
bf_checkpoint = json.loads((OUT / "scoring_review/backflip_checkpoint_identity.json").read_text(encoding="utf-8"))
static_audit = json.loads((OUT / "backflip_review/static_input_audit.json").read_text(encoding="utf-8"))
assert static_audit["status"] == "PASS" and static_audit["passed_records"] == static_audit["records"] == 3063
assert static_audit["failures"] == []
assert audit["records"] == len(names) == 2203


def p(dataset, method):
    frame = profile[(profile.dataset == dataset) & (profile.method == method)]
    assert len(frame) == 1
    return frame.iloc[0]


atlas_full, atlas_bf = p("atlas", "Full"), p("atlas", "backflip-2.1 (chains)")
dyn_full, dyn_bf = p("dynarepo", "Full"), p("dynarepo", "backflip-2.1 (chains)")
share = contribution[(contribution.comparator == "backflip-2.1 (chains)") & (contribution.metric == "spearman") & (contribution.dataset == "atlas")].iloc[0].share_of_total_signed_gap
core = regions[(regions.dataset == "atlas") & (regions.region == "remaining_chain_core")].set_index("method")
folds = attributes[attributes.stratifier == "original_fold"]
lines = ["# 蛋白组与BackFlip-2.1差距调查", "",
    "调查结论：当前保存的预测中，FlexPPI在蛋白组的RMSF排序性能确实低于BackFlip-2.1。全量复核没有发现能撤销这一结果的样本错配、输出列、有效位点、单位或汇总错误。已经定位到差距主要发生在哪一批蛋白、误差表现是什么；训练条件为何导致这种差异仍需受控训练实验，不能直接归因于某一个因素或声称外部模型泄漏。", "",
    f"这2,203条实际测试输入包含1,929条ATLAS蛋白（87.56%）及274条DynaRepo蛋白。按用户确认的生物学规则，2,193条为单链输入、10条为同源多聚体或同一成熟蛋白的拆分链输入。原始PDB可以包含其他伙伴，但本次类别针对实际模拟/模型输入；本表中的ATLAS全部是单条链。完整清单按差距从大到小保存在all_2203_test_proteins_with_names.csv，包含原PDB链ID、长度、测试折、官方链名称、双方ρ/MAE和分数差。", "",
    "| 蛋白来源 | 条数 | FlexPPI ρ | BackFlip-2.1 ρ | FlexPPI MAE (Å) | BackFlip-2.1 MAE (Å) |",
    "|---|---:|---:|---:|---:|---:|",
    f"| ATLAS | 1,929 | {atlas_full.mean_rho:.3f} | {atlas_bf.mean_rho:.3f} | {atlas_full.mean_mae:.3f} | {atlas_bf.mean_mae:.3f} |",
    f"| DynaRepo蛋白 | 274 | {dyn_full.mean_rho:.3f} | {dyn_bf.mean_rho:.3f} | {dyn_full.mean_mae:.3f} | {dyn_bf.mean_mae:.3f} |", "",
    f"以上主对照为BackFlip-2.1逐链版，沿用主图口径。ATLAS占总ρ差距的{100*share:.2f}%，DynaRepo蛋白的平均ρ差约{dyn_bf.mean_rho-dyn_full.mean_rho:.3f}。因此这次调查首先应针对ATLAS，而不是把全部蛋白混在一起寻找一个通用解释。DynaRepo结果也不是已证明无训练暴露的外部测试。", "",
    f"差距分布很广：BackFlip-2.1在{audit['bf2_rho_wins']:,}/2,203条（{100*audit['bf2_rho_wins']/2203:.2f}%）蛋白上ρ更高，FlexPPI在{audit['full_rho_wins']}条更高；逐条ρ差值中位数（Full减BackFlip）为{audit['delta_rho_quantiles']['0.5']:.3f}。五个原始测试折的平均差值均在{folds.delta_rho.min():.3f}至{folds.delta_rho.max():.3f}之间，没有单个测试折独自崩坏。即使仅作为离群影响诊断，移除差距最大的100条后，两者均值仍相差{abs(influence.iloc[-1].delta_rho):.3f}。这不是删样本后的替代成绩；原基准仍保留全部记录。", "",
    f"逐残基误差也有明确特征。在ATLAS中，FlexPPI平均预测比参考低{abs(atlas_full.mean_bias):.3f} Å；同一蛋白不同残基之间的预测RMSF标准差与参考RMSF标准差之比，中位数为{atlas_full.median_std_ratio:.3f}，BackFlip为{atlas_bf.median_std_ratio:.3f}。这个无量纲比值描述残基之间的起伏，不是模型预测不确定性或轨迹的时间方差。直观地说，FlexPPI往往把柔性峰压低，残基之间的起伏也压小。这个现象影响幅度误差；但ρ只看残基排序，统一加一个常数或乘一个正数不会改善ρ，所以只做单位/标定修正不能解决排序落后。", "",
    f"为了判断是否只是柔性末端造成差距，额外将每条输入链两端各10个残基与剩余链内部位置分开计算。ATLAS去掉这些末端后的平均ρ仍为FlexPPI {core.loc['Full','rho']:.3f}、BackFlip {core.loc['backflip-2.1 (chains)','rho']:.3f}；MAE分别{core.loc['Full','mae']:.3f}与{core.loc['backflip-2.1 (chains)','mae']:.3f} Å。这里的链内部位置只按两端定义，不表示结构埋藏核心；该区域仅用于诊断，未修改主基准掩码。说明问题也发生在链内部的柔性分布，不局限于少数末端尖峰。", "",
    "以下是差距最大的几条具体测试输入。它们是按误差排序选出的诊断案例，不代表全部测试蛋白的生物学分布。官方名称按ATLAS原PDB作者链映射，而不是把整个PDB条目的所有蛋白名称混作这一条链。", "",
    "| 实际测试输入 | 长度 | 原测试折 | FlexPPI ρ | BackFlip-2.1 ρ |",
    "|---|---:|---:|---:|---:|"]
for row in names.head(6).to_dict("records"):
    lines.append(f"| {row['sample_id']} | {row['n_residues']} | {row['original_fold']} | {row['spearman::Full']:.3f} | {row['spearman::backflip-2.1 (chains)']:.3f} |")
lines += ["",
    "例如1NG6_A的参考平均RMSF为2.393 Å，FlexPPI为1.474 Å，BackFlip为2.422 Å；FlexPPI的残基起伏幅度只有参考的约54%，但BackFlip约为102%。4KE2_C则不是简单低估：FlexPPI平均幅度偏高，排序ρ仍接近零，说明也存在柔性区域位置判断错误。不能把所有失败案例都解释成同一种缩放偏差。具体每个残基的参考值及7种模型预测保存在worst20_residue_profiles_diagnostic_only.csv。", "",
    "全量身份复核覆盖Full的3,063个原NPZ和BackFlip两种输入的6,126个NPZ。Full确实取mean[:,0]的RMSF均值列，不是Hχ或预测方差；原标签、有效位点、氨基酸、链、作者残基号与共同冻结参考一致。9,189条记录的独立ρ、Pearson和MAE复算与已发布表差值全部为0。Full五折实际checkpoint哈希、seed7、epoch40和输出列顺序也与原导出清单逐项吻合。不是把五个模型训练集预测当测试，也不是对测试样本做五模型平均。", "",
    "输入线索核查：BackFlip官方loader虽然读取PDB字段，但实际网络入口只使用骨架坐标/刚体和残基身份；本机规范化PDB把B-factor固定为0.00，因而没有通过B-factor读入MD答案。Full前向也不读取batch中的标签字段。两者都使用同一冻结结构来源，骨架静态文件的全量对应检查及官方源码身份核验见backflip_review。", "",
    "最值得优先验证的是训练设定。Full从头训练，同时学习RMSF和Hχ，并给ATLAS、DynaRepo、ITrace三个来源各1/3的抽样概率；它不是只针对ATLAS训练的模型。各折ATLAS样本在40轮中的预期平均访问次数约21.2次。Full固定采用第40轮，没有验证集选点；BackFlip本次权重是官方ATLAS/monomeric版本，实际checkpoint为第221轮，并按验证损失选点。这个训练目标和抽样差异与ATLAS上预测偏保守的现象一致，但只是机制候选。epoch数、batch大小和每轮抽样不同，不能直接把221/40称为计算量倍率，更不能据此断定单纯加轮数就能追平。", "",
    "另一项需要诚实保留的限制是训练暴露：Full是记录/生物学组留出，序列家族未全部隔离；protein中274条仍有同家族成员进入对应训练折。BackFlip-2.1的配置指向alphaflow train/val/test CSV，本地没有可绑定实际checkpoint的名单，因此它对这些ATLAS记录究竟是训练、验证还是未见，仍未知。不能拿BackFlip-1或DynaProt的名单代替BackFlip-2，也不能把“在ATLAS上特别好”当成已经证明数据泄漏。", "",
    "标签定义也并非完全等同：本模型ATLAS参考由各重复轨迹预计算RMSF取平均；BackFlip数据生成源码则支持合并重复轨迹计算协方差后开方，并存在不同拟合/参考选项。实际发布训练NPZ和全部生成参数未绑定，所以这只是目标迁移限制。它不会推翻双方在当前同一参考上的真实分数，也不能单独解释这次差距。原论文全3,063条结果已复现：Full总体ρ约0.796、BackFlip-2.1逐链约0.848；新蛋白子集不应与总体平均分混比。", "",
    "后续最有判别力的实验，应先固定一个独立验证集合，对Full做ATLAS专门抽样／训练、RMSF单任务、训练收敛与不确定性损失权重的受控比较，再在一次性未见测试集上评价。每次只改变一个因素，并记录验证曲线，才能区分来源平衡、双任务取舍与优化目标的影响。不能根据本次测试蛋白挑参数、删失败样本或在测试集拟合缩放。本轮只完成调查和可复核数据产物，没有启动重训练、改分或重画主图。", "",
    "完整证据分为split_review（测试成员、训练折、来源贡献与原论文口径）、scoring_review（全部原始预测、checkpoint及评分一致性）和backflip_review（输入/特征/静态结构核验）。所有诊断脚本、逐记录清单、逐残基案例及文件哈希与本报告同目录。"]
text = "\n".join(lines) + "\n"
text = text.replace("实际checkpoint为第221轮，并按验证损失选点", "checkpoint内记录epoch=221，并按验证损失保留/选点")
text = text.replace("原论文全3,063条结果已复现", "完整3,063条队列的结果已复核")
text = text.replace("官方链名称、双方ρ/MAE和分数差。", "官方链名称、双方ρ/MAE和分数差。ATLAS中的7uw7_A未在现有官方缓存找到链实体名称，名称栏留空并明确标记；它的预测与评分仍完整保留。")
text = text.replace("骨架静态文件的全量对应检查及官方源码身份核验见backflip_review。", "3,063条输入的静态结构核验全部通过，8,036个规范化结构文件的7,680,628个原子中非零B-factor数量为0；骨架坐标、来源哈希和官方源码身份核验见backflip_review。")
(OUT / "REPORT_zh.md").write_text(text, encoding="utf-8")
source_paths = [OUT / name for name in ("prediction_profile_by_dataset.csv", "gap_by_record_attributes.csv", "region_summary.csv",
    "all_2203_test_proteins_with_names.csv", "outlier_influence_diagnostic_only.csv", "profile_audit.json",
    "split_review/protein_gap_source_contribution.csv", "split_review/REPORT_zh.md", "scoring_review/SCORING_REVIEW_zh.md",
    "scoring_review/backflip_checkpoint_identity.json", "backflip_review/static_input_audit.json", "backflip_review/source_feature_audit.json")]
report_audit = {"created_utc": datetime.now(timezone.utc).isoformat(), "report_script_sha256": sha(Path(__file__)),
    "input_sha256": {str(path.relative_to(OUT)): sha(path) for path in source_paths},
    "report_sha256": sha(OUT / "REPORT_zh.md"), "benchmark_modified": False, "new_inference_or_training": False}
(OUT / "report_audit.json").write_text(json.dumps(report_audit, indent=2), encoding="utf-8")
print(json.dumps({"report": str(OUT / "REPORT_zh.md"), "records": len(names), "atlas_fraction": 1929/2203, "atlas_gap_contribution": share}, indent=2))
