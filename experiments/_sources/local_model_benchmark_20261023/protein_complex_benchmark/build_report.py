"""Factual two-class report and manuscript text from completed frozen scores."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "evaluation"
FIG = ROOT / "figures"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    audit = json.loads((OUT / "scoring_audit.json").read_text(encoding="utf-8"))
    summary = pd.read_csv(OUT / "group_summary.csv", float_precision="round_trip")
    paired = pd.read_csv(OUT / "paired_differences.csv", float_precision="round_trip")
    exposure_path = ROOT / "references/training_exposure_protocol_sources_sha256.json"
    exposure = json.loads(exposure_path.read_text(encoding="utf-8"))
    exposure_members = ROOT / "references/training_exposure_membership_frozen.csv"
    assert sha(exposure_members) == exposure["output_hashes"][exposure_members.name]
    members = pd.read_csv(ROOT / "classification/membership_frozen.csv", encoding="utf-8-sig")
    assert sha(ROOT / "classification/membership_frozen.csv") == audit["membership_sha256"]
    methods = audit["complete_methods"]
    groups = ("protein", "protein_complex")
    chinese = {"protein": "蛋白", "protein_complex": "蛋白复合物"}
    english = {"protein": "proteins", "protein_complex": "protein complexes"}
    counts = audit["group_counts"]

    def row(method, group, metric):
        result = summary[(summary.method == method) & (summary.group == group) & (summary.metric == metric)]
        assert len(result) == 1
        assert int(result.iloc[0].n_records) == counts[group]
        return result.iloc[0]

    def mean(method, group, metric):
        return float(row(method, group, metric)["mean"])

    def interval(method, group, metric):
        result = row(method, group, metric)
        return f"{result['mean']:.3f} [{result.ci_low:.3f}, {result.ci_high:.3f}]"

    own = exposure["own_Full_Intra_OOF"]
    lines = ["# 本机模型横评：蛋白／蛋白复合物", "",
        f"较快阶段已完成14个方法／配置，在相同的全部3,063条记录、961,321个冻结有效RMSF位点上覆盖完整。14项包含本文消融及同一BackFlip版本的两种输入方式，不等于14个独立外部模型。PLANET、RMSF-net全量运行按用户选择留待后补。", "",
        f"本版只按生物学身份分成两组：蛋白{counts['protein']:,}条，蛋白复合物{counts['protein_complex']:,}条。同源多聚体与同一成熟蛋白的拆分结构域归蛋白；不同蛋白或独立异源肽伙伴结合归蛋白复合物。判断针对实际输入，逐条以链序列匹配官方实体，不能由链数或旧source标签直接决定。所有记录均已分类，未决为0，没有丢弃记录。244条TCR–pMHC包含在蛋白复合物中。", "",
        "分组政策和逐条证据先冻结，之后才读取两组的性能结果。沿用已完成的逐记录预测与评分，未重跑模型，未修改目标、掩码、训练折或序列家族。PDB转置/撤销替换等来源冲突按实际序列修复注释对应关系，原始标签和冲突说明保留。", "",
        "每条记录先在同一批有效残基上评分，再对记录等权平均。ρ衡量残基柔性排序；AP衡量筛出参考最高10%柔性残基的能力（保留阈值并列）；MAE为物理RMSF幅度误差。ANM、GNM、FastProtFlex只输出相对柔性，不给Å误差，也不在测试集拟合缩放。区间按序列家族成组重抽样2,000次，再按记录数加权计算，避免把近亲记录当成完全独立样本。", ""]
    draft = ["Local completed-cohort comparison of proteins and protein complexes.",
        f"Fourteen method/configuration entries were evaluated on identical frozen valid residues across 3,063 records (961,321 RMSF positions), comprising {counts['protein']:,} proteins and {counts['protein_complex']:,} protein complexes. Protein identity was assigned from sequence-matched official molecular annotations before binary-group performance aggregation. Homomultimers and split domains of the same mature protein were classified as proteins; assemblies containing different proteins or independent heterologous peptide partners were classified as protein complexes. All 244 TCR–pMHC records are included in the latter group. Existing predictions, targets, masks and per-record scores were retained unchanged."]
    for group in groups:
        comparator = "backflip-2.1 (chains)"
        delta = paired[(paired.group == group) & (paired.comparator == comparator) & (paired.metric == "spearman")].iloc[0]
        external = [m for m in methods if m not in ("Full", "Intra", "Static", "AA")]
        best_rho = max(external, key=lambda m: mean(m, group, "spearman"))
        lines += [f"## {chinese[group]}：n={counts[group]:,}", "",
            f"FlexPPI的平均ρ为{mean('Full', group, 'spearman'):.3f}，BackFlip-2.1逐链为{mean(comparator, group, 'spearman'):.3f}；两者逐记录配对差值为{delta.full_minus_comparator:+.3f}，95%家族成组区间[{delta.ci_low:.3f}, {delta.ci_high:.3f}]。本组最高外部平均ρ来自{best_rho}（{mean(best_rho, group, 'spearman'):.3f}）。AP和幅度误差独立列出，不合成任意总分。", "",
            "| 方法／配置 | ρ [95%区间] ↑ | AP ↑ | MAE (Å) ↓ | 标准化MAE ↓ |",
            "|---|---:|---:|---:|---:|"]
        for method in methods:
            physical = audit["methods"][method]["absolute_angstrom"]
            lines.append(f"| {method} | {interval(method, group, 'spearman')} | {mean(method, group, 'top10_ap'):.3f} | " +
                (f"{mean(method, group, 'mae'):.3f}" if physical else "不适用（相对分数）") + f" | {mean(method, group, 'z_mae'):.3f} |")
        lines.append("")
        draft.append(f"For {english[group]}, FlexPPI achieved mean Spearman correlation {mean('Full', group, 'spearman'):.3f}, high-RMSF average precision {mean('Full', group, 'top10_ap'):.3f} and mean absolute error {mean('Full', group, 'mae'):.3f} Angstroms. Corresponding chain-wise BackFlip-2.1 results were {mean(comparator, group, 'spearman'):.3f}, {mean(comparator, group, 'top10_ap'):.3f} and {mean(comparator, group, 'mae'):.3f} Angstroms. The paired Spearman difference (FlexPPI minus chain-wise BackFlip-2.1) was {delta.full_minus_comparator:.3f} (95% family-cluster bootstrap interval {delta.ci_low:.3f} to {delta.ci_high:.3f}). The highest external mean Spearman correlation in this group was {mean(best_rho, group, 'spearman'):.3f} from {best_rho}. Full results for all completed entries and metrics are provided without an aggregate ranking score.")

    # A rule-boundary sensitivity analysis uses the existing complete scores only.
    # It does not change the primary frozen assignment or select a favorable group.
    edge = members[members.pdb_id.astype(str).str.upper() == "3F4Y"]
    if len(edge):
        per = pd.read_csv(OUT / "per_record_metrics.csv", float_precision="round_trip")
        assert set(edge.benchmark_group) == {"protein"}
        altered = per.copy()
        altered.loc[altered.sample_id.isin(edge.sample_id), "group"] = "protein_complex"
        sensitivity = []
        for group in groups:
            for method in methods:
                for metric in ("spearman", "top10_ap", "mae"):
                    values = altered[(altered.group == group) & (altered.method == method)][metric]
                    if values.notna().all():
                        primary = mean(method, group, metric)
                        sensitivity.append(dict(group=group, method=method, metric=metric,
                            primary_mean=primary, alternative_mean=float(values.mean()),
                            change=float(values.mean()) - primary, alternative_n_records=len(values),
                            moved_records=";".join(edge.sample_id), policy="3F4Y gp41-derived CHR designed peptide treated as independent partner instead of same-protein domain variant; sensitivity only"))
        sensitivity = pd.DataFrame(sensitivity)
        sensitivity.to_csv(OUT / "boundary_sensitivity_3f4y.csv", index=False)
        max_changes = {metric: float(sensitivity.loc[sensitivity.metric == metric, "change"].abs().max())
            for metric in ("spearman", "top10_ap", "mae")}
        lines += [f"边界敏感性：3F4Y的gp41同源派生CHR设计肽按同一成熟蛋白域变体归蛋白，并保留其约50%位点修改和抑制剂用途的事实。若另将该{len(edge)}条记录归为独立肽伙伴复合物，全部方法及两组中，平均ρ的最大绝对变化为{max_changes['spearman']:.6f}，平均AP为{max_changes['top10_ap']:.6f}，物理MAE为{max_changes['mae']:.6f} Å；具体数值见evaluation/boundary_sensitivity_3f4y.csv。该备选检查未改变主分组。", ""]

    lines += ["## 解释与复现", "",
        f"Full/Intra复用冻结五折OOF预测：被评记录本身及同一生物学组均不进入对应训练折，但{own['families_crossing_original_folds']}个序列家族跨原折，{own['per_variant']['full']['records_with_same_family_in_train']}条测试记录仍有同家族成员在训练中。244条TCR–pMHC属于同一依赖家族，分散在原五折；归入复合物不会使它们成为244个独立家族。图中家族区间不能被解读为全新序列家族泛化证明。", "",
        "外部训练名单的精确ID匹配亦不相同：DynaProt-M在ATLAS队列有1,262条训练ID匹配；BackFlip-1/seq各有1,108条Flexpert候选训练名单匹配，原checkpoint CSV未获得；BackFlip-2训练分组未知。公开ID测试子集的同记录敏感性比较保留在上级目录TRAINING_EXPOSURE_SENSITIVITY_zh.md，不在新分组中重新选择有利子集。来源ID修复不回写已冻结训练暴露清单，非ATLAS来源的精确ID检查因此仍有限制。", "",
        "参考RMSF分别去除每条链的刚体平移和转动；外部模型原训练目标可能采用不同对齐方式，这里评价的是向共同参考的迁移表现。整个输入结构和逐链输入分别保留，不能由“能读取多链坐标”推出“经过复合物专门训练”。外部模型未提供同一定义的侧链联合状态熵Hχ，因此本次横评集中于RMSF。", "",
        "FastProtFlex完整覆盖3,063条：3,055条保留官方原版缓存，8条使用已验证只改变中间矩阵释放时机的版本；三个可承受原版对照均数值完全相同。最大体系资源探测与正式输出的身份/有限值检查有记录，不能将没有原版参考的最大体系声称为逐位原版对照。ANM/GNM采用文档化ProDy方程；GNM接触边界修正和简并模态处理原始审计均保留。", "",
        "figures/benchmark_current.pdf为183 mm矢量主图，benchmark_distributions_current.pdf为全部记录ρ的补充分布图；均提供PNG预览和可编辑SVG。evaluation/保留逐记录分数、两组汇总、完整覆盖和配对差值；classification/保留本次冻结成员、政策及证据哈希。旧三组产物仍在上级目录，避免混淆时应使用本目录。QA_REPORT.md只认实际已检查图像的哈希。当前稿未自动替换入Word。", ""]
    (ROOT / "RESULTS_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    draft.append("These results compare released weights and frozen out-of-fold predictions under a common reference, with unequal and partly unresolved training exposure. Full/Intra exclude each evaluated record and its biological group from the corresponding fold, but 904 evaluated records have related-family members in training. The 244 TCR–pMHC records form one dependency family within the protein-complex group, so these results do not establish generalization to unseen sequence families. Intervals resample families while retaining record-weighted means. Relative flexibility scores are not calibrated to Angstroms or substituted for the joint side-chain state entropy Hχ. Full-cohort PLANET and RMSF-net runs remain deferred.")
    FIG.mkdir(exist_ok=True)
    (FIG / "RESULTS_en_draft.txt").write_text("\n\n".join(draft) + "\n", encoding="utf-8")
    caption = (f"Figure. Flexibility prediction for proteins and protein complexes. (a–c) Proteins (n = {counts['protein']:,}); (d–f) protein complexes (n = {counts['protein_complex']:,}). Columns show record-mean Spearman rank correlation, average precision for high-RMSF residues (reference RMSF at or above each record's 90th percentile, retaining ties), and RMSF mean absolute error in Angstroms. All 14 completed method/configuration entries use the same 3,063 records and 961,321 frozen valid RMSF positions, with equal record weights. Homomultimers and split domains of one mature protein are proteins; distinct proteins or independent heterologous peptide partners are protein complexes. Classification was frozen from sequence-matched official annotations before binary-group performance aggregation. Blue, ochre and teal distinguish FlexPPI configurations, published prediction models and simple/graph controls, respectively; marker shapes and direct labels identify configurations. Intervals are 95% record-weighted family-cluster bootstrap intervals (2,000 resamples). The 244 TCR–pMHC records remain one dependency family within the protein-complex group. ANM, GNM and FastProtFlex produce relative scores and are not assigned physical MAE values. Full, Intra, Static and AA reuse frozen five-fold out-of-fold predictions; training exposure differs between methods and related-family exposure remains. The reference removes chain-wise rigid-body motion. Only completed full-cohort entries are plotted, and heavier models are deferred. The companion empirical cumulative distributions include all records without downsampling.\n")
    (FIG / "CAPTION_en.txt").write_text(caption, encoding="utf-8")
    sources = [OUT / "scoring_audit.json", OUT / "group_summary.csv", OUT / "paired_differences.csv", OUT / "per_record_metrics.csv", ROOT / "classification/membership_frozen.csv", exposure_path, exposure_members]
    outputs = [ROOT / "RESULTS_zh.md", FIG / "RESULTS_en_draft.txt", FIG / "CAPTION_en.txt"]
    if (OUT / "boundary_sensitivity_3f4y.csv").exists():
        outputs.append(OUT / "boundary_sensitivity_3f4y.csv")
    report_audit = {"created_utc": datetime.now(timezone.utc).isoformat(), "report_script_sha256": sha(__file__),
        "input_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "outputs": {str(path.relative_to(ROOT)): sha(path) for path in outputs},
        "no_inference": True, "primary_groups": counts, "training_exposure_equalized": False}
    (ROOT / "report_audit.json").write_text(json.dumps(report_audit, indent=2), encoding="utf-8")
    print(json.dumps({"status": "report_written", "primary_groups": counts}), flush=True)


if __name__ == "__main__":
    main()
