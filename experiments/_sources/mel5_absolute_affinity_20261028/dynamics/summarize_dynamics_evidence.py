"""Transparent summary and full-chain sequence audit for the mutation analyses."""
from pathlib import Path
import hashlib
import json
import sys

OUT = Path(__file__).resolve().parent
EXP = OUT.parent.parent
sys.path.insert(0, str(EXP / "label_sensitivity/runtime"))
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

STATES = ("3hg1", "4jff", "4jfd", "4jfe")
PRED = EXP / "tcr_case_v40_20261017/family_predictions/full"
ROLES = {"A": "HLA", "B": "beta2M", "C": "peptide", "D": "TCR_alpha", "E": "TCR_beta"}


def main():
    predictions = {}
    sequences = []
    for pdb in STATES:
        with np.load(PRED / f"itrace_{pdb}.npz", allow_pickle=False) as z:
            p = {k: z[k].copy() for k in z.files}
        predictions[pdb] = p
        for chain in dict.fromkeys(p["chain"].tolist()):
            mask = p["chain"] == chain
            sequences.append(dict(pdb=pdb, chain=chain, role=ROLES[chain],
                                  n_residues=int(mask.sum()),
                                  model_residue_ids=",".join(p["author_residue_id"][mask].astype(str)),
                                  sequence="".join(p["amino_acid"][mask])))
    sequence_frame = pd.DataFrame(sequences)
    sequence_frame.to_csv(OUT / "all_four_states_complete_chain_sequences.csv", index=False)
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-2,
                              open_gap_score=-6, extend_gap_score=-0.5)
    sequence_audit, alignment_rows = [], []
    for b in ("3hg1", "4jfd", "4jfe"):
        for chain, role in ROLES.items():
            pa, pb = predictions["4jff"], predictions[b]
            ia = np.flatnonzero(pa["chain"] == chain)
            ib = np.flatnonzero(pb["chain"] == chain)
            sa = "".join(pa["amino_acid"][ia])
            sb = "".join(pb["amino_acid"][ib])
            aln = aligner.align(sa, sb)[0]
            seen_a, seen_b = set(), set()
            mismatches = []
            for (start_a, stop_a), (start_b, stop_b) in zip(*aln.aligned):
                for a_idx, b_idx in zip(range(start_a, stop_a), range(start_b, stop_b)):
                    a_pos, b_pos = ia[a_idx], ib[b_idx]
                    seen_a.add(a_idx)
                    seen_b.add(b_idx)
                    row = dict(pdb_reference="4jff", pdb_comparison=b, chain=chain, role=role,
                               reference_id=str(pa["author_residue_id"][a_pos]),
                               comparison_id=str(pb["author_residue_id"][b_pos]),
                               reference_aa=sa[a_idx], comparison_aa=sb[b_idx],
                               same_amino_acid=sa[a_idx] == sb[b_idx])
                    alignment_rows.append(row)
                    if not row["same_amino_acid"]:
                        mismatches.append(row)
            missing_a = [dict(id=str(pa["author_residue_id"][ia[i]]), aa=sa[i]) for i in range(len(sa)) if i not in seen_a]
            missing_b = [dict(id=str(pb["author_residue_id"][ib[i]]), aa=sb[i]) for i in range(len(sb)) if i not in seen_b]
            sequence_audit.append(dict(pdb_reference="4jff", pdb_comparison=b, chain=chain, role=role,
                                       reference_length=len(sa), comparison_length=len(sb),
                                       identical_sequence=sa == sb, mismatches=mismatches,
                                       reference_unaligned=missing_a, comparison_unaligned=missing_b,
                                       alignment=str(aln)))
    pd.DataFrame(alignment_rows).to_csv(OUT / "full_chain_alignment_to_4jff.csv", index=False)
    (OUT / "full_chain_sequence_audit.json").write_text(json.dumps(sequence_audit, ensure_ascii=False, indent=2), encoding="utf-8")

    f = pd.read_csv(OUT / "peptide_mutations_complete_per_position.csv")
    s = pd.read_csv(OUT / "peptide_mutations_predefined_region_change_metrics.csv")
    main_scopes = ("all45_CDR", "all10_peptide", "all30_pMHC_interface", "HLA_interface", "CDR2α", "CDR3α", "CDR3β")
    main_summary = s[s.scope.isin(main_scopes) & (s.metric.eq("rmsf") | s.scope.isin(("HLA_interface", "all30_pMHC_interface")))]
    main_summary.to_csv(OUT / "manuscript_relevant_summary_all_prespecified_scopes.csv", index=False)
    # Full eligibility, agreement and replica-separation fields stay visible;
    # sorting these complete rows is descriptive and not an additional benchmark.
    for metric in ("rmsf", "hchi"):
        q = f[f[metric + "_eligible"]].copy()
        q["model_absolute_delta"] = q[metric + "_model_delta"].abs()
        q["MD_absolute_delta"] = q[metric + "_MD_delta"].abs()
        q.sort_values(["contrast", "set", "MD_absolute_delta"], ascending=[True, True, False]).to_csv(OUT / f"all_eligible_positions_sorted_by_MD_{metric}_change.csv", index=False)
    # Existing v125 claims specifically highlighted these HLA residues. Show all,
    # including failures; do not replace them only with newly favorable sites.
    highlights = f[f.set.eq("pMHC_interface") & f.role.eq("mhc") & f.id_before.astype(str).isin(("65", "66", "68", "72", "75"))]
    highlights.to_csv(OUT / "HLA_previously_highlighted_and_R65_validation.csv", index=False)
    loops = f[f.set.eq("CDR") & f.region.isin(("CDR2α", "CDR3β"))]
    loops.to_csv(OUT / "complete_CDR2alpha_and_CDR3beta_mutation_validation.csv", index=False)

    text = """# 肽替换动态证据复核

本分析复用了与 v125 文稿和 v124 图4一致的家族留出完整模型预测，没有重训、重新推理或运行新 MD。四态共12条独立200 ns轨迹均重新读取全部保存帧，逐蛋白链独立拟合Cα后复算RMSF，并以每个χ角固定3个120°区间计算联合侧链状态熵。全部有效位置的均值和样本方差均在3×10⁻⁶内复现原预测包参考标签；图源模型读数逐值在10⁻¹²内一致。每态三个独立模拟才是重复单位，轨迹帧不是独立重复。

I7A 的证据较强。完整45个CDR位置的RMSF变化方向在37位与MD相同，变化排序Spearman相关为0.498。现有文稿已关注的CDR2α保留全部6位，模型和MD的变化方向6/6相同，变化排序相关为0.943；MD区域均值由1.153升至1.427 Å，模型由0.878升至0.989 Å。六位中只有两位的MD均值变化超过其描述性2SE阈值，整个环的三个重复范围仍有重叠，因此不把六位方向一致写成三个重复均显著增高。完整CDR3β十三位的MD区域均值由0.926升至1.194 Å，三个突变态重复的区域均值均高于三个参考态重复；模型只由0.648升至0.671 Å，明显低估变化幅度，且逐位置排序相关只有0.291。

I7A 在同型有效HLA界面18位的Hχ变化方向有14位一致，变化排序相关为0.364。此前文稿突出的位置必须完整核查：Q72模型下降而MD升高；K68均值方向一致但重复间波动大；R75均值下降方向一致。HLA R65提供较强的局部侧链例证，其MD均值由1.146升至2.139 nats，模型由1.412升至1.492 nats，三条I7A轨迹的均值均高于三条参考轨迹；这是完整界面核查后的描述性位点，不可代替完整界面成绩。

G4A 的整体恢复明显有限。完整45个CDR位置的RMSF变化方向有29位相同，但变化排序相关为−0.095；完整肽10位为5/10同向、相关−0.152。CDR2α仍有5/6主链方向一致、相关0.600，但CDR3α的变化排序相关为−0.883。HLA界面18位Hχ只有10/18同向，相关−0.146；此前突出K68与Q72的模型变化方向均与MD均值相反。这些结果应保留，不能把G4A和I7A一并写成成功恢复全部肽替换响应。G4A的七个同型有效肽Hχ变化排序相关为0.821，但其平均变化方向模型与MD相反，不能只报告这个相关系数。

这些是已有不同复合物结构之间的回顾性比较，不是所有输入只改变一个肽残基的严格受控实验。全链核查发现4JFD的TCRα链较4JFF缺少末端F，4JFE的TCRα链首两位为QK而4JFF为KQ；完整差异、序列及对齐见全链审计。CDR和pMHC比较沿用文稿冻结的逐对映射，所有相同残基位置仍可作结构状态比较，但不将所有响应因果归于单个肽替换。处理后轨迹没有显式水分子，不能由本分析验证动态水桥或溶剂重排。

建议主文明确给出G4A/I7A的实验亲和力差异，然后以I7A的完整CDR2α方向和全45位结果作为肽替换动态验证，结合接触全候选分析提供解释线索。G4A与逐位失败情况放补充材料并在主文保留“并非全部肽替换响应均得到一致恢复”的限制。结合态RMSF、Hχ均不能直接等同总结合熵或证明亲和力巨变的因果来源。

主要文件：complete_three_contrasts_per_position.csv保留全部255条位置对照，peptide_mutations_predefined_region_change_metrics.csv给出两种肽替换的所有预定区域成绩，predefined_region_means_per_independent_replica.csv给出每条独立轨迹的区域均值，analysis_audit.json保留原始轨迹、模型、源表哈希和原始标签核对结果。abs(delta)>2SE是描述性过滤，不是显著性检验；3×3跨重复差值不是九条独立轨迹。
"""
    (OUT / "动态证据复核说明.md").write_text(text, encoding="utf-8")
    manifest = {str(p.relative_to(OUT)): dict(bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in OUT.rglob("*") if p.is_file() and p.name != "output_manifest.json"}
    (OUT / "output_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps([x for x in sequence_audit if x["pdb_comparison"] != "3hg1" and not x["identical_sequence"]], ensure_ascii=True, indent=2))
    print("Summary saved", OUT)


if __name__ == "__main__":
    main()
