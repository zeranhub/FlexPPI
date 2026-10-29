"""Freeze the scientific proposal without recording user approval or calling Origin."""
from pathlib import Path
import json
import hashlib

ROOT = Path(__file__).resolve().parent
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

panels = [
    dict(panel="b", purpose="检验完整复合物环境是否恢复CaMKI连续调控核心的低RMSF。",
         primary="全部25位点的原生编号294–318、MD/Full/Intra RMSF；303–316为预先定义的结构核心。",
         elements="三条全序列绝对值曲线；可附同25位点MD/Full细节，两个视图均不裁剪观测。",
         source="cam_all25.csv", n=25),
    dict(panel="c", purpose="检验泛素全部有效侧链位置的预测幅值与MD是否一致。",
         primary="MD Hχ、Full Hχ；66有效侧链，原72位置清单和有效掩码保留。",
         elements="等物理比例的一致性散点、数学y=x线、MAE和Spearman；保留合法零值，标Q40/Q41/Q62。",
         source="ubiquitin_72_position_profile_with_original_mask.csv", n=66),
    dict(panel="d", purpose="检验同为Gln的六个位置能否恢复不同的局部侧链多样性。",
         primary="全部6个Gln的MD/Full联合Hχ及起始复合物侧链SASA。",
         elements="按MD Hχ排序的六组成对点，全部位点具名、右侧对齐SASA数值；不作SASA回归。",
         source="ubiquitin_all6_Gln_Hchi_SASA_plot.csv", n=6),
    dict(panel="e", purpose="检验同一个联合Hχ预测是否保留实际MD单角占据差异。",
         primary="原表已经按氨基酸类型去均值的联合Full Hχ及逐χ角细箱熵。",
         elements="四幅散点保留χ1/2/3/4的66/52/23/10点，零参考线、ρ和n；Q41/Q2特殊点形，不作拟合。",
         source="defined_angles_long.csv", n=[66,52,23,10]),
    dict(panel="f", purpose="检验侧链Hχ相对主链RMSF的占据关联优势是否跨记录存在。",
         primary="83记录内Full Hχ与Full RMSF对两个MD指标的同类型相关。",
         elements="两指标各83个配对相关差Δρ点、零线、记录均值及来源分层重抽10000次的95%区间；来源用点形区分。",
         source="retained_absolute_correlations_all498.csv", n=83,
         explicit_analysis="对同一记录作Full Hχ相关减Full RMSF相关；均值/CI保持记录等权，seed=20261026，来源ATLAS/DynaRepo/ITrace=30/23/30。"),
]
proposals = []
for file in sorted((ROOT / "origin_inputs").glob("*_understanding.json")):
    obj=json.loads(file.read_text(encoding="utf-8"))
    assert obj["confirmation_gate"]["can_confirm_now"]
    proposals.append(dict(path=str(file), file_sha256=digest(file),
        proposal_hash=obj["understanding"]["proposal_hash"],
        source_sha256=obj["source"]["sha256"],
        column_decisions=obj["understanding"]["column_decisions"],
        status="awaiting_user_confirmation"))
contract = dict(status="awaiting_user_confirmation", purpose="用完整源数据重设计图3全部非结构子图，区分模型准确性、局部环境差异与MD占据关联。",
    panels=panels, editaplot_semantic_proposals=proposals,
    source_hashes={f.name:digest(f) for f in sorted((ROOT/"source_data").glob("*.csv"))},
    visual_contract=dict(font="Arial", title_pt=10, panel_letter_pt=11, axis_and_tick_pt=9, axis_width_pt=0.6,
                         Full=dict(color_hex="#2166AC", line="solid"), Intra=dict(color_hex="#9A633D", line="dashed"),
                         MD="black open symbols; dotted line where residue order is meaningful"),
    origin_route=dict(kind="task_specific_native_Origin_objects_using_EditaPlot_session_and_verification",
        reason="The public CLI exact visual-style override currently covers XPS; this task requires ordinary line/scatter typography and line styles.",
        claim="Only report an applied field after native object readback and exported visual QA pass."),
    support_and_retained="原始结构、掩码、编号、分组和来源仅用于身份/验证/标注；其余原始列完整保存，不增加曲线。Intra Hχ跨记录完整数据保留，不补造泛素Intra。",
    exclusions=["No smoothing, regression, outlier deletion or new experimental measurements.",
                "Full/Intra are trained models using the same bound input, not measured bound/free states.",
                "No MD resimulation, model refitting, generated conformations or binding free-energy prediction.",
                "Original DOCX and structure panels are not modified by this figure design task."],
    requested_outputs=["editable Origin OPJU", "vector PDF", "PNG", "TIF", "source data and native object QA"],
    approval_record=None)
out=ROOT/"scientific_design_proposal.json"
out.write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding="utf-8")
(ROOT/"scientific_design_proposal.sha256").write_text(digest(out)+"\n",encoding="ascii")
print(json.dumps(dict(proposal_sha256=digest(out), panels=[x['panel'] for x in panels], semantic_proposals=len(proposals)),ensure_ascii=False))
