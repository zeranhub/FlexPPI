"""Create exact, source-bound manuscript replacements after metrics exist.

This script only writes a JSON proposal. It never edits the DOCX.
"""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
from lxml import etree

ROOT = Path("E:/Project/FlexPPI/experiments/_sources/predictive_figure_v44_20261018")
SOURCE = ROOT.parent.parent / "paper_v43_20261018_中文修订_能力与构象合并图_三线表.docx"
SOURCE_SHA = "5d32c2af90a6aa8268ac0ef50ba46480370fc798bcc890d91fa71123dfb51c5a"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def main():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    with ZipFile(SOURCE) as z:
        tree = etree.fromstring(z.read("word/document.xml"))
    paragraphs = ["".join(p.xpath(".//w:t/text()", namespaces=NS)) for p in tree.xpath("//w:body/w:p", namespaces=NS)]
    summary_path = ROOT / "analysis/summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    metrics = summary["metrics"]
    rmsf_key = next(k for k in metrics if k.lower() == "rmsf")
    hchi_key = next(k for k in metrics if k.lower() in {"hchi", "hχ", "chi_entropy"})
    r, h = metrics[rmsf_key], metrics[hchi_key]
    full_key = next(k for k in r if k.lower() == "full")
    intra_key = next(k for k in r if k.lower() in {"intra", "intrachain"})
    static_key = next(k for k in r if k.lower() == "static")
    aa_key = next(k for k in r if k.lower() in {"aa", "aa_only"})

    def locate(prefix):
        found = [p for p in paragraphs if p.startswith(prefix)]
        assert len(found) == 1, (prefix, len(found))
        return found[0]

    replacements = []

    def replace(prefix, new):
        replacements.append({"prefix": prefix, "old": locate(prefix), "new": new})

    old = locate("在3,063条记录的五折体系留出评价中，FlexPPI")
    new = old.replace("（表1，图2）。", "（表1）。预测与MD参考的逐位置对应关系见图2a、b。")
    replace("在3,063条记录的五折体系留出评价中，FlexPPI", new)

    old = locate("完整结构编码相对仅链内编码")
    removed = "图2e、f进一步按实际输入链数展示同一记录的误差差值，以区分单链输入和具有伙伴链的结构。"
    assert removed in old
    replace("完整结构编码相对仅链内编码", old.replace(removed, ""))

    old = locate("同类型残基比较检验模型")
    replace("同类型残基比较检验模型", old.replace("图2c、d", "图2e、f"))
    old = locate("氨基酸类型均值对照的整体Hχ相关已达")
    replace("氨基酸类型均值对照的整体Hχ相关已达", old.replace("图2d", "图2f"))
    old = locate("RMSF方面，完整模型在20种氨基酸类型")
    replace("RMSF方面，完整模型在20种氨基酸类型", old.replace("平均误差。", "平均误差（图2e）。"))

    method = (
        "高动态位点筛选以每条记录MD参考的第90百分位及以上为阳性，四种方法在相同有效位置计分；并列参考全部纳入，并列预测作为同一阈值组。"
        "平均精确率（AP）按召回率增量乘对应精确率求和；ROC曲线下面积（AUROC）为所有阳性–阴性位点对中，阳性预测值更高的比例，并列计半分。两指标逐记录等权平均。"
        "精确率–召回率（PR）曲线在间隔0.01的共同召回率上，对首次达到该召回率的阈值组精确率取记录均值，曲线面积不代替平均AP。"
        "单类记录不计AP与AUROC，保留覆盖统计；两指标的95%区间按来源分层重抽记录10,000次。"
        "Hχ另在同记录同类型至少4个位点内重定第90百分位，先对有效类型、再对记录等权平均。"
        "预测–参考密度及分箱预测均值对记录等权，各记录内位置等权。"
    )
    assert len(method) <= 350, len(method)
    pr = (
        "进一步将连续预测用于高动态位点筛选，以每条记录MD参考第90百分位及以上的位置为目标集合。"
        f"完整模型对RMSF与Hχ的平均AP分别为{r[full_key]['mean_ap']:.3f}"
        f"（95%区间{r[full_key]['ap_ci'][0]:.3f}–{r[full_key]['ap_ci'][1]:.3f}）和"
        f"{h[full_key]['mean_ap']:.3f}（{h[full_key]['ap_ci'][0]:.3f}–{h[full_key]['ap_ci'][1]:.3f}），"
        f"对应AUROC为{r[full_key]['mean_auroc']:.3f}和{h[full_key]['mean_auroc']:.3f}（图2c、d）。"
        f"链内模型的两项AP为{r[intra_key]['mean_ap']:.3f}和{h[intra_key]['mean_ap']:.3f}，"
        f"静态特征回归为{r[static_key]['mean_ap']:.3f}和{h[static_key]['mean_ap']:.3f}，"
        f"类型均值为{r[aa_key]['mean_ap']:.3f}和{h[aa_key]['mean_ap']:.3f}。"
        "完整模型在两个目标上的AP均高于三项对照，支持按模型读数优先筛选高动态位置。"
    )
    insertions = [
        {"anchor_prefix": "评价先在每条记录的有效残基上计算", "position": "after", "text": method},
        {"anchor_prefix": "在3,063条记录的五折体系留出评价中，FlexPPI", "position": "after", "text": pr},
    ]
    type_metrics = summary["within_type"]
    within_type = (
        "在同一记录、同一种氨基酸内重新定义高熵位置后，"
        f"完整模型的平均AP为{type_metrics[full_key]['mean_ap']:.3f}，"
        f"高于链内模型的{type_metrics[intra_key]['mean_ap']:.3f}、"
        f"静态特征回归的{type_metrics[static_key]['mean_ap']:.3f}及"
        f"类型均值对照的{type_metrics[aa_key]['mean_ap']:.3f}；"
        f"完整模型的AUROC为{type_metrics[full_key]['mean_auroc']:.3f}。"
        f"该评价覆盖{type_metrics[full_key]['n_records']:,}条记录中的"
        f"{type_metrics[full_key]['n_eligible_record_types']:,}个记录–类型组，"
        "说明结构信息能够在固定化学种类后进一步区分高熵位置。"
    )
    insertions.append({"anchor_prefix": "氨基酸类型均值对照的整体Hχ相关已达", "position": "after", "text": within_type})
    proposal = {
        "source_sha256": SOURCE_SHA,
        "metrics_summary": str(summary_path),
        "replacements": replacements,
        "insertions": insertions,
        "notes": [
            "Caption is owned by the figure worker; do not consume from this proposal.",
            "Do not change any original table, equation, reference, comment, other figure, or biological background.",
            "No AP added to abstract: the existing abstract already covers the main quantitative evidence.",
            "Within-type auxiliary AP/AUROC is explicitly a separate within-record/type q90 task; no plot call claims it is the whole-record PR curve.",
        ],
    }
    out = ROOT / "text_review/manuscript_proposal.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(proposal, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"path": str(out), "method_characters": len(method), "result_paragraph": pr}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
