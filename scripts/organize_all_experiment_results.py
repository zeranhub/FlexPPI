#!/usr/bin/env python3
"""Collect experiment-generated results into one auditable package.

This script intentionally excludes benchmark input datasets, trajectories,
prepared labels, smoke fixtures, and model checkpoints.  It preserves every
file currently stored below artifacts/evaluation and artifacts/paper, and it
adds the small JSON training histories/configurations needed for provenance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DESTINATION = ROOT / "output" / "all_experiment_results"


@dataclass(frozen=True)
class CopySpec:
    source: Path
    destination: Path
    collection: str


COPY_SPECS = (
    CopySpec(
        ROOT / "artifacts" / "evaluation",
        Path("01_raw_evaluation_outputs"),
        "evaluation_outputs",
    ),
    CopySpec(
        ROOT / "artifacts" / "paper",
        Path("03_paper_tables_figures"),
        "paper_outputs",
    ),
)

TRAINING_RECORDS = (
    (ROOT / "artifacts" / "pretrain" / "history.json", Path("02_training_records/pretrain/history.json")),
    (
        ROOT / "artifacts" / "pretrain" / "resolved_config.json",
        Path("02_training_records/pretrain/resolved_config.json"),
    ),
    (ROOT / "artifacts" / "finetune" / "history.json", Path("02_training_records/finetune/history.json")),
    (
        ROOT / "artifacts" / "finetune" / "resolved_config.json",
        Path("02_training_records/finetune/resolved_config.json"),
    ),
)


EXPERIMENT_DEFINITIONS = (
    {
        "experiment": "aggregate_model_validation",
        "source": "artifacts/evaluation (root files only)",
        "package": "01_raw_evaluation_outputs (root files only)",
        "content": "FlexPPI state/pair validation summaries and FastProtFlex smoke outputs",
        "basis": "ATLAS/DynaRepo evaluation labels where applicable",
        "status": "complete_as_stored",
        "caveat": "Aggregate FlexPPI validation JSON files do not contain per-residue predictions.",
    },
    {
        "experiment": "7L1B_7L1C_flexppi",
        "source": "artifacts/evaluation/7L1B_7L1C_flexppi",
        "package": "01_raw_evaluation_outputs/7L1B_7L1C_flexppi",
        "content": "Downloaded wild-type and mutant PDB structures",
        "basis": "No experimental prediction result was generated locally",
        "status": "input_only_incomplete",
        "caveat": "Only 7L1B.pdb and 7L1C.pdb exist; predictions, colored structures and PyMOL scripts are absent.",
    },
    {
        "experiment": "7L1D_7RRG_tcr_pmhc_flexppi",
        "source": "artifacts/evaluation/7L1D_7RRG_tcr_pmhc_flexppi",
        "package": "01_raw_evaluation_outputs/7L1D_7RRG_tcr_pmhc_flexppi",
        "content": "FlexPPI whole/subcomplex predictions, residue tables, colored PDBs, analyses and PyMOL scripts",
        "basis": "Within-case structural comparison; not MD RMSF ground truth",
        "status": "complete_as_stored",
        "caveat": "PyMOL .pse session files are not present.",
    },
    {
        "experiment": "fastprotflex_compare_smoke",
        "source": "artifacts/evaluation/fastprotflex_compare_smoke",
        "package": "01_raw_evaluation_outputs/fastprotflex_compare_smoke",
        "content": "FastProtFlex comparison smoke-test predictions and timings",
        "basis": "Pipeline smoke test",
        "status": "diagnostic_only",
        "caveat": "Not suitable as paper evidence.",
    },
    {
        "experiment": "fastprotflex_comparison",
        "source": "artifacts/evaluation/fastprotflex_comparison",
        "package": "01_raw_evaluation_outputs/fastprotflex_comparison",
        "content": "ATLAS, DynaRepo, chain-context and TCR3d comparison outputs",
        "basis": "ATLAS/DynaRepo labels; TCR3d structural proxies",
        "status": "complete_with_known_gap",
        "caveat": "chain_context has summary metrics but lacks a separate independent FlexPPI per-residue export.",
    },
    {
        "experiment": "fastprotflex_tcr3d_smoke",
        "source": "artifacts/evaluation/fastprotflex_tcr3d_smoke",
        "package": "01_raw_evaluation_outputs/fastprotflex_tcr3d_smoke",
        "content": "TCR3d comparison smoke-test predictions and timings",
        "basis": "Pipeline smoke test with structural proxies",
        "status": "diagnostic_only",
        "caveat": "TCR3d has no MD trajectory truth.",
    },
    {
        "experiment": "tcr3d_full",
        "source": "artifacts/evaluation/tcr3d_full",
        "package": "01_raw_evaluation_outputs/tcr3d_full",
        "content": "Finetuned-checkpoint TCR3d per-residue predictions and report",
        "basis": "B-factor/CDR structural proxy only",
        "status": "complete_proxy_evaluation",
        "caveat": "Must not be described as MD RMSF validation.",
    },
    {
        "experiment": "tcr3d_smoke",
        "source": "artifacts/evaluation/tcr3d_smoke",
        "package": "01_raw_evaluation_outputs/tcr3d_smoke",
        "content": "FlexPPI TCR3d smoke-test output",
        "basis": "Pipeline smoke test with structural proxies",
        "status": "diagnostic_only",
        "caveat": "Not suitable as paper evidence.",
    },
    {
        "experiment": "tcr3d_static",
        "source": "artifacts/evaluation/tcr3d_static",
        "package": "01_raw_evaluation_outputs/tcr3d_static",
        "content": "Static/pretrain-checkpoint TCR3d per-residue predictions and report",
        "basis": "B-factor/CDR structural proxy only",
        "status": "complete_proxy_evaluation",
        "caveat": "Some legacy column names contain 'flexppi_finetune' although the report identifies the pretrain checkpoint.",
    },
    {
        "experiment": "tcr_pmhc_pymol",
        "source": "artifacts/evaluation/tcr_pmhc_pymol",
        "package": "01_raw_evaluation_outputs/tcr_pmhc_pymol",
        "content": "FlexPPI/FastProtFlex residue predictions, colored PDBs and PyMOL script for 8PJG",
        "basis": "Qualitative structure-coloring comparison",
        "status": "complete_as_stored",
        "caveat": "No MD ground truth and no PyMOL .pse session file.",
    },
    {
        "experiment": "pretrain_training_record",
        "source": "artifacts/pretrain/history.json + resolved_config.json",
        "package": "02_training_records/pretrain",
        "content": "Epoch-level training history and resolved configuration",
        "basis": "Training run logs",
        "status": "complete_as_stored",
        "caveat": "Checkpoint .pt files are intentionally excluded.",
    },
    {
        "experiment": "finetune_training_record",
        "source": "artifacts/finetune/history.json + resolved_config.json",
        "package": "02_training_records/finetune",
        "content": "Epoch-level training history and resolved configuration",
        "basis": "Training run logs",
        "status": "complete_as_stored",
        "caveat": "Checkpoint .pt files are intentionally excluded; FlexPPI-FT is not part of the revised main table.",
    },
    {
        "experiment": "paper_tables_figures",
        "source": "artifacts/paper",
        "package": "03_paper_tables_figures",
        "content": "Publication tables, figures, CSV/JSON/TEX sources and exports",
        "basis": "Derived from evaluation outputs",
        "status": "complete_as_stored",
        "caveat": "The loose files are canonical; the older ZIP may not reflect the latest revised Table 1.",
    },
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify(relative_path: Path) -> str:
    text = relative_path.as_posix().lower()
    name = relative_path.name.lower()
    suffix = relative_path.suffix.lower()
    tags: list[str] = []

    if "smoke" in text:
        tags.append("smoke_test")
    if "7l1b_7l1c_flexppi" in text:
        tags.append("input_only_case")
    if "history.json" in text:
        tags.append("training_history")
    if "resolved_config.json" in text:
        tags.append("run_config")
    if any(token in name for token in ("prediction", "predictions", "residue_predictions")):
        tags.append("prediction")
    if any(token in name for token in ("per_protein", "per_chain", "per_system", "metrics")):
        tags.append("target_level_metrics")
    if "common_residue" in name or "residue_analysis" in name:
        tags.append("aligned_residue_table")
    if "timing" in name:
        tags.append("runtime_measurement")
    if "manifest" in name:
        tags.append("run_manifest")
    if "report" in name or "validation" in name:
        tags.append("summary_report")
    if suffix in {".pdb", ".cif", ".mmcif"}:
        tags.append("structure")
    if suffix in {".pml", ".pse"}:
        tags.append("pymol_visualization")
    if text.startswith("03_paper_tables_figures"):
        tags.append("paper_derived_output")
    if suffix in {".png", ".svg", ".pdf", ".tex"}:
        tags.append("figure_or_table")
    if suffix == ".zip":
        tags.append("archive")
    return ";".join(dict.fromkeys(tags)) if tags else "other_experiment_output"


def source_for_package(relative_path: Path) -> Path:
    parts = relative_path.parts
    if parts[0] == "01_raw_evaluation_outputs":
        return ROOT / "artifacts" / "evaluation" / Path(*parts[1:])
    if parts[0] == "03_paper_tables_figures":
        return ROOT / "artifacts" / "paper" / Path(*parts[1:])
    if parts[0] == "02_training_records":
        return ROOT / "artifacts" / parts[1] / Path(*parts[2:])
    raise ValueError(f"Unknown packaged path: {relative_path}")


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def copy_sources(staging: Path) -> None:
    for spec in COPY_SPECS:
        if not spec.source.is_dir():
            raise FileNotFoundError(spec.source)
        shutil.copytree(spec.source, staging / spec.destination, copy_function=shutil.copy2)
    for source, destination in TRAINING_RECORDS:
        if not source.is_file():
            raise FileNotFoundError(source)
        target = staging / destination
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def build_file_manifest(staging: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(staging.rglob("*")):
        if not path.is_file() or path.parts[-2:-1] == ("00_metadata",):
            continue
        relative = path.relative_to(staging)
        source = source_for_package(relative)
        if not source.is_file():
            raise FileNotFoundError(source)
        if path.stat().st_size != source.stat().st_size:
            raise RuntimeError(f"Size mismatch after copy: {relative}")
        packaged_hash = sha256(path)
        source_hash = sha256(source)
        if packaged_hash != source_hash:
            raise RuntimeError(f"SHA-256 mismatch after copy: {relative}")
        rows.append(
            {
                "package_path": relative.as_posix(),
                "source_path": source.relative_to(ROOT).as_posix(),
                "data_role": classify(relative),
                "bytes": path.stat().st_size,
                "sha256": packaged_hash,
            }
        )
    return rows


def group_stats(manifest: list[dict[str, object]]) -> dict[str, tuple[int, int]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    evaluation_root: list[dict[str, object]] = []
    for row in manifest:
        path = Path(str(row["package_path"]))
        parts = path.parts
        if parts[0] == "01_raw_evaluation_outputs":
            if len(parts) == 2:
                evaluation_root.append(row)
            else:
                grouped[parts[1]].append(row)
        elif parts[0] == "02_training_records":
            grouped[f"{parts[1]}_training_record"].append(row)
        elif parts[0] == "03_paper_tables_figures":
            grouped["paper_tables_figures"].append(row)
    grouped["aggregate_model_validation"] = evaluation_root
    return {
        key: (len(rows), sum(int(row["bytes"]) for row in rows))
        for key, rows in grouped.items()
    }


def write_metadata(staging: Path, manifest: list[dict[str, object]]) -> None:
    metadata = staging / "00_metadata"
    metadata.mkdir(parents=True, exist_ok=True)

    write_csv(
        metadata / "FILE_MANIFEST.csv",
        manifest,
        ["package_path", "source_path", "data_role", "bytes", "sha256"],
    )
    with (metadata / "SHA256SUMS.txt").open("w", encoding="utf-8", newline="\n") as handle:
        for row in manifest:
            handle.write(f"{row['sha256']}  {row['package_path']}\n")

    stats = group_stats(manifest)
    index_rows: list[dict[str, object]] = []
    for item in EXPERIMENT_DEFINITIONS:
        count, size = stats.get(str(item["experiment"]), (0, 0))
        index_rows.append({**item, "file_count": count, "bytes": size})
    write_csv(
        metadata / "EXPERIMENT_INDEX.csv",
        index_rows,
        [
            "experiment",
            "status",
            "file_count",
            "bytes",
            "content",
            "basis",
            "caveat",
            "source",
            "package",
        ],
    )

    role_counts: Counter[str] = Counter()
    role_bytes: Counter[str] = Counter()
    for row in manifest:
        for role in str(row["data_role"]).split(";"):
            role_counts[role] += 1
            role_bytes[role] += int(row["bytes"])
    role_rows = [
        {"data_role": role, "file_count": role_counts[role], "bytes": role_bytes[role]}
        for role in sorted(role_counts)
    ]
    write_csv(metadata / "DATA_ROLE_INDEX.csv", role_rows, ["data_role", "file_count", "bytes"])

    package_metadata = {
        "package_name": "all_experiment_results",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "project_root": str(ROOT),
        "scope": "Experiment-generated results and derived paper outputs only",
        "data_file_count": len(manifest),
        "data_bytes": sum(int(row["bytes"]) for row in manifest),
        "included_sources": [
            "artifacts/evaluation/**",
            "artifacts/pretrain/history.json",
            "artifacts/pretrain/resolved_config.json",
            "artifacts/finetune/history.json",
            "artifacts/finetune/resolved_config.json",
            "artifacts/paper/**",
        ],
        "excluded": [
            "ATLAS, DynaRepo, TCR3d and other source datasets",
            "MD trajectories and prepared labels/manifests in artifacts/data",
            "model checkpoint files (*.pt)",
            "smoke input fixtures in artifacts/smoke_data",
            "duplicate render scratch files in tmp and output/pdf",
            "source code and virtual environments",
        ],
    }
    (metadata / "PACKAGE_METADATA.json").write_text(
        json.dumps(package_metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    total_mib = package_metadata["data_bytes"] / 1024**2
    readme = f"""# FlexPPI 全部实验结果归档

本目录整理的是**实验运行产生的结果及其必要的单案例结构附件**，不包含 ATLAS、DynaRepo、TCR3d 等成套输入数据库，也不包含 MD 轨迹、预处理标签、模型权重或虚拟环境。为保证 PyMOL 案例能够直接打开，原评测目录中随结果存放的少量 PDB 被一并保留；`7L1B_7L1C_flexppi` 的两个 PDB 还用于如实记录该案例目前只有输入、没有预测结果的状态。当前共收录 {len(manifest)} 个数据/结果文件，合计 {total_mib:.2f} MiB。所有源文件均按字节原样复制，并逐一与源文件进行 SHA-256 校验。

## 目录说明

- `00_metadata/FILE_MANIFEST.csv`：每个文件的归档位置、原始位置、数据角色、字节数和 SHA-256。
- `00_metadata/EXPERIMENT_INDEX.csv`：每组实验的文件数、用途、验证依据、完成状态和已知缺口。
- `00_metadata/DATA_ROLE_INDEX.csv`：按预测表、指标表、运行时间、结构可视化等角色统计。
- `00_metadata/SHA256SUMS.txt`：可供外部校验工具直接使用的哈希清单。
- `00_metadata/VERIFICATION_REPORT.json`：本次独立复核的文件数、字节数和哈希检查结论。
- `01_raw_evaluation_outputs/`：所有评测运行的逐残基预测、对齐结果、逐目标指标、报告、耗时和结构可视化文件；保持 `artifacts/evaluation` 原目录结构。
- `02_training_records/`：预训练和微调的 `history.json` 与 `resolved_config.json`，不含 `.pt` 权重。
- `03_paper_tables_figures/`：由实验结果派生的论文表格、图、CSV/JSON/TEX 和导出文件。

## 解释结果时必须保留的限制

1. TCR3d 不含 MD 轨迹，因此 `tcr3d_full`、`tcr3d_static` 和相应 FastProtFlex 比较使用的是 B-factor/CDR 等结构代理指标，不能写成“用 MD RMSF 真值验证”。
2. `tcr3d_static` 的个别列名仍带有历史遗留前缀 `flexppi_finetune`，但该目录报告中记录的实际模型是静态/预训练检查点；分析时应以报告中的 checkpoint 为准。
3. `7L1B_7L1C_flexppi` 当前只有下载的两个静态 PDB，尚未生成 FlexPPI 预测、着色 PDB、逐残基表或 PyMOL 脚本，因此它是“输入已准备、实验结果缺失”，不能算已完成实验。
4. `fastprotflex_comparison/chain_context` 已有汇总指标和 FastProtFlex 逐残基输出，但没有单独导出的 FlexPPI independent-mode 逐残基表；若需对所有方法做逐位点复算，应先补齐该文件。
5. 结构案例中已保存着色 PDB 和 `.pml` 脚本，但当前没有生成 `.pse` 会话文件。
6. 名称含 `smoke` 的文件只用于检查流程是否能运行，不应进入论文主结果或显著性结论。
7. `03_paper_tables_figures` 中散装文件是当前可核验版本；其中旧 ZIP 可能早于最新版 Table 1，引用时应优先使用 `benchmark_figures/` 下的文件。

## 数据层级

“原始实验输出”指模型运行直接得到的 prediction/residue CSV、JSON、耗时、运行清单和结构着色文件；逐蛋白/逐链指标与报告是这些输出的汇总；论文图表属于二次派生结果。归档保留了这三个层级，但没有把评测所依赖的输入数据库再次复制进来。
"""
    (staging / "README.md").write_text(readme, encoding="utf-8")


def verify_manifest(staging: Path, manifest: list[dict[str, object]]) -> None:
    for row in manifest:
        path = staging / str(row["package_path"])
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"Verification size mismatch: {path}")
        if sha256(path) != row["sha256"]:
            raise RuntimeError(f"Verification SHA-256 mismatch: {path}")


def write_verification_report(staging: Path, manifest: list[dict[str, object]]) -> None:
    total_bytes = sum(int(row["bytes"]) for row in manifest)
    report = {
        "verified_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "verification_scope": f"{len(manifest)} packaged data/result files listed in FILE_MANIFEST.csv",
        "source_file_count": len(manifest),
        "manifest_row_count": len(manifest),
        "source_bytes": total_bytes,
        "manifest_bytes": total_bytes,
        "experiment_index_file_sum": len(manifest),
        "experiment_index_byte_sum": total_bytes,
        "sha256_errors": 0,
        "size_errors": 0,
        "result": "PASS",
        "note": "Package metadata files themselves are intentionally outside the data-file checksum manifest to avoid recursive self-hashing.",
    }
    path = staging / "00_metadata" / "VERIFICATION_REPORT.json"
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    destination = args.destination.resolve()
    staging = destination.with_name(destination.name + ".building")
    if destination.exists() or staging.exists():
        raise FileExistsError(
            f"Destination or staging directory already exists: {destination} / {staging}. "
            "Remove it explicitly after reviewing its contents."
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=True)
    try:
        copy_sources(staging)
        manifest = build_file_manifest(staging)
        write_metadata(staging, manifest)
        verify_manifest(staging, manifest)
        write_verification_report(staging, manifest)
        # destination is guaranteed not to exist, so a plain rename avoids
        # Windows' stricter replace-existing permission path for directories.
        staging.rename(destination)
    except Exception:
        # Keep staging on failure so partial work is inspectable and no data is hidden.
        raise
    print(
        json.dumps(
            {
                "destination": str(destination),
                "data_file_count": len(manifest),
                "data_bytes": sum(int(row["bytes"]) for row in manifest),
                "verified": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
