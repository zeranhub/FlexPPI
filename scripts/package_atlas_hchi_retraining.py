# -*- coding: utf-8 -*-
"""Package validated ATLAS Hchi labels and unchanged v2 fold membership for Linux."""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import sys
import tarfile

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from flexppi.config import load_config
from flexppi.data.records import remap_path

RELEASE = "artifacts/data/atlas_rmsf_hchi_v1"
SPLIT = "data/splits/state_rmsf_hchi_v2_atlas_hchi_v1"
RUN = "artifacts/performance_first_v2_atlas_hchi_v1"
CONFIG = "configs/performance_first_v2_atlas_hchi_v1_fold"
LOCAL_MAP = {
    "/ranze/repos/FlexPPI": str(ROOT),
    "/ranze/data/ATLAS_extracted": "E:/data/protein_md_datasets/ATLAS_extracted",
    "/ranze/data/dynarepo/dynarepo_non_nucleic_full": "D:/data/dynarepo_non_nucleic_full",
    "/ranze/data/ITrace": "E:/data/ITrace",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_sha(path):
    # Python treats CRLF and LF identically; compare source across Windows/Linux
    # after universal-newline conversion, while data hashes remain byte-exact.
    return hashlib.sha256(Path(path).read_text(encoding="utf-8").encode("utf-8")).hexdigest()


def write_text(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(value)


def write_json(path, value):
    write_text(path, json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def read_records(path):
    return [json.loads(line) for line in Path(path).read_text(encoding="utf-8").splitlines() if line]


def key(row):
    return row["dataset"], row["sample_id"]


def remote_path(value, repo):
    normal = str(value).replace("\\", "/")
    mappings = {
        str(ROOT).replace("\\", "/"): repo,
        "/ranze/repos/FlexPPI": repo,
        "E:/data/protein_md_datasets/ATLAS_extracted": "/ranze/data/ATLAS_extracted",
        "D:/data/dynarepo_non_nucleic_full": "/ranze/data/dynarepo/dynarepo_non_nucleic_full",
        "E:/data/ITrace": "/ranze/data/ITrace",
    }
    for prefix, replacement in sorted(mappings.items(), key=lambda item: len(item[0]), reverse=True):
        if normal == prefix or normal.startswith(prefix + "/"):
            normal = replacement + normal[len(prefix):]
            break
    require(normal.startswith("/") and ":" not in normal, "Unmapped Windows path: " + normal)
    return normal


PREFLIGHT = r'''"""Read-only verification of the atlas_hchi_v1 training payload and server inputs."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--path-map", action="append", default=[], metavar="OLD=NEW")
    parser.add_argument("--fold", type=int, choices=range(5))
    parser.add_argument("--check-new-run", action="store_true")
    args = parser.parse_args()
    root = args.repo.resolve()
    mappings = [value.split("=", 1) for value in args.path_map]
    def resolve(text):
        normal = text.replace("\\", "/")
        for old, new in sorted(mappings, key=lambda pair: len(pair[0]), reverse=True):
            if normal == old or normal.startswith(old.rstrip("/") + "/"):
                return Path(new) / normal[len(old):].lstrip("/")
        return Path(text)
    errors = []
    checksum_file = root / "artifacts/data/atlas_rmsf_hchi_v1/PAYLOAD_SHA256SUMS"
    for line in checksum_file.read_text().splitlines():
        expected, name = line.split("  ", 1)
        path = root / name
        if not path.is_file() or digest(path) != expected:
            errors.append("Payload missing or changed: " + name)
    inventory = json.loads((root / "artifacts/data/atlas_rmsf_hchi_v1/training_inputs.json").read_text())
    code_hashes = json.loads((root / "artifacts/data/atlas_rmsf_hchi_v1/source_code_hashes.json").read_text())
    for name, expected in code_hashes.items():
        path = root / name
        actual = hashlib.sha256(path.read_text(encoding="utf-8").encode("utf-8")).hexdigest() if path.is_file() else None
        if actual != expected:
            errors.append("Training source missing or changed; synchronize verified code: " + name)
    for row in inventory:
        for field in ("structure", "label"):
            path = resolve(row[field + "_path"])
            if not path.is_file():
                errors.append("Missing " + field + ": " + str(path))
            elif digest(path) != row[field + "_sha256"]:
                errors.append("Changed " + field + ": " + str(path))
    sys.path.insert(0, str(root / "src"))
    try:
        from flexppi.config import load_config
        for fold in range(5):
            config = load_config(str(root / f"configs/performance_first_v2_atlas_hchi_v1_fold{fold}.yaml"))
            if config.model.metrics != ["rmsf", "chi_entropy"] or config.init_checkpoint is not None:
                errors.append("Unexpected task or checkpoint initialization")
            if config.train.epochs != 40 or config.train.val_fraction != 0:
                errors.append("Unexpected training/validation schedule")
            if config.data.path_mappings:
                errors.append("Unexpected path mappings")
        if args.check_new_run:
            if args.fold is None:
                errors.append("--check-new-run needs --fold")
            else:
                out = root / f"artifacts/performance_first_v2_atlas_hchi_v1/fold_{args.fold}"
                if out.exists() and any(out.iterdir()):
                    errors.append("Run output already populated: " + str(out))
    except Exception as error:
        errors.append(type(error).__name__ + ": " + str(error))
    print(json.dumps({"passed": not errors, "systems": len(inventory),
                      "datasets": dict(Counter(row["dataset"] for row in inventory)),
                      "errors": errors}, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
'''


def package(args):
    release = args.label_release_dir.resolve()
    summary = json.loads((release / "summary.json").read_text(encoding="utf-8"))
    audit = json.loads((release / "release_audit.json").read_text(encoding="utf-8"))
    require(summary["scope"] == "full" and not summary["failures"], "Full extraction is not successful")
    require(summary["requested_systems"] == summary["done"] == summary["complete"] == 1938, "Incomplete ATLAS extraction")
    require(audit["passed"] and audit["summary_sha256"] == sha(release / "summary.json"), "Missing or stale independent audit")
    atlas = {row["sample_id"]: row for row in read_records(summary["manifest"])}
    require(len(atlas) == 1938, "ATLAS source manifest is incomplete")
    repo = args.remote_repo.rstrip("/")
    require(repo.startswith("/") and ":" not in repo, "Use an absolute Linux repository path")
    destination = args.output_dir.resolve()
    payload = destination / "payload"
    require(not payload.exists(), "Output payload exists; choose a new --output-dir to avoid overwriting")
    source = args.source_split_dir
    source_summary = json.loads((source / "summary.json").read_text(encoding="utf-8"))
    require(source_summary["version"] == "state_rmsf_hchi_v2_performance_first",
            "This package preserves the v2 splits used in this task; another split version needs its own run/config identity")
    original = read_records(source / "all_ready.jsonl")
    require(len(original) == len({key(row) for row in original}) == 3068, "Original pool differs from frozen v2")
    require(Counter(row["dataset"] for row in original) == {"atlas": 1929, "dynarepo": 895, "itrace": 244},
            "Unexpected source composition")
    reports = {}
    for sid, row in atlas.items():
        report = json.loads((release / "systems" / f"{sid}.json").read_text(encoding="utf-8"))
        expected_fingerprint = summary.get("quality_mask_fingerprints", {}).get(sid, summary["fingerprint"])
        require(report["status"] == "complete" and report["fingerprint"] == expected_fingerprint, "Incomplete/mixed extraction: " + sid)
        require(sha(row["label_path"]) == report["output_sha256"], "Changed ATLAS label: " + sid)
        reports[sid] = report
    updated, inventory, totals = {}, [], Counter()
    for row in original:
        new = copy.deepcopy(row)
        sid = row["sample_id"]
        new["metadata"]["source_split_version"] = row["metadata"]["split_version"]
        new["metadata"]["split_version"] = PurePosixPath(SPLIT).name
        new["structure_path"] = remote_path(row["structure_path"], repo)
        new["trajectory_paths"] = [remote_path(path, repo) for path in row["trajectory_paths"]]
        if row["dataset"] == "atlas":
            label = Path(atlas[sid]["label_path"])
            new["label_path"] = f"{repo}/{RELEASE}/labels/{sid}.npz"
            new["metadata"].update({
                "hchi_labels": reports[sid]["hchi_labels"], "hchi_replicates": 3,
                "hchi_frames_per_replica": reports[sid]["frames_per_replica"],
                "label_protocol": summary["protocol"], "label_sha256": reports[sid]["output_sha256"],
                "original_label_path": remote_path(row["label_path"], repo),
                "original_label_sha256": reports[sid]["source"]["old_label_sha256"],
            })
            if sid in summary.get("quality_mask_fingerprints", {}):
                new["metadata"]["quality_masked_residues"] = summary["quality_mask_policy"]["exclusions"][sid]
                new["metadata"]["quality_mask_policy"] = summary["quality_mask_policy"]["policy_id"]
        else:
            label = Path(remap_path(row["label_path"], LOCAL_MAP))
            new["label_path"] = remote_path(row["label_path"], repo)
        structure = Path(remap_path(row["structure_path"], LOCAL_MAP))
        require(label.is_file() and structure.is_file(), "Missing local training input: " + sid)
        label_hash = sha(label)
        with np.load(label, allow_pickle=False) as data:
            for index, name in enumerate(data["metric_names"]):
                if name in ("rmsf", "chi_entropy"):
                    totals[name] += int(data["mask"][:, index].sum())
        inventory.append({"sample_id": sid, "dataset": row["dataset"], "fold": row["metadata"]["fold"],
                          "structure_path": new["structure_path"], "structure_sha256": sha(structure),
                          "label_path": new["label_path"], "label_sha256": label_hash})
        updated[key(new)] = new
    for fold in range(5):
        train = read_records(source / f"fold_{fold}_train.jsonl")
        test = read_records(source / f"fold_{fold}_test.jsonl")
        train_ids, test_ids = {key(row) for row in train}, {key(row) for row in test}
        require(len(train) == len(train_ids) and len(test) == len(test_ids), "Duplicate samples in original fold")
        require(not train_ids & test_ids and train_ids | test_ids == set(updated), "Fold membership mismatch")
        require(all(row["metadata"]["fold"] == fold for row in test), "Wrong test-fold identity")
        require(test_ids == {key(row) for row in original if row["metadata"]["fold"] == fold},
                "Test manifest and all-ready assignments disagree")
        for field in ("cluster_id", "biological_system", "complex_signature"):
            require(not ({row["metadata"][field] for row in train} & {row["metadata"][field] for row in test}),
                    f"Original fold leakage: {field}, fold {fold}")
        for role, records in (("train", train), ("test", test)):
            write_text(payload / SPLIT / f"fold_{fold}_{role}.jsonl",
                       "".join(json.dumps(updated[key(row)], sort_keys=True, ensure_ascii=False) + "\n" for row in records))
        config = yaml.safe_load((ROOT / f"configs/performance_first_v2_fold{fold}.yaml").read_text())
        expected_train = copy.deepcopy(config["train"])
        expected_model = copy.deepcopy(config["model"])
        config["output_dir"] = f"{RUN}/fold_{fold}"
        config["init_checkpoint"] = None
        config["data"]["state_manifest"] = f"{SPLIT}/fold_{fold}_train.jsonl"
        config["data"]["test_state_manifest"] = f"{SPLIT}/fold_{fold}_test.jsonl"
        config["data"]["path_mappings"] = {}
        config_file = payload / f"{CONFIG}{fold}.yaml"
        write_text(config_file, yaml.safe_dump(config, sort_keys=False))
        loaded = load_config(str(config_file))
        require(loaded.train.epochs == 40 and loaded.init_checkpoint is None, "Unexpected training config")
        require(config["train"] == expected_train and config["model"] == expected_model,
                "Model or training hyperparameters changed")
    write_text(payload / SPLIT / "all_ready.jsonl", "".join(json.dumps(updated[key(row)], sort_keys=True) + "\n" for row in original))
    release_dir = payload / RELEASE
    remote_atlas = []
    for sid, row in atlas.items():
        converted = copy.deepcopy(row)
        converted["label_path"] = f"{repo}/{RELEASE}/labels/{sid}.npz"
        converted["structure_path"] = remote_path(row["structure_path"], repo)
        converted["trajectory_paths"] = [remote_path(path, repo) for path in row["trajectory_paths"]]
        converted["metadata"]["original_label_path"] = remote_path(row["metadata"]["original_label_path"], repo)
        remote_atlas.append(converted)
    write_text(release_dir / "atlas_states_hchi.jsonl",
               "".join(json.dumps(row, sort_keys=True, ensure_ascii=False) + "\n" for row in remote_atlas))
    for sid, row in atlas.items():
        target = release_dir / "labels" / f"{sid}.npz"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(row["label_path"], target)
        target_report = release_dir / "systems" / f"{sid}.json"
        target_report.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(release / "systems" / f"{sid}.json", target_report)
    for name in ("summary.json", "run_config.json", "release_audit.json", "audit_per_system.csv",
                 "model_loading_smoke_test.json", "strict_extraction_summary.json", "quality_mask_policy.json",
                 "geometry_diagnostic.json", "quality_mask_equivalence_test.json"):
        shutil.copyfile(release / name, release_dir / name)
    write_json(release_dir / "training_inputs.json", inventory)
    new_summary = {"version": PurePosixPath(SPLIT).name, "source_split": source.name, "fold_membership_unchanged": True,
                   "systems": len(original), "sources": dict(Counter(row["dataset"] for row in original)),
                   "all_atlas_labels_in_package": len(atlas), "training_atlas_systems": 1929,
                   "metric_mask_counts": dict(totals), "training_result_status": "not_trained",
                   "quality_masked_residues": summary.get("quality_mask_policy", {}).get("exclusions", {}),
                   "label_protocol": summary["protocol"], "source_split_summary_sha256": sha(source / "summary.json")}
    write_json(payload / SPLIT / "summary.json", new_summary)
    write_text(payload / "scripts/preflight_atlas_hchi_retraining.py", PREFLIGHT)
    launch = f'''#!/usr/bin/env bash
set -euo pipefail
task_fold="@@{{1:?Usage: bash scripts/train_atlas_hchi_fold.sh FOLD_INDEX (0..4)}}"
case "$task_fold" in 0|1|2|3|4) ;; *) echo "Fold must be 0..4" >&2; exit 2;; esac
task_repo={json.dumps(repo)}
cd "$task_repo"
task_python="@@{{FLEXPPI_PYTHON:-$task_repo/.venv-linux/bin/python}}"
"$task_python" scripts/preflight_atlas_hchi_retraining.py --fold "$task_fold" --check-new-run
export PYTHONPATH="$task_repo/src@@{{PYTHONPATH:+:$PYTHONPATH}}"
exec "$task_python" -u -m flexppi.cli train --config "{CONFIG}@@{{task_fold}}.yaml"
'''.replace("@@", "$")
    write_text(payload / "scripts/train_atlas_hchi_fold.sh", launch)
    excluded = sorted(set(atlas) - {row["sample_id"] for row in original if row["dataset"] == "atlas"})
    write_json(release_dir / "atlas_outside_frozen_training_pool.json", {
        "sample_ids": excluded, "count": len(excluded),
        "reason": "Not in the original v2 ready pool; extraction does not change sample eligibility."})
    for name in ("backfill_atlas_hchi.py", "backfill_atlas_hchi_buffered.py", "extract_itrace_hchi.py",
                 "audit_atlas_hchi_release.py", "package_atlas_hchi_retraining.py",
                 "complete_atlas_hchi_quality_masks.py", "diagnose_atlas_hchi_geometry.py"):
        target = release_dir / "provenance_code" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "scripts" / name, target)
    write_json(release_dir / "source_code_hashes.json", {
        str(path.relative_to(ROOT)).replace("\\", "/"): source_sha(path)
        for path in (ROOT / "src").rglob("*.py")})
    write_json(release_dir / "source_code_raw_hashes.json", {
        str(path.relative_to(ROOT)).replace("\\", "/"): sha(path)
        for path in (ROOT / "src").rglob("*.py")})
    readme = f"""# ATLAS Hχ 标签及远程五折重训包

本包补齐全部1,938个ATLAS体系的Hχ标签，其中1,929个属于原v2冻结训练池。其余9个标签一并保存，但不改变原来的入选规则。训练池仍为3,068个体系，包含ATLAS 1,929、DynaRepo 895、ITrace 244；五折成员、分组、seed7、40 epoch、模型及优化参数沿用原v2。新训练输出为 {RUN}/fold_N，不会复用旧checkpoint。当前包只准备数据，尚未重新训练。

每个残基按所有定义侧链χ角的三个120°区间组成联合状态，使用自然对数计算每条100 ns轨迹的Shannon熵，再对三条重复取均值及样本方差。每条轨迹读取1,001帧。RMSF/SASA的原值、方差和mask已逐数组核对保持不变。ALA/GLY等无定义χ角的位置使用NaN及false mask，不能填零作为监督。

质量例外仅有两处：4g6d_B 与 4g94_B 的空链、作者残基编号148、LYS。两者CE–NZ距离在源PDB及三条轨迹中异常，最小镜像不能修复。仅屏蔽这两个残基的Hχ（NaN及false mask），不修改坐标、不调整阈值、不删除蛋白、不改变原RMSF/SASA。两体系仍在原fold2。geometry_diagnostic.json记录六条轨迹的逐对距离诊断，quality_mask_policy.json记录处理规则；strict_extraction_summary.json保留原始1,936成功、2个被拦截的结果，summary.json记录补处理后的完整1,938体系。原1,936个成功NPZ逐字节保持不变。原run_config对应首次严格提取，补处理代码哈希另行保存于quality_mask_policy及两体系fingerprint。

五折训练池当前标签mask计数：RMSF={totals['rmsf']:,}，Hχ={totals['chi_entropy']:,}。这些是NPZ标签mask计数，模型加载时仍会结合几何有效性mask。

## 服务器使用

将 atlas_hchi_retraining_v1.tar.gz 上传到服务器，然后在现有仓库运行：

    cd {repo}
    tar -xzf /实际上传位置/atlas_hchi_retraining_v1.tar.gz
    .venv-linux/bin/python scripts/preflight_atlas_hchi_retraining.py
    bash scripts/train_atlas_hchi_fold.sh 0

将最后的0依次替换为1、2、3、4。也可以按顺序运行：

    for fold in 0 1 2 3 4; do bash scripts/train_atlas_hchi_fold.sh "$fold" || break; done

推荐使用已有tmux会话保持训练进程；本包不会自动登录服务器或启动训练。若Python路径不同，设置 FLEXPPI_PYTHON。若仓库不是 {repo}，请使用打包脚本的 --remote-repo 重新生成，以确保所有绝对路径匹配；外部数据根路径如需变化，在清单中明确改写后重新验收。

本包是现有服务器仓库的增量数据包，不包含约91GB原始轨迹，也不包含已有DynaRepo/ITrace标签和结构。服务器应已有 /ranze/data/ATLAS_extracted、/ranze/data/dynarepo/dynarepo_non_nucleic_full、/ranze/data/ITrace 以及原仓库artifacts/data下DynaRepo/ITrace的NPZ。训练只读取PDB与NPZ，不读取轨迹。预检会对3,068个体系的输入结构和标签逐一核对SHA-256，缺失或变化即返回失败；请先解决预检列出的路径或文件问题。

原ITrace清单中的Windows路径已经转为Linux。提取过程的原始summary/system报告保留本地路径作为来源记录，训练入口使用新五折清单；atlas_states_hchi.jsonl另提供全部1,938个体系的Linux清单。provenance_code保存提取和打包代码供追溯，不覆盖服务器现有src实现；source_code_hashes.json按统一LF换行计算训练源码哈希，预检逐个核对，避免Windows/Linux换行差异导致误报，实际源码不一致时须先同步已经验收的实现。source_code_raw_hashes.json同时保留未经换行转换的原始字节哈希。结构与NPZ标签始终按原始字节核验，不做这种换行转换。

## 文件

{RELEASE}/labels：1,938个增强NPZ；其中values/variances/mask兼容现有加载器，额外保存三重复Hχ、χ角数和残基身份。

{SPLIT}：all_ready及五折train/test清单、更新后的统计。

{CONFIG}0..4.yaml：五个新配置。scripts/preflight_atlas_hchi_retraining.py为只读预检，scripts/train_atlas_hchi_fold.sh启动单折。

PAYLOAD_SHA256SUMS校验包内所有其他文件。archive_sha256.txt校验压缩包本身。本次只修改标签和路径，不把旧v2性能值当成此次重训结果。补入ATLAS Hχ后，每体系两任务平均损失中的任务构成会发生预期变化，应重新完成五折测试评价。
"""
    write_text(payload / "ATLAS_HCHI_RETRAINING_README.md", readme)
    files = sorted(path for path in payload.rglob("*") if path.is_file())
    checksums = "".join(f"{sha(path)}  {path.relative_to(payload).as_posix()}\n" for path in files)
    write_text(release_dir / "PAYLOAD_SHA256SUMS", checksums)
    archive = destination / "atlas_hchi_retraining_v1.tar.gz"
    with tarfile.open(archive, "w:gz", compresslevel=6) as tar:
        for path in sorted(payload.rglob("*")):
            if path.is_file():
                tar.add(path, arcname=path.relative_to(payload).as_posix(), recursive=False)
    write_text(destination / "archive_sha256.txt", f"{sha(archive)}  {archive.name}\n")
    shutil.copyfile(payload / "ATLAS_HCHI_RETRAINING_README.md", destination / "README.md")
    print(json.dumps({"archive": str(archive), "bytes": archive.stat().st_size, **new_summary}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label-release-dir", type=Path, default=ROOT / RELEASE)
    parser.add_argument("--source-split-dir", type=Path, default=ROOT / "data/splits/state_rmsf_hchi_v2_performance_first")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/atlas_hchi_retraining_v1_linux")
    parser.add_argument("--remote-repo", default="/ranze/repos/FlexPPI")
    package(parser.parse_args())


if __name__ == "__main__":
    main()
