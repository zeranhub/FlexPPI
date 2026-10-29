"""Create ten matched single-target configurations without changing frozen inputs/src."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path, PurePosixPath

import yaml

from completion_common import ARM_METRICS, jsonl, sha256, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--guard", type=Path, default=Path(__file__).with_name("frozen_source_guard.json"))
    parser.add_argument("--runtime-training-root", help="Optional absolute Linux root when generating a transfer package on Windows")
    parser.add_argument("--runtime-output-root", help="Optional absolute Linux completion root when generating a transfer package on Windows")
    args = parser.parse_args()
    training = args.training_root.resolve()
    output = args.output_root.resolve()
    if bool(args.runtime_training_root) != bool(args.runtime_output_root):
        raise ValueError("Both runtime roots must be provided together")
    runtime_training = PurePosixPath(args.runtime_training_root) if args.runtime_training_root else training
    runtime_output = PurePosixPath(args.runtime_output_root) if args.runtime_output_root else output
    if not runtime_training.is_absolute() or not runtime_output.is_absolute():
        raise ValueError("Runtime roots must be absolute")
    if output == training or training in output.parents:
        raise ValueError("Choose a separate completion directory outside the frozen training directory")
    guard = json.loads(args.guard.read_text(encoding="utf-8"))
    for relative, expected in guard["source_hashes"].items():
        if sha256(training / relative) != expected:
            raise ValueError("Frozen source differs: " + relative)
    freeze = json.loads((training / "training_freeze.json").read_text(encoding="utf-8"))
    assert freeze["records"] == 3063 and freeze["valid_counts"] == [961321, 824988]
    assert sha256(training / "training_freeze.json") == guard["freeze_sha256"]
    for relative, expected in freeze["input_hashes"].items():
        normalized = relative.replace("\\", "/")
        if sha256(training / normalized) != expected:
            raise ValueError("Frozen label/structure differs: " + normalized)
    seen = set()
    fold_counts = {}
    file_hashes = {"training_freeze.json": sha256(training / "training_freeze.json")}
    jobs = []
    for fold in range(5):
        train_path = training / "splits" / ("fold_{}_train.jsonl".format(fold))
        test_path = training / "splits" / ("fold_{}_test.jsonl".format(fold))
        train_records, test_records = jsonl(train_path), jsonl(test_path)
        train_ids = {row["sample_id"] for row in train_records}
        test_ids = {row["sample_id"] for row in test_records}
        assert len(train_ids) == len(train_records) and len(test_ids) == len(test_records)
        assert not train_ids & test_ids and len(train_ids | test_ids) == 3063
        assert not seen & test_ids
        seen.update(test_ids)
        fold_counts[str(fold)] = dict(train=len(train_ids), test=len(test_ids))
        for path in [train_path, test_path]:
            file_hashes[path.relative_to(training).as_posix()] = sha256(path)
        base_path = training / "configs" / ("full_fold{}.yaml".format(fold))
        base = yaml.safe_load(base_path.read_text(encoding="utf-8"))
        assert base["seed"] == 7 and base["init_checkpoint"] is None
        assert base["train"]["epochs"] == 40 and base["train"]["sampling_strategy"] == "source_system_balanced"
        assert base["model"]["metrics"] == ARM_METRICS["full"]
        assert base["model"]["information_ablation"] == "full"
        assert base["model"]["sequence_backend"] == "learned"
        assert not base["model"]["use_condition_features"] and not base["model"]["enable_response_head"]
        assert base["train"]["val_fraction"] == 0 and base["data"]["validation_state_manifest"] is None
        file_hashes[base_path.relative_to(training).as_posix()] = sha256(base_path)
        for arm in ["rmsf_only", "chi_entropy_only"]:
            config = copy.deepcopy(base)
            config["model"]["metrics"] = ARM_METRICS[arm]
            config["output_dir"] = str(runtime_output / "runs" / arm / ("fold_{}".format(fold)))
            config["data"]["state_manifest"] = str(runtime_training / "splits" / train_path.name)
            config["data"]["test_state_manifest"] = str(runtime_training / "splits" / test_path.name)
            # Retain every other value, including seed, resume, all original metric weights,
            # data sampling and the from-scratch initialization. Unused weights are ignored.
            path = output / "configs" / ("{}_fold{}.yaml".format(arm, fold))
            path.parent.mkdir(parents=True, exist_ok=True)
            text = yaml.safe_dump(config, allow_unicode=True, sort_keys=False)
            if path.exists() and path.read_text(encoding="utf-8") != text:
                raise ValueError("Existing completion config differs: " + str(path))
            path.write_text(text, encoding="utf-8")
            jobs.append(dict(arm=arm, fold=fold, config=str(runtime_output / "configs" / path.name), config_sha256=sha256(path),
                             output_dir=config["output_dir"], final_checkpoint=str(runtime_output / "runs" / arm / ("fold_{}".format(fold)) / "epoch_040.pt")))
    assert len(seen) == 3063
    # Alternate targets within each fold. Preserve all jobs even if their result is poor.
    jobs.sort(key=lambda row: (row["fold"], row["arm"]))
    plan = dict(protocol="chemical_identity_corrected_v1_single_target_completion_v1", training_root=str(runtime_training),
                output_root=str(runtime_output), frozen_records=3063, seed=7, epochs=40, jobs=jobs,
                fold_counts=fold_counts, input_file_hashes=file_hashes, source_hashes=guard["source_hashes"],
                selection="fixed epoch040 for every arm/fold; no best.pt selection; no test-based exclusion",
                changes="metrics and independent output directory only; manifest paths made absolute without changing bytes")
    write_json(output / "single_target_plan.json", plan)
    print(json.dumps(dict(jobs=len(jobs), output_root=str(output), fold_counts=fold_counts), indent=2))


if __name__ == "__main__":
    main()
