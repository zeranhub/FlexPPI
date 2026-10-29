"""Final-epoch OOF inference with explicit output/label metric mappings."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

from completion_common import ARM_METRICS, LABEL_NAMES, fold_summary, jsonl, sha256, statistics, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--arms", nargs="+", choices=list(ARM_METRICS), default=["rmsf_only", "chi_entropy_only"])
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--dual-guard", type=Path, default=Path(__file__).with_name("dual_checkpoint_guard.json"))
    parser.add_argument("--skip-unfinished", action="store_true", help="Export only complete fits; audit remains partial")
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding="utf-8"))
    training, output = Path(plan["training_root"]), Path(plan["output_root"])
    for relative, expected in {**plan["source_hashes"], **plan["input_file_hashes"]}.items():
        if sha256(training / relative) != expected:
            raise ValueError("Frozen input/source differs: " + relative)
    sys.path.insert(0, str(training / "src"))
    import torch
    from flexppi.config import ModelConfig
    from flexppi.constants import METRIC_NAMES
    from flexppi.data.dataset import StateDataset, collate_states
    from flexppi.data.structures import parse_pdb
    from flexppi.model.information_ablation import build_information_ablation
    from flexppi.trainer import move_to_device
    assert tuple(METRIC_NAMES) == LABEL_NAMES
    torch.set_num_threads(2)
    evaluation = output / "evaluation"
    evaluation.mkdir(parents=True, exist_ok=True)
    rows, predictions, missing, totals = [], [], [], {}
    device = torch.device(args.device)
    dual_guard = json.loads(args.dual_guard.read_text(encoding="utf-8"))["weights"]
    job_index = {(row["arm"], row["fold"]): row for row in plan["jobs"]}
    for arm in args.arms:
        totals[arm] = dict(records=0, rmsf_valid=0, chi_entropy_valid=0, completed_folds=[])
        for fold in range(5):
            path = training / "runs" / "full" / ("fold_{}".format(fold)) / "epoch_040.pt" if arm == "full" else Path(job_index[arm, fold]["final_checkpoint"])
            if not path.exists():
                missing.append(dict(arm=arm, fold=fold, path=str(path)))
                if args.skip_unfinished:
                    continue
                raise FileNotFoundError(str(path))
            checkpoint = torch.load(path, map_location="cpu", weights_only=False)
            config = checkpoint["config"]
            assert checkpoint["epoch"] == 40 and config["seed"] == 7 and config["init_checkpoint"] is None
            assert config["model"]["metrics"] == ARM_METRICS[arm]
            assert config["model"]["information_ablation"] == "full"
            assert config["train"]["sampling_strategy"] == "source_system_balanced"
            model = build_information_ablation(ModelConfig(**config["model"]), "full").to(device).eval()
            model.load_state_dict(checkpoint["model"], strict=True)
            del checkpoint
            checkpoint_hash = sha256(path)
            if arm == "full":
                assert checkpoint_hash == dual_guard["trained_checkpoints/full/fold_{}/epoch_040.pt".format(fold)]
            manifest_path = training / "splits" / ("fold_{}_test.jsonl".format(fold))
            dataset = StateDataset(str(manifest_path), max_residues=1024)
            expected = jsonl(manifest_path)
            assert len(dataset) == len(expected)
            assert [row.sample_id for row in dataset.records] == [row["sample_id"] for row in expected]
            folder = evaluation / arm
            folder.mkdir(exist_ok=True)
            for index, record in enumerate(dataset.records):
                example = dataset[index]
                structure = parse_pdb(record.structure_path)
                batch = move_to_device(collate_states([example]), device)
                with torch.inference_mode():
                    prediction = model.forward_state(batch)
                mean = prediction["mean"][0].cpu().numpy()
                variance = prediction["log_variance"][0].exp().cpu().numpy()
                target = example["label_values"].numpy()
                valid = (example["label_mask"] & example["target_mask"][:, None] & prediction["mask"][0].cpu()[:, None]).numpy()
                assert mean.shape == (len(target), len(ARM_METRICS[arm]))
                assert np.isfinite(mean).all() and np.isfinite(variance).all() and (variance > 0).all()
                cache = folder / (record.sample_id + ".npz")
                np.savez_compressed(cache, mean=mean, model_variance=variance, target=target, valid=valid,
                                    label_variance=example["label_variances"].numpy(),
                                    metric_names=np.asarray(ARM_METRICS[arm]), label_names=np.asarray(LABEL_NAMES),
                                    amino_acid=np.asarray(list(structure.sequence)), chain=np.asarray(structure.chain_ids),
                                    author_residue_id=np.asarray(structure.author_residue_ids), aatype=structure.aatype)
                predictions.append(dict(variant=arm, fold=fold, sample_id=record.sample_id, dataset=record.dataset,
                                        cache_path=str(cache), checkpoint_sha256=checkpoint_hash,
                                        input_sha256=sha256(record.structure_path), label_sha256=sha256(record.label_path),
                                        record=record.to_dict(), metric_names=ARM_METRICS[arm]))
                for output_column, metric in enumerate(ARM_METRICS[arm]):
                    label_column = LABEL_NAMES.index(metric)
                    mask = valid[:, label_column]
                    if mask.any():
                        rows.append(dict(variant=arm, sample_id=record.sample_id, dataset=record.dataset, fold=fold,
                                         metric=metric, **statistics(target[mask, label_column], mean[mask, output_column])))
                    totals[arm][metric + "_valid"] += int(mask.sum())
                totals[arm]["records"] += 1
                if (index + 1) % 50 == 0:
                    print(arm, fold, index + 1, len(dataset), flush=True)
            totals[arm]["completed_folds"].append(fold)
            del model
            if device.type == "cuda":
                torch.cuda.empty_cache()
    pd.DataFrame(rows).to_csv(evaluation / "per_record_metrics.csv", index=False)
    (evaluation / "prediction_manifest.jsonl").write_text("".join(json.dumps(row, ensure_ascii=False) + "\n" for row in predictions), encoding="utf-8")
    frame = pd.DataFrame(rows)
    if len(frame):
        folds, summary = fold_summary(frame, ["mae", "rmse", "bias", "spearman", "pearson"], ["variant", "metric"])
        folds.to_csv(evaluation / "source_fold_macro_metrics.csv", index=False)
        summary.to_csv(evaluation / "five_fold_mean_sample_sd.csv", index=False)
    for arm, counts in totals.items():
        if len(counts["completed_folds"]) == 5:
            assert counts["records"] == 3063
            for metric in ARM_METRICS[arm]:
                expected = 961321 if metric == "rmsf" else 824988
                assert counts[metric + "_valid"] == expected
    write_json(evaluation / "evaluation_audit.json", dict(status="partial" if missing else "complete", arms=args.arms,
               counts=totals, missing_fits=missing, seed=7, epoch=40, plan_sha256=sha256(args.plan),
               export_code_sha256=sha256(__file__), label_mapping="Explicit metric name -> original METRIC_NAMES index",
               completed_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())))


if __name__ == "__main__":
    main()
