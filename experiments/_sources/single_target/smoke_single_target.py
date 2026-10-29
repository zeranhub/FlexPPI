"""One finite forward/backward per target; no optimizer step or checkpoint write."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from completion_common import ARM_METRICS, sha256, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding="utf-8"))
    training = Path(plan["training_root"])
    for relative, expected in plan["source_hashes"].items():
        assert sha256(training / relative) == expected, relative
    sys.path.insert(0, str(training / "src"))
    import torch
    import numpy as np
    from flexppi.config import load_config
    from flexppi.constants import METRIC_NAMES
    from flexppi.data.dataset import StateDataset, collate_states
    from flexppi.losses import state_loss, heteroscedastic_nll
    from flexppi.model.information_ablation import build_information_ablation
    from flexppi.trainer import seed_everything, move_to_device
    torch.set_num_threads(2)
    device = torch.device(args.device)
    audit = dict(device=str(device), torch=torch.__version__, numpy=np.__version__, cuda_runtime=torch.version.cuda,
                 operation="one forward/backward per target; no optimizer step/checkpoint write", targets=[])
    for arm in ["rmsf_only", "chi_entropy_only"]:
        job = next(row for row in plan["jobs"] if row["arm"] == arm and row["fold"] == 0)
        config = load_config(job["config"])
        dataset = StateDataset(config.data.state_manifest, 1024, config.data.path_mappings)
        # A short, deterministically selected real training record exercises the actual
        # parser and masks while keeping the smoke independent of model performance.
        candidates = sorted((int(row.metadata.get("residue_count", 0)), row.sample_id, index)
                            for index, row in enumerate(dataset.records)
                            if 32 <= int(row.metadata.get("residue_count", 0)) <= 256)
        metric = ARM_METRICS[arm][0]
        label_column = METRIC_NAMES.index(metric)
        example = None
        for _, _, index in candidates:
            candidate = dataset[index]
            if int((candidate["label_mask"][:, label_column] & candidate["target_mask"]).sum()) >= 4:
                example = candidate
                break
        if example is None:
            raise ValueError("No eligible short labelled record for " + metric)
        seed_everything(7)
        if device.type == "cuda":
            torch.cuda.reset_peak_memory_stats(device)
        model = build_information_ablation(config.model, "full").to(device).train()
        batch = move_to_device(collate_states([example]), device)
        output = model(batch)
        assert output["mean"].shape[-1] == 1 and torch.isfinite(output["mean"]).all()
        loss = state_loss(output, batch, [1.0], [metric])
        selected_values = batch["label_values"][..., [label_column]]
        selected_variance = batch["label_variances"][..., [label_column]]
        mask = batch["label_mask"][..., [label_column]] & batch["target_mask"][..., None] & output["mask"][..., None]
        reference = heteroscedastic_nll(output["mean"], output["log_variance"], selected_values, selected_variance, mask, [1.0])
        torch.testing.assert_close(loss, reference)
        assert torch.isfinite(loss)
        loss.backward()
        gradients = [parameter.grad for parameter in model.parameters() if parameter.grad is not None]
        assert gradients and all(torch.isfinite(gradient).all() for gradient in gradients)
        assert any(torch.count_nonzero(gradient) for gradient in gradients)
        record = dict(arm=arm, metric=metric, original_label_column=label_column, output_column=0,
                      sample_id=example["sample_id"], residues=int(example["aatype"].shape[0]), valid_labels=int(mask.sum()),
                      loss=float(loss.detach().cpu()), parameters_with_gradient=len(gradients),
                      parameter_count=sum(parameter.numel() for parameter in model.parameters()))
        if device.type == "cuda":
            record["peak_allocated_bytes"] = int(torch.cuda.max_memory_allocated(device))
            audit["gpu"] = torch.cuda.get_device_name(device)
        audit["targets"].append(record)
        print(json.dumps(record), flush=True)
        del model, batch, output, loss, reference, gradients
        if device.type == "cuda":
            torch.cuda.empty_cache()
    write_json(Path(plan["output_root"]) / "smoke_audit.json", audit)


if __name__ == "__main__":
    main()
