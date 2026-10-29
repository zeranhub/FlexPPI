"""Read-only CPU audit of completed checkpoint tensors and a fresh backward pass."""
from pathlib import Path
import hashlib
import json
import math
import sys

TRAIN = Path("/root/flexppi_revision_20261002/corrected_training_v1")
COMPLETION = Path("/root/flexppi_completion_20261010/single_target")
sys.path.insert(0, str(TRAIN / "src"))

import torch
torch.set_num_threads(2)
from flexppi.config import ModelConfig
from flexppi.data.dataset import StateDataset, collate_states
from flexppi.model.information_ablation import build_information_ablation
from flexppi.losses import state_loss


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


guard = json.loads((COMPLETION / "frozen_source_guard.json").read_text())
source_checks = {name: sha(TRAIN / name) == value for name, value in guard["source_hashes"].items()}
plan = json.loads((COMPLETION / "single_target_plan.json").read_text())
split_checks = {name: sha(TRAIN / name) == value for name, value in plan["input_file_hashes"].items()}
result = dict(source_files_checked=len(source_checks), all_source_hashes_match=all(source_checks.values()),
              source_mismatches=[name for name, ok in source_checks.items() if not ok],
              frozen_input_hashes_match=all(split_checks.values()), input_mismatches=[name for name, ok in split_checks.items() if not ok],
              device_used="cpu", no_optimizer_step=True, no_checkpoint_write=True, no_gpu_use=True, checkpoints=[])
folder = COMPLETION / "runs" / "chi_entropy_only" / "fold_0"
latest = None
for epoch in range(1, 6):
    path = folder / ("epoch_{:03d}.pt".format(epoch))
    checkpoint = torch.load(path, map_location="cpu", weights_only=False)
    config = checkpoint["config"]
    weights = checkpoint["model"]
    weight_tensor_count = len(weights)
    nonfinite_weights = [name for name, tensor in weights.items() if tensor.is_floating_point() and not bool(torch.isfinite(tensor).all())]
    optimizer_tensors = []
    for state in checkpoint["optimizer"]["state"].values():
        optimizer_tensors += [value for value in state.values() if torch.is_tensor(value) and value.is_floating_point()]
    optimizer_finite = all(bool(torch.isfinite(tensor).all()) for tensor in optimizer_tensors)
    result["checkpoints"].append(dict(epoch=checkpoint["epoch"], global_step=checkpoint["global_step"],
      path=str(path), sha256=sha(path), train_loss=float(checkpoint["train_loss"]),
      train_loss_finite=math.isfinite(checkpoint["train_loss"]), val_loss_is_nan=math.isnan(checkpoint["val_loss"]),
      model_tensors=weight_tensor_count, nonfinite_model_tensor_names=nonfinite_weights,
      optimizer_float_tensors=len(optimizer_tensors), optimizer_tensors_all_finite=optimizer_finite,
      all_history_train_losses_finite=all(math.isfinite(row["train_loss"]) for row in checkpoint["history"])))
    latest = checkpoint
assert latest is not None
config = latest["config"]
result["configuration"] = dict(epochs=config["train"]["epochs"], val_fraction=config["train"]["val_fraction"],
 validation_state_manifest=config["data"]["validation_state_manifest"], metrics=config["model"]["metrics"],
 state_manifest=config["data"]["state_manifest"], test_state_manifest=config["data"]["test_state_manifest"], seed=config["seed"],
 init_checkpoint=config["init_checkpoint"], sampling_strategy=config["train"]["sampling_strategy"])
dataset = StateDataset(config["data"]["state_manifest"], 1024)
sample_index = next(index for index, row in enumerate(dataset.records) if row.sample_id == "3tvj_I")
example = dataset[sample_index]
model = build_information_ablation(ModelConfig(**config["model"]), "full").cpu().train()
model.load_state_dict(latest["model"], strict=True)
batch = collate_states([example])
output = model(batch)
loss = state_loss(output, batch, [1.0], ["chi_entropy"])
loss.backward()
gradients = [parameter.grad for parameter in model.parameters() if parameter.grad is not None]
result["latest_checkpoint_cpu_forward_backward"] = dict(epoch=latest["epoch"], sample_id=example["sample_id"],
 residues=int(example["aatype"].shape[0]), mean_all_finite=bool(torch.isfinite(output["mean"]).all()),
 log_variance_all_finite=bool(torch.isfinite(output["log_variance"]).all()), loss=float(loss.detach()),
 loss_finite=bool(torch.isfinite(loss)), gradient_tensors=len(gradients),
 all_gradients_finite=all(bool(torch.isfinite(gradient).all()) for gradient in gradients),
 any_gradient_nonzero=any(bool(torch.count_nonzero(gradient)) for gradient in gradients))
result["prior_gpu_smoke"] = json.loads((COMPLETION / "smoke_audit.json").read_text())
print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
