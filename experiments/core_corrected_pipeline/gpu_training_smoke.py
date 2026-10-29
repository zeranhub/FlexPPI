"""One real labelled training step, discarded; measure memory without a checkpoint."""
from pathlib import Path
import sys
import json
import time
import hashlib
import torch
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from flexppi.config import load_config
from flexppi.data.dataset import StateDataset, collate_states
from flexppi.trainer import move_to_device, _batch_loss, seed_everything
from flexppi.model.information_ablation import build_information_ablation

root = Path(__file__).resolve().parent
config = load_config(str(root / "smoke_config.yaml"))
seed_everything(config.seed)
torch.set_num_threads(2)
device = torch.device("cuda")
dataset = StateDataset(config.data.state_manifest, max_residues=1024)
assert len(dataset) == 1
example = dataset[0]
batch = move_to_device(collate_states([example]), device)
model = build_information_ablation(config.model, "full").to(device).train()
optimizer = torch.optim.AdamW(model.parameters(), lr=config.train.learning_rate, weight_decay=config.train.weight_decay)
torch.cuda.reset_peak_memory_stats()
started = time.monotonic()
loss = _batch_loss(config, model, batch)["loss"]
assert torch.isfinite(loss)
loss.backward()
torch.nn.utils.clip_grad_norm_(model.parameters(), config.train.grad_clip)
optimizer.step()
torch.cuda.synchronize()
report = dict(device=torch.cuda.get_device_name(), torch_version=torch.__version__, length=len(example["aatype"]),
              sample_id=example["sample_id"], real_labelled_sample=True, dtype="float32", optimizer_steps=1,
              discarded_model=True, checkpoints_written=0, finite_loss=float(loss.detach()),
              max_memory_allocated_bytes=torch.cuda.max_memory_allocated(),
              max_memory_reserved_bytes=torch.cuda.max_memory_reserved(), seconds=time.monotonic()-started,
              purpose="Resource feasibility only. This model is discarded, not a scientific training run.")
(root / "gpu_training_smoke.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print(json.dumps(report,indent=2))
