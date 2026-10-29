import json
from pathlib import Path

import pytest
import torch

from flexppi.config import ExperimentConfig
from flexppi.data.records import remap_path


def _tiny_config(tmp_path: Path, epochs: int) -> ExperimentConfig:
    config = ExperimentConfig()
    config.stage = "pretrain"
    config.seed = 7
    config.device = "cpu"
    config.output_dir = str(tmp_path)
    config.data.max_residues = 1024
    config.model.metrics = ["rmsf", "chi_entropy"]
    config.model.enable_response_head = False
    config.train.epochs = epochs
    config.train.val_fraction = 0.0
    config.train.sampling_strategy = "uniform_records"
    config.train.log_every = 1
    config.train.checkpoint_every = 1
    return config


def _write_smoke_manifest(tmp_path: Path, count: int = 2):
    import numpy as np

    from flexppi.constants import METRIC_NAMES
    from flexppi.data.labels import save_label_file

    rows = []
    for index in range(count):
        sample = "smoke_{}".format(index)
        structure_dir = tmp_path / "struct" / sample
        structure_dir.mkdir(parents=True, exist_ok=True)
        pdb = structure_dir / (sample + ".pdb")
        lines = []
        for residue in range(24):
            serial = 1 + residue * 4
            lines.append(
                "ATOM  {:5d}  N   GLY A{:3d}     {:8.3f}{:8.3f}{:8.3f}".format(
                    serial, residue + 1, 0.0, float(residue), 0.0
                )
            )
            lines.append(
                "ATOM  {:5d}  CA  GLY A{:3d}     {:8.3f}{:8.3f}{:8.3f}".format(
                    serial + 1, residue + 1, 1.4, float(residue), 0.0
                )
            )
            lines.append(
                "ATOM  {:5d}  C   GLY A{:3d}     {:8.3f}{:8.3f}{:8.3f}".format(
                    serial + 2, residue + 1, 2.5, float(residue), 0.0
                )
            )
            lines.append(
                "ATOM  {:5d}  O   GLY A{:3d}     {:8.3f}{:8.3f}{:8.3f}".format(
                    serial + 3, residue + 1, 3.6, float(residue), 0.0
                )
            )
        lines.append("END")
        pdb.write_text("\n".join(lines), encoding="utf-8")
        label_path = tmp_path / "labels" / (sample + ".npz")
        values = np.full((24, len(METRIC_NAMES)), 1.0, dtype=np.float32)
        values[:, 0] = np.linspace(0.1, 2.0, 24)
        variances = np.zeros_like(values)
        mask = np.ones_like(values, dtype=bool)
        mask[:, 2] = False
        save_label_file(str(label_path), values, variances, mask)
        rows.append(
            {
                "sample_id": sample,
                "dataset": "smoke",
                "structure_path": str(pdb),
                "label_path": str(label_path),
                "trajectory_paths": [],
                "temperature": 300.0,
                "window_ns": 100.0,
                "target_chains": [],
                "source_group": "smoke_absolute",
                "metadata": {
                    "biological_system": "sys_{}".format(index),
                    "cluster_id": "group_{}".format(index),
                    "pdb_id": "TEST{}".format(index),
                    "residue_count": 24,
                },
            }
        )
    manifest = tmp_path / "smoke_train.jsonl"
    with manifest.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")
    test_manifest = tmp_path / "smoke_test.jsonl"
    with test_manifest.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(rows[0]) + "\n")
    return manifest, test_manifest


def test_resume_restores_optimizer_epoch_and_history(tmp_path):
    from flexppi.trainer import train

    manifest, _ = _write_smoke_manifest(tmp_path)
    config = _tiny_config(tmp_path, epochs=2)
    config.data.state_manifest = str(manifest)
    first = train(config)
    assert first["train_samples"] == 2
    history_path = tmp_path / "history.json"
    history = json.loads(history_path.read_text(encoding="utf-8"))
    assert [row["epoch"] for row in history] == [1, 2]
    epoch_one = torch.load(tmp_path / "epoch_001.pt", weights_only=False)
    epoch_two = torch.load(tmp_path / "epoch_002.pt", weights_only=False)
    assert epoch_two["epoch"] == 2
    assert "optimizer" in epoch_two
    assert "rng" in epoch_two

    config.train.epochs = 4
    resumed = train(config)
    history = json.loads(history_path.read_text(encoding="utf-8"))
    assert [row["epoch"] for row in history] == [1, 2, 3, 4]
    assert resumed["train_samples"] == 2
    assert (tmp_path / "epoch_004.pt").exists()
    assert epoch_one["global_step"] == 2


def test_train_refuses_to_resume_finished_run(tmp_path):
    from flexppi.trainer import train

    manifest, _ = _write_smoke_manifest(tmp_path)
    config = _tiny_config(tmp_path, epochs=1)
    config.data.state_manifest = str(manifest)
    train(config)
    with pytest.raises(RuntimeError, match="refusing to resume"):
        train(config)


def test_train_rejects_silent_record_drops(tmp_path):
    from flexppi.trainer import train

    manifest, _ = _write_smoke_manifest(tmp_path)
    extra = tmp_path / "extra.jsonl"
    with manifest.open("r", encoding="utf-8") as source, extra.open("w", encoding="utf-8") as out:
        for line in source:
            out.write(line)
    import json as _json

    broken = _json.loads(extra.read_text(encoding="utf-8").strip().split("\n")[0])
    broken["sample_id"] = "dropped_sample"
    broken["structure_path"] = str(tmp_path / "missing.pdb")
    with extra.open("a", encoding="utf-8") as handle:
        handle.write(_json.dumps(broken) + "\n")
    config = _tiny_config(tmp_path, epochs=1)
    config.data.state_manifest = str(extra)
    with pytest.raises(RuntimeError, match="dropped"):
        train(config)


def _norm(path: str) -> str:
    return path.replace("\\", "/").rstrip("/")


def test_remap_path_is_bidirectional_and_prefers_existing(tmp_path):
    real = tmp_path / "real" / "sub" / "file.txt"
    real.parent.mkdir(parents=True)
    real.write_text("x", encoding="utf-8")
    mappings = {str(tmp_path / "linux"): r"E:\windows"}
    assert remap_path(str(real), mappings) == str(real)
    linux_missing = str(tmp_path / "linux" / "sub" / "file.txt")
    assert _norm(remap_path(linux_missing, mappings)) == "E:/windows/sub/file.txt"
    assert _norm(remap_path(r"E:\windows\sub\file.txt", mappings)) == _norm(linux_missing)
    mappings_with_existing = {str(tmp_path / "linux"): r"E:\windows", r"E:\data": str(tmp_path)}
    assert remap_path(r"E:\data\sub\file.txt", mappings_with_existing) == str(
        tmp_path / "sub" / "file.txt"
    )


def test_evaluate_split_test_without_validation_manifest(tmp_path):
    from flexppi.evaluation import evaluate
    from flexppi.trainer import train

    manifest, test_manifest = _write_smoke_manifest(tmp_path)
    config = _tiny_config(tmp_path, epochs=1)
    config.data.state_manifest = str(manifest)
    config.data.validation_state_manifest = None
    config.data.test_state_manifest = str(test_manifest)
    train(config)
    result = evaluate(config, str(tmp_path / "epoch_001.pt"), split="test")
    assert result["split"]["name"] == "test"
    assert result["split"]["samples_evaluated"] == 1
    assert result["split"]["validation_samples"] == 0
    assert "rmsf" in result["metrics"]
