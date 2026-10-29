"""Export frozen epoch-40, held-out predictions for the v3 balanced split.

The v2 exporter supplies the unchanged manifest, input, mask, and atomic-cache
checks. This adapter binds those checks to v3 and validates the frozen split
checksums and the checkpoint's training/test manifest identity before inference.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import torch

import export_v2_oof_predictions as shared


ROOT = Path(__file__).resolve().parents[1]
SPLIT_DIRECTORY = "data/splits/state_rmsf_hchi_v3_difficulty_balanced"
CHECKPOINT_DIRECTORY = "artifacts/performance_first_v3_difficulty_balanced"
PROTOCOL = {
    **shared.PROTOCOL,
    "name": "v3_difficulty_balanced_final_epoch40_oof",
    "split_version": "state_rmsf_hchi_v3_difficulty_balanced",
    "checkpoint_directory": CHECKPOINT_DIRECTORY,
    "split_design": "frozen v3; label-distribution balancing after inspection of v2",
}


def audit_frozen_checksums(split_dir: Path) -> None:
    """Refuse altered frozen manifests or audit files before checking membership."""
    required = {"all_ready.jsonl", "leakage_audit.json", "summary.json"}
    required |= {f"fold_{fold}_{role}.jsonl" for fold in range(5)
                 for role in ("train", "test")}
    checked = set()
    for line in (split_dir / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, name = line.split(maxsplit=1)
        name = name.lstrip("*")
        path = (split_dir / name).resolve()
        if path.parent != split_dir.resolve():
            raise ValueError(f"Checksum entry is outside the frozen split: {name}")
        if name in checked or shared.sha256(path) != expected:
            raise ValueError(f"Frozen split checksum mismatch/duplicate: {name}")
        checked.add(name)
    if not required.issubset(checked):
        raise ValueError(f"Missing frozen checksum entries: {sorted(required - checked)}")
    leakage = json.loads((split_dir / "leakage_audit.json").read_text(encoding="utf-8"))
    if not leakage.get("passed") or any(leakage.get("violations", {}).values()):
        raise ValueError("The frozen v3 leakage audit does not pass")


def parse_path_mappings(values: list[str]) -> dict[str, str]:
    mappings = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--path-map must have the form remote_prefix=local_prefix")
        source, target = value.split("=", 1)
        if not source or not target:
            raise ValueError("--path-map requires nonempty remote and local prefixes")
        if source in mappings and mappings[source] != target:
            raise ValueError(f"Conflicting path mappings for {source}")
        mappings[source] = target
    return mappings


def implementation_hash() -> str:
    paths = [Path(__file__).resolve(), Path(shared.__file__).resolve()]
    paths += sorted((ROOT / "src" / "flexppi").rglob("*.py"))
    return shared.digest_json({str(path.relative_to(ROOT)): shared.sha256(path)
                               for path in paths})


def load_fold(root: Path, fold: int, device: torch.device, infer: bool,
              path_mappings: dict[str, str]) -> tuple:
    folder = root / CHECKPOINT_DIRECTORY / f"fold_{fold}"
    checkpoint_path = folder / "epoch_040.pt"
    resolved = json.loads((folder / "resolved_config.json").read_text(encoding="utf-8"))
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if checkpoint.get("epoch") != 40:
        raise ValueError(f"Fold {fold} checkpoint is not epoch 40")
    if checkpoint["config"] != resolved:
        raise ValueError(f"Fold {fold} checkpoint and resolved configs disagree")
    data = resolved["data"]
    for field, role in (("state_manifest", "train"), ("test_state_manifest", "test")):
        expected = f"{SPLIT_DIRECTORY}/fold_{fold}_{role}.jsonl"
        actual = str(data.get(field, "")).replace("\\", "/")
        if actual != expected and not actual.endswith("/" + expected):
            raise ValueError(f"Fold {fold} checkpoint uses an unexpected {field}: {actual}")
    if (data.get("validation_state_manifest") is not None
            or resolved["train"]["val_fraction"] != 0.0
            or resolved["train"]["epochs"] != 40):
        raise ValueError(f"Fold {fold} does not match the frozen 40-epoch, no-validation protocol")
    config = shared.ModelConfig(**checkpoint["config"]["model"])
    if (list(config.metrics) != list(shared.METRICS) or config.use_condition_features
            or config.enable_response_head or config.sequence_backend != "learned"):
        raise ValueError(f"Fold {fold} is not the frozen structure-only two-output model")
    model = None
    if infer:
        model = shared.FlexPPI(config)
        model.load_state_dict(checkpoint["model"], strict=True)
        model.to(device).eval()
    # Mapping only changes where identical inputs are found; original record and
    # input file hashes remain in every cache identity.
    data["path_mappings"] = {**data.get("path_mappings", {}), **path_mappings}
    return model, resolved, checkpoint_path, shared.sha256(checkpoint_path)


def main(argv: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--output-dir", default=str(
        ROOT / "artifacts/evaluation_v3/oof"))
    parser.add_argument("--fold", type=int, choices=range(5), action="append")
    parser.add_argument("--limit", type=int, help="Smoke only: maximum systems per fold")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--path-map", action="append", default=[],
                        help="Repeatable remote_prefix=local_prefix mapping")
    args = parser.parse_args(argv)
    if args.threads < 1 or (args.limit is not None and args.limit < 1):
        parser.error("threads and limit must be positive")
    try:
        mappings = parse_path_mappings(args.path_map)
    except ValueError as error:
        parser.error(str(error))
    audit_frozen_checksums(Path(args.root).resolve() / SPLIT_DIRECTORY)
    shared.SPLIT_DIRECTORY = SPLIT_DIRECTORY
    shared.PROTOCOL = PROTOCOL
    shared.implementation_hash = implementation_hash
    shared.load_fold = lambda root, fold, device, infer: load_fold(
        root, fold, device, infer, mappings)
    shared.run(args)


if __name__ == "__main__":
    main()
