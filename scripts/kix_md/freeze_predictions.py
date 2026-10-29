"""Freeze predictions on the exact four MD constructs before seeing new MD labels."""
from pathlib import Path
import csv
import datetime as dt
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "src")]
import numpy as np
import torch
from run_apo_holo_exploration import state, infer, digest
from export_atlas_hchi_v1_oof_predictions import load_fold, implementation_hash


def main():
    base = ROOT / "artifacts/kix_md_20260914"
    output = base / "frozen_predictions"
    output.mkdir(exist_ok=True)
    if (output / "FROZEN.json").exists():
        raise FileExistsError("Frozen baseline already exists; do not overwrite after MD starts")
    torch.set_num_threads(6)
    cases = {name: state(dict(id=name, path=str(base / "inputs" / f"{name}.pdb"), chains=list(chains)))
             for name, chains in {"apo": "B", "myb": "AB", "mll": "BC", "both": "ABC"}.items()}
    implementation = implementation_hash()
    rows, checks = [], []
    for fold in range(5):
        start = time.monotonic()
        model, config, path, checksum = load_fold(ROOT, fold, torch.device("cpu"), True)
        checks.append(dict(fold_zero_based=fold, checkpoint=str(path), sha256=checksum))
        for name, case in cases.items():
            values = infer(model, case)
            assert values.shape == (len(case["aa"]), 2) and np.isfinite(values).all()
            for index in range(len(values)):
                rows.append(dict(state=name, fold_zero_based=fold, chain=case["chain"][index],
                                 author_residue_id=case["author"][index], amino_acid=case["aa"][index],
                                 predicted_rmsf_angstrom=float(values[index, 0]),
                                 predicted_hchi_nats=float(values[index, 1])))
        del model
        print(f"Frozen fold {fold}: four states, {time.monotonic()-start:.1f} s", flush=True)
    with (output / "residue_predictions.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    info = dict(created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                checkpoints=checks, implementation_sha256=implementation,
                input_sha256={name: digest(case["spec"]["path"]) for name, case in cases.items()},
                prediction_csv_sha256=digest(output / "residue_predictions.csv"),
                protocol="2AGH model1, MLL839-857, no new MD labels, no refitting",
                model_conditions_used=False, nominal_inference_temperature_k=300,
                new_md_temperature_k=298, inference_duration_ns=100,
                uncertainty="Five checkpoint spread is not independent MD repeat uncertainty.",
                caveat="Near-identical KIX sequence in training; not independent family-held-out validation. Preserve all folds, do not select by new MD agreement.")
    (output / "FROZEN.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(dict(rows=len(rows), output=str(output))), flush=True)


if __name__ == "__main__":
    main()
