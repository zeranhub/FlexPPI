"""Reproduce the CaM–CaMKI case numbers and the partner-removal ablation.

Inputs are the frozen cohort record A001K (1MXE), the main fold-4 checkpoint,
and the existing OOF NPZ. No training and no label changes are performed.

Run from this directory with the project Python:
    python reproduce_case.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from flexppi.inference import predict

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
CORE = Path(r"E:\Project\FlexPPI\experiments\core_corrected_pipeline")
COHORT = Path(r"E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\cohort\records.csv")


def subset_pdb(src: Path, chains: set[str], path: Path) -> Path:
    lines = src.read_text().splitlines()
    keep = [l for l in lines if (not l.startswith(("ATOM  ", "HETATM"))) or l[21:22] in chains]
    path.write_text("\n".join(keep) + "\n")
    return path


def shifted_pdb(src: Path, chains: set[str], dx: float, path: Path) -> Path:
    out = []
    for line in src.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[21:22] in chains:
            x = float(line[30:38]) + dx
            line = line[:30] + f"{x:8.3f}" + line[38:]
        out.append(line)
    path.write_text("\n".join(out) + "\n")
    return path


def main() -> None:
    r = pd.read_csv(COHORT).query("sample_id == 'A001K'").iloc[0]
    src = Path(r.structure_path)
    ckpt = CORE / "trained_checkpoints" / "full" / f"fold_{int(r.original_fold)}" / "epoch_040.pt"
    pep_only = subset_pdb(src, {"B"}, WORK / "A001K_peptide_only.pdb")
    cam_far = shifted_pdb(src, {"A"}, 100.0, WORK / "A001K_cam_shifted_100A.pdb")

    runs = {
        "full_complex": (src, ["A", "B"]),
        "full_peptide_only": (pep_only, ["B"]),
        "full_cam_shifted_100A": (cam_far, ["A", "B"]),
    }
    frames = []
    for tag, (pdb, chains) in runs.items():
        out = predict(checkpoint_path=str(ckpt), apo_pdb=str(pdb), target_chains=chains,
                      temperature=310.0, window_ns=100.0, device_name="cpu")
        d = pd.DataFrame(out["rows"])
        d["variant"] = tag
        frames.append(d)
    pred = pd.concat(frames, ignore_index=True)
    pred["chain"] = pred.chain.astype(str)
    pred["author_residue_id"] = pred.author_residue_id.astype(str)
    pred.to_csv(WORK / "ablation_predictions.csv", index=False)

    npz = np.load(CORE / "neural_exports" / "evaluation" / "full" / "A001K.npz", allow_pickle=True)
    target = pd.DataFrame({
        "chain": npz["chain"].astype(str),
        "author_residue_id": npz["author_residue_id"].astype(str),
        "amino_acid": npz["amino_acid"].astype(str),
        "target_rmsf": npz["target"][:, 0],
        "valid_rmsf": npz["valid"][:, 0],
    })
    pivot = pred.pivot_table(index=["chain", "author_residue_id"], columns="variant", values="rmsf", aggfunc="first").reset_index()
    q = target.merge(pivot, on=["chain", "author_residue_id"], how="left")
    q.to_csv(WORK / "ablation_aligned.csv", index=False)

    peptide = q[(q.chain == "B") & q.valid_rmsf]
    rows = []
    for label, col in [("MD target", "target_rmsf"), ("Full complex", "full_complex"),
                       ("Full without partner", "full_peptide_only"),
                       ("Full, partner shifted 100 A", "full_cam_shifted_100A")]:
        d = peptide[peptide[col].notna()]
        rows.append(dict(condition=label, n=len(d), mean_rmsf=float(d[col].mean()),
                         mae_vs_md=float(np.mean(abs(d[col] - d.target_rmsf))),
                         spearman_vs_md=float(spearmanr(d.target_rmsf, d[col]).correlation)))
    summary = pd.DataFrame(rows)
    summary.to_csv(WORK / "peptide_ablation_summary.csv", index=False)
    print(summary.to_string(index=False))
    print("wrote", WORK)


if __name__ == "__main__":
    main()
