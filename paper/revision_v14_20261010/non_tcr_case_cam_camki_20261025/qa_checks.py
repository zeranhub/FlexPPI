"""Independent QA checks for the CaM–CaMKI case tables and figure values."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
data = pd.read_csv(HERE / "cam_camki_per_residue.csv", dtype={"author_residue_id": str})
metrics = pd.read_csv(HERE / "case_metrics.csv")
ablation = pd.read_csv(HERE / "work/ablation_aligned.csv", dtype={"author_residue_id": str})
checks = {}

checks["rows_169"] = len(data) == 169
checks["unique_chain_residue"] = not data.duplicated(["chain", "author_residue_id"]).any()
checks["valid_rmsf_169"] = int(data.valid_rmsf.sum()) == 169
checks["valid_chi_145"] = int(data.valid_chi.sum()) == 145
checks["interface_rmsf_43"] = int(data.interface_ca8.fillna(False).sum()) == 43
checks["peptide_25"] = int((data.chain == "B").sum()) == 25
checks["cam_144"] = int((data.chain == "A").sum()) == 144
checks["interface_chi_38"] = int((data.interface_ca8.fillna(False) & data.valid_chi).sum()) == 38

def recompute(subset, metric):
    if subset == "all":
        mask = np.ones(len(data), bool)
    elif subset == "interface_ca8":
        mask = data.interface_ca8.fillna(False).to_numpy(bool)
    elif subset == "outside_interface":
        mask = ~data.interface_ca8.fillna(False).to_numpy(bool)
    elif subset == "peptide_chain_B":
        mask = (data.chain == "B").to_numpy(bool)
    elif subset == "cam_chain_A":
        mask = (data.chain == "A").to_numpy(bool)
    else:
        raise KeyError(subset)
    y = data["target_" + ("rmsf" if metric == "rmsf" else "chi")].to_numpy(float)
    p = data["full_" + ("rmsf" if metric == "rmsf" else "chi")].to_numpy(float)
    q = data["intra_" + ("rmsf" if metric == "rmsf" else "chi")].to_numpy(float)
    good = mask & np.isfinite(y) & np.isfinite(p) & np.isfinite(q)
    y, p, q = y[good], p[good], q[good]
    return dict(n=int(good.sum()), full_mae=float(np.mean(abs(p - y))), intra_mae=float(np.mean(abs(q - y))),
                full_rho=float(spearmanr(y, p).correlation), intra_rho=float(spearmanr(y, q).correlation))

errors = []
for _, row in metrics[metrics.model == "main_oof"].iterrows():
    rec = recompute(row.subset, row.metric)
    for key in ["n", "full_mae", "intra_mae", "full_rho", "intra_rho"]:
        if abs(rec[key] - row[key]) > 1e-6:
            errors.append((row.subset, row.metric, key, rec[key], row[key]))
checks["metrics_recompute_exact"] = not errors

pep = data[data.chain == "B"]
vals = {
    "MD": float(pep.target_rmsf.mean()),
    "Full": float(pep.full_rmsf.mean()),
    "Full_no_partner": float(pep.full_no_partner_rmsf.mean()),
    "Intra": float(pep.intra_rmsf.mean()),
}
figure_values = {"MD": 0.789, "Full": 0.856, "Full_no_partner": 4.219, "Intra": 2.003}
checks["figure_values_match"] = all(abs(vals[k] - figure_values[k]) < 0.001 for k in vals)

pep_ab = pd.read_csv(HERE / "peptide_ablation_summary.csv")
checks["ablation_rows_4"] = len(pep_ab) == 4
checks["ablation_full_mae_0.099"] = abs(float(pep_ab.loc[pep_ab.condition == "Full complex", "mae_vs_md"].iloc[0]) - 0.09857) < 1e-3
checks["ablation_without_mae_3.43"] = abs(float(pep_ab.loc[pep_ab.condition == "Full without partner", "mae_vs_md"].iloc[0]) - 3.42996) < 1e-3

print(pd.Series(checks).to_string())
print("\nvalues", vals)
print("errors", errors[:5])
print("all checks", all(checks.values()))
(HERE / "QA_report.md").write_text(
    "# QA report\n\n"
    + "\n".join(f"- `{k}`: {v}" for k, v in checks.items())
    + "\n\nFigure peptide means: " + ", ".join(f"{k}={v:.3f}" for k, v in vals.items())
    + "\n\nAll checks passed: " + str(all(checks.values())) + "\n",
    encoding="utf-8")
