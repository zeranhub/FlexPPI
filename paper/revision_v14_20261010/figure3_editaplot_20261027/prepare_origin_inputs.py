"""Create traceable display tables; never modify the frozen manuscript source CSVs."""
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "source_data"
DEST = ROOT / "origin_inputs"
DEST.mkdir(exist_ok=True)
items = []

def emit(name, frame, source, purpose, template, transformations):
    path = DEST / f"{name}.csv"
    frame.to_csv(path, index=False, encoding="utf-8-sig")
    items.append(dict(id=name, path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                      purpose=purpose, template_id=template, rows=len(frame), columns=list(frame.columns),
                      transformations=transformations))

source = SOURCE / "cam_all25.csv"
cam = pd.read_csv(source)
assert len(cam) == 25 and cam.CaMKI_native_position.tolist() == list(range(294, 319))
emit("b_camki_rmsf", cam[["CaMKI_native_position", "target_rmsf", "pred_rmsf", "intra_rmsf"]].rename(
    columns={"CaMKI_native_position": "CaMKI residue", "target_rmsf": "MD RMSF (A)",
             "pred_rmsf": "Full RMSF (A)", "intra_rmsf": "Intra RMSF (A)"}), source,
    "Test recovery of the continuous W303-M316 constrained CaMKI core, retaining all 25 residues.",
    "trend", ["Rename columns only; source order and all 25 residue values retained.",
              "303-316 shading is defined by structural anchors, not by observed values."])

source = SOURCE / "ubiquitin_72_position_profile_with_original_mask.csv"
ub = pd.read_csv(source)
valid = ub[ub.valid_hchi].copy()
assert len(ub) == 72 and len(valid) == 66
emit("c_ubiquitin_agreement", valid[["MD_Hchi_nats", "Full_Hchi_nats"]].rename(
    columns={"MD_Hchi_nats": "MD Hchi (nats)", "Full_Hchi_nats": "Full Hchi (nats)"}), source,
    "Test agreement of Full and MD Hchi over all 66 valid ubiquitin positions.", "scatter",
    ["Apply the supplied valid_hchi mask: six undefined positions excluded from numeric scatter, retained in frozen source.",
     "Valid zero values (including prolines) retained.", "Equal coordinate ranges and a mathematical y=x guide; no fit."])

source = SOURCE / "ubiquitin_all6_Gln_Hchi_SASA_plot.csv"
gln = pd.read_csv(source).sort_values("MD_Hchi_nats", kind="stable")
emit("d_gln_pairs", pd.DataFrame({"MD Hchi (nats)": gln.MD_Hchi_nats.to_numpy(),
                                 "Full Hchi (nats)": gln.Full_Hchi_nats.to_numpy()}), source,
    "Compare all six same-type sidechains; pair MD and Full and annotate source SASA without regression.",
    "scatter", ["Sort six entire residue rows by MD Hchi, preserve paired residue identity and all values.",
                "Final view uses six ordinal lanes, residue names and aligned sidechain SASA annotations from the frozen source."])

source = SOURCE / "defined_angles_long.csv"
angles = pd.read_csv(source)
for chi, n in [("chi1", 66), ("chi2", 52), ("chi3", 23), ("chi4", 10)]:
    data = angles[angles.chi == chi]
    assert len(data) == n
    emit(f"e_{chi}_association", data[["predicted_Hchi_nats_AA_centered", "fine_entropy_four_origins_mean_nats_AA_centered"]].rename(
        columns={"predicted_Hchi_nats_AA_centered": "AA-centered joint Full Hchi (nats)",
                 "fine_entropy_four_origins_mean_nats_AA_centered": "AA-centered MD chi entropy (nats)"}), source,
        f"Check whether the same joint Hchi prediction retains {chi} MD occupancy differences within residue type.",
        "scatter", [f"Select the predefined {chi} layer; all {n} valid source rows retained.",
                    "Use the already supplied AA-centered coordinates; no new normalization, fitting or outlier removal.",
                    "Q41 and Q2 highlighted by source identity where chi is defined; zero-valued singleton groups retained."])

(DEST / "input_manifest.json").write_text(json.dumps({"status": "awaiting_semantic_confirmation", "items": items,
    "style": {"font": "Arial", "title_pt": 10, "panel_letter_pt": 11, "axis_and_tick_pt": 9,
              "axis_width_pt": 0.6, "Full": {"hex": "#2166AC", "line": "solid"},
              "Intra": {"hex": "#9A633D", "line": "dashed"}},
    "scope": "Non-structural Figure 3 panels only; manuscript and structures remain source context."},
    ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"prepared_tables": len(items), "rows": {x['id']: x['rows'] for x in items}}, ensure_ascii=False))
