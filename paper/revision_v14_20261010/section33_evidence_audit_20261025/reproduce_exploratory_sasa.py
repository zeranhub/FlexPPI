"""Exploratory one-frame SASA analysis; does not modify manuscript/source files.

Runtime: E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe
All results are descriptive, not a calibrated model-vs-baseline evaluation.
"""
import csv
import json
from collections import Counter
from pathlib import Path

import mdtraj as md
import numpy as np
from scipy.stats import spearmanr

PDB = Path('E:/Project/FlexPPI/paper/figures_v59/03_KIX与泛素局部动态验证/source_archive/ubiquitin/inputs/A003R_input_structure.pdb')
SOURCE = Path('E:/Project/FlexPPI/experiments/_case_validation_sources/data/ubiquitin_fold0_all_residues.csv')
OUTPUT = Path(__file__).resolve().parent / 'reproduced'
OUTPUT.mkdir(exist_ok=True)

trajectory = md.load_pdb(str(PDB))
trajectory = trajectory.atom_slice([atom.index for atom in trajectory.topology.atoms if atom.element.symbol != 'H'])
assert trajectory.n_atoms == 935
# MDTraj coordinates and probe radius are in nm; output nm^2 is converted to A^2.
areas = md.shrake_rupley(trajectory, probe_radius=0.14, n_sphere_points=960, mode='atom')[0] * 100
sidechain_areas = {}
for atom, area in zip(trajectory.topology.atoms, areas):
    key = (atom.residue.chain.chain_id, atom.residue.resSeq)
    if atom.name not in ('N', 'CA', 'C', 'O', 'OXT'):
        sidechain_areas[key] = sidechain_areas.get(key, 0.0) + float(area)

with SOURCE.open(encoding='utf-8-sig') as stream:
    records = [record for record in csv.DictReader(stream) if record['chain'] == 'A']
assert len(records) == 72
records = [record for record in records if record['valid_hchi'] == 'True']
assert len(records) == 66
rows = []
for record in records:
    rows.append(dict(
        residue=record['amino_acid'] + record['residue_id'],
        sidechain_SASA_A2=sidechain_areas.get((record['chain'], int(record['residue_id'])), 0.0),
        md_hchi=float(record['md_hchi']),
        pred_hchi=float(record['pred_hchi']),
        md_hchi_minus_type_mean=float(record['md_hchi_minus_type_mean']),
        pred_hchi_minus_type_mean=float(record['pred_hchi_minus_type_mean']),
    ))
with (OUTPUT / 'all66_exploratory_initial_sidechain_sasa.csv').open('w', newline='', encoding='utf-8') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

amino_acids = np.array([row['residue'][0] for row in rows])
counts = Counter(amino_acids)
sasa = np.array([row['sidechain_SASA_A2'] for row in rows])
truth = np.array([row['md_hchi'] for row in rows])
prediction = np.array([row['pred_hchi'] for row in rows])
results = []
for name, keep in (
    ('all66_raw', np.ones(len(rows), bool)),
    ('all66_within_type_centered', np.ones(len(rows), bool)),
    ('type_n_ge_2_within_type_centered', np.array([counts[aa] >= 2 for aa in amino_acids])),
    ('type_n_ge_4_within_type_centered', np.array([counts[aa] >= 4 for aa in amino_acids])),
    ('all6_glutamines', amino_acids == 'Q'),
):
    x, y, z = sasa[keep].copy(), truth[keep].copy(), prediction[keep].copy()
    groups = amino_acids[keep]
    if 'within_type_centered' in name:
        for aa in np.unique(groups):
            selected = groups == aa
            x[selected] -= x[selected].mean()
            y[selected] -= y[selected].mean()
            z[selected] -= z[selected].mean()
    results.append(dict(subset=name, n=len(x), amino_acids=sorted(set(groups)),
                        rho_sasa_md_hchi=float(spearmanr(x, y).statistic),
                        rho_prediction_md_hchi=float(spearmanr(z, y).statistic)))
report = dict(
    method='Exploratory MDTraj Shrake-Rupley, complete 935-heavy-atom complex; 1.4 A probe; 960 points/atom; starting structure; sidechain excludes N,CA,C,O,OXT; ubiquitin chain A only',
    counts_by_amino_acid=dict(counts), results=results,
    limitation='No independent baseline calibration, no significance inference, one protein and one starting coordinate model; centered values use this protein type means.',
)
(OUTPUT / 'reproducible_summary_with_sensitivity.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))

