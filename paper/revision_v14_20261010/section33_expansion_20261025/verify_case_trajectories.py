"""Recompute all-residue case dynamics from raw XTC with canonical atom names.

New outputs only. Saving interval and physical duration come from XTC, not
metadata. Blocks are time segments of one trajectory, never independent trials.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.formats.libmdaxdr import XTCFile
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
CORE = Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
sys.path.insert(0, str(CORE/'training_payload/src'))
from flexppi.constants import CHI_ATOMS, THREE_TO_ONE
from flexppi.data.labels import _dihedral_degrees, _entropy, _chainwise_rmsf


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def circular_metrics(angles):
    # Return means over chi angles and four origins of 5-degree histograms.
    values = []
    for offset in (0., 1.25, 2.5, 3.75):
        for col in range(angles.shape[1]):
            codes = np.floor(((angles[:, col]+180-offset) % 360)/5).astype(int)
            counts = np.bincount(codes, minlength=72)
            probabilities = counts[counts > 0]/len(codes)
            ent = -np.sum(probabilities*np.log(probabilities))
            width = float((counts >= counts.max()/np.e).sum()*5)
            values.append((ent, width))
    return np.array(values).mean(axis=0)


def joint_entropy(angles, offset=0.):
    if offset == 0.:
        # Exact frozen-label rule: +180 belongs to the last bin, not wrapped
        # into the first. Keep this endpoint convention during reproduction.
        bins = np.clip(np.floor((angles+180.)/120.).astype(int), 0, 2)
    else:
        bins = np.floor(((angles+180-offset) % 360)/120).astype(int)
    codes = (bins*np.power(3, np.arange(angles.shape[1]))).sum(axis=1)
    return _entropy(codes)


def main(sid, do_rmsf=False):
    out = ROOT/'verified_case_trajectories'/sid
    out.mkdir(parents=True, exist_ok=True)
    pdb = CORE/f'training_payload/dataset/structures/{sid}.pdb'
    xtc = Path('D:/data/dynarepo_non_nucleic_full')/sid/'trajectory.xtc'
    cache = CORE/f'neural_exports/evaluation/full/{sid}.npz'
    universe = mda.Universe(str(pdb))
    residues = [r for r in universe.residues if r.resname in THREE_TO_ONE and 'CA' in r.atoms.names]
    with np.load(cache) as z:
        truth, valid = z['target'][:, :2], z['valid'][:, :2]
        aa, chains, authors = z['amino_acid'], z['chain'], z['author_residue_id']
        predictions = z['mean']
    assert len(residues) == len(truth)
    assert [(str(r.atoms.chainIDs[0]).strip() or '_', str(r.resid)+str(r.icode).strip(), THREE_TO_ONE[r.resname])
            for r in residues] == list(zip(chains, authors, aa))
    quads, cols, ca = [], {}, []
    for i, residue in enumerate(residues):
        atoms = {atom.name: int(atom.index) for atom in residue.atoms}
        ca.append(atoms['CA'])
        definitions = CHI_ATOMS.get(residue.resname, [])
        if definitions and all(all(name in atoms for name in quad) for quad in definitions):
            cols[i] = list(range(len(quads), len(quads)+len(definitions)))
            quads.extend([[atoms[name] for name in quad] for quad in definitions])
    quad = np.asarray(quads)
    angles, coords, times = [], [], []
    max_bond, min_bond = 0., np.inf
    with XTCFile(str(xtc), 'r') as stream:
        for frame in stream:
            xyz = np.asarray(frame.x, dtype=np.float32)*10
            assert len(xyz) == len(universe.atoms), (len(xyz), len(universe.atoms))
            points = xyz[quad]
            bonds = np.linalg.norm(np.diff(points, axis=1), axis=2)
            max_bond = max(max_bond, float(bonds.max()))
            min_bond = min(min_bond, float(bonds.min()))
            if max_bond > 2.25 or min_bond < .95:
                raise ValueError(f'Broken chi covalent geometry; no silent PBC repair {sid}: {min_bond}, {max_bond}')
            angles.append(_dihedral_degrees(points))
            if do_rmsf:
                coords.append(xyz[ca].copy())
            times.append(float(frame.time))
    times, angles = np.asarray(times), np.asarray(angles, dtype=np.float32)
    dt = np.diff(times)
    assert np.allclose(dt, dt[0], rtol=1e-5, atol=1e-3)
    nblocks = len(times)//1000
    if nblocks == 0:
        raise ValueError('No original 1000-frame window')
    used = nblocks*1000
    blocks = [slice(i*1000, (i+1)*1000) for i in range(nblocks)]
    np.savez_compressed(out/'all_chi_angles.npz', angles_degrees=angles[:used], time_ps=times[:used],
                        quadruple_atoms=quad, residue_index=np.array(list(cols)),
                        chi_start=np.array([c[0] for c in cols.values()]),
                        chi_count=np.array([len(c) for c in cols.values()]))
    values = np.full((nblocks, len(truth)), np.nan)
    side_rows, chi_rows, window_rows = [], [], []
    for i, localcols in cols.items():
        theta = angles[:used, localcols]
        block_hchi = np.array([joint_entropy(theta[s]) for s in blocks])
        values[:, i] = block_hchi
        fine = circular_metrics(theta)
        block_fine = np.array([circular_metrics(theta[s]) for s in blocks])
        shifted = np.array([joint_entropy(theta[s], 60.) for s in blocks])
        side_rows.append(dict(sample_id=sid, residue_index=i, chain=chains[i], author_residue_id=authors[i],
                              amino_acid=aa[i], n_chi=len(localcols), frozen_hchi=float(truth[i, 1]),
                              model_hchi=float(predictions[i, 1]), recomputed_hchi=float(block_hchi.mean()),
                              hchi_block_min=float(block_hchi.min()), hchi_block_max=float(block_hchi.max()),
                              shifted60_joint_entropy_mean=float(shifted.mean()),
                              fine_entropy_mean_per_chi=float(fine[0]), width_mean_per_chi=float(fine[1]),
                              fine_entropy_window_mean=float(block_fine[:, 0].mean()),
                              width_window_mean=float(block_fine[:, 1].mean())))
        for b, (h, (e, w), sh) in enumerate(zip(block_hchi, block_fine, shifted)):
            window_rows.append(dict(sample_id=sid, residue_index=i, chain=chains[i], author_residue_id=authors[i],
                                    amino_acid=aa[i], block=b+1, start_ps=float(times[b*1000]),
                                    end_ps=float(times[(b+1)*1000-1]), hchi=float(h),
                                    fine_entropy_mean_per_chi=float(e), width_mean_per_chi=float(w),
                                    shifted60_joint_entropy=float(sh)))
        for chi, col in enumerate(localcols):
            e, w = circular_metrics(angles[:used, [col]])
            chi_rows.append(dict(sample_id=sid, residue_index=i, chain=chains[i], author_residue_id=authors[i],
                                 amino_acid=aa[i], chi=chi+1, fine_entropy=float(e), angular_width=float(w)))
    s = pd.DataFrame(side_rows)
    s.to_csv(out/'all_sidechains_occupancy.csv', index=False)
    pd.DataFrame(window_rows).to_csv(out/'all_sidechains_per_window.csv', index=False)
    pd.DataFrame(chi_rows).to_csv(out/'all_chi_occupancy.csv', index=False)
    good = valid[:, 1] & np.isfinite(values.mean(axis=0))
    discrepancy = float(np.max(np.abs(values.mean(axis=0)[good]-truth[good, 1])))
    assert discrepancy < 1e-5, (sid, discrepancy)
    rmsf_error = None
    if do_rmsf:
        coords = np.asarray(coords, dtype=np.float32)
        ci = np.array([list(dict.fromkeys(chains)).index(ch) for ch in chains])
        rb = np.array([_chainwise_rmsf(coords[segment], ci) for segment in blocks])
        rmsf_error = float(np.max(np.abs(rb.mean(0)[valid[:, 0]]-truth[valid[:, 0], 0])))
        assert rmsf_error < 1e-5, (sid, rmsf_error)
        rr = []
        for i in range(len(truth)):
            for b in range(nblocks):
                rr.append(dict(sample_id=sid, residue_index=i, chain=chains[i], author_residue_id=authors[i],
                               amino_acid=aa[i], block=b+1, rmsf=float(rb[b, i]),
                               frozen_rmsf=float(truth[i, 0]), model_rmsf=float(predictions[i, 0])))
        pd.DataFrame(rr).to_csv(out/'all_residues_rmsf_per_window.csv', index=False)
        np.savez_compressed(out/'ca_coordinates.npz', coords_A=coords[:used], time_ps=times[:used],
                            chain=chains, author_residue_id=authors, amino_acid=aa)
    audit = dict(status='passed', sample_id=sid, n_frames=len(times), n_frames_used=used,
                 first_time_ps=float(times[0]), last_time_ps=float(times[-1]), dt_ps=float(dt[0]),
                 duration_ns=float((times[-1]-times[0])/1000), n_original_windows=nblocks,
                 original_window_duration_ns=float(1000*dt[0]/1000),
                 n_sidechains=int(good.sum()), hchi_max_recomputation_error=discrepancy,
                 rmsf_max_recomputation_error=rmsf_error, min_chi_bond_A=min_bond, max_chi_bond_A=max_bond,
                 source_files=[dict(path=str(p), sha256=sha(p)) for p in [pdb, xtc, cache]],
                 interpretation='All sidechains; original bins revalidated; 5deg angular histograms use four origins. '
                                'Offset60deg joint bins center on conventional rotamers and are only a sensitivity. '
                                'Time windows of one trajectory, not independent repeats.')
    (out/'audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('samples', nargs='+')
    ap.add_argument('--rmsf', action='store_true')
    args = ap.parse_args()
    for sample in args.samples:
        main(sample, args.rmsf)
