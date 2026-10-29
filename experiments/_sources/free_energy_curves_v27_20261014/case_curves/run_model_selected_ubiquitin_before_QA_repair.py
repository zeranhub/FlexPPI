"""Prediction-first ubiquitin GLN cases; existing raw trajectory is read-only.

The primary curve grid (501 observations, 0-100 ns at 200 ps) matches the
population analysis. Original 5000 observations remain a separate label audit.
Only this case_curves directory is written. No model inference or new MD.
"""
from __future__ import annotations
import sys
sys.dont_write_bytecode = True
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
from pathlib import Path
import ast
import argparse
import hashlib
import json
import shutil
import tempfile
from datetime import datetime, timezone
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
CORRECTED = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
CACHE = CORRECTED / 'neural_exports/evaluation/full/A003R.npz'
LABEL = CORRECTED / 'training_payload/dataset/labels/A003R.npz'
LABEL_CODE = CORRECTED / 'training_payload/src/flexppi/data/labels.py'
PROVENANCE = Path(r'E:\Project\FlexPPI\experiments\_case_validation_sources\ubiquitin\provenance.json')
Q41_OLD = BASE / 'energy_landscape_v21_20261012/case_analysis/Q41_angles_used.csv'
RUNTIME = BASE / 'label_sensitivity/runtime'
ORIGINS = (0., 1.25, 2.5, 3.75)
CHI_ATOMS = (('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'OE1'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf8')


def frozen_plan():
    HERE.mkdir(parents=True, exist_ok=True)
    prov = json.loads(PROVENANCE.read_text(encoding='utf8'))
    sources = {str(PROVENANCE): sha(PROVENANCE), str(CACHE): 'b8a764739993e3bce18dedb3cff493c50a1a52613a08baac16f1803c4edabb5a',
               str(LABEL): 'a528abc00e9df39370da04bf76957acb471b9b900ca66e53ccdde955e7d409c2', str(LABEL_CODE): sha(LABEL_CODE),
               str(Q41_OLD): 'd17cf3115073b20a9126aca7f8d186650d9e63ad9445a45f830890b341842b15',
               prov['source_structure']: prov['source_structure_sha256'], prov['source_trajectory']: prov['source_trajectory_sha256']}
    for path, digest in sources.items(): assert sha(path) == digest, path
    cache = np.load(CACHE, allow_pickle=False)
    indices = np.flatnonzero((cache['chain'] == 'A') & (cache['amino_acid'] == 'Q') & cache['valid'][:, 1] & np.isfinite(cache['mean'][:, 1]))
    inventory = pd.DataFrame([dict(residue='Q' + str(cache['author_residue_id'][i]), chain='A', author_residue_id=str(cache['author_residue_id'][i]),
        residue_index=int(i), predicted_joint_Hchi_nats=float(cache['mean'][i, 1])) for i in indices])
    assert inventory.residue.tolist() == ['Q2', 'Q31', 'Q40', 'Q41', 'Q49', 'Q62']
    ordered = inventory.sort_values(['predicted_joint_Hchi_nats', 'residue_index'], kind='mergesort')
    low, high = str(ordered.residue.iloc[0]), str(ordered.residue.iloc[-1])
    inventory['prediction_rank_ascending'] = inventory.residue.map({r: j + 1 for j, r in enumerate(ordered.residue)})
    inventory['selection'] = inventory.residue.map({low: 'lowest_predicted_Hchi', high: 'highest_predicted_Hchi'}).fillna('not_selected')
    plan = dict(protocol='OOF_prediction_extreme_Gln_ubiquitin_MD_curves_v1', created_utc=datetime.now(timezone.utc).isoformat(),
        analysis_code_sha256=sha(__file__), source_hash_files=sources,
        model='Original frozen corrected Full fold0, epoch40 seed7; output1=joint Hchi; no training/inference',
        selection='Among every valid chainA GLN in A003R, choose minimum and maximum predicted joint Hchi; stable tie by canonical residue_index; selection reads only model scores and valid mask, never target/MD angle results',
        prediction_inventory=inventory.to_dict('records'), low_predicted_residue=low, high_predicted_residue=high,
        primary_curve_grid=dict(frame_indices=list(range(0, 1001, 2)), frames=501, first_ns=0., last_ns=100., step_ps=200.,
            rationale='Fixed common time/sampling grid of population analysis; adopted before new case-angle/width computation; no trajectory-window search'),
        original_label_grid=dict(frame_indices=list(range(5000)), frames=5000, first_ns=0., last_ns=499.9, step_ps=100.,
            windows='Five consecutive nonoverlapping 1000-frame windows; window arithmetic mean/variance reproduces original joint entropy target/label_variance'),
        chi_definitions=[list(x) for x in CHI_ATOMS], angle_code='Original unmodified labels.py _dihedral_degrees using float32 XTC coordinates; independent float64 MDAnalysis.calc_dihedrals check',
        primary_histogram='72 equal 5-degree periodic bins with origin0: [-180,-175),...,[175,180); wrap+180 to-180; all bins retained including empty',
        Fstar='-ln(p_bin/p_max), p=count/numberframes; empty bins missing; no smoothing/interpolation/pseudocount; marginal MD observations, not model-predicted free energy',
        width='W(1): union measure in degrees of occupied bins with Fstar<=1, including disconnected bins; per-angle mean over four origins [0,1.25,2.5,3.75] then mean across three GLN chi angles; identical to v26 population definition',
        source_observations='Extract every six GLN to verify identity/targets; only prediction-extreme pair are primary plotted cases; no case replacement or angle exclusion',
        all_original_results_preserved='Primary501 and original5000 have separate angle/curve/summary files; all six inventory outcomes retained; no result-based omissions')
    destination = HERE / 'analysis_plan.json'
    if destination.exists():
        previous = json.loads(destination.read_text(encoding='utf8'))
        assert previous['analysis_code_sha256'] == sha(__file__)
        assert previous['source_hash_files'] == sources
        assert previous['prediction_inventory'] == plan['prediction_inventory']
        plan = previous
    else:
        write_json(destination, plan)
        inventory.to_csv(HERE / 'selection_inventory_prediction_only.csv', index=False)
    return plan


def original_angle_function():
    tree = ast.parse(LABEL_CODE.read_text(encoding='utf8'))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == '_dihedral_degrees')
    namespace = {'np': np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(LABEL_CODE), 'exec'), namespace)
    return namespace['_dihedral_degrees']


def histogram(angles, origin=0.):
    bins = np.floor(((angles.astype(float) + 180 - origin) % 360) / 5).astype(np.int64)
    assert bins.min() >= 0 and bins.max() < 72
    counts = np.bincount(bins, minlength=72)
    assert int(counts.sum()) == len(angles)
    fstar = np.full(72, np.nan, float)
    nonzero = counts > 0
    fstar[nonzero] = -np.log(counts[nonzero] / counts.max())
    # Alternative logarithmic definition, independently evaluated.
    np.testing.assert_allclose(fstar[nonzero], np.log(counts.max()) - np.log(counts[nonzero]), atol=1e-14, rtol=0)
    return counts, fstar


def widths(angles, site, grid):
    rows = []
    for qi in range(3):
        for origin in ORIGINS:
            count, fstar = histogram(angles[:, qi], origin)
            mask = (count > 0) & (fstar <= 1.)
            equivalent = (count > 0) & (count >= count.max() * np.exp(-1.))
            np.testing.assert_array_equal(mask, equivalent)
            rows.append(dict(grid=grid, residue=site, chi=f'chi{qi + 1}', origin_degrees=origin,
                n_frames=len(angles), low_Fstar_threshold=1., width_degrees=float(5 * mask.sum()),
                maximum_bin_count=int(count.max()), occupied_bins=int((count > 0).sum()), qualifying_bins=int(mask.sum())))
    return rows


def curves(angles, site, grid, window):
    rows = []
    for qi in range(3):
        count, fstar = histogram(angles[:, qi])
        for b in range(72):
            rows.append(dict(grid=grid, residue=site, window=window, chi=f'chi{qi + 1}', angle_lower_degrees=-180 + 5 * b,
                angle_upper_degrees=-175 + 5 * b, angle_center_degrees=-177.5 + 5 * b,
                n_frames=len(angles), count=int(count[b]), probability=float(count[b] / len(angles)),
                density_per_degree=float(count[b] / len(angles) / 5), F_star=float(fstar[b]),
                low_Fstar_le1=bool(count[b] > 0 and fstar[b] <= 1.)))
    return rows


def joint_entropy(angles):
    # Exact original float32 endpoint rule. +180 belongs to final original sector.
    bins = np.clip(np.floor((angles + np.float32(180)) / np.float32(120)).astype(np.int16), 0, 2)
    codes = (bins * (3 ** np.arange(3))).sum(axis=1)
    count = np.bincount(codes, minlength=27)
    p = count[count > 0].astype(float) / len(angles)
    value = float(-np.dot(p, np.log(p)))
    _, independent_counts = np.unique(codes, return_counts=True)
    independent_p = independent_counts.astype(float) / independent_counts.sum()
    assert abs(value - float(-np.sum(independent_p * np.log(independent_p)))) <= 1e-14
    return value


def analyze(plan):
    sys.path.insert(0, str(RUNTIME))
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import calc_dihedrals
    prov = json.loads(PROVENANCE.read_text(encoding='utf8'))
    angle_function = original_angle_function()
    inventory = pd.DataFrame(plan['prediction_inventory'])
    sites = inventory.residue.tolist()
    selected = (plan['low_predicted_residue'], plan['high_predicted_residue'])
    source_stat = {p: (Path(p).stat().st_size, Path(p).stat().st_mtime_ns) for p in plan['source_hash_files']}
    stage = HERE / '_staging'
    stage.mkdir(exist_ok=True)
    angles = np.empty((5000, len(sites), 3), np.float32)
    times = np.empty(5000, float)
    maximum_independent_error = 0.
    maximum_bond_length = 0.
    atom_mappings = []
    with tempfile.TemporaryDirectory(prefix='existing_md_', dir=stage) as folder:
        folder = Path(folder)
        assert folder.resolve().is_relative_to(stage.resolve())
        top, traj = folder / 'structure.pdb', folder / 'trajectory.xtc'
        shutil.copy2(prov['source_structure'], top)
        shutil.copy2(prov['source_trajectory'], traj)
        assert sha(top) == prov['source_structure_sha256'] and sha(traj) == prov['source_trajectory_sha256']
        universe = mda.Universe(str(top), str(traj))
        try:
            assert len(universe.trajectory) == 5001
            quads = []
            for row in inventory.to_dict('records'):
                residue = universe.select_atoms(f"chainID A and resid {row['author_residue_id']}").residues
                assert len(residue) == 1 and residue[0].resname == 'GLN'
                assert residue[0].ix == row['residue_index']
                lookup = {a.name: int(a.index) for a in residue[0].atoms}
                assert len(lookup) == len(residue[0].atoms)
                mapped = [[lookup[a] for a in quad] for quad in CHI_ATOMS]
                atom_mappings.append(dict(residue=row['residue'], cache_index=row['residue_index'], mda_resindex=int(residue[0].ix), atom_indices=mapped, atom_names=CHI_ATOMS))
                quads.extend(mapped)
            quads = np.array(quads, np.int64)
            for frame, ts in enumerate(universe.trajectory[:5000]):
                assert abs(ts.time - frame * 100) <= .05
                times[frame] = ts.time
                points = ts.positions[quads]
                values = angle_function(points).reshape(len(sites), 3)
                assert np.isfinite(values).all()
                angles[frame] = values
                points64 = points.astype(np.float64)
                reference = np.degrees(calc_dihedrals(points64[:, 0], points64[:, 1], points64[:, 2], points64[:, 3])).reshape(len(sites), 3)
                errors = abs((values.astype(float) - reference + 180) % 360 - 180)
                maximum_independent_error = max(maximum_independent_error, float(errors.max()))
                maximum_bond_length = max(maximum_bond_length, float(np.linalg.norm(np.diff(points, axis=1), axis=-1).max()))
        finally:
            universe.trajectory.close()
    assert maximum_independent_error < 1e-3
    assert maximum_bond_length < 5.
    old = pd.read_csv(Q41_OLD)
    assert len(old) == 5000
    old_values = old[[f'chi{q}_degrees' for q in (1, 2, 3)]].to_numpy(float)
    q41_values = angles[:, sites.index('Q41'), :].astype(float)
    q41_errors = abs((q41_values - old_values + 180) % 360 - 180)
    q41_max_error = float(q41_errors.max())
    np.testing.assert_allclose(q41_values, old_values, atol=1e-5, rtol=0)
    cache = np.load(CACHE, allow_pickle=False)
    label = np.load(LABEL, allow_pickle=False)
    np.testing.assert_array_equal(cache['valid'], label['mask'])
    np.testing.assert_allclose(cache['target'], label['values'], rtol=0, atol=0, equal_nan=True)
    angle_rows_5000, angle_rows_501 = [], []
    curve_rows_5000, curve_rows_501, width_rows, summaries, joint_rows = [], [], [], [], []
    for si, site in enumerate(sites):
        metadata = inventory.iloc[si].to_dict()
        values = angles[:, si, :]
        short = values[::2][:501]
        assert len(short) == 501
        for frame in range(5000):
            row = dict(residue=site, selection=metadata['selection'], frame_zero_based=frame, time_ns=float(times[frame] / 1000),
                **{f'chi{q + 1}_degrees': float(values[frame, q]) for q in range(3)})
            angle_rows_5000.append(row)
            if frame <= 1000 and frame % 2 == 0: angle_rows_501.append(row)
        means = []
        for wi in range(5):
            value = joint_entropy(values[wi * 1000:(wi + 1) * 1000])
            means.append(value)
            joint_rows.append(dict(residue=site, window=wi + 1, start_frame=wi * 1000, end_frame_exclusive=(wi + 1) * 1000, Hchi_nats=value))
        index = int(metadata['residue_index'])
        md_h = float(np.mean(means))
        md_var = float(np.var(means, ddof=0))
        np.testing.assert_allclose(md_h, cache['target'][index, 1], atol=2e-6, rtol=0)
        np.testing.assert_allclose(md_var, cache['label_variance'][index, 1], atol=2e-6, rtol=0)
        own_width_501 = widths(short, site, 'primary_501frames')
        own_width_5000 = widths(values, site, 'label_source_5000frames')
        width_rows += own_width_501 + own_width_5000
        summaries.append(metadata | dict(MD_joint_Hchi_mean5windows_nats=md_h, frozen_MD_joint_Hchi_target_nats=float(cache['target'][index, 1]),
            prediction_absolute_error_nats=abs(float(metadata['predicted_joint_Hchi_nats']) - md_h),
            original_joint_label_variance=md_var,
            MD_width_501frames_deg=float(np.mean([r['width_degrees'] for r in own_width_501])),
            MD_width_5000frames_deg=float(np.mean([r['width_degrees'] for r in own_width_5000])),
            **{f'MD_width_501_chi{q + 1}_deg': float(np.mean([r['width_degrees'] for r in own_width_501 if r['chi'] == f'chi{q + 1}'])) for q in range(3)}))
        if site in selected:
            curve_rows_501 += curves(short, site, 'primary_501frames', 'pooled')
            curve_rows_5000 += curves(values, site, 'label_source_5000frames', 0)
            for wi in range(5): curve_rows_5000 += curves(values[wi * 1000:(wi + 1) * 1000], site, 'label_source_5000frames', wi + 1)
    outputs = {}
    for filename, rows in [('angles_501frames_all6GLN.csv', angle_rows_501), ('angles_5000frames_all6GLN.csv', angle_rows_5000),
        ('histograms_501frames.csv', curve_rows_501), ('histograms_5000frames.csv', curve_rows_5000),
        ('width_per_angle_origin.csv', width_rows), ('selection_inventory_with_MD_outcomes.csv', summaries),
        ('joint_entropy_original_five_windows.csv', joint_rows)]:
        destination = HERE / filename
        pd.DataFrame(rows).to_csv(destination, index=False)
        outputs[filename] = dict(sha256=sha(destination), rows=len(rows))
    assert outputs['histograms_501frames.csv']['rows'] == 432
    assert outputs['histograms_5000frames.csv']['rows'] == 2592
    assert outputs['angles_501frames_all6GLN.csv']['rows'] == 3006
    assert outputs['angles_5000frames_all6GLN.csv']['rows'] == 30000
    for path, digest in plan['source_hash_files'].items():
        assert sha(path) == digest
        assert (Path(path).stat().st_size, Path(path).stat().st_mtime_ns) == source_stat[path]
    audit = dict(status='complete', plan_sha256=sha(HERE / 'analysis_plan.json'), source_files_sha256_unchanged=True,
        source_files_size_and_mtime_unchanged=True, selection_uses_only_prediction=True,
        selected_low=selected[0], selected_high=selected[1], source_frames=5001, analyzed_original_frames=5000,
        primary_frames=501, primary_time_range_ns=[0., 100.], primary_step_ps=200.,
        raw_MDA_access='Only verified staging copy, never source trajectory; staging auto-deleted after close',
        MDAnalysis_version=mda.__version__, numpy_version=np.__version__,
        independent_float64_dihedral_max_error_degrees=maximum_independent_error, maximum_adjacent_atom_distance_angstrom=maximum_bond_length,
        Q41_original_archive_max_angle_error_degrees=q41_max_error,
        Q41_matches_original_archive_atol_degrees=1e-5, original_joint_Hchi_and_variance_verified_all6GLN=True,
        atom_mappings=atom_mappings, outputs=outputs, all6_summaries=summaries)
    write_json(HERE / 'extraction_audit.json', audit)
    print(json.dumps(clean({'selected_low': selected[0], 'selected_high': selected[1], 'all6_summaries': summaries,
        'Q41_max_error': q41_max_error, 'independent_dihedral_max_error': maximum_independent_error}), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze-only', action='store_true')
    args = parser.parse_args()
    plan = frozen_plan()
    if args.freeze_only:
        print(json.dumps({'selected_low': plan['low_predicted_residue'], 'selected_high': plan['high_predicted_residue'],
            'inventory': plan['prediction_inventory'], 'plan_sha256': sha(HERE / 'analysis_plan.json')}, ensure_ascii=False, indent=2))
    else:
        analyze(plan)
