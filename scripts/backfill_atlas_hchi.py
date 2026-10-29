"""Add trajectory-derived H-chi to ATLAS in a separate, resumable label release.

The existing RMSF/SASA columns and the old manifests are preserved. Extraction
reuses the ITrace full-replica protocol: CHI_ATOMS, three 120-degree bins, joint
rotamer states, natural-log Shannon entropy, then mean/sample variance across
the three replicas. No model training or split changes are performed here.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Dict, List

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
import numpy as np
from extract_itrace_hchi import dihedral_degrees, entropy
from flexppi.constants import CHI_ATOMS, METRIC_NAMES, THREE_TO_ONE
from flexppi.data.labels import load_label_file
from flexppi.data.records import load_state_records, remap_state_record, write_jsonl
from flexppi.data.structures import parse_pdb

PROTOCOL = 'atlas_hchi_full_replica_120deg_joint_3rep_v1'


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf-8')
    os.replace(temporary, path)


def trajectory_stamp(path):
    p = Path(path); s = p.stat()
    return {'path': str(p.resolve()), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}


def extract_replica(pdb_path, trajectory_path):
    """Use the shared entropy definition with explicit ATLAS atom-name handling."""
    import MDAnalysis as mda
    universe = mda.Universe(str(pdb_path), str(trajectory_path))
    quadruples = []; residue_info = []; alias_count = 0
    for residue in universe.residues:
        name = str(residue.resname).upper()
        definitions = CHI_ATOMS.get(name, [])
        if not definitions or name not in THREE_TO_ONE: continue
        lookup = {}
        for atom in residue.atoms:
            if str(getattr(atom, 'altLoc', '')).strip() not in ('', 'A', '1'): continue
            lookup.setdefault(str(atom.name).upper(), int(atom.index))
        # CHARMM-style ATLAS ILE CD and standard PDB ILE CD1 are the same
        # terminal atom attached to CG1, as checked by consecutive bond lengths.
        if name == 'ILE' and 'CD1' not in lookup and 'CD' in lookup:
            lookup['CD1'] = lookup['CD']; alias_count += 1
        if not all(atom in lookup for definition in definitions for atom in definition): continue
        local = []
        for definition in definitions:
            local.append(len(quadruples)); quadruples.append([lookup[a] for a in definition])
        ca = universe.atoms[lookup['CA']]
        key = {'chain': str(getattr(ca, 'chainID', '')).strip() or '_',
               'residue_number': str(ca.resid), 'insertion_code': str(getattr(ca, 'icode', '')).strip(),
               'amino_acid_3': name, 'n_chi_angles': len(local)}
        residue_info.append((key, local))
    if not quadruples: raise ValueError('No complete side-chain dihedrals')
    quad = np.asarray(quadruples, dtype=np.int64)
    frames = []; min_bond = float('inf'); max_bond = 0.; first_time = None; last_time = None
    dt_ps = float(universe.trajectory.dt)
    try:
        for ts in universe.trajectory:
            positions = ts.positions[quad]
            bonds = np.linalg.norm(np.diff(positions, axis=1), axis=-1)
            min_bond = min(min_bond, float(bonds.min())); max_bond = max(max_bond, float(bonds.max()))
            if min_bond < .5 or max_bond > 3.0:
                raise ValueError('Broken/periodic side-chain geometry; bond range {:.3f}–{:.3f} Å'.format(min_bond, max_bond))
            angles = dihedral_degrees(positions)
            if not np.isfinite(angles).all(): raise ValueError('Nonfinite side-chain angles')
            frames.append(np.clip(np.floor((angles + 180.) / 120.).astype(np.int8), 0, 2))
            if first_time is None: first_time = float(ts.time)
            last_time = float(ts.time)
    finally:
        universe.trajectory.close()
    states = np.stack(frames, axis=0)
    rows = []
    for key, local in residue_info:
        codes = np.zeros(len(states), dtype=np.int32)
        for power, column in enumerate(local): codes += states[:, column].astype(np.int32) * (3 ** power)
        rows.append(dict(key, chi_entropy=entropy(codes), n_frames=len(states)))
    qa = {'frames': len(states), 'dt_ps': dt_ps, 'start_ps': first_time, 'end_ps': last_time,
          'duration_ns': (last_time - first_time) / 1000., 'ile_cd_aliases': alias_count,
          'min_chi_bond_angstrom': min_bond, 'max_chi_bond_angstrom': max_bond}
    return rows, qa


def merge_replica_rows(structure, replica_rows, old_values, old_variances, old_mask):
    """Map by chain/author ID/name; require all replicas for each valid target."""
    keys = list(zip(structure.chain_ids, structure.author_residue_ids, structure.residue_names))
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate structure residue keys')
    index = {key: i for i, key in enumerate(keys)}
    per_replica = np.full((3, len(keys)), np.nan, dtype=np.float64)
    chi_count = np.zeros(len(keys), dtype=np.int16)
    frame_counts = []
    if len(replica_rows) != 3:
        raise ValueError('Exactly three replicas are required')
    for rep, rows in enumerate(replica_rows):
        seen = set(); frames = set()
        # ATLAS PDBs often have a blank chain column. MDAnalysis then reports
        # the segment name (e.g. SYSTEM) through the ITrace helper. A unique
        # blank PDB chain is reconciled by its full author-ID/residue-name key.
        blank_chain = set(structure.chain_ids) == {'_'}
        if blank_chain and len({r['chain'] for r in rows}) != 1:
            raise ValueError('Ambiguous segment-to-blank-chain mapping')
        for row in rows:
            chain = '_' if blank_chain else row['chain']
            key = (chain, str(row['residue_number']) + row['insertion_code'], row['amino_acid_3'])
            if key not in index or key in seen:
                raise ValueError('Unmapped or duplicate trajectory residue: {}'.format(key))
            seen.add(key); pos = index[key]
            count = int(row['n_chi_angles']); value = float(row['chi_entropy'])
            if count != len(CHI_ATOMS.get(key[2], [])) or count == 0:
                raise ValueError('Incomplete/undefined chi angles: {}'.format(key))
            if not np.isfinite(value) or value < -1e-7 or value > count * np.log(3) + 1e-6:
                raise ValueError('Invalid H-chi value: {}'.format(key))
            if chi_count[pos] not in (0, count):
                raise ValueError('Replica chi definitions disagree')
            chi_count[pos] = count; per_replica[rep, pos] = value
            frames.add(int(row['n_frames']))
        if len(frames) != 1 or next(iter(frames)) < 2:
            raise ValueError('Empty replica or inconsistent frame counts')
        frame_counts.append(next(iter(frames)))
    finite = np.isfinite(per_replica)
    if np.any(finite.any(axis=0) != finite.all(axis=0)):
        raise ValueError('Replica residue coverage differs')
    valid = finite.all(axis=0)
    if not valid.any():
        raise ValueError('No valid H-chi labels')
    values, variances, mask = old_values.copy(), old_variances.copy(), old_mask.copy()
    h = METRIC_NAMES.index('chi_entropy')
    values[:, h] = np.nan; variances[:, h] = 0; mask[:, h] = valid
    # Match existing ITrace consensus: float64 aggregation, float32 storage.
    data = per_replica[:, valid].astype(np.float64)
    values[valid, h] = data.mean(axis=0)
    variances[valid, h] = data.var(axis=0, ddof=1)
    for col in range(len(METRIC_NAMES)):
        if col != h:
            for before, after in [(old_values, values), (old_variances, variances), (old_mask, mask)]:
                if not np.array_equal(before[:, col], after[:, col], equal_nan=True):
                    raise AssertionError('Non-H-chi labels changed')
    return values, variances, mask, per_replica, chi_count, np.asarray(frame_counts)


def process_system(job):
    row, out_text, fingerprint = job
    started = time.monotonic(); out = Path(out_text)
    record_id = row['sample_id']
    label = out / 'labels' / (record_id + '.npz')
    report = out / 'systems' / (record_id + '.json')
    source_label = Path(row['label_path']); pdb = Path(row['structure_path'])
    if source_label.resolve() == label.resolve():
        raise ValueError('Output must be separate from original labels')
    trajectories = [Path(p) for p in row['trajectory_paths']]
    if len(trajectories) != 3 or len(set(trajectories)) != 3:
        raise ValueError('Expected three distinct ATLAS trajectories')
    source = {'old_label_sha256': sha256(source_label), 'structure_sha256': sha256(pdb),
              'trajectories': [trajectory_stamp(p) for p in trajectories]}
    if report.exists() and label.exists():
        cached = json.loads(report.read_text(encoding='utf-8'))
        if (cached.get('source') == source and cached.get('fingerprint') == fingerprint
                and cached.get('output_sha256') == sha256(label)):
            cached['resumed'] = True
            return cached
    structure = parse_pdb(str(pdb))
    old = load_label_file(str(source_label))
    if old[0].shape[0] != structure.num_residues:
        raise ValueError('Structure and original label lengths differ')
    rows = []; replica_qa = []
    for replica, trajectory in enumerate(trajectories, 1):
        result, qa = extract_replica(pdb, trajectory)
        rows.append(result); replica_qa.append(qa)
    values, variances, mask, replicas, nchi, frames = merge_replica_rows(structure, rows, *old)
    if source['trajectories'] != [trajectory_stamp(p) for p in trajectories]:
        raise ValueError('Trajectory changed during extraction')
    if source['old_label_sha256'] != sha256(source_label) or source['structure_sha256'] != sha256(pdb):
        raise ValueError('Source label or structure changed during extraction')
    temporary = label.with_suffix('.npz.tmp')
    with temporary.open('wb') as handle:
        np.savez_compressed(handle, values=values, variances=variances, mask=mask,
                            metric_names=np.asarray(METRIC_NAMES), hchi_replica_values=replicas,
                            hchi_replica_frames=frames, hchi_chi_count=nchi,
                            chain_ids=np.asarray(structure.chain_ids),
                            author_residue_ids=np.asarray(structure.author_residue_ids),
                            residue_names=np.asarray(structure.residue_names),
                            hchi_protocol=np.asarray(PROTOCOL))
    os.replace(temporary, label)
    h = METRIC_NAMES.index('chi_entropy')
    missing = {}
    for name, is_valid in zip(structure.residue_names, mask[:, h]):
        if not is_valid: missing[name] = missing.get(name, 0) + 1
    updated = dict(row); updated['label_path'] = str(label.resolve())
    updated['target_chains'] = structure.unique_chains()
    updated['metadata'] = dict(row.get('metadata', {}), hchi_labels=int(mask[:, h].sum()),
                               label_protocol=PROTOCOL, hchi_replicates=3,
                               hchi_frames_per_replica=frames.tolist(),
                               original_label_path=str(source_label.resolve()),
                               original_label_sha256=source['old_label_sha256'])
    result = {'sample_id': record_id, 'status': 'complete', 'fingerprint': fingerprint,
              'source': source, 'output_sha256': sha256(label), 'record': updated,
              'residues': structure.num_residues, 'rmsf_labels': int(mask[:, 0].sum()),
              'hchi_labels': int(mask[:, h].sum()), 'missing_hchi_by_resname': missing,
              'frames_per_replica': frames.tolist(), 'non_hchi_columns_unchanged': True,
              'replica_qa': replica_qa,
              'residue_mapping': 'exact chain/author-ID/name; unique blank PDB chain reconciles segment names',
              'hchi_min': float(values[mask[:, h], h].min()),
              'hchi_max': float(values[mask[:, h], h].max()),
              'seconds': round(time.monotonic() - started, 3)}
    atomic_json(report, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=ROOT / 'artifacts/data/atlas_states.jsonl')
    parser.add_argument('--atlas-root', type=Path, default=Path(r'E:\data\protein_md_datasets\ATLAS_extracted'))
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'artifacts/data/atlas_rmsf_hchi_v1')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    if args.workers < 1: parser.error('--workers must be positive')
    out = args.output_dir.resolve()
    for directory in [out, out / 'labels', out / 'systems']: directory.mkdir(parents=True, exist_ok=True)
    mappings = {'/ranze/data/ATLAS_extracted': str(args.atlas_root.resolve()), '/ranze/repos/FlexPPI': str(ROOT)}
    records = [remap_state_record(r, mappings).to_dict() for r in load_state_records(str(args.manifest))]
    if args.limit: records = records[:args.limit]
    if any(r['dataset'] != 'atlas' for r in records): raise ValueError('Manifest must contain only ATLAS')
    if len({r['sample_id'] for r in records}) != len(records): raise ValueError('Duplicate sample IDs')
    fingerprint = {'protocol': PROTOCOL, 'extractor_sha256': sha256(ROOT / 'scripts/extract_itrace_hchi.py'),
                   'constants_sha256': sha256(ROOT / 'src/flexppi/constants.py'),
                   'backfill_sha256': sha256(Path(__file__))}
    run = {'protocol': PROTOCOL, 'started_utc': datetime.now(timezone.utc).isoformat(),
           'source_manifest': str(args.manifest.resolve()), 'source_manifest_sha256': sha256(args.manifest),
           'fingerprint': fingerprint, 'requested_systems': len(records), 'workers': args.workers,
           'scope': 'pilot' if args.limit else 'full', 'original_labels_preserved': True}
    atomic_json(out / 'run_config.json', run)
    results = {}; failures = {}; started = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_system, (r, str(out), fingerprint)): r['sample_id'] for r in records}
        for future in as_completed(futures):
            sid = futures[future]
            try:
                result = future.result(); results[sid] = result
            except Exception as error:
                failures[sid] = '{}: {}'.format(type(error).__name__, error)
            done = len(results) + len(failures)
            status = {'done': done, 'total': len(records), 'complete': len(results),
                      'failures': failures, 'elapsed_seconds': round(time.monotonic() - started, 1)}
            atomic_json(out / 'progress.json', status)
            if done <= 3 or done % 25 == 0 or done == len(records) or sid in failures:
                print(json.dumps(dict(status, latest=sid)), flush=True)
    rows = [results[r['sample_id']]['record'] for r in records if r['sample_id'] in results]
    name = 'atlas_states_hchi.jsonl' if not failures else 'atlas_states_hchi.partial.jsonl'
    write_jsonl(str(out / name), rows)
    missing = {}
    for v in results.values():
        for aa, n in v['missing_hchi_by_resname'].items(): missing[aa] = missing.get(aa, 0) + n
    summary = dict(run, **status, manifest=str(out / name),
                   rmsf_labels=sum(v['rmsf_labels'] for v in results.values()),
                   hchi_labels=sum(v['hchi_labels'] for v in results.values()),
                   missing_hchi_by_resname=missing,
                   frame_count_histogram={str(n): sum(v['frames_per_replica'].count(n) for v in results.values())
                                          for n in sorted({n for v in results.values() for n in v['frames_per_replica']})})
    atomic_json(out / 'summary.json', summary)
    print(json.dumps(summary, ensure_ascii=False), flush=True)
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
