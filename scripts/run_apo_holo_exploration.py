"""Audited five-checkpoint absolute-state inference for the exploratory cohort.

No weights are fitted. Intermediate predictions retain both raw residue numbers
and common-residue correspondence. See the adjacent recorded protocol.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts/apo_holo_exploration_20260913"
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
from export_atlas_hchi_v1_oof_predictions import load_fold, implementation_hash
from flexppi.inference import _structure_example
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb
from flexppi.data.pairing import global_sequence_alignment
from flexppi.data.labels import load_label_file, _sphere_points, _shrake_rupley_residue_sasa
from flexppi.constants import THREE_TO_ONE

KEYS = ('aatype', 'residue_index', 'chain_index', 'backbone', 'backbone_mask',
        'geometry_mask', 'target_mask', 'label_values', 'label_variances', 'label_mask')


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def subset(example, indices):
    result = dict(example)
    for key in KEYS:
        result[key] = example[key][indices].clone()
    chain = result['chain_index'].numpy()
    mapping = {old: new for new, old in enumerate(dict.fromkeys(chain.tolist()))}
    result['chain_index'] = torch.tensor([mapping[x] for x in chain], dtype=torch.long)
    return result


def state(spec):
    structure = parse_pdb(spec['path'])
    chains = spec.get('chains') or structure.unique_chains()
    if not set(chains).issubset(structure.unique_chains()):
        raise ValueError(f"Missing chains {chains}; actual {structure.unique_chains()} in {spec['path']}")
    ix = np.asarray([i for i, c in enumerate(structure.chain_ids) if c in chains], dtype=int)
    if len(ix) > 1024 or len(ix) < 2:
        raise ValueError(f"Input size {len(ix)} outside 2..1024; no truncation")
    example = subset(_structure_example(spec['path'], chains, 300., 100., 0.), ix)
    return dict(spec=spec, structure=structure, indices=ix, example=example,
                chain=np.array(structure.chain_ids)[ix], aa=np.array(list(structure.sequence))[ix],
                author=np.array(structure.author_residue_ids)[ix])


def align_states(apo, holo, target):
    achains = list(dict.fromkeys(apo['chain']))
    if len(achains) < len(target):
        raise ValueError("Too few apo chains")
    candidates = []
    for assignment in itertools.permutations(achains, len(target)):
        pa, ph, quality = [], [], []
        for ac, hc in zip(assignment, target):
            ai = np.flatnonzero(apo['chain'] == ac)
            hi = np.flatnonzero(holo['chain'] == hc)
            if not len(ai) or not len(hi):
                break
            aln = global_sequence_alignment(''.join(apo['aa'][ai]), ''.join(holo['aa'][hi]))
            cov = min(len(aln.indices_a)/len(ai), len(aln.indices_b)/len(hi))
            if aln.identity < .95 or cov < .85:
                break
            pa.extend(ai[aln.indices_a]); ph.extend(hi[aln.indices_b])
            quality.append(dict(apo_chain=str(ac), holo_chain=str(hc), identity=aln.identity, coverage=cov))
        if len(quality) == len(target):
            candidates.append((len(pa), pa, ph, quality, list(assignment)))
    if not candidates:
        raise ValueError("No chain mapping passing identity >=0.95 and bidirectional exact coverage >=0.85")
    return max(candidates, key=lambda x: x[0])


def heavy_atoms(st):
    structure = st['structure']
    lookup = {(structure.chain_ids[i], structure.author_residue_ids[i]): k for k, i in enumerate(st['indices'])}
    positions, residues, radii, seen = [], [], [], set()
    for line in Path(st['spec']['path']).read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if not line.startswith(('ATOM  ', 'HETATM')) or len(line) < 54:
            continue
        if line[16] not in ' A1' or line[17:20].strip() not in THREE_TO_ONE:
            continue
        atom = line[12:16].strip()
        element = line[76:78].strip() or atom.lstrip('0123456789')[0]
        if element in ('H', 'D'):
            continue
        key = (line[21].strip() or '_', line[22:26].strip() + line[26].strip())
        if key not in lookup or (key, atom) in seen:
            continue
        seen.add((key, atom))
        positions.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        residues.append(lookup[key])
        radii.append(dict(C=1.7, N=1.55, O=1.52, S=1.8, SE=1.9).get(element, 1.7) + 1.4)
    return np.asarray(positions), np.asarray(residues), np.asarray(radii)


def geometry(holo, target_mask):
    xyz, res, radii = heavy_atoms(holo)
    is_target = target_mask[res]
    if not np.any(~is_target):
        raise ValueError("No protein partner atoms")
    atom_distance = cKDTree(xyz[~is_target]).query(xyz[is_target])[0]
    distances = np.full(len(target_mask), np.inf)
    np.minimum.at(distances, res[is_target], atom_distance)
    if distances.min() > 5:
        raise ValueError("Selected components have no <=5 Angstrom heavy-atom contact")
    sphere = _sphere_points(96)
    sasa_full = _shrake_rupley_residue_sasa(xyz, radii, res, len(target_mask), sphere)
    sasa_iso = _shrake_rupley_residue_sasa(xyz[is_target], radii[is_target], res[is_target], len(target_mask), sphere)
    ca = holo['example']['backbone'][:, 1].numpy()
    ca_valid = holo['example']['backbone_mask'][:, 1].numpy()
    # Missing C-alpha atoms are never counted at the parser's zero placeholder.
    counts = np.full(len(ca), np.nan)
    counts[ca_valid] = [len(x)-1 for x in cKDTree(ca[ca_valid]).query_ball_point(ca[ca_valid], 8.)]
    return distances, sasa_iso, sasa_iso-sasa_full, counts


def prepare(case):
    holo = state(case['holo'])
    target = case['target_chains']
    trials, failures = [], []
    for spec in case.get('apo_candidates', [case.get('apo')]):
        try:
            apo = state(spec)
            mapped = align_states(apo, holo, target)
            # Apo entries with unspecified chains may contain ASU copies. Keep
            # only the matched component, not arbitrary crystallographic mates.
            chosen = dict(spec, chains=mapped[4])
            apo = state(chosen)
            mapped = align_states(apo, holo, target)
            trials.append((mapped[0], apo, mapped))
        except Exception as exc:
            failures.append(str(exc))
    if not trials:
        raise ValueError('; '.join(failures))
    _, apo, mapped = max(trials, key=lambda x: x[0])
    ia, ih = np.asarray(mapped[1]), np.asarray(mapped[2])
    if 'exact_indices' in case:
        actual = set(zip(apo['indices'][ia], holo['indices'][ih]))
        expected = set(zip(*case['exact_indices']))
        if actual != expected:
            raise ValueError("Recomputed mapping differs from local MD manifest")
    target_mask = np.isin(holo['chain'], target)
    selected = np.flatnonzero(target_mask)
    iso = dict(holo, example=subset(holo['example'], selected), indices=holo['indices'][selected],
               spec=dict(holo['spec'], chains=target))
    # Infer B in target-only context but preserve original within-chain indices.
    local_iso = {v: k for k, v in enumerate(selected)}
    ii = np.array([local_iso[v] for v in ih])
    valid = (apo['example']['geometry_mask'][ia].numpy() & holo['example']['geometry_mask'][ih].numpy())
    ia, ih, ii = ia[valid], ih[valid], ii[valid]
    if len(ia) < 10 or not np.array_equal(apo['aa'][ia], holo['aa'][ih]):
        raise ValueError("Insufficient valid exact residue correspondences")
    distance, sasa, burial, contacts = geometry(holo, target_mask)
    metadata = pd.DataFrame(dict(case_id=case['case_id'], direction=case['direction'], group=case['group'],
        amino_acid=holo['aa'][ih], chain=holo['chain'][ih], author_id=holo['author'][ih],
        apo_chain=apo['chain'][ia], apo_author_id=apo['author'][ia],
        apo_index=apo['indices'][ia], holo_index=holo['indices'][ih],
        partner_distance=distance[ih], isolated_sasa=sasa[ih], burial=burial[ih], ca_contacts8=contacts[ih]))
    metadata['region'] = np.where(metadata.partner_distance <= 5, 'interface', np.where(metadata.partner_distance <= 12, 'proximal', 'distal'))
    for metric_index, metric in enumerate(('rmsf', 'chi_entropy')):
        for role, st, ix in [('apo', apo, ia), ('holo', holo, ih)]:
            values = np.full(len(ix), np.nan)
            if st['spec'].get('label_path'):
                y, _, mask = load_label_file(st['spec']['label_path'])
                if len(y) != st['structure'].num_residues:
                    raise ValueError('Label/structure residue count mismatch')
                original = st['indices'][ix]
                values = np.where(mask[original, metric_index], y[original, metric_index], np.nan)
            metadata[f'true_{role}_{metric}'] = values
        metadata['true_delta_' + metric] = metadata['true_holo_' + metric] - metadata['true_apo_' + metric]
    provenance = dict(case=case, selected_apo=apo['spec'], mapping=mapped[3], valid_residues=len(ia),
        apo_residues=len(apo['indices']), holo_residues=len(holo['indices']),
        apo_sha256=digest(apo['spec']['path']), holo_sha256=digest(holo['spec']['path']),
        label_sha256={r:digest(s['spec']['label_path']) for r,s in [('apo',apo),('holo',holo)] if s['spec'].get('label_path')})
    return dict(apo=apo, iso=iso, holo=holo, ia=ia, ii=ii, ih=ih, metadata=metadata, provenance=provenance)


def training_audit(prepared):
    cache = OUT / 'training_sequence_index.json'
    if cache.exists():
        records = json.loads(cache.read_text())
    else:
        records = []
        manifest = ROOT / 'artifacts/evaluation_atlas_hchi_v1/oof/prediction_manifest.jsonl'
        for line in manifest.read_text().splitlines():
            row = json.loads(line)
            z = np.load(row['cache_path'], allow_pickle=False)
            meta = row['record'].get('metadata', {})
            records.append(dict(sample_id=row['sample_id'], held_out_fold=row['fold'],
                pdb_ids=meta.get('pdb_ids', [meta.get('pdb_id', '')]),
                sequences=[''.join(z['amino_acid'][z['chain']==c]) for c in dict.fromkeys(z['chain'].tolist())]))
        write_json(cache, records)
    results = []
    for cid, data in prepared.items():
        pdbs = {data[r]['spec'].get('pdb_id', '').upper() for r in ('apo', 'holo')}
        seqs = {''.join(data[r]['aa'][data[r]['chain']==c]) for r in ('apo','holo') for c in dict.fromkeys(data[r]['chain'].tolist())}
        matches = []
        for rec in records:
            by_pdb = bool(pdbs.intersection(x.upper() for x in rec['pdb_ids'] if x))
            by_seq = bool(seqs.intersection(rec['sequences']))
            if by_pdb or by_seq:
                matches.append(dict(sample_id=rec['sample_id'], held_out_fold=rec['held_out_fold'], exact_pdb=by_pdb, exact_sequence=by_seq))
        clean = [fold for fold in range(5) if all(x['held_out_fold']==fold for x in matches)]
        results.append(dict(case_id=cid, exact_matches=matches, folds_without_detected_exact_overlap=clean,
                            limitation='Exact PDB/chain sequence only; homology audit pending. Five-fold ensemble is descriptive, not OOF.'))
    write_json(OUT / 'training_overlap.json', results)


def infer(model, st):
    with torch.inference_mode():
        prediction = model.forward_state(collate_states([st['example']]))
    values = prediction['mean'][0].cpu().numpy()
    values[~prediction['mask'][0].cpu().numpy()] = np.nan
    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', type=int)
    parser.add_argument('--folds', default='0,1,2,3,4')
    args = parser.parse_args()
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    cases = json.loads((OUT / 'cohorts.json').read_text(encoding='utf-8'))
    if args.limit:
        cases = cases[:args.limit]
    prepared, statuses = {}, []
    for case in cases:
        try:
            if case.get('status') == 'unavailable':
                raise ValueError(case['error'])
            prepared[case['case_id']] = prepare(case)
            write_json(OUT / 'mapping' / (case['case_id'] + '.json'), prepared[case['case_id']]['provenance'])
            statuses.append(dict(case_id=case['case_id'], direction=case['direction'], status='prepared'))
        except Exception as exc:
            statuses.append(dict(case_id=case['case_id'], direction=case['direction'], status='failed', error=str(exc)))
        print('prepare', statuses[-1], flush=True)
    write_json(OUT / 'run_status.json', statuses)
    training_audit(prepared)
    code_hash = hashlib.sha256((implementation_hash()+digest(__file__)).encode()).hexdigest()
    frames, checkpoints = [], []
    pred_dir = OUT / 'state_predictions'
    pred_dir.mkdir(exist_ok=True)
    for fold in map(int, args.folds.split(',')):
        model, config, checkpoint, checksum = load_fold(ROOT, fold, torch.device('cpu'), True)
        checkpoints.append(dict(fold=fold, path=str(checkpoint), sha256=checksum, implementation_sha256=code_hash))
        for cid, data in prepared.items():
            row = data['metadata'].copy()
            row['fold'] = fold
            values = {}
            for role, ix in [('apo','ia'), ('iso','ii'), ('holo','ih')]:
                st = data[role]
                identity = dict(path_sha256=digest(st['spec']['path']), original_indices=st['indices'].tolist(),
                                checkpoint=checksum, implementation=code_hash)
                key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
                path = pred_dir / (key + '.npz')
                if path.exists():
                    prediction = np.load(path, allow_pickle=False)['mean']
                else:
                    prediction = infer(model, st)
                    np.savez_compressed(path, mean=prediction, identity_json=json.dumps(identity, sort_keys=True))
                values[role] = prediction[data[ix]]
                for j, metric in enumerate(('rmsf','chi_entropy')):
                    row[role+'_'+metric] = values[role][:,j]
            for metric in ('rmsf','chi_entropy'):
                row['delta_'+metric] = row['holo_'+metric]-row['apo_'+metric]
                row['context_'+metric] = row['holo_'+metric]-row['iso_'+metric]
                row['geometry_'+metric] = row['iso_'+metric]-row['apo_'+metric]
                assert np.allclose(row['delta_'+metric], row['context_'+metric]+row['geometry_'+metric], atol=1e-6, equal_nan=True)
            row.loc[row.amino_acid.isin(['A','G']), [c for c in row if 'chi_entropy' in c]] = np.nan
            frames.append(row)
            print(f'fold {fold}: {cid} ({len(row)} paired residues)', flush=True)
        del model
    pd.concat(frames, ignore_index=True).to_csv(OUT / 'residue_predictions.csv', index=False)
    write_json(OUT / 'checkpoints.json', checkpoints)
    write_json(OUT / 'run_status.json', [dict(x, status='complete') if x['status']=='prepared' else x for x in statuses])
    print('DONE', len(prepared), 'cases,', len(frames), 'case-fold evaluations', flush=True)


if __name__ == '__main__':
    main()
