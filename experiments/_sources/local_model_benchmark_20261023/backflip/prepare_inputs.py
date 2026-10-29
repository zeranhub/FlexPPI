"""Freeze identity-preserving BackFlip whole-complex and independent-chain inputs.

Only reads the existing corrected cohort.  No MD, training, or label adjustment.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
from collections import Counter

import numpy as np

ROOT = Path(__file__).resolve().parent
CORE = Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
sys.path.insert(0, str(CORE / 'training_payload/src'))
from flexppi.data.structures import parse_pdb, target_mask


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_records():
    original = CORE / 'training_payload/splits/all_ready.jsonl'
    rows = [json.loads(x) for x in original.read_text(encoding='utf-8').splitlines()]
    for row in rows:
        row['role'] = 'training'
        for key, kind in [('structure_path', 'structures'), ('label_path', 'labels')]:
            row[key] = str(CORE / 'training_payload/dataset' / kind / Path(row[key]).name)
    return rows


def freeze(row):
    sample = row['sample_id']
    folder = ROOT / 'inputs' / sample
    folder.mkdir(parents=True, exist_ok=True)
    structure = parse_pdb(row['structure_path'])
    chains = structure.unique_chains()
    assert len(chains) <= 26
    chain_mapping = {chain: chr(65+i) for i, chain in enumerate(chains)}
    mapping = []
    keep = np.flatnonzero(structure.frame_mask)
    rank = {int(idx): pos for pos, idx in enumerate(keep)}
    for i in range(structure.num_residues):
        chain = structure.chain_ids[i]
        mapping.append(dict(record_index=i, input_index=rank.get(i, -1),
            original_chain=chain, original_author_residue_id=structure.author_residue_ids[i],
            amino_acid=structure.sequence[i], residue_name=structure.residue_names[i],
            canonical_chain=chain_mapping[chain],
            canonical_residue_id=int(structure.residue_index[i])+1,
            valid_frame=bool(structure.frame_mask[i])))
    with (folder / 'residue_mapping.csv').open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(mapping[0]))
        w.writeheader()
        w.writerows(mapping)
    contexts = []
    for name, chain in [('whole_complex', None)] + [('chain_'+chain_mapping[c], c) for c in chains]:
        indices = keep if chain is None else keep[np.array(structure.chain_ids)[keep] == chain]
        if not len(indices):
            continue
        text = []
        last_chain = None
        for i in indices:
            info = mapping[int(i)]
            ch = info['canonical_chain'] if chain is None else 'A'
            if last_chain is not None and ch != last_chain:
                text.append('TER\n')
            last_chain = ch
            for atom_index, atom in enumerate(['N', 'CA', 'C', 'O']):
                if not structure.backbone_mask[i, atom_index]:
                    continue
                x,y,z = structure.backbone[i, atom_index]
                element = atom[0]
                text.append(f'ATOM  {len(text)+1:5d} {atom:^4s} {info["residue_name"]:3s} {ch}{info["canonical_residue_id"]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n')
        path = folder / (name+'.pdb')
        path.write_text(''.join(text)+'TER\nEND\n', encoding='ascii')
        parsed = parse_pdb(str(path))
        assert parsed.sequence == ''.join(structure.sequence[int(i)] for i in indices)
        assert parsed.frame_mask.all()
        assert np.array_equal(parsed.backbone_mask, structure.backbone_mask[indices])
        assert np.array_equal(parsed.backbone, structure.backbone[indices])
        contexts.append(dict(context=name, original_chain=chain,
            pdb_path=str(path), sha256=sha(path), record_indices=indices.tolist()))
    with np.load(row['label_path']) as z:
        y = z['values'][:,0]
        valid = z['mask'][:,0] & structure.frame_mask & target_mask(structure, row['target_chains'])
        assert len(y) == structure.num_residues
    np.savez_compressed(folder/'reference.npz', rmsf=y, valid=valid,
        frame_mask=structure.frame_mask, sequence=np.array(list(structure.sequence)),
        chain=np.array(structure.chain_ids), author_residue_id=np.array(structure.author_residue_ids))
    return dict(sample_id=sample, role=row['role'], dataset=row['dataset'],
        fold=row['metadata'].get('fold',-1),
        biological_system=row['metadata'].get('biological_system',row.get('source_group')),
        n_residues=structure.num_residues, n_valid_rmsf=int(valid.sum()),
        n_chains=len(chains), original_structure_path=row['structure_path'],
        original_structure_sha256=sha(row['structure_path']), label_sha256=sha(row['label_path']),
        reference_path=str(folder/'reference.npz'), mapping_path=str(folder/'residue_mapping.csv'),
        mapping_sha256=sha(folder/'residue_mapping.csv'),
        insertion_code_residues=sum(any(c.isalpha() for c in a) for a in structure.author_residue_ids),
        contexts=contexts)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--smoke-only', action='store_true')
    args = p.parse_args()
    rows = load_records()
    if args.smoke_only:
        rows = [r for r in rows if r['sample_id'] == 'itrace_1ao7']
    manifest = []
    for i,row in enumerate(rows):
        manifest.append(freeze(row))
        if (i+1)%100 == 0:
            print('Frozen',i+1,len(rows),flush=True)
    name = 'smoke_manifest.json' if args.smoke_only else 'manifest.json'
    (ROOT/name).write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    audit = dict(records=len(manifest), datasets=dict(Counter(x['dataset'] for x in manifest)),
        chains=sum(x['n_chains'] for x in manifest), valid_rmsf=sum(x['n_valid_rmsf'] for x in manifest),
        label_definition='Existing corrected RMSF labels, per-chain Kabsch; no trajectory recomputation.',
        input_definition='All valid N/CA/C frames; chain names remapped A-Z and insertion codes retained by temporary integer renumbering with original identity map.',
        sampling='No selection by prediction errors. Smoke record fixed before inference: itrace_1ao7.',
        code_sha256=sha(__file__), manifest_sha256=sha(ROOT/name))
    (ROOT/('smoke_input_audit.json' if args.smoke_only else 'input_audit.json')).write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2),flush=True)


if __name__ == '__main__':
    main()
