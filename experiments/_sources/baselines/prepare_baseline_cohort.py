"""Create identity-checked complete-structure and per-chain baseline inputs.

The canonical PDB keeps all protein heavy atoms from the frozen first model.
Only atom serialization/altloc duplication changes; original coordinates,
residue order, chain identities and author numbering are retained. Truth is
the exact original valid-mask intersection, independent of baseline outputs.
"""
import argparse
import csv
import gzip
import hashlib
import json
import sys
import importlib.util
import types
from pathlib import Path
import numpy as np


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical_atoms(path, structure):
    residue_keys = {(c, a): i for i, (c, a) in enumerate(zip(structure.chain_ids, structure.author_residue_ids))}
    atoms, seen = [[] for _ in structure.sequence], set()
    for line in Path(path).read_text(encoding='utf-8').splitlines(True):
        if line.startswith('ENDMDL'):
            break
        if not line.startswith(('ATOM  ', 'HETATM')) or len(line) < 54 or line[16] not in (' ', 'A', '1'):
            continue
        chain = line[21].strip() or '_'
        author = line[22:26].strip() + line[26].strip()
        key = (chain, author)
        if key not in residue_keys:
            continue
        index = residue_keys[key]
        name = line[12:16].strip()
        if (index, name) in seen:
            continue
        # Hydrogens are excluded equally from all relative-score inputs.
        element = line[76:78].strip().upper() if len(line) >= 78 else ''
        if element in ('H', 'D') or (not element and name.lstrip('0123456789').startswith('H')):
            continue
        seen.add((index, name))
        atoms[index].append(line)
    assert all(any(line[12:16].strip() == 'CA' for line in group) for group in atoms)
    return atoms


def write_pdb(path, atoms, indices, chain_override=None):
    lines = []
    for j, index in enumerate(indices):
        for line in atoms[index]:
            output = 'ATOM  ' + f'{len(lines)+1:5d}' + line[11:16] + ' ' + line[17:]
            if chain_override is not None:
                output = output[:21] + chain_override + f'{j+1:4d}' + ' ' + output[27:]
            lines.append(output.rstrip('\r\n') + '\n')
    path.write_text(''.join(lines) + 'END\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--training-root', required=True)
    parser.add_argument('--source-root', required=True, help='Directory containing the frozen flexppi package')
    parser.add_argument('--output', required=True)
    parser.add_argument('--runtime-output-prefix', help='Optional Linux output prefix for manifests generated on Windows')
    args = parser.parse_args()
    root, output = Path(args.training_root), Path(args.output)
    if (output / 'cohort_freeze.json').exists():
        raise RuntimeError('Frozen cohort already exists')
    # Import the frozen dependency-free parser without executing flexppi's
    # package initializer, which imports the GPU model and therefore PyTorch.
    source = Path(args.source_root)
    for name in ['flexppi', 'flexppi.data']:
        package = types.ModuleType(name)
        package.__path__ = []
        sys.modules[name] = package
    for name, filename in [('flexppi.constants', 'flexppi/constants.py'), ('flexppi.data.structures', 'flexppi/data/structures.py')]:
        specification = importlib.util.spec_from_file_location(name, source / filename)
        module = importlib.util.module_from_spec(specification)
        sys.modules[name] = module
        specification.loader.exec_module(module)
    parse_pdb = sys.modules['flexppi.data.structures'].parse_pdb
    target_mask = sys.modules['flexppi.data.structures'].target_mask
    records = [json.loads(line) for line in (root / 'splits/all_ready.jsonl').read_text(encoding='utf-8').splitlines()]
    training_freeze = json.loads((root / 'training_freeze.json').read_text(encoding='utf-8'))
    frozen_hashes = {key.replace('\\', '/'): value for key, value in training_freeze['input_hashes'].items()}
    assert len(records) == 3063 and len({r['sample_id'] for r in records}) == 3063
    for folder in ['structures', 'chains']:
        (output / folder).mkdir(parents=True, exist_ok=True)
    runtime = args.runtime_output_prefix.rstrip('/') if args.runtime_output_prefix else output.as_posix()
    manifest, chains, record_index = [], [], []
    truth_fields = ['sample_id', 'dataset', 'biological_system', 'fold', 'residue_index', 'amino_acid', 'chain', 'author_residue_id', 'input_id', 'input_index', 'target_rmsf', 'valid_rmsf', 'target_chi_entropy', 'valid_chi_entropy']
    counts, offset = np.zeros(2, dtype=np.int64), 0
    with gzip.open(output / 'truth_residues.csv.gz', 'wt', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=truth_fields)
        writer.writeheader()
        for position, record in enumerate(records):
            sample = record['sample_id']
            structure_path = root / 'dataset/structures' / (sample + '.pdb')
            label_filename = record['label_path'].replace('\\', '/').split('/')[-1]
            label_path = root / 'dataset/labels' / label_filename
            structure = parse_pdb(str(structure_path))
            atoms = canonical_atoms(structure_path, structure)
            key = 'record_' + hashlib.sha256(sample.encode()).hexdigest()[:24]
            pdb = output / 'structures' / (key + '.pdb')
            write_pdb(pdb, atoms, range(structure.num_residues))
            canonical = parse_pdb(str(pdb))
            assert canonical.sequence == structure.sequence and canonical.author_residue_ids == structure.author_residue_ids
            assert canonical.chain_ids == structure.chain_ids
            assert np.array_equal(canonical.backbone_mask, structure.backbone_mask)
            assert np.array_equal(canonical.backbone, structure.backbone)
            with np.load(label_path) as label:
                target = label['values'][:, :2].copy()
                valid = label['mask'][:, :2] & structure.frame_mask[:, None] & target_mask(structure, record['target_chains'])[:, None]
            assert target.shape == valid.shape == (structure.num_residues, 2)
            assert np.isfinite(target[valid]).all()
            assert sha(label_path) == frozen_hashes['dataset/labels/' + label_filename]
            assert sha(structure_path) == frozen_hashes['dataset/structures/' + structure_path.name]
            counts += valid.sum(axis=0)
            manifest.append(dict(sample_id=sample, cohort=record['dataset'], biological_system=record['metadata']['biological_system'],
                fold=record['metadata']['fold'], container_pdb_path=runtime + '/structures/' + pdb.name,
                n_residues=structure.num_residues, n_chains=len(structure.unique_chains()), input_sha256=sha(pdb),
                original_input_sha256=sha(structure_path), label_sha256=sha(label_path)))
            record_index.append(dict(sample_id=sample, dataset=record['dataset'], biological_system=record['metadata']['biological_system'],
                fold=record['metadata']['fold'], start=offset, stop=offset+structure.num_residues))
            chain_mapping = {}
            for chain in structure.unique_chains():
                indices = structure.indices_for_chain(chain)
                input_id = 'chain_' + hashlib.sha256((sample + '|' + chain).encode()).hexdigest()[:24]
                chain_pdb = output / 'chains' / (input_id + '.pdb')
                write_pdb(chain_pdb, atoms, indices, chain_override='A')
                isolated = parse_pdb(str(chain_pdb))
                assert isolated.sequence == structure.chain_sequence(chain)
                assert np.array_equal(isolated.backbone, structure.backbone[indices])
                chains.append(dict(input_id=input_id, sample_id=sample, chain=chain, n_residues=len(indices),
                    sequence=isolated.sequence, pdb_file=runtime + '/chains/' + chain_pdb.name,
                    input_sha256=sha(chain_pdb)))
                for input_index, index in enumerate(indices):
                    chain_mapping[int(index)] = (input_id, input_index)
            for index in range(structure.num_residues):
                input_id, input_index = chain_mapping[index]
                writer.writerow(dict(sample_id=sample, dataset=record['dataset'], biological_system=record['metadata']['biological_system'],
                    fold=record['metadata']['fold'], residue_index=index, amino_acid=structure.sequence[index],
                    chain=structure.chain_ids[index], author_residue_id=structure.author_residue_ids[index],
                    input_id=input_id, input_index=input_index, target_rmsf=float(target[index, 0]), valid_rmsf=int(valid[index, 0]),
                    target_chi_entropy=float(target[index, 1]), valid_chi_entropy=int(valid[index, 1])))
            offset += structure.num_residues
            if (position + 1) % 250 == 0:
                print('Prepared', position + 1, len(records), flush=True)
    for filename, rows, delimiter in [('structure_manifest.tsv', manifest, '\t'), ('chain_manifest.csv', chains, ',')]:
        with (output / filename).open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
    (output / 'record_index.json').write_text(json.dumps(record_index, indent=2), encoding='utf-8')
    assert counts.tolist() == [961321, 824988]
    assert len(chains) == 4973
    report = dict(status='frozen', records=3063, chains=len(chains), total_residues=offset,
        valid_rmsf=int(counts[0]), valid_chi_entropy=int(counts[1]), corrected_labels=True,
        protocol='Complete original protein assembly for ProDy/RMSF-net; independent isolated chain for PLANET-MD; record-level identity and common valid residues preserved.',
        canonicalization='First model, known amino-acid residues, unique heavy atoms, accepted altlocs blank/A/1; exact sequence, chain, author numbering and backbone coordinate validation.',
        source_manifest_sha256=sha(root / 'splits/all_ready.jsonl'), source_freeze_sha256=sha(root / 'training_freeze.json'),
        code_sha256=sha(__file__), runtime_prefix=runtime,
        files={filename: sha(output / filename) for filename in ['structure_manifest.tsv', 'chain_manifest.csv', 'truth_residues.csv.gz', 'record_index.json']})
    (output / 'cohort_freeze.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
