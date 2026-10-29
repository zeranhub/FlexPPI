"""PLANET-MD released full/sequence checkpoints, independently per frozen chain.

ESM3 is loaded once. Structure feature stage is read from the released
checkpoint rather than assumed from CLI defaults. No training or calibration.
This file is prepared for root-agent GPU dispatch; it does not launch itself.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import numpy as np


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--source-root', required=True)
    parser.add_argument('--weights', required=True)
    parser.add_argument('--variant', choices=['full', 'seq'], required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--retry-failed', action='store_true')
    args = parser.parse_args()
    os.environ.setdefault('HF_HOME', str(Path(args.source_root).parent / 'hf_cache'))
    os.environ.setdefault('HF_HUB_OFFLINE', '1')
    sys.path.insert(0, args.source_root)
    import torch
    from biotite.structure.io import pdb
    from planet_md import PlanetMD
    from planet_md.esm3 import get_model, get_structure_vae, get_tokenizers
    from planet_md.features import esm3_sequence, esm3_vqvae
    from planet_md.structure.protein_chain import ProteinChain
    cohort, output = Path(args.cohort), Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    freeze = json.loads((cohort / 'cohort_freeze.json').read_text())
    manifest_path = cohort / 'chain_manifest.csv'
    assert sha(manifest_path) == freeze['files']['chain_manifest.csv']
    checkpoint = torch.load(args.weights, map_location='cpu', weights_only=True)
    parameters = checkpoint['hyper_parameters']
    assert bool(parameters['struct_features']) == (args.variant == 'full')
    stage = parameters.get('struct_stage', 'quantized')
    metadata = dict(manifest_sha256=sha(manifest_path), weights_sha256=sha(args.weights),
        code_sha256=sha(__file__), variant=args.variant, structure_feature_stage=stage, device=args.device,
        input='Independent chain from corrected frozen first static structure; exact original residue mapping',
        temperature_policy='Released default300K for every input, preserving the static sequence/structure-only baseline protocol',
        esm_loading_policy='Load the exact released float32 weights on CPU before transfer to inference device, avoiding simultaneous duplicate GPU state dictionaries; no dtype or weight changes',
        rmsf_native_unit='nm from unscaled MDTraj md.rmsf labels in released preprocessing', rmsf_conversion='10*native_nm to Angstrom, no fitted calibration',
        limitations='Released weights may overlap ATLAS/mdCATH; independent-chain input does not encode partners; SHP is not the chi-state entropy target.')
    metadata_path = output / 'metadata.json'
    if metadata_path.exists():
        assert json.loads(metadata_path.read_text()) == metadata
    else:
        metadata_path.write_text(json.dumps(metadata, indent=2))
    device = torch.device(args.device)
    model = PlanetMD.load_from_checkpoint(str(Path(args.weights))).to(device).eval()
    del checkpoint
    esm_model = get_model(device='cpu').to(device).eval()
    tokenizers = get_tokenizers()
    structure_model = get_structure_vae(device=device).to(device).eval() if args.variant == 'full' else None
    with manifest_path.open(encoding='utf-8', newline='') as handle:
        entries = list(csv.DictReader(handle))
    if args.limit is not None:
        entries = entries[:args.limit]
    for position, entry in enumerate(entries):
        input_id = entry['input_id']
        path, status_path = output / (input_id+'.npz'), output / (input_id+'.json')
        if path.exists() and status_path.exists():
            old = json.loads(status_path.read_text())
            if old.get('status') == 'complete':
                assert old['input_sha256'] == entry['input_sha256'] and sha(path) == old['output_sha256']
                with np.load(path) as saved:
                    assert saved['rmsf_angstrom'].shape == (int(entry['n_residues']),)
                    assert np.isfinite(saved['rmsf_angstrom']).all()
                continue
        if status_path.exists() and json.loads(status_path.read_text()).get('status') == 'failed' and not args.retry_failed:
            continue
        started = time.perf_counter()
        status = dict(input_id=input_id, sample_id=entry['sample_id'], chain=entry['chain'], input_sha256=entry['input_sha256'], status='running')
        try:
            assert sha(entry['pdb_file']) == entry['input_sha256']
            sequence = entry['sequence']
            assert len(sequence) == int(entry['n_residues'])
            with torch.inference_mode():
                sequence_features = esm3_sequence(sequence, esm_model, tokenizers)[:, 1:-1, :]
                assert sequence_features.shape[1] == len(sequence)
                structure_features = None
                if structure_model is not None:
                    structure = pdb.PDBFile.read(entry['pdb_file']).get_structure(model=1)
                    chain = ProteinChain.from_atomarray(structure)
                    assert chain.sequence == sequence
                    assert stage in ['encoded', 'pre-quantized', 'quantized'], 'Unsupported released feature stage'
                    structure_features = esm3_vqvae(chain, structure_model, stage=stage)
                    assert structure_features.shape[1] == len(sequence)
                prediction = model({'seq_feats': sequence_features, 'struct_feats': structure_features})
                native = prediction['rmsf'].reshape(-1).cpu().numpy()
            assert native.shape == (len(sequence),) and np.isfinite(native).all()
            temporary = output / (input_id+'.tmp.npz')
            np.savez_compressed(temporary, rmsf_native_nm=native, rmsf_angstrom=10.0*native, sequence=np.asarray(sequence))
            os.replace(temporary, path)
            status.update(status='complete', residues=len(sequence), output_sha256=sha(path))
        except Exception as error:
            status.update(status='failed', error=type(error).__name__ + ': ' + str(error))
            if device.type == 'cuda':
                torch.cuda.empty_cache()
        status['seconds'] = time.perf_counter()-started
        status_path.write_text(json.dumps(status, indent=2))
        print(position+1, len(entries), input_id, status['status'], f"{status['seconds']:.3f}s", flush=True)


if __name__ == '__main__':
    main()
