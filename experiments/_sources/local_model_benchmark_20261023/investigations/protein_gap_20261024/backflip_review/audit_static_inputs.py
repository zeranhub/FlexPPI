"""Read-only static-input/source audit. No network, torch, model, or performance reads."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import importlib.util
import json
import sys
import types
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CORE = Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
OUT = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024**2),b''):
            digest.update(block)
    return digest.hexdigest()


def csv_rows(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        return list(csv.DictReader(handle))


def main():
    assert ROOT.name=='local_model_benchmark_20261023'
    source = CORE/'training_payload/src'
    for name in ('flexppi','flexppi.data'):
        module = types.ModuleType(name)
        module.__path__ = []
        sys.modules[name] = module
    for name,relative in (('flexppi.constants','flexppi/constants.py'),
                           ('flexppi.data.structures','flexppi/data/structures.py')):
        spec = importlib.util.spec_from_file_location(name,source/relative)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    parse_pdb = sys.modules['flexppi.data.structures'].parse_pdb
    records = csv_rows(ROOT/'cohort/records.csv')
    bf = json.loads((ROOT/'backflip/manifest.json').read_text(encoding='utf-8-sig'))
    bf_by = {r['sample_id']:r for r in bf}
    assert len(records)==len(bf)==len(bf_by)==3063
    manifest_path = CORE/'neural_exports/evaluation/prediction_manifest.jsonl'
    full = {}
    for line in manifest_path.open(encoding='utf-8-sig'):
        row = json.loads(line)
        if row['variant']=='full':
            assert row['sample_id'] not in full
            full[row['sample_id']] = row
    assert set(full)==set(bf_by)=={r['sample_id'] for r in records}
    freeze_path = CORE/'training_payload/training_freeze.json'
    freeze = json.loads(freeze_path.read_text(encoding='utf-8-sig'))
    frozen = {k.replace('\\','/'):v for k,v in freeze['input_hashes'].items()}
    metadata_path = ROOT/'backflip/results/backflip-2.1/run_metadata.json'
    metadata = json.loads(metadata_path.read_text(encoding='utf-8-sig'))
    sources = [ROOT/'cohort/records.csv',ROOT/'backflip/manifest.json',ROOT/'backflip/input_audit.json',
        ROOT/'backflip/prepare_inputs.py',ROOT/'backflip/run_benchmark.py',metadata_path,
        ROOT/'backflip/weights/backflip-2.1/config.yaml',ROOT/'backflip/weights/backflip-2.1/backflip-2.1.ckpt',
        source/'flexppi/constants.py',source/'flexppi/data/structures.py',source/'flexppi/data/dataset.py',
        source/'flexppi/model/network.py',freeze_path,manifest_path,
        ROOT/'backflip/upstream_bf2_887c6a7.zip',ROOT/'backflip/dependency_patch_audit.json']
    original_hashes = {str(p):sha(p) for p in sources}
    assert original_hashes[str(sources[1])]==metadata['input_manifest_sha256']
    assert original_hashes[str(ROOT/'backflip/run_benchmark.py')]==metadata['runner_sha256']
    assert original_hashes[str(ROOT/'backflip/weights/backflip-2.1/config.yaml')]==metadata['checkpoint_config_sha256']
    assert original_hashes[str(ROOT/'backflip/weights/backflip-2.1/backflip-2.1.ckpt')]==metadata['checkpoint_sha256']
    numeric_source = Path(metadata['source_path'])
    names = ['backflip/deployment/utils.py','backflip/deployment/inference_class.py',
        'backflip/models/protein_embedder.py','backflip/models/node_embedder.py',
        'backflip/models/edge_embedder.py','backflip/models/flexibility_model.py']
    official_rows = []
    with zipfile.ZipFile(ROOT/'backflip/upstream_bf2_887c6a7.zip') as archive:
        for relative in names:
            path = numeric_source/relative
            if not path.exists():
                continue
            matches = [n for n in archive.namelist() if n.endswith('/'+relative)]
            assert len(matches)==1
            current = sha(path)
            assert hashlib.sha256(archive.read(matches[0])).hexdigest()==current,relative
            original_hashes[str(path)] = current
            official_rows.append(dict(path=str(path),sha256=current,downloaded_official_archive_bytes_equal=True))
    results,failures = [],[]
    totals = dict(canonical_context_files=0,canonical_atoms=0,canonical_nonzero_bfactor_atoms=0,
        canonical_nonunit_occupancy_atoms=0,source_backbone_atoms=0,source_nonzero_backbone_bfactor_atoms=0,
        complete_geometry_positions=0,canonical_coordinate_positions=0)
    for position,record in enumerate(records,1):
        sample = record['sample_id']
        row = bf_by[sample]
        report = dict(sample_id=sample,dataset=record['dataset'],n_chains=int(record['n_chains']),
            n_residues=int(record['n_residues']),n_valid_rmsf=int(record['n_valid_rmsf']))
        try:
            original = Path(row['original_structure_path'])
            actual_original_sha = sha(original)
            expected = {row['original_structure_sha256'],record['original_input_sha256'],full[sample]['input_sha256'],
                        frozen['dataset/structures/'+original.name]}
            assert expected=={actual_original_sha},(sample,'different original structure hashes')
            report.update(original_input_sha256=actual_original_sha,
                full_and_backflip_frozen_source_SHA_identical=True,source_file=str(original))
            protein = parse_pdb(str(original))
            assert protein.num_residues==row['n_residues']==int(record['n_residues'])
            baseline_path = Path(record['structure_path'])
            assert sha(baseline_path)==record['structure_sha256']
            baseline = parse_pdb(str(baseline_path))
            assert baseline.sequence==protein.sequence and baseline.chain_ids==protein.chain_ids and baseline.author_residue_ids==protein.author_residue_ids
            assert np.array_equal(baseline.backbone_mask,protein.backbone_mask)
            assert np.array_equal(baseline.backbone,protein.backbone)
            report['common_benchmark_backbone_exactly_equal_to_full_source'] = True
            mapping = csv_rows(row['mapping_path'])
            assert sha(row['mapping_path'])==row['mapping_sha256'] and len(mapping)==protein.num_residues
            for i,m in enumerate(mapping):
                assert int(m['record_index'])==i and m['original_chain']==protein.chain_ids[i]
                assert m['original_author_residue_id']==protein.author_residue_ids[i] and m['amino_acid']==protein.sequence[i]
                assert (m['valid_frame']=='True')==bool(protein.frame_mask[i])
            complete = int(protein.frame_mask.sum())
            totals['complete_geometry_positions'] += complete
            for chunk in row['contexts']:
                path = Path(chunk['pdb_path'])
                assert sha(path)==chunk['sha256']
                indices = np.asarray(chunk['record_indices'],dtype=int)
                canonical = parse_pdb(str(path))
                assert canonical.sequence==''.join(protein.sequence[i] for i in indices)
                assert canonical.frame_mask.all()
                assert np.array_equal(canonical.backbone_mask,protein.backbone_mask[indices])
                assert np.array_equal(canonical.backbone,protein.backbone[indices])
                for line in path.read_text(encoding='utf-8').splitlines():
                    if line.startswith(('ATOM  ','HETATM')):
                        totals['canonical_atoms'] += 1
                        totals['canonical_nonzero_bfactor_atoms'] += int(float(line[60:66])!=0.)
                        totals['canonical_nonunit_occupancy_atoms'] += int(float(line[54:60])!=1.)
                totals['canonical_context_files'] += 1
                totals['canonical_coordinate_positions'] += len(indices)
            for line in original.read_text(encoding='utf-8',errors='replace').splitlines():
                if line.startswith(('ATOM  ','HETATM')) and len(line)>=66 and line[12:16].strip() in ('N','CA','C','O'):
                    totals['source_backbone_atoms'] += 1
                    totals['source_nonzero_backbone_bfactor_atoms'] += int(float(line[60:66])!=0.)
            report.update(backflip_all_canonical_context_SHA_and_coordinates_exact=True,
                canonical_context_files=len(row['contexts']),complete_geometry_positions=complete,
                canonical_bfactors_all_zero=True,status='pass')
        except Exception as error:
            report.update(status='failed',error=repr(error))
            failures.append(dict(sample_id=sample,error=repr(error)))
        results.append(report)
        if position%250==0:
            print(json.dumps(dict(progress=position,records=len(records),failures=len(failures))),flush=True)
    assert all(sha(Path(p))==h for p,h in original_hashes.items()),'Evidence source changed during read-only audit'
    assert not (OUT/'static_input_per_record_audit.csv').exists()
    fields = list(dict.fromkeys(k for r in results for k in r))
    with (OUT/'static_input_per_record_audit.csv').open('x',encoding='utf-8-sig',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(results)
    output = dict(created_utc=datetime.now(timezone.utc).isoformat(),status='PASS' if not failures else 'FAIL',
        scope='Read-only all3063 frozen PDB/source-hash/coordinate/Bfactor audit. No model or feature forward was executed.',
        records=len(records),passed_records=sum(r['status']=='pass' for r in results),failures=failures,
        counts=totals,BackFlip_2_1=dict(checkpoint_sha256=metadata['checkpoint_sha256'],
            config_sha256=metadata['checkpoint_config_sha256'],source_commit=metadata['source_commit'],
            source_archive_checks=official_rows,expected_input_features=['N/CA/C-derived translations','N/CA/C-derived rotations',
                'amino-acid identity','all-ones inference residue mask','default concatenated residue ordinal index'],
            excluded_from_model_input=['PDB Bfactor','occupancy','MD label values','reference NPZ arrays','source dataset label','chain identity']),
        no_network=True,no_torch_import=True,no_inference=True,no_original_files_modified=True,
        source_sha256=[dict(path=p,sha256=h) for p,h in sorted(original_hashes.items())],
        output_per_record_sha256=sha(OUT/'static_input_per_record_audit.csv'))
    with (OUT/'static_input_audit.json').open('x',encoding='utf-8') as handle:
        json.dump(output,handle,indent=2,ensure_ascii=False)
    print(json.dumps(dict(status=output['status'],records=len(records),passed=output['passed_records'],
        failures=failures,counts=totals,output=str(OUT))),flush=True)


if __name__=='__main__':
    main()
