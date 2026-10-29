"""Verify the local prepared archive against the committed remote audit hashes."""
from pathlib import Path
import datetime
import hashlib
import json

BATCH = Path(r'E:\Project\FlexPPI\experiments\_sources\kix_md_triplicate_rebuild_20261010')


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    prepared = BATCH / 'prepared'
    ready = json.loads((prepared / 'READY.json').read_text())
    audit = json.loads((prepared / 'input_audit.json').read_text())
    repair = json.loads((prepared / 'repair_audit.json').read_text())
    assert audit['status'] == 'all_four_rebuilt_inputs_passed'
    for state in ('apo', 'myb', 'mll', 'both'):
        folder = prepared / state
        meta = json.loads((folder / 'metadata.json').read_text())
        assert sha(folder / 'system.xml') == ready['states'][state]['system_sha256'] == audit['states'][state]['system_sha256'] == meta['system_sha256']
        assert sha(folder / 'protein.pdb') == ready['states'][state]['protein_pdb_sha256'] == audit['states'][state]['protein_sha256'] == meta['protein_pdb_sha256']
        assert sha(folder / 'solvated.pdb') == audit['states'][state]['solvated_pdb_sha256']
        assert audit['states'][state]['passed']
    files = {str(path.relative_to(BATCH)).replace('\\', '/'):
             {'bytes': path.stat().st_size, 'sha256': sha(path)}
             for path in sorted(prepared.rglob('*')) if path.is_file()}
    assert len(files) == 20, f'Unexpected prepared archive file count: {len(files)}'
    report = dict(status='prepared_cpu_rebuilt_and_audited_local_copy_verified',
                  created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  batch='KIX_rebuilt_triplicate_20261010', old_rep01_included=False,
                  physical_protocol=ready['parameters'], openmm_version=ready['openmm_version'],
                  environment_interpreter='/root/flexppi_completion_20261010/md_conda/bin/python',
                  remote_prepared='/root/flexppi_completion_20261010/kix_md_triplicate_rebuild_20261010/prepared',
                  cpu_prepare_environment={'CUDA_VISIBLE_DEVICES': '', 'OPENMM_CPU_THREADS': '4'},
                  this_agent_created_gpu_context=False,
                  gpu_smoke_and_production='Root agent owns subsequent GPU smoke and production; this is a CPU preparation handoff only.',
                  existing_heavy_atom_repair_displacement_nm=repair['existing_atom_max_displacement_nm'],
                  state_input_audit=audit['states'], prepared_files=files)
    (BATCH / 'PREPARED_HANDOFF.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'files': len(files),
                      'total_bytes': sum(row['bytes'] for row in files.values()),
                      'max_preserved_heavy_atom_deviation_angstrom': max(row['max_coordinate_deviation_after_translation_angstrom'] for row in audit['states'].values())}, indent=2))


if __name__ == '__main__':
    main()
