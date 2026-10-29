"""Read-only frozen payload/unit audit for the isolated PLANET CPU smoke."""
from pathlib import Path
import csv
import datetime
import hashlib
import importlib.metadata as metadata
import json
import os
import sys
import torch

ROOT = Path('/root/flexppi_completion_20261010/baselines')
OUT = ROOT / 'planet_cpu_smoke_audit_20261010'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def package_version(name):
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return 'not installed (preprocessing-only dependency not required for this inference adapter)'


def main():
    OUT.mkdir(exist_ok=True)
    torch.set_num_threads(2)
    frozen_path = ROOT / 'planet_payload_freeze.json'
    frozen = json.loads(frozen_path.read_text())
    verified = []
    for entry in frozen['files']:
        path = ROOT / entry['path']
        actual = sha(path)
        assert path.stat().st_size == entry['bytes'] and actual == entry['sha256'], str(path)
        verified.append(dict(path=entry['path'], bytes=path.stat().st_size, sha256=actual))
    payload = ROOT / 'planet_payload.tar'
    assert sha(payload) == frozen['payload_sha256']
    cohort = ROOT / 'cohort'
    cohort_freeze = json.loads((cohort / 'cohort_freeze.json').read_text())
    for name, expected in cohort_freeze['files'].items():
        assert sha(cohort / name) == expected, name
    with (cohort / 'chain_manifest.csv').open(newline='') as stream:
        chains = list(csv.DictReader(stream))
    assert len(chains) == cohort_freeze['chains']
    for chain in chains:
        assert sha(Path(chain['pdb_file'])) == chain['input_sha256']
        assert len(chain['sequence']) == int(chain['n_residues'])
    sys.path.insert(0, str(ROOT / 'PLANET-MD'))
    checkpoints = {}
    for variant, filename in [('full', 'rocketshp_pretrained_20250505_v1.ckpt'),
                              ('seq', 'rocketshp_pretrained_20250507_v1_seq.ckpt')]:
        path = ROOT / 'weights' / filename
        checkpoint = torch.load(path, map_location='cpu', weights_only=True)
        hp = checkpoint['hyper_parameters']
        assert bool(hp['struct_features']) == (variant == 'full')
        state = checkpoint['state_dict']
        checkpoints[variant] = dict(path=str(path), sha256=sha(path),
                                   struct_features=bool(hp['struct_features']),
                                   structure_stage=hp.get('struct_stage', 'quantized'),
                                   hyper_parameters={k: str(v) for k, v in hp.items()},
                                   float_dtypes=sorted({str(v.dtype) for v in state.values() if torch.is_floating_point(v)}),
                                   tensor_count=len(state))
        default_temperature = next((float(value.item()) for key, value in state.items()
                                    if key.endswith('default_temp')), None)
        checkpoints[variant]['default_temp_from_state_dict'] = default_temperature
        assert default_temperature == 300.0
        del state, checkpoint
    adapter = ROOT / 'run_planet_md_3063.py'
    text = adapter.read_text()
    assert "get_model(device='cpu').to(device).eval()" in text
    assert 'rmsf_angstrom=10.0*native' in text
    assert not any(term in text for term in ('.half(', '.bfloat16(', 'autocast(', 'set_default_dtype('))
    trajectory = ROOT / 'PLANET-MD' / 'planet_md' / 'trajectory.py'
    trajectory_text = trajectory.read_text()
    assert 'return md.rmsf(traj, traj, 0, atom_indices=atom_indices)' in trajectory_text
    report = dict(status='all_frozen_planet_payload_and_cohort_files_passed_before_cpu_smoke',
                  created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  python=sys.executable, upstream_commit=frozen['upstream_planet_commit'],
                  esm_snapshot=frozen['esm_snapshot'], rocket_snapshot=frozen['rocket_snapshot'],
                  frozen_manifest_sha256=sha(frozen_path), payload_sha256=sha(payload),
                  verified_payload_files=verified, total_payload_bytes=sum(row['bytes'] for row in verified),
                  cohort_chains_verified=len(chains), cohort_records=cohort_freeze['records'],
                  smoke_inputs=chains[:2], checkpoints=checkpoints,
                  adapter_sha256=sha(adapter), no_dtype_conversion_in_adapter=True,
                  released_weights_float_dtypes=[checkpoints[v]['float_dtypes'] for v in ['full', 'seq']],
                  native_rmsf_unit='nm: frozen trajectory.py calls MDTraj md.rmsf directly; preprocessing and dataset labels perform no scale conversion',
                  conversion='rmsf_angstrom = 10.0 * rmsf_native_nm; no fitted calibration',
                  default_temperature='300 K: upstream PlanetMD forward uses default_temp when no temp key is provided',
                  hf_home=os.environ.get('HF_HOME'), hf_hub_offline=os.environ.get('HF_HUB_OFFLINE'),
                  cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'),
                  torch_num_threads=torch.get_num_threads(),
                  cpu_only=True, smoke_predictions_started=False,
                  runtime_versions={name:package_version(name) for name in ['torch','numpy','esm','mdtraj','biotite','huggingface-hub']})
    (OUT / 'payload_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key:report[key] for key in ['status','cohort_chains_verified','total_payload_bytes','adapter_sha256','native_rmsf_unit','released_weights_float_dtypes']}, indent=2))


if __name__ == '__main__':
    main()
