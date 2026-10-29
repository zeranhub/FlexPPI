"""CPU-only environment audit; imports and platform enumeration create no Context."""
from pathlib import Path
import datetime
import hashlib
import importlib.metadata as metadata
import json
import platform
import sys
import tarfile

TASK = Path('/root/flexppi_completion_20261010')
BATCH = TASK / 'kix_md_triplicate_rebuild_20261010'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    import numpy
    import scipy
    import MDAnalysis
    import pdbfixer
    import openmm
    from openmm import app

    expected_versions = {'numpy': '2.2.6', 'scipy': '1.15.3',
                         'MDAnalysis': '2.9.0', 'pdbfixer': '1.12.0', 'OpenMM': '8.4.0'}
    actual_versions = {name: metadata.version(name) for name in expected_versions}
    versions_match = actual_versions == expected_versions
    historical_path = BATCH / 'historical_audit' / 'prepared_READY.json'
    historical = json.loads(historical_path.read_text())
    expected_ff = historical['forcefield_sha256']
    data_folder = Path(app.__file__).parent / 'data'
    actual_ff = {name: sha(data_folder / name) for name in expected_ff}
    ff_match = actual_ff == expected_ff
    comparisons = {name: {'historical': expected_ff[name], 'actual': actual_ff[name],
                          'match': expected_ff[name] == actual_ff[name]}
                   for name in expected_ff}
    platforms = [openmm.Platform.getPlatform(i).getName()
                 for i in range(openmm.Platform.getNumPlatforms())]
    plugin_failures = list(openmm.Platform.getPluginLoadFailures())
    conda_file = TASK / 'md_conda_pkgs' / 'openmm-8.4.0-py311h1c77f12_0.conda'
    official_conda_sha = '8542b87f2b86a567322db41a2f28616f9a57227b0e7a8f98326a8f6cb9417b7f'
    conda_sha = sha(conda_file)
    source_file = TASK / 'md_install_sources' / 'pdbfixer_94cfa4c0ca551cdc5f13320f9a658efd59f2b881.tar.gz'
    with tarfile.open(source_file) as archive:
        setup_entry = next(member for member in archive.getmembers() if member.name.endswith('/setup.py'))
        setup_text = archive.extractfile(setup_entry).read().decode()
    pdbfixer_source_version_verified = '1.12.0' in setup_text
    report = dict(
        status='environment_ready_cpu_verified' if versions_match and ff_match and
               conda_sha == official_conda_sha and 'CUDA' in platforms and
               pdbfixer_source_version_verified else 'environment_audit_failed',
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        python=sys.executable, python_version=sys.version, libc=platform.libc_ver(),
        gpu_context_created=False, simulation_run=False,
        gpu_smoke_required_before_production=True,
        expected_versions=expected_versions, actual_versions=actual_versions,
        versions_match=versions_match,
        historical_openmm='8.4.0.post2 (PyPI; Linux wheel requires glibc >= 2.34)',
        current_openmm='8.4.0 py311h1c77f12_0 (conda-forge; glibc >= 2.17)',
        version_difference_explicit=True,
        openmm_runtime_full_version=openmm.version.full_version,
        platforms=platforms, plugin_load_failures=plugin_failures,
        forcefield_data_directory=str(data_folder),
        historical_forcefield_manifest=str(historical_path),
        historical_forcefield_manifest_sha256=sha(historical_path),
        forcefield_sha256_match=ff_match, forcefield_files=comparisons,
        openmm_conda_archive=dict(path=str(conda_file), sha256=conda_sha,
                                 official_sha256=official_conda_sha,
                                 official_source='https://api.anaconda.org/release/conda-forge/openmm/8.4.0'),
        pdbfixer_source=dict(commit='94cfa4c0ca551cdc5f13320f9a658efd59f2b881',
                             official_tag='v1.12', path=str(source_file), sha256=sha(source_file),
                             version_verified=pdbfixer_source_version_verified),
        driver_version_text=Path('/proc/driver/nvidia/version').read_text(),
        caveat='CUDA 12.6 NVRTC may require newer PTX support than driver 550.78; plugin enumeration is not a GPU smoke test.'
    )
    report_path = TASK / 'md_environment_audit.json'
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    if report['status'] != 'environment_ready_cpu_verified':
        raise RuntimeError(json.dumps(report, indent=2))
    (TASK / 'MD_ENV_READY').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in
                      ['status', 'python', 'actual_versions', 'platforms',
                       'plugin_load_failures', 'forcefield_sha256_match', 'gpu_context_created']}, indent=2))


if __name__ == '__main__':
    main()
