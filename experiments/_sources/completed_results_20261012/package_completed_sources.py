"""Read-only scientific assets; archive only 20 final weights and source tables."""
from pathlib import Path
import tarfile, hashlib, json
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'completed_results_20261012'
OUT.mkdir(exist_ok=True)
files = []
for folder in ['single_target', 'family_folds']:
    files.extend(sorted((ROOT/folder/'runs').glob('*/fold_*/epoch_040.pt')))
assert len(files) == 20
for rel in [
    'family_folds/evaluation/oof_per_residue.csv.gz',
    'family_folds/external_2430_evaluation/per_residue.csv.gz',
    'baselines/additional_original_oof/common_residue_predictions.csv.gz',
    'baselines/additional_family_oof/common_residue_predictions.csv.gz',
    'single_target/single_target_plan.json', 'single_target/single_target_queue_state.json',
    'single_target/analysis/analysis_audit.json', 'single_target/evaluation/evaluation_audit.json',
    'family_folds/training_freeze.json', 'family_folds/record_family_folds.csv',
    'family_folds/evaluation/evaluation_audit.json', 'family_folds/external_2430_evaluation/evaluation_audit.json',
    'pipeline_state.json', 'pipeline_control.json',
    'baselines/rmsf_net/run_audit.json', 'baselines/planet_full/run_audit.json', 'baselines/planet_seq/run_audit.json',
    'baselines/additional_original_oof/scoring_audit.json', 'baselines/additional_family_oof/scoring_audit.json',
]:
    files.append(ROOT/rel)
files.extend(sorted((ROOT/'single_target/configs').glob('*.yaml')))
files.extend(sorted((ROOT/'family_folds/configs').glob('*.yaml')))
assert all(p.is_file() and p.resolve().is_relative_to(ROOT.resolve()) for p in files)
def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()
manifest = {str(p.relative_to(ROOT)).replace('\\','/'): {'bytes':p.stat().st_size, 'sha256':sha(p)} for p in files}
package = OUT/'completed_scientific_sources_20261012.tar.gz'
with tarfile.open(package, 'w:gz', compresslevel=1) as archive:
    for p in files:
        archive.add(p, arcname=str(p.relative_to(ROOT)), recursive=False)
report = {'archive':str(package), 'archive_bytes':package.stat().st_size, 'archive_sha256':sha(package),
          'final_checkpoints':20, 'members':manifest, 'historical_checkpoints_excluded':True}
(OUT/'completed_scientific_sources_manifest.json').write_text(json.dumps(report, indent=2), encoding='utf8')
print(json.dumps({k:report[k] for k in ['archive','archive_bytes','archive_sha256','final_checkpoints']}))
