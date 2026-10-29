"""Freeze v110 inputs; validate exact values. Does not plan or render a layout."""
from pathlib import Path
import hashlib
import json
import shutil
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'revision_v110_compact_five_panels'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare():
    manifest = []
    directories = ['figure3/source_data', 'pymol', 'source_inputs', 'sasa_audit']
    filenames = ['case_profiles.py', 'plot_overlay.py', 'case_profile_manifest.json',
                 'case_profiles_all97_with_masks.csv', 'case_SASA_all97_starting_complex.csv',
                 'standalone_all368_plotted_values.csv', 'same_residue_MD_and_Full_all824988.npz',
                 'direct_relationship_statistics.json', 'cohort_all3063_records.csv',
                 'figure3/publication_style.py', 'figure3/main_original_panels_a_to_e_values.csv']
    inputs = [SOURCE / name for name in filenames]
    for directory in directories:
        inputs.extend(p for p in (SOURCE / directory).rglob('*') if p.is_file())
    for source in sorted(set(inputs)):
        relative = source.relative_to(SOURCE)
        target = ROOT / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        assert digest(source) == digest(target)
        manifest.append({'relative_path': relative.as_posix(), 'sha256': digest(source)})

    # Change only dependency paths in copied helpers, never the plotting/data logic.
    case_path = ROOT / 'case_profiles.py'
    text = case_path.read_text(encoding='utf-8')
    text = text.replace("ub_source=ROOT.parent/'revision_v93_dual_case_entropy/ubiquitin_audit/all119_input_environment.csv'",
                        "ub_source=INPUT/'all119_input_environment.csv'")
    text = text.replace("shutil.copy2(ub_source,INPUT/ub_source.name)",
                        "assert ub_source.exists(), 'Frozen solvent-exposure input is missing'")
    case_path.write_text(text, encoding='utf-8')

    overlay_path = ROOT / 'plot_overlay.py'
    text = overlay_path.read_text(encoding='utf-8')
    text = text.replace("if not (ROOT/name).exists():shutil.copy2(OLD/name,ROOT/name)\n        assert sha256((ROOT/name).read_bytes()).digest()==sha256((OLD/name).read_bytes()).digest()",
                        "assert (ROOT/name).exists(), 'Run prepare_frozen_inputs.py first'\n        frozen=json.loads((ROOT/'frozen_input_manifest.json').read_text(encoding='utf-8'))\n        expected={r['relative_path']:r['sha256'] for r in frozen['copied_inputs']}\n        assert sha256((ROOT/name).read_bytes()).hexdigest()==expected[name]")
    overlay_path.write_text(text, encoding='utf-8')

    profiles = pd.read_csv(ROOT / 'case_profiles_all97_with_masks.csv', float_precision='round_trip')
    curves = pd.read_csv(ROOT / 'standalone_all368_plotted_values.csv', float_precision='round_trip')
    sasa = pd.read_csv(ROOT / 'case_SASA_all97_starting_complex.csv', float_precision='round_trip')
    old = pd.read_csv(ROOT / 'figure3/main_original_panels_a_to_e_values.csv', float_precision='round_trip')
    assert len(profiles) == 97 and len(curves) == 368 and len(sasa) == 97 and len(old) == 767
    assert old.groupby('panel').size().to_dict() == {'b': 117, 'c': 348, 'e': 302}
    p = profiles[profiles.case.eq('Ubiquitin') & profiles.amino_acid.eq('P')]
    assert p.residue_id.tolist() == [19, 37, 38] and p.valid_hchi.all() and p.MD_Hchi_nats.eq(0).all()
    assert profiles[profiles.case.eq('CaMKI')].valid_hchi.sum() == 21
    assert profiles[profiles.case.eq('Ubiquitin')].valid_hchi.sum() == 66

    # This native row ledger is the exact old b/c/e plotted ledger, not regenerated values.
    old.to_csv(ROOT / 'frozen_all767_plotted_values.csv', index=False)
    with np.load(ROOT / 'same_residue_MD_and_Full_all824988.npz', allow_pickle=False) as z:
        for key in ['MD_RMSF', 'MD_Hchi', 'Full_RMSF', 'Full_Hchi']:
            assert len(z[key]) == 824988 and np.isfinite(z[key]).all()
    report = {'source': str(SOURCE), 'copied_inputs': manifest,
              'native_curves': 368, 'static_SASA': 97, 'chi_scatter_coordinates': 302,
              'total_plotted_scalar_values': 767, 'common_density_positions': 824988,
              'case_positions': 97, 'invalid_Hchi_gaps_retained': True, 'legal_Pro_zeros_retained': True,
              'layout_generated': False, 'paperbanana_layout_available': False,
              'helper_changes': ['case_profiles.load_sasa: use frozen local Ubiquitin SASA',
                                 'plot_overlay.load_data: use frozen local SHA256 manifest'],
              'local_helper_sha256': {p.name: digest(p) for p in [case_path, overlay_path]}}
    (ROOT / 'frozen_input_manifest.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    return report


if __name__ == '__main__':
    report = prepare()
    print(json.dumps({k: v for k, v in report.items() if k != 'copied_inputs'}, ensure_ascii=False, indent=2))
