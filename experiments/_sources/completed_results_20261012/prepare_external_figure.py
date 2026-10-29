from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
source = ROOT.parent / 'figures/redraw_quantitative.py'
out = Path(__file__).resolve().parent / 'figures'
out.mkdir(exist_ok=True)
text = source.read_text(encoding='utf8')
text = text.replace("EXT = BASE / 'figures/corrected_external_results_preparation/completed/20261003T163520Z_03a3b308'", "EXT = Path(r'" + str(ROOT / 'family_folds/external_2430_evaluation') + "')")
text = text.replace('2,431条', '2,430条')
text = text.replace("str(EXT/'all_32_per_record_metrics.csv')", "str(EXT/'per_record_metrics.csv')")
text = text[:text.index('def main():')] + '''def main():
    path = EXT / 'per_record_metrics.csv'
    audit_path = EXT / 'evaluation_audit.json'
    audit = json.loads(audit_path.read_text(encoding='utf8'))
    assert audit['status'] == 'complete'
    d = read(path)
    assert len(d) == 32 * 8
    assert d.groupby(['variant', 'metric']).sample_id.nunique().eq(32).all()
    assert np.isfinite(d[['mae', 'rmse', 'bias', 'spearman', 'pearson']]).all().all()
    d['protocol'] = 'strict_single_seed7'
    d['family_components'] = 1
    fig4(d)
    MANIFEST['training_records'] = 2430
    MANIFEST['excluded_records'] = 633
    MANIFEST['new_training'] = 0
    MANIFEST['new_inference'] = 0
    MANIFEST['visual_review'] = 'pending'
    MANIFEST['audit_sha256'] = sha(audit_path)
    (OUT/'external_2430_figure_manifest.json').write_text(json.dumps(MANIFEST, ensure_ascii=False, indent=2), encoding='utf8')
    print(json.dumps({'figures': list(MANIFEST['figures']), 'records':32, 'training_records':2430}, ensure_ascii=False))

if __name__ == '__main__':
    main()
'''
(out / 'redraw_external_2430.py').write_text(text, encoding='utf8')
print(out / 'redraw_external_2430.py')
