"""Scale-free sensitivity analysis requested after the frozen pilot readout.

No fitting to MD, no model/checkpoint selection, and no raw-data modification.
State-specific normalization changes the estimand from motion amplitude to
within-state relative flexibility; it cannot establish absolute rigidification.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'artifacts/kix_md_20260914/pilot_results_20260916'

def transform(values, mode):
    values = np.asarray(values, float)
    assert np.isfinite(values).all() and len(values) > 1
    if mode == 'raw':
        return values
    if mode == 'percentile':
        return 100 * (rankdata(values, method='average') - 1) / (len(values) - 1)
    if mode == 'standardized':
        return (values - values.mean()) / values.std(ddof=0)
    if mode == 'relative_to_mean':
        return values / values.mean()
    raise ValueError(mode)

def summarize(x, y):
    x, y = np.asarray(x), np.asarray(y)
    variable = np.ptp(x) > 1e-12 and np.ptp(y) > 1e-12
    return dict(n=len(x), pearson=float(np.corrcoef(x, y)[0, 1]) if variable else np.nan,
                spearman=float(spearmanr(x, y)[0]) if variable else np.nan, mae=float(np.abs(x-y).mean()),
                md_mean=float(x.mean()), model_mean=float(y.mean()))

def main():
    source = pd.read_csv(BASE / 'comparison/matched_residues.csv', dtype={'model_fold': str})
    rows = []
    for (state, metric, fold), group in source.groupby(['state', 'metric', 'model_fold']):
        group = group.sort_values('author_residue_number')
        assert len(group) == (87 if metric == 'rmsf' else 74)
        for mode in ['raw', 'percentile', 'standardized', 'relative_to_mean']:
            item = group.copy()
            item['normalization'] = mode
            item['md_value'] = transform(group.md_value, mode)
            item['model_value'] = transform(group.model_value, mode)
            rows.append(item)
    transformed = pd.concat(rows, ignore_index=True)
    transformed.to_csv(BASE / 'relative_residues.csv', index=False)
    regional, metrics, deltas, top_sets = [], [], [], []
    for keys, group in transformed.groupby(['state', 'metric', 'model_fold', 'normalization']):
        state, metric, fold, mode = keys
        header = dict(state=state, metric=metric, model_fold=fold, normalization=mode)
        metrics.append(dict(**header, **summarize(group.md_value, group.model_value)))
        region = group[group.author_residue_number.between(614, 621)]
        regional.append(dict(**header, n=len(region), md_mean=region.md_value.mean(), model_mean=region.model_value.mean()))
        if mode == 'percentile':
            a = set(group.loc[group.md_value >= 80, 'author_residue_number'])
            b = set(group.loc[group.model_value >= 80, 'author_residue_number'])
            top_sets.append(dict(**header, md_top_count=len(a), model_top_count=len(b),
                                 overlap=len(a & b), union=len(a | b), jaccard=len(a & b)/len(a | b),
                                 md_residues=';'.join(map(str, sorted(a))), model_residues=';'.join(map(str, sorted(b)))))
    delta_stats = []
    for before, after in [('apo', 'myb'), ('apo', 'mll'), ('apo', 'both'), ('myb', 'mll')]:
        a = transformed[transformed.state == after]
        b = transformed[transformed.state == before]
        keys = ['metric', 'model_fold', 'normalization', 'author_residue_number', 'amino_acid']
        pair = a.merge(b, on=keys, suffixes=('_after', '_before'), validate='one_to_one')
        table = pair[keys].copy()
        table['contrast'] = after + '_minus_' + before
        table['md_delta'] = pair.md_value_after - pair.md_value_before
        table['model_delta'] = pair.model_value_after - pair.model_value_before
        deltas.append(table)
        for (metric, fold, mode), group in table.groupby(['metric', 'model_fold', 'normalization']):
            for region, subset in [('KIX', group), ('L12', group[group.author_residue_number.between(614, 621)])]:
                delta_stats.append(dict(contrast=after+'_minus_'+before, metric=metric, model_fold=fold,
                                       normalization=mode, region=region,
                                       zero_response_mae=float(np.abs(subset.md_delta).mean()),
                                       mae_improvement_over_zero=float(np.abs(subset.md_delta).mean()-np.abs(subset.md_delta-subset.model_delta).mean()),
                                       sign_agreement_fraction=float((np.sign(subset.md_delta) == np.sign(subset.model_delta)).mean()),
                                       **summarize(subset.md_delta, subset.model_delta)))
    for name, data in [('relative_regions', regional), ('relative_metrics', metrics), ('top_quintile_overlap', top_sets), ('relative_delta_metrics', delta_stats)]:
        pd.DataFrame(data).to_csv(BASE / (name+'.csv'), index=False)
    pd.concat(deltas, ignore_index=True).to_csv(BASE / 'relative_deltas.csv', index=False)
    report = dict(status='complete', requested_after_MD_readout=True, not_preregistered_primary_endpoint=True,
                  normalization_reference='Same KIX residues within each state separately for MD and each model/checkpoint; RMSF87, Hchi74 excluding AGP',
                  percentile='100*(average rank-1)/(n-1): 0=least flexible;100=most flexible within the state',
                  standardized='(value-state mean)/state population SD', relative_to_mean='value/state arithmetic mean',
                  ensemble_rule='Average original five checkpoint predictions first, then normalize. Individual checkpoints retained separately.',
                  interpretations=['High per-state rank agreement supports local flexibility mapping, not necessarily response prediction.',
                                   'Rank shifts refer to relative position among residues, not Angstrom/nat amplitude changes.',
                                   'State-wise normalization removes global shifts and can obscure real global changes.',
                                   'No MD-fitted calibration, no optimized threshold, no selected model or state.',
                                   'One MD repeat and near-identical KIX training sequence limit generalization.'],
                  source_sha256=hashlib.sha256((BASE/'comparison/matched_residues.csv').read_bytes()).hexdigest())
    (BASE/'relative_analysis_contract.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(pd.DataFrame(regional).query("model_fold == 'ensemble_mean' and metric == 'rmsf' and normalization in ['raw','percentile']").to_string(index=False))
    print(pd.DataFrame(delta_stats).query("model_fold == 'ensemble_mean' and region == 'KIX' and normalization == 'percentile'").to_string(index=False))

if __name__ == '__main__':
    main()
