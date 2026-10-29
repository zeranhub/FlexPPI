"""Audit and compact the frozen Figure 5 readouts; no training, MD or inference.

The newly calculated confidence intervals resample whole records within each
source, using identical paired resamples for every method and readout.  These
are descriptive record-level intervals, not family-level intervals.
"""
from pathlib import Path
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
BASE = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OLD = OUT.parent.parent / 'revision_v80_redesign' / 'figure5'
ENT = Path(r'E:\Project\FlexPPI\experiments\_sources\model_argument_v24_20261013\experiment_predictive_occupancy')
WID = Path(r'E:\Project\FlexPPI\experiments\_sources\free_energy_model_v26_20261013\experiment')
SEED = 20261026
B = 10000
PREDICTORS = ['full_hchi', 'intrachain_hchi', 'full_rmsf']
LABELS = dict(full_hchi='Full Hχ', intrachain_hchi='Intra Hχ', full_rmsf='Full RMSF')
read = lambda p: pd.read_csv(p, float_precision='round_trip')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

def filt(frame, entropy=False):
    frame = frame[frame.scope.eq('all')].copy()
    if 'threshold' in frame:
        frame = frame[np.isclose(frame.threshold, 1)]
    if entropy:
        frame = frame[frame.metric.eq('fine_entropy_mean_per_chi')]
    return frame

OUT.mkdir(parents=True, exist_ok=True)
inputs = []
corr_parts, sel_parts, old_contrasts, old_selections = [], [], [], []
for key, folder, original, is_entropy, unit in [
    ('fine_entropy', OLD/'source/entropy', ENT, True, 'nats/χ'),
    ('lowF_width', OLD/'source/width', WID, False, 'degrees/χ'),
]:
    path = folder/'replica1_per_record_correlations.csv'
    inputs.append(path)
    corr = filt(read(path), is_entropy)
    check = filt(read(original/path.name), is_entropy)
    join = corr.merge(check, on=['sample_id','predictor'], suffixes=('_fig5','_source'), validate='one_to_one')
    np.testing.assert_allclose(join.rho_type_centered_fig5, join.rho_type_centered_source, rtol=0, atol=1e-14)
    assert len(corr) == 249 and corr.sample_id.nunique() == 83
    assert corr.groupby('predictor').positions.sum().to_dict() == {p:33142 for p in PREDICTORS}
    corr = corr[['sample_id','dataset','fold','predictor','positions','AA_types','rho_type_centered']]
    corr['metric_key'] = key
    corr['method_label'] = corr.predictor.map(LABELS)
    corr_parts.append(corr)

    path = folder/'replica1_high_low_per_record.csv'
    inputs.append(path)
    selected = filt(read(path), is_entropy)
    # The old entropy record table omitted absolute means; reconstruct them
    # from all rankable, equally weighted AA groups, without revising halves.
    aa_path = original/'replica1_high_low_per_AA.csv'
    inputs.append(aa_path)
    groups = filt(read(aa_path), is_entropy)
    groups = groups[groups.rankable].copy()
    abs_means = groups.groupby(['sample_id','predictor'], as_index=False).agg(
        low_MD_mean=('low_MD_mean','mean'), high_MD_mean=('high_MD_mean','mean'))
    if 'low_MD_mean' in selected:
        joined = selected.merge(abs_means, on=['sample_id','predictor'], suffixes=('_frozen','_AA'), validate='one_to_one')
        for col in ['low_MD_mean','high_MD_mean']:
            np.testing.assert_allclose(joined[col+'_frozen'], joined[col+'_AA'], rtol=0, atol=1e-12)
        selected = selected.drop(columns=['low_MD_mean','high_MD_mean'])
    selected = selected.merge(abs_means, on=['sample_id','predictor'], validate='one_to_one')
    np.testing.assert_allclose(selected.high_MD_mean-selected.low_MD_mean, selected.high_minus_low_MD, rtol=0, atol=1e-12)
    selected = selected[['sample_id','dataset','fold','predictor','AA_groups','positions','low_positions','high_positions','omitted_middle_positions','low_MD_mean','high_MD_mean','high_minus_low_MD']]
    selected['metric_key'] = key
    selected['unit'] = unit
    selected['method_label'] = selected.predictor.map(LABELS)
    assert len(selected) == 249
    sel_parts.append(selected)

    path = folder/'replica1_paired_correlation_difference_summary.csv'
    inputs.append(path)
    contrasts = filt(read(path), is_entropy)
    contrasts = contrasts[contrasts.dataset.eq('all')].copy()
    contrasts['metric_key'] = key
    old_contrasts.append(contrasts)
    path = folder/'replica1_high_low_summary.csv'
    inputs.append(path)
    ss = filt(read(path), is_entropy)
    ss = ss[ss.dataset.eq('all')].copy()
    ss['metric_key'] = key
    old_selections.append(ss)

correlations = pd.concat(corr_parts, ignore_index=True)
selections = pd.concat(sel_parts, ignore_index=True)
pd.concat(old_contrasts).to_csv(OUT/'frozen_paired_correlation_statistics.csv', index=False, float_format='%.17g')
pd.concat(old_selections).to_csv(OUT/'frozen_selection_statistics.csv', index=False, float_format='%.17g')
correlations.to_csv(OUT/'retained_absolute_correlations_all498.csv', index=False, float_format='%.17g')
selections.to_csv(OUT/'retained_prediction_groups_all498.csv', index=False, float_format='%.17g')

meta = correlations[['sample_id','dataset']].drop_duplicates().sort_values('sample_id').reset_index(drop=True)
assert meta.dataset.value_counts().to_dict() == {'atlas':30,'itrace':30,'dynarepo':23}
record_ids = meta.sample_id.to_list()
rng = np.random.default_rng(SEED)
parts = []
for source in ['atlas','dynarepo','itrace']:
    loc = np.flatnonzero(meta.dataset.to_numpy() == source)
    parts.append(loc[rng.integers(0,len(loc),size=(B,len(loc)))])
indices = np.concatenate(parts, axis=1)
np.savez_compressed(OUT/'paired_source_stratified_bootstrap_indices.npz',
    record_ids=meta.sample_id.to_numpy(str), record_source=meta.dataset.to_numpy(str),
    record_indices=indices, seed=SEED, resamples=B)

summary = []
comparisons = []
for kind, table, values in [
    ('correlation', correlations, ['rho_type_centered']),
    ('prediction_selection', selections, ['low_MD_mean','high_MD_mean','high_minus_low_MD']),
]:
    for metric, sub in table.groupby('metric_key'):
        arrays = {}
        for pred in PREDICTORS:
            group = sub[sub.predictor.eq(pred)].set_index('sample_id').loc[record_ids]
            for col in values:
                arr = group[col].to_numpy(float)
                assert np.isfinite(arr).all() and len(arr) == 83
                boot = arr[indices].mean(axis=1)
                lo, hi = np.quantile(boot,[.025,.975])
                summary.append(dict(panel_kind=kind, metric_key=metric, predictor=pred,
                    method_label=LABELS[pred], value_key=col, n_records=83,
                    mean=arr.mean(), bootstrap_95_low=lo, bootstrap_95_high=hi,
                    median=np.median(arr), min=arr.min(), max=arr.max(),
                    positive_records=int((arr>0).sum()),negative_records=int((arr<0).sum())))
                arrays[pred,col] = arr
        value = 'rho_type_centered' if kind == 'correlation' else 'high_minus_low_MD'
        for comparator in ['intrachain_hchi','full_rmsf']:
            arr = arrays['full_hchi',value] - arrays[comparator,value]
            boot = arr[indices].mean(axis=1)
            lo, hi = np.quantile(boot,[.025,.975])
            comparisons.append(dict(panel_kind=kind,metric_key=metric,value_key=value,
                contrast='full_hchi_minus_'+comparator,n_records=83,mean=arr.mean(),
                bootstrap_95_low=lo,bootstrap_95_high=hi,
                positive_records=int((arr>0).sum()),negative_records=int((arr<0).sum())))

summary = pd.DataFrame(summary)
comparisons = pd.DataFrame(comparisons)
summary.to_csv(OUT/'retained_panel_statistics.csv', index=False, float_format='%.17g')
comparisons.to_csv(OUT/'retained_paired_comparison_statistics.csv', index=False, float_format='%.17g')

# Save relevant existing prose for exact section-position inspection.
docx = BASE/'paper_v90_20261026_中文修订_图3数据呈现与机制推断.docx'
ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with zipfile.ZipFile(docx) as z:
    root = ET.fromstring(z.read('word/document.xml'))
paras = [''.join(p.itertext()) for p in []] # text extraction below uses w:t only
paras = [''.join(t.text or '' for t in p.findall('.//w:t',ns)) for p in root.findall('.//w:body/w:p',ns)]
selected = []
active = False
for i, text in enumerate(paras):
    if text.startswith('3.5'):
        active = True
    if active:
        selected.append(dict(paragraph_index=i,text=text))
    if active and text.startswith('4 '):
        break
for i, text in enumerate(paras):
    if text.startswith('2.7'):
        selected += [dict(paragraph_index=j,text=paras[j]) for j in range(i,min(i+12,len(paras)))]
        break
(OUT/'existing_sections_2_7_and_3_5.json').write_text(json.dumps(selected,ensure_ascii=False,indent=2),encoding='utf8')

audit = dict(status='complete',no_model_retraining=True,no_new_MD=True,no_inference=True,
    frozen_source_records=83,source_counts=meta.dataset.value_counts().to_dict(),eligible_positions=33142,
    retained_correlations=len(correlations),retained_selection_rows=len(selections),
    all_three_predictors_and_both_readouts_present=True,source_absolute_rho_values_exact=True,
    selection_absolute_MD_means_reconstructed_from_original_AA_groups=True,
    bootstrap=dict(seed=SEED,resamples=B,unit='whole record within source; same paired indices for all outputs',
        family_cluster_resampling=False),
    unit_boundary='Fine-bin entropy nats/chi; low-Fstar width degrees/chi; Fstar dimensionless log-count ratio, no temperature or kcal/mol conversion',
    model_boundary='Original pipeline frozen OOF, not the new family-separated analysis; descriptors from the same existing MD ensemble, not experimental thermodynamics or causal validation',
    source_sha256={str(p):sha(p) for p in inputs},
    output_sha256={p.name:sha(p) for p in OUT.glob('*.csv')})
(OUT/'entropy_panel_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(summary.to_string(index=False))
print(comparisons.to_string(index=False))
