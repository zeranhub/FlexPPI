"""Descriptive partner-specificity audit; no raw between-state subtraction.

All states, both fixed initial contact surfaces, L12 and all eligible residues
are retained. Time windows and model checkpoints are NOT MD repeats.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / 'output/flexppi_kix_statewise_showcase_20260916'
QUALITY = ROOT / 'artifacts/kix_md_20260914/pilot_results_20260916/quality'
OUT = FIGURE / 'partner_specificity_audit'
STATES = ['apo', 'myb', 'mll', 'both']


def main():
    paths = dict(residues=FIGURE/'source_data/statewise_residues.csv',
                 folds=FIGURE/'source_data/all_five_checkpoints_and_ensemble.csv',
                 contacts=QUALITY/'native_contact_occupancy.csv',
                 windows=QUALITY/'residue_window_metrics.csv')
    source, folds, contacts, windows = [pd.read_csv(paths[k]) for k in paths]
    assert len(source) == 644 and len(folds) == 7728 and len(windows) == 5232
    assert not source.duplicated(['state','metric','author_residue_number']).any()
    assert set(source.state) == set(STATES)
    fixed = contacts[contacts.state.eq('both')]
    regions = {name: sorted(fixed.loc[fixed.partner_chain.eq(chain),
                                    'kix_author_residue'].unique().tolist())
               for name, chain in [('cMyb_contact_surface','A'),('MLL_contact_surface','C')]}
    regions['L12'] = list(range(614,622))
    assert [len(regions[x]) for x in regions] == [17,15,8]
    assert not set(regions['cMyb_contact_surface']) & set(regions['MLL_contact_surface'])
    paired, within, summaries = [], [], []
    for metric,n in [('rmsf',87),('hchi',74)]:
        data = source[source.metric.eq(metric)]
        md = data.pivot(index='author_residue_number',columns='state',values='md_percentile')[STATES]
        model = data.pivot(index='author_residue_number',columns='state',values='flexppi_percentile')[STATES]
        assert len(md)==len(model)==n and not md.isna().any().any()
        match = np.array([[spearmanr(md[a],model[b]).statistic for b in STATES] for a in STATES])
        for i,a in enumerate(STATES):
            for j,b in enumerate(STATES):
                paired.append(dict(metric=metric,md_state=a,model_state=b,n_residues=n,
                                   spearman=match[i,j],correct_state=a==b))
        for name,table in [('MD',md),('FlexPPI',model)]:
            for a,b in combinations(STATES,2):
                within.append(dict(metric=metric,source=name,state_a=a,state_b=b,
                                   spearman=spearmanr(table[a],table[b]).statistic))
        summaries.append(dict(metric=metric,correct_pair_mean=float(np.trace(match)/4),
                              mismatched_pair_mean=float(match[~np.eye(4,dtype=bool)].mean()),
                              best_match_correct_count=int((match.argmax(axis=1)==np.arange(4)).sum()),
                              best_matching_model_state={a:STATES[match[i].argmax()] for i,a in enumerate(STATES)}))

    region_rows = []
    for region,positions in regions.items():
        for (state,metric),group in source[source.author_residue_number.isin(positions)].groupby(['state','metric']):
            for name,col in [('MD','md_percentile'),('FlexPPI','flexppi_percentile')]:
                region_rows.append(dict(region=region,state=state,metric=metric,source=name,
                                        n_residues=len(group),mean_percentile=float(group[col].mean())))
    region_table = pd.DataFrame(region_rows)
    # Preserve all model checkpoints; no fold is selected using MD performance.
    fold_rows = []
    pct = folds[folds.normalization.eq('percentile')].copy()
    for region,positions in regions.items():
        for (state,metric,fold),group in pct[pct.author_residue_number.isin(positions)].groupby(['state','metric','model_fold']):
            fold_rows.append(dict(region=region,state=state,metric=metric,model_fold=fold,
                                  n_residues=len(group),mean_percentile=float(group.model_value.mean())))
    fold_table = pd.DataFrame(fold_rows)

    # Rank each whole-KIX window before computing any local region summary.
    kix = windows[windows.chain.eq('B')].copy()
    assert len(kix)==4176
    window_keys = ['state','window_kind','start_ns_exclusive','end_ns_inclusive']
    rank_rows, window_region_rows, half_agreement = [], [], []
    for keys,group in kix.groupby(window_keys):
        identifiers = dict(zip(window_keys,keys))
        for metric,col in [('rmsf','rmsf_angstrom'),('hchi','chi_entropy_nats')]:
            eligible = group if metric=='rmsf' else group[~group.amino_acid.isin(list('AGP'))]
            eligible = eligible.sort_values('author_residue_number').copy()
            assert len(eligible)==(87 if metric=='rmsf' else 74)
            assert np.isfinite(eligible[col]).all()
            eligible['percentile'] = 100*(rankdata(eligible[col])-1)/(len(eligible)-1)
            for row in eligible.itertuples():
                rank_rows.append(dict(**identifiers,metric=metric,author_residue_number=row.author_residue_number,
                                      amino_acid=row.amino_acid,percentile=row.percentile))
            if keys[1]=='cumulative' and keys[3]==100:
                expected = source[source.state.eq(keys[0]) & source.metric.eq(metric)].sort_values('author_residue_number')
                np.testing.assert_allclose(eligible.percentile,expected.md_percentile,atol=1e-10)
            for region,positions in regions.items():
                local=eligible[eligible.author_residue_number.isin(positions)]
                window_region_rows.append(dict(**identifiers,region=region,metric=metric,
                                                n_residues=len(local),mean_percentile=float(local.percentile.mean())))
    rank_table = pd.DataFrame(rank_rows)
    assert len(rank_table)==4*12*(87+74)
    window_regions = pd.DataFrame(window_region_rows)
    for (state,metric),group in rank_table[rank_table.window_kind.eq('disjoint_half')].groupby(['state','metric']):
        halves=group.pivot(index='author_residue_number',columns='start_ns_exclusive',values='percentile')
        half_agreement.append(dict(state=state,metric=metric,n_residues=len(halves),
                                   spearman=spearmanr(halves[0.],halves[50.]).statistic))
    pattern_rows=[]
    for metric in ['rmsf','hchi']:
        g=fold_table[(fold_table.metric.eq(metric)) & (~fold_table.model_fold.astype(str).eq('ensemble_mean'))]
        for fold in sorted(g.model_fold.astype(str).unique()):
            cg=g[g.model_fold.astype(str).eq(fold) & g.region.eq('cMyb_contact_surface')].set_index('state').mean_percentile
            mg=g[g.model_fold.astype(str).eq(fold) & g.region.eq('MLL_contact_surface')].set_index('state').mean_percentile
            pattern_rows.append(dict(metric=metric,model_fold=fold,
                                     cmyb_surface_complete_order=bool(cg['myb']<cg['both']<cg['apo']<cg['mll']),
                                     mll_surface_cmyb_uniquely_highest=bool(mg['myb']>mg.drop('myb').max())))

    OUT.mkdir(exist_ok=True)
    exports={'cross_state_matching':pd.DataFrame(paired),'within_source_state_similarity':pd.DataFrame(within),
             'fixed_region_percentiles':region_table,'checkpoint_region_percentiles':fold_table,
             'md_window_region_percentiles':window_regions,'md_window_residue_percentiles':rank_table,
             'md_half_profile_agreement':pd.DataFrame(half_agreement),
             'checkpoint_pattern_support':pd.DataFrame(pattern_rows)}
    for name,table in exports.items():
        table.to_csv(OUT/(name+'.csv'),index=False)
    manifest=dict(scope='Exploratory partner-dependent within-state ranks; no raw state subtraction',
                  states=STATES,regions=regions,
                  region_definition='Both-state initial <=4.5 Angstrom heavy-atom contact residue sets from existing quality/native_contact_occupancy.csv; all states use the same sets; no contact-occupancy filtering',
                  local_selection='Both contact surfaces retained regardless of outcome; L12 retained as prior region',
                  hchi_exclusion='A/G/P only, following existing evaluation; cMyb/MLL/L12 eligible n=16/14/4',
                  n_md_trajectories_per_state=1,production_ns_per_state=100,
                  no_p_values=True,no_independent_repeat_intervals=True,
                  caveat='Within-state percentile is relative placement, not absolute flexibility change or proof of allosteric causality. Time windows are not independent repeats. Near-identical KIX exists in training.',
                  global_matching=summaries,
                  input_sha256={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in paths.items()},
                  output_rows={k:len(v) for k,v in exports.items()},
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'audit_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(region_table.to_string(index=False))
    print(json.dumps(summaries,indent=2))
    print(pd.DataFrame(pattern_rows).to_string(index=False))
    print(OUT)


if __name__=='__main__':
    main()
