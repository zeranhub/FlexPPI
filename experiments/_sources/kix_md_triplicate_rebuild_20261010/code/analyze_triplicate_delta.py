"""Frozen, exhaustive triplicate KIX analysis; no fitting, selection, CI or p.

Requires all 12 new-batch repeats complete. Historical pilot is never pooled.
Predictions CSV: state,method,chain,author_residue_number,amino_acid,
predicted_rmsf_angstrom,predicted_hchi_nats. Methods full/static mandatory;
intra/aa optional. SHA and original input hashes must be in FROZEN.json.
Only numpy is needed in addition to Python's standard library.
"""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import numpy as np

STATES = ('apo', 'myb', 'mll', 'both')
REPEATS = ('rep01', 'rep02', 'rep03')
METRICS = {'rmsf': ('rmsf_angstrom', 'predicted_rmsf_angstrom'),
           'hchi': ('chi_entropy_nats', 'predicted_hchi_nats')}
REGIONS = ('KIX', 'L12')
NORMALIZATIONS = ('raw', 'state_percentile')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        raise ValueError('Empty output: ' + str(path))
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def ranks(x):
    x = np.asarray(x, dtype=float)
    _, inverse, counts = np.unique(x, return_inverse=True, return_counts=True)
    ends = np.cumsum(counts)
    starts = ends - counts + 1
    return ((starts + ends) / 2)[inverse]


def percentile(x):
    return 100 * (ranks(x) - 1) / (len(x) - 1)


def corr(x, y):
    if np.ptp(x) < 1e-12 or np.ptp(y) < 1e-12:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def score(reference, prediction):
    x, y = np.asarray(reference), np.asarray(prediction)
    if x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Unmatched or non-finite scoring data')
    return dict(n_residues=len(x), mae=float(np.mean(np.abs(y-x))),
                rmse=float(np.sqrt(np.mean((y-x)**2))), bias=float(np.mean(y-x)),
                mean_md=float(x.mean()), mean_model=float(y.mean()),
                pearson=corr(x,y), spearman=corr(ranks(x),ranks(y)))


def load_completed(root):
    md, sources = {}, []
    all_rows, identities = [], None
    for state, repeat in itertools.product(STATES, REPEATS):
        job = root / 'runs' / state / repeat
        done = json.loads((job / 'DONE.json').read_text())
        summary_path = job / 'analysis' / 'summary.json'
        summary = json.loads(summary_path.read_text())
        metrics_path = job / 'analysis' / 'residue_metrics.csv'
        if (done['state'] != state or done['repeat'] != int(repeat[3:]) or
            done['production_ns'] != 100 or done['equilibration_ns'] != 10 or
            done['final_step'] != 55000000 or summary['n_frames'] != 100000 or
            summary['replicate'] != repeat or summary['state'] != state or
            len(summary['segments']) != 100 or summary['KIX_CA_count'] != 87):
            raise ValueError('Incomplete or mismatched new-batch repeat: ' + str(job))
        if sha(metrics_path) != summary['output_csv_sha256']:
            raise ValueError('Per-repeat CSV checksum mismatch: ' + str(job))
        if sha(root / 'prepared' / state / 'protein.pdb') != summary['protein_pdb_sha256']:
            raise ValueError('Prepared topology mismatch: ' + str(job))
        equilibration = sorted((job / 'equilibration').glob('segment_*.json'))
        if [p.name for p in equilibration] != [f'segment_{i:04d}.json' for i in range(10)]:
            raise ValueError('Equilibration does not contain all ten committed segments')
        for marker in equilibration:
            commit = json.loads(marker.read_text())
            checkpoint = marker.parent / commit['checkpoint']
            if sha(checkpoint) != commit['checkpoint_sha256']:
                raise ValueError('Equilibration checkpoint checksum mismatch')
        names = [s['filename'] for s in summary['segments']]
        if names != [f'segment_{i:04d}.dcd' for i in range(100)]:
            raise ValueError('Production segments are not exactly 0000..0099')
        for segment in summary['segments']:
            trajectory = job / 'production' / segment['filename']
            marker = job / 'production' / (trajectory.stem + '.json')
            commit = json.loads(marker.read_text())
            if trajectory.stat().st_size != segment['bytes'] or commit['expected_frames'] != 1000:
                raise ValueError('Invalid production segment: ' + str(trajectory))
            if sha(trajectory) != segment['sha256']:
                raise ValueError('Production DCD checksum changed after per-repeat analysis')
            checkpoint = marker.parent / commit['checkpoint']
            if sha(checkpoint) != commit['checkpoint_sha256']:
                raise ValueError('Production checkpoint checksum mismatch')
        rows = sorted([r for r in read_csv(metrics_path) if r['chain'] == 'B'],
                      key=lambda r: int(r['author_residue_number']))
        keys = [(int(r['author_residue_number']), r['amino_acid'], int(r['n_chi_angles'])) for r in rows]
        if [k[0] for k in keys] != list(range(586, 673)):
            raise ValueError('KIX author residue numbering mismatch')
        if identities is None:
            identities = keys
        if identities != keys:
            raise ValueError('Residue identity/mask differs across states or repeats')
        if any(r['state'] != state or r['replicate'] != repeat for r in rows):
            raise ValueError('MD residue rows identify a different job')
        for metric, (column, _) in METRICS.items():
            selected = [r for r in rows if metric == 'rmsf' or int(r['n_chi_angles']) > 0]
            md[state, repeat, metric] = np.array([float(r[column]) for r in selected])
            if not np.isfinite(md[state, repeat, metric]).all():
                raise ValueError('Non-finite MD readout')
        all_rows += rows
        sources.append(dict(state=state,repeat=repeat,path=str(job),metrics_sha256=sha(metrics_path),
                            summary_sha256=sha(summary_path),frames=100000))
    return md, identities, sources, all_rows


def load_frozen(root, identities):
    folder = root / 'frozen_predictions'
    frozen = json.loads((folder / 'FROZEN.json').read_text())
    path = folder / 'residue_predictions.csv'
    if sha(path) != frozen['prediction_csv_sha256']:
        raise ValueError('Frozen prediction CSV changed')
    for state in STATES:
        if sha(root / 'inputs' / f'{state}.pdb') != frozen['input_sha256'][state]:
            raise ValueError('Frozen predictor input differs from this MD construct')
    rows = [r for r in read_csv(path) if r['chain'] == 'B']
    methods = sorted({r['method'] for r in rows})
    if not {'full','static'}.issubset(methods) or not set(methods).issubset({'full','static','intra','aa'}):
        raise ValueError('Full and Static predictions are required; unknown method is not allowed')
    pred = {}
    for state, method in itertools.product(STATES, methods):
        selected = sorted([r for r in rows if r['state']==state and r['method']==method],
                          key=lambda r:int(r['author_residue_number']))
        keys = [(int(r['author_residue_number']), r['amino_acid']) for r in selected]
        if keys != [(r,a) for r,a,n in identities]:
            raise ValueError('Every predictor must cover all 87 KIX residues once')
        for metric, (_, column) in METRICS.items():
            values = [float(row[column]) for row, (_,_,n) in zip(selected, identities)
                      if metric=='rmsf' or n>0]
            pred[state, method, metric] = np.asarray(values)
            if not np.isfinite(pred[state, method, metric]).all():
                raise ValueError('Non-finite frozen prediction')
    return pred, methods, frozen


def analyze(root, output):
    policy = json.loads((root/'SCHEDULE_POLICY.json').read_text())
    if policy['old_rep01_included'] or policy['allowed_repeats'] != [1,2,3]:
        raise ValueError('This analysis requires three fresh repeats, with historical pilot excluded')
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Output must be empty to preserve previous analysis')
    md, identities, sources, all_rows = load_completed(root)
    pred, methods, frozen = load_frozen(root, identities)
    mdv, pv, numbers = {}, {}, {}
    for metric in METRICS:
        numbers[metric] = np.array([r for r,a,n in identities if metric=='rmsf' or n>0])
        for normalization in NORMALIZATIONS:
            transform = (lambda x:x) if normalization=='raw' else percentile
            for state, repeat in itertools.product(STATES, REPEATS):
                mdv[state, repeat, metric, normalization] = transform(md[state,repeat,metric])
            for state, method in itertools.product(STATES, methods):
                pv[state, method, metric, normalization] = transform(pred[state,method,metric])
    means = {key:np.mean([mdv[key[0],r,key[1],key[2]] for r in REPEATS],axis=0)
             for key in itertools.product(STATES,METRICS,NORMALIZATIONS)}
    absolute, mean_absolute, regional, delta, pair_scores, pair_rows, mean_rows, response_gain = [],[],[],[],[],[],[],[]
    for state,metric,normalization in itertools.product(STATES,METRICS,NORMALIZATIONS):
        ids = numbers[metric]
        for region in REGIONS:
            mask = np.ones(len(ids),dtype=bool) if region=='KIX' else ((ids>=614)&(ids<=621))
            for repeat in REPEATS:
                x=mdv[state,repeat,metric,normalization][mask]
                regional.append(dict(state=state,repeat=repeat,metric=metric,normalization=normalization,
                    region=region,n_residues=len(x),md_mean=float(x.mean()),
                    md_sqrt_mean_squared_rmsf=float(np.sqrt(np.mean(x*x))) if metric=='rmsf' and normalization=='raw' else ''))
                for method in methods:
                    absolute.append(dict(state=state,repeat=repeat,metric=metric,normalization=normalization,
                        region=region,method=method,**score(x,pv[state,method,metric,normalization][mask])))
            for method in methods:
                mean_absolute.append(dict(state=state,metric=metric,normalization=normalization,region=region,
                    method=method,md_repeats=3,**score(means[state,metric,normalization][mask],pv[state,method,metric,normalization][mask])))
        for i,residue in enumerate(ids):
            values=np.array([mdv[state,r,metric,normalization][i] for r in REPEATS])
            mean_rows.append(dict(state=state,author_residue_number=int(residue),metric=metric,normalization=normalization,
                md_mean=float(values.mean()),md_sd_three_repeats=float(values.std(ddof=1)),md_min=float(values.min()),md_max=float(values.max())))
    contrasts=[(b,a) for b,a in itertools.combinations(STATES,2)]
    for before,after in contrasts:
        contrast=after+'_minus_'+before
        for metric,normalization in itertools.product(METRICS,NORMALIZATIONS):
            ids=numbers[metric]
            mean_delta=means[after,metric,normalization]-means[before,metric,normalization]
            pdeltas={m:pv[after,m,metric,normalization]-pv[before,m,metric,normalization] for m in methods}
            for region in REGIONS:
                mask=np.ones(len(ids),dtype=bool) if region=='KIX' else ((ids>=614)&(ids<=621))
                zero=float(np.abs(mean_delta[mask]).mean())
                static=float(np.abs(pdeltas['static'][mask]-mean_delta[mask]).mean())
                pair_region=[]
                for rb,ra in itertools.product(REPEATS,REPEATS):
                    values=mdv[after,ra,metric,normalization]-mdv[before,rb,metric,normalization]
                    pair_region.append(float(values[mask].mean()))
                    for method in methods:
                        stats=score(values[mask],pdeltas[method][mask])
                        pair_scores.append(dict(contrast=contrast,before_repeat=rb,after_repeat=ra,region=region,
                            metric=metric,normalization=normalization,method=method,independent_repeats_per_state=3,
                            zero_response_mae=float(np.abs(values[mask]).mean()),**stats))
                for method in methods:
                    stats=score(mean_delta[mask],pdeltas[method][mask])
                    delta.append(dict(contrast=contrast,region=region,metric=metric,normalization=normalization,
                        method=method,independent_repeats_per_state=3,reference_combinations=9,
                        zero_response_mae=zero,mae_gain_over_zero=zero-stats['mae'],
                        static_delta_mae=static,mae_gain_over_static=static-stats['mae'],
                        sign_agreement_fraction=float(np.mean(np.sign(mean_delta[mask])==np.sign(pdeltas[method][mask]))),
                        regional_delta_min_over_9=float(min(pair_region)),regional_delta_max_over_9=float(max(pair_region)),
                        all_9_regional_deltas_same_sign=bool(min(pair_region)>0 or max(pair_region)<0),**stats))
            for rb,ra in itertools.product(REPEATS,REPEATS):
                values=mdv[after,ra,metric,normalization]-mdv[before,rb,metric,normalization]
                for i,residue in enumerate(ids):
                    pair_rows.append(dict(contrast=contrast,before_repeat=rb,after_repeat=ra,metric=metric,
                        normalization=normalization,author_residue_number=int(residue),md_delta=float(values[i]),
                        full_delta=float(pdeltas['full'][i]),static_delta=float(pdeltas['static'][i]),
                        independent_repeats_per_state=3,combinations_are_not_independent=True))
    interactions=[]
    for metric,normalization,region in itertools.product(METRICS,NORMALIZATIONS,REGIONS):
        ids=numbers[metric]; mask=np.ones(len(ids),dtype=bool) if region=='KIX' else ((ids>=614)&(ids<=621))
        x=means['both',metric,normalization]-means['myb',metric,normalization]-means['mll',metric,normalization]+means['apo',metric,normalization]
        for method in methods:
            y=pv['both',method,metric,normalization]-pv['myb',method,metric,normalization]-pv['mll',method,metric,normalization]+pv['apo',method,metric,normalization]
            interactions.append(dict(descriptor='both-myb-mll+apo',region=region,metric=metric,normalization=normalization,
                method=method,interpretation='Nonadditive simulated readout; not binding free-energy cooperativity',**score(x[mask],y[mask])))
    output.mkdir(parents=True,exist_ok=True)
    for name,rows in [('md_all_KIX_residues',all_rows),('absolute_per_repeat',absolute),('absolute_vs_mean_MD',mean_absolute),
        ('region_per_repeat',regional),('MD_per_residue_repeat_variation',mean_rows),('delta_vs_mean_MD',delta),
        ('delta_all_9_reference_combination_scores',pair_scores),('delta_all_9_reference_combination_residues',pair_rows),
        ('nonadditive_readout_descriptor',interactions)]:
        write_csv(output/(name+'.csv'),rows)
    shutil.copy2(root/'frozen_predictions/residue_predictions.csv',output/'frozen_predictions_source.csv')
    report=dict(status='all_12_new_batch_repeats_analyzed',batch=policy['batch'],historical_rep01_included=False,
        MD_repeats_per_state=3,MD_repeats_total=12,KIX_rmsf_sites=len(numbers['rmsf']),KIX_hchi_sites=len(numbers['hchi']),
        L12_rmsf_sites=int(np.sum((numbers['rmsf']>=614)&(numbers['rmsf']<=621))),
        L12_hchi_sites=int(np.sum((numbers['hchi']>=614)&(numbers['hchi']<=621))),
        hchi_mask='All residues with n_chi_angles>0; Ala/Gly unavailable, Pro included',
        raw_reference='First compute each repeat separately, then equally average the three per-residue readouts',
        percentile_reference='Within each state/repeat/metric use 100*(average_rank-1)/(n-1), then equally average three repeats; model ranks separately within each state/method',
        region_definition='KIX586..672; L12 614..621; arithmetic mean primary, sqrt(mean(RMSF^2)) separately named',
        delta_definition='Difference of state means; all 9 cross-repeat combinations retained for descriptive sensitivity; NOT 9 independent replicates',
        comparator='Zero response and frozen static model, same residues and no fit on new MD',
        uncertainty='Three independent starts per state; SD/min/max descriptive. No bootstrap CI or p value.',
        sources=sources,frozen_predictions=frozen,analysis_script_sha256=sha(Path(__file__)),
        output_sha256={p.name:sha(p) for p in output.glob('*.csv')})
    (output/'summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    report=analyze(args.root.resolve(strict=True),(args.output_dir or args.root/'triplicate_delta_analysis').resolve())
    print(json.dumps({k:report[k] for k in ['status','MD_repeats_total','KIX_rmsf_sites','KIX_hchi_sites','L12_hchi_sites']},indent=2))
