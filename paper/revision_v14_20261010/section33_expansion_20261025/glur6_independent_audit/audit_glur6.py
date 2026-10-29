"""Read-only audit of the 3H6G monomer/dimer candidate.

All differences use dimer minus monomer. Continuous trajectory blocks are
sensitivity checks, never independent simulation replicates.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
IDS = ['A00C7', 'A00DH']
MODELS = {
    'Full': 'pred_full_rmsf', 'Intra': 'pred_intra_rmsf',
    'AA_SASA_contacts': 'pred_aa_sasa_contacts_rmsf',
    'Static47_SASA': 'pred_static47_sasa_rmsf',
}

def correlation(a, b, method='spearman'):
    if len(a) < 3 or np.ptp(a) < 1e-15 or np.ptp(b) < 1e-15:
        return np.nan
    if method == 'spearman':
        a = pd.Series(a).rank(method='average').to_numpy()
        b = pd.Series(b).rank(method='average').to_numpy()
    return float(np.corrcoef(a,b)[0,1])

def score(d, label, scope, mask):
    v = d.loc[mask]
    y = v.delta_md.to_numpy()
    rows = []
    for model in list(MODELS) + ['Zero_change']:
        p = np.zeros(len(v)) if model == 'Zero_change' else v[f'delta_{model}'].to_numpy()
        err = np.abs(y-p)
        row = dict(dimer_target=label, scope=scope, model=model, n=len(v),
                   delta_mae=float(err.mean()) if len(v) else np.nan,
                   delta_median_abs_error=float(np.median(err)) if len(v) else np.nan,
                   delta_rmse=float(np.sqrt(np.mean((y-p)**2))) if len(v) else np.nan,
                   delta_spearman=correlation(y,p), delta_pearson=correlation(y,p,'pearson'),
                   mean_md_delta=float(y.mean()) if len(v) else np.nan,
                   median_md_delta=float(np.median(y)) if len(v) else np.nan,
                   md_rigidified_fraction=float(np.mean(y < 0)) if len(v) else np.nan,
                   mean_pred_delta=float(p.mean()) if len(v) else np.nan)
        for cutoff in [0, .1, .2, .3, .5]:
            keep = np.abs(y) >= cutoff
            row[f'effect_n_ge_{cutoff:g}'] = int(keep.sum())
            row[f'direction_accuracy_ge_{cutoff:g}'] = float(np.mean(np.sign(p[keep])==np.sign(y[keep]))) if keep.any() else np.nan
        rows.append(row)
    return rows

pred = pd.read_parquet(ROOT/'feature_baseline/oof_predictions_all.parquet', filters=[('sample_id','in',IDS)])
featcols = ['sample_id','residue_index','chain','author_residue_id','amino_acid','fold',
            'interface_ca8','min_partner_ca_distance_A','partner_ca_contacts_8A',
            'distance_to_chain_end','sidechain_sasa_complex_A2','sidechain_sasa_isolated_A2',
            'residue_sasa_complex_A2','delta_residue_sasa_A2']
feat = pd.read_parquet(ROOT/'feature_baseline/features_all.parquet', columns=featcols, filters=[('sample_id','in',IDS)])
meta_cols = [c for c in featcols if c not in pred or c in ['sample_id','residue_index']]
base = pred.merge(feat[meta_cols], on=['sample_id','residue_index'], validate='one_to_one')
assert base.fold.nunique()==1 and base.fold.iloc[0]==2
source_check=[]
for sid in IDS:
    g=base[base.sample_id==sid].sort_values('residue_index')
    source=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline/neural_exports/evaluation/full')/f'{sid}.npz'
    with np.load(source,allow_pickle=False) as z:
        assert len(g)==len(z['mean'])
        assert np.array_equal(g.residue_index.to_numpy(),np.arange(len(g)))
        model_error=float(np.max(np.abs(g.pred_full_rmsf.to_numpy()-z['mean'][:,0])))
        target_error=float(np.max(np.abs(g.md_rmsf.to_numpy()-z['target'][:,0])))
        assert model_error<1e-7 and target_error<1e-7
        source_check.append(dict(sample_id=sid,npz=str(source),full_prediction_max_error=model_error,target_max_error=target_error))
window = pd.concat([pd.read_csv(ROOT/f'verified_case_trajectories/{s}/all_residues_rmsf_per_window.csv') for s in IDS])
assert set(window.block)==set(range(1,6))
assert not window.duplicated(['sample_id','residue_index','block']).any()
assert np.max(np.abs(window.groupby(['sample_id','residue_index']).rmsf.mean().to_numpy() - window.groupby(['sample_id','residue_index']).frozen_rmsf.first().to_numpy())) < 2e-6
mono = base[(base.sample_id=='A00DH')&(base.chain=='B')].sort_values('residue_index').reset_index(drop=True)
assert len(mono)==383
mw = window[window.sample_id=='A00DH'].pivot(index='residue_index',columns='block',values='rmsf').sort_index().to_numpy()
tables = []
for chain in ['A','B']:
    bound = base[(base.sample_id=='A00C7')&(base.chain==chain)].sort_values('residue_index').reset_index(drop=True)
    assert len(bound)==383
    assert bound.amino_acid.tolist()==mono.amino_acid.tolist()
    # Prepared files number chains consecutively: chain B in the dimer starts
    # at 384, while the extracted monomer starts at 1. Align exact sequences.
    bw = window[(window.sample_id=='A00C7')&(window.chain==chain)].pivot(index='residue_index',columns='block',values='rmsf').sort_index().to_numpy()
    d = pd.DataFrame({'dimer_target':chain, 'sequence_position':np.arange(1,384),
                      'author_residue_id':bound.author_residue_id,
                      'monomer_author_residue_id':mono.author_residue_id,'amino_acid':bound.amino_acid,
                      'monomer_md':mono.md_rmsf,'dimer_md':bound.md_rmsf,
                      'delta_md':bound.md_rmsf-mono.md_rmsf,
                      'distance_to_chain_end':bound.distance_to_chain_end,
                      'interface_ca8':bound.interface_ca8,
                      'min_partner_ca_distance_A':bound.min_partner_ca_distance_A,
                      'partner_ca_contacts_8A':bound.partner_ca_contacts_8A,
                      'dimer_delta_residue_sasa_A2':bound.delta_residue_sasa_A2,
                      'dimer_residue_sasa_complex_A2':bound.residue_sasa_complex_A2,
                      'monomer_residue_sasa_complex_A2':mono.residue_sasa_complex_A2})
    for model, col in MODELS.items():
        d[f'monomer_{model}'] = mono[col]
        d[f'dimer_{model}'] = bound[col]
        d[f'delta_{model}'] = bound[col]-mono[col]
    for i in range(5):
        d[f'monomer_block_{i+1}'] = mw[:,i]
        d[f'dimer_block_{i+1}'] = bw[:,i]
    tables.append(d)

avg = tables[0].copy()
avg['dimer_target'] = 'AB_mean'
numeric = [c for c in avg if c.startswith(('dimer_','delta_')) and c!='dimer_target']
numeric += ['min_partner_ca_distance_A','partner_ca_contacts_8A']
for c in numeric:
    avg[c] = (tables[0][c].to_numpy()+tables[1][c].to_numpy())/2
avg.interface_ca8 = tables[0].interface_ca8 | tables[1].interface_ca8
tables.append(avg)

allmetrics, blockmetrics, regional, stability, centered = [], [], [], [], []
for d in tables:
    label = d.dimer_target.iloc[0]
    delta_blocks = np.stack([d[f'dimer_block_{i}']-d[f'monomer_block_{i}'] for i in range(1,6)],axis=1)
    # Cross all bound/monomer windows so an arbitrary same-index pairing is not decisive.
    cross = np.stack([d[f'dimer_block_{i}']-d[f'monomer_block_{j}'] for i in range(1,6) for j in range(1,6)],axis=1)
    d['paired_blocks_negative_n'] = (delta_blocks<0).sum(axis=1)
    d['paired_blocks_direction_consistent_n'] = (np.sign(delta_blocks)==np.sign(d.delta_md.to_numpy()[:,None])).sum(axis=1)
    d['cross25_negative_n'] = (cross<0).sum(axis=1)
    d['cross25_direction_consistent_n'] = (np.sign(cross)==np.sign(d.delta_md.to_numpy()[:,None])).sum(axis=1)
    d['paired_block_min_delta'] = delta_blocks.min(axis=1)
    d['paired_block_max_delta'] = delta_blocks.max(axis=1)
    d['cross25_min_delta'] = cross.min(axis=1)
    d['cross25_max_delta'] = cross.max(axis=1)
    masks = {'all':np.ones(len(d),bool)}
    for end in [3,10,20]:
        masks[f'exclude_terminal_{end}'] = d.distance_to_chain_end>=end
    for frac in [.05,.1,.2]:
        # Ranking on truth, applied equally to every model.
        rank = np.argsort(-np.abs(d.delta_md.to_numpy()))
        keep = np.ones(len(d),bool)
        keep[rank[:int(np.ceil(frac*len(d)))]] = False
        masks[f'exclude_top_abs_delta_{frac:g}'] = keep
    dist = d.min_partner_ca_distance_A
    masks.update(interface_ca8=d.interface_ca8,
                 noninterface=~d.interface_ca8,
                 partner_distance_le8=dist<=8,
                 partner_distance_8to12=(dist>8)&(dist<=12),
                 partner_distance_12to20=(dist>12)&(dist<=20),
                 partner_distance_gt20=dist>20)
    for scope, mask in masks.items():
        allmetrics += score(d,label,scope,mask)
    for block in range(1,6):
        db = d.copy()
        db.delta_md = delta_blocks[:,block-1]
        for scope in ['all','exclude_terminal_10','interface_ca8','noninterface']:
            for row in score(db,label,scope,masks[scope]):
                row['block'] = block
                blockmetrics.append(row)
    for start in range(1,384,25):
        end = min(start+24,383)
        mask = d.sequence_position.between(start,end)
        for row in score(d,label,f'positions_{start}_{end}',mask):
            row.update(start=start,end=end)
            regional.append(row)
    absdelta = np.abs(d.delta_md.to_numpy())
    order = np.sort(absdelta)[::-1]
    ss = dict(dimer_target=label,n=383,
              mean_delta=d.delta_md.mean(),median_delta=d.delta_md.median(),
              negative_n=int((d.delta_md<0).sum()),positive_n=int((d.delta_md>0).sum()),
              paired_all5_consistent_n=int((d.paired_blocks_direction_consistent_n==5).sum()),
              paired_ge4_consistent_n=int((d.paired_blocks_direction_consistent_n>=4).sum()),
              cross25_all_consistent_n=int((d.cross25_direction_consistent_n==25).sum()),
              cross25_ge20_consistent_n=int((d.cross25_direction_consistent_n>=20).sum()),
              top5percent_abs_effect_share=float(order[:int(np.ceil(.05*383))].sum()/absdelta.sum()),
              top10percent_abs_effect_share=float(order[:int(np.ceil(.1*383))].sum()/absdelta.sum()),
              terminal10_abs_effect_share=float(absdelta[d.distance_to_chain_end<10].sum()/absdelta.sum()))
    for cutoff in [.1,.2,.3,.5]:
        mask=absdelta>=cutoff
        ss[f'n_abs_delta_ge_{cutoff:g}']=int(mask.sum())
        ss[f'paired_all5_consistent_ge_{cutoff:g}']=int(((d.paired_blocks_direction_consistent_n==5)&mask).sum())
        ss[f'cross25_ge20_consistent_ge_{cutoff:g}']=int(((d.cross25_direction_consistent_n>=20)&mask).sum())
    stability.append(ss)
    y = d.delta_md.to_numpy()
    for model in list(MODELS)+['Zero_change','Posthoc_oracle_median_constant']:
        if model == 'Zero_change':
            p=np.zeros(len(d))
        elif model == 'Posthoc_oracle_median_constant':
            p=np.full(len(d), np.median(y))
        else:
            p=d[f'delta_{model}'].to_numpy()
        cy=y-y.mean()
        cp=p-p.mean()
        row=dict(dimer_target=label, model=model, n=len(d), md_mean=y.mean(),
                 pred_mean=p.mean(), md_median=np.median(y),
                 delta_mae=float(np.mean(np.abs(y-p))),
                 mean_centered_mae=float(np.mean(np.abs(cy-cp))),
                 mean_centered_rmse=float(np.sqrt(np.mean((cy-cp)**2))),
                 delta_spearman=correlation(y,p),
                 note='Posthoc descriptor using this case truth, not a trained or fair predictive baseline.' if model.startswith('Posthoc') else 'Both vectors separately centered on their own all-residue mean; diagnostic only.')
        for cutoff in [0,.1,.2,.3,.5]:
            mask=np.abs(y)>=cutoff
            row[f'n_ge_{cutoff:g}']=int(mask.sum())
            row[f'direction_accuracy_ge_{cutoff:g}']=float(np.mean(np.sign(p[mask])==np.sign(y[mask]))) if mask.any() else np.nan
            row[f'always_decrease_accuracy_ge_{cutoff:g}']=float(np.mean(y[mask]<0)) if mask.any() else np.nan
        centered.append(row)

pd.concat(tables).to_csv(OUT/'all_383_residue_deltas_and_windows.csv',index=False)
pd.DataFrame(allmetrics).to_csv(OUT/'delta_prediction_metrics.csv',index=False)
pd.DataFrame(blockmetrics).to_csv(OUT/'timeblock_prediction_metrics.csv',index=False)
pd.DataFrame(regional).to_csv(OUT/'prespecified_25_residue_regions.csv',index=False)
pd.DataFrame(stability).to_csv(OUT/'effect_spread_and_window_stability.csv',index=False)
pd.DataFrame(centered).to_csv(OUT/'centered_shape_and_posthoc_constant_diagnostics.csv',index=False)
fam=pd.read_csv('E:/Project/FlexPPI/experiments/_sources/family_folds/record_family_folds.csv')
rows=fam[fam.sample_id.isin(IDS)]
components=rows.component_id.unique()
members=fam[fam.component_id.isin(components)]
members.to_csv(OUT/'family_component_members.csv',index=False)
audit=dict(samples=IDS,fold=2,sequence_length=383,all_target_sequences_identical=True,
           direction='dimer minus monomer',AB_mean='Arithmetic mean of A/B labels and predictions, then subtract the same monomer.',
           AB_interface='Either A or B has another-chain CA within 8 Angstrom.',
           AB_distance='Arithmetic mean of the A/B minimum partner CA distances.',
           family_clean=bool((members.original_fold==2).all()),family_component_members=len(members),
           full_export_source_check=source_check,
           max_block_mean_label_error=float(np.max(np.abs(window.groupby(['sample_id','residue_index']).rmsf.mean().to_numpy()-window.groupby(['sample_id','residue_index']).frozen_rmsf.first().to_numpy()))),
           uncertainty='Five consecutive 100 ns blocks from one 500 ns trajectory per state, not independent simulation replicas. No inferential P values or block bootstrap confidence intervals.',
           selection='Exploratory selected case; no population-generalization claim from this example alone.')
(OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(pd.DataFrame(allmetrics).query("scope == 'all'").to_string(index=False))
print(pd.DataFrame(stability).to_string(index=False))
