"""Compare existing 5IVX MD labels with the held-out checkpoint, not apo/holo MD.

This does not run or create MD. It checks original-file hashes and the frozen
fold-3 split, and maps the MD topology to deposited author IDs by exact full
chain sequence equality. Both the topology-input and deposited-input predictions
are retained. Nothing here validates the predicted binding-induced delta.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from analyze_predictions import ROOT,BASE,OUT


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    cache=ROOT/'artifacts/evaluation_atlas_hchi_v1/oof/fold_3/itrace_5ivx_050ba0aa59cc.npz'
    with np.load(cache,allow_pickle=False) as data:
        z={k:data[k] for k in data.files}
    meta=json.loads(str(z['metadata_json']))
    assert meta['fold']==3 and meta['sample_id']=='itrace_5ivx'
    assert meta['protocol']['epoch']==40 and meta['protocol']['name']=='atlas_hchi_v1_final_epoch40_oof'
    checkpoint=ROOT/'artifacts/performance_first_v2_atlas_hchi_v1/fold_3/epoch_040.pt'
    for k,path in [('checkpoint',checkpoint),('input',meta['structure_path']),('label',meta['label_path'])]:
        assert sha(path)==meta['identity'][k+'_sha256'], 'Original input/checkpoint/label changed'
    payload=ROOT/'output/atlas_hchi_retraining_v1_linux/payload/data/splits/state_rmsf_hchi_v2_atlas_hchi_v1'
    for role in ['train','test']:
        records=[json.loads(x) for x in (payload/f'fold_3_{role}.jsonl').read_text().splitlines() if x.strip()]
        found=any(r['sample_id']=='itrace_5ivx' for r in records)
        assert found == (role=='test')
    valid=z['label_mask'] & z['target_mask'][:,None] & z['geometry_mask'][:,None] & z['output_mask'][:,None]
    assert np.array_equal(valid,z['valid'])
    assert np.isfinite(z['mean'][valid]).all() and np.isfinite(z['target'][valid]).all()
    spec=next(s for s in json.loads((BASE/'states.json').read_text()) if s['id']=='b4_5ivx_full_full')
    with np.load(BASE/'prediction_cache/b4_5ivx_full_full_fold3.npz',allow_pickle=False) as new:
        predictions={(str(c),str(a)):v for c,a,v in zip(new['chain'],new['author'],new['prediction'])}
    mapping,rows=[],[]
    for ch in np.unique(z['chain']):
        indices=np.where(z['chain']==ch)[0]
        seq=''.join(z['amino_acid'][indices])
        matches=[c for c,r in spec['residues'].items() if r['sequence']==seq]
        assert len(matches)==1, 'Need an unambiguous complete-sequence chain match'
        newchain=matches[0]
        mapping.append(dict(topology_chain=str(ch),deposited_chain=newchain,n_residues=len(indices),mapping='Exact full-chain sequence equality, no alignment gaps.'))
        for idx,author in zip(indices,spec['residues'][newchain]['author_ids']):
            for mi,metric in enumerate(['rmsf','hchi']):
                rows.append(dict(chain=newchain,author_id=author,topology_chain=ch,topology_author_id=str(z['author_residue_id'][idx]),
                    amino_acid=str(z['amino_acid'][idx]),metric=metric,valid=bool(valid[idx,mi]),
                    md_reference=float(z['target'][idx,mi]),predicted_deposited_fold3=float(predictions[(newchain,author)][mi]),
                    predicted_topology_oof=float(z['mean'][idx,mi]),label_variance=float(z['label_variance'][idx,mi]),
                    scope='holo only; existing MD aggregate; fold 3, not independent family holdout'))
    df=pd.DataFrame(rows);df.to_csv(OUT/'b423_existing_holo_md_comparison.csv',index=False)
    summaries=[]
    for scope,chs in [('TCR alpha beta','EF'),('TCR beta','F'),('Full protein context','EFABP')]:
        for metric in ['rmsf','hchi']:
            q=df[df.chain.isin(list(chs))&df.metric.eq(metric)&df.valid]
            if metric=='hchi':q=q[~q.amino_acid.isin(list('AGP'))]
            for predictor in ['predicted_deposited_fold3','predicted_topology_oof']:
                err=q[predictor]-q.md_reference
                pred_resid=q[predictor]-q.groupby('amino_acid')[predictor].transform('mean')
                md_resid=q.md_reference-q.groupby('amino_acid').md_reference.transform('mean')
                summaries.append(dict(scope=scope,metric=metric,predictor=predictor,n_residues=len(q),
                    spearman=float(spearmanr(q[predictor],q.md_reference).statistic),mae=float(abs(err).mean()),
                    rmse=float(np.sqrt((err**2).mean())),within_case_AA_centered_spearman=float(spearmanr(pred_resid,md_resid).statistic)))
    pd.DataFrame(summaries).to_csv(OUT/'b423_holo_md_metrics.csv',index=False)
    report=dict(cache=str(cache),cache_sha256=sha(cache),checkpoint_sha256=sha(checkpoint),chain_mapping=mapping,
        source='Existing ITrace 5IVX labels, protocol rmsf_existing_chainfit+hchi_120deg_joint_rotamer_3replica; not new simulations.',
        n_replicates_from_frozen_manifest=3,scope='Absolute holo-state validation within the existing combination-level split; no matched apo MD.',
        caveats=['Same system had historical diagnostic use; this is not untouched external validation.',
            'Fold 3 excludes itrace_5ivx but has homologous training sequences.',
            'MD topology-input prediction and deposited-input prediction are distinguished.',
            'AA-centered correlation is a descriptive within-case calculation, not a training-only baseline or held-out adjusted model.',
            'Aggregate label variance is stored but not converted to an experimental error bar or confidence interval.'])
    (OUT/'b423_holo_md_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)


if __name__=='__main__':
    main()
