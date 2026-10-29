"""Evaluate true apo/holo differences from corrected chemical-identity Full OOF.

No training, response heads, rescaling, or trajectory modification.  The family
clean subset is audited against the existing conservative all-chain components.
"""
import json
import hashlib
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(r'E:\Project\FlexPPI')
OUT = Path(__file__).resolve().parent
CORE = ROOT / 'experiments/core_corrected_pipeline'
EXPORT = CORE / 'neural_exports/evaluation/full'
READY = CORE / 'training_payload/splits/all_ready.jsonl'
PAIR = ROOT / 'artifacts/data/dynarepo_pairs.jsonl'
FAMILY = ROOT / 'experiments/_sources/family_folds/record_family_folds.csv'
sys.path.insert(0,str(ROOT/'src'))
from flexppi.data.pairing import global_sequence_alignment

def read_jsonl(path):
    return [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def rank_corr(a,b):
    return float(spearmanr(a,b).statistic) if len(a)>2 and np.std(a)>1e-10 and np.std(b)>1e-10 else np.nan

def metrics(g, truth='md_delta', pred='model_delta'):
    y,p = g[truth].to_numpy(float),g[pred].to_numpy(float)
    strong=np.abs(y)>=0.1
    n=max(1,int(np.ceil(.2*len(y))))
    return dict(residues=len(y),mae_model=float(np.mean(np.abs(p-y))),
                mae_zero=float(np.mean(np.abs(y))),
                mae_gain_over_zero=float(np.mean(np.abs(y))-np.mean(np.abs(p-y))),
                spearman=rank_corr(p,y),bias=float(np.mean(p-y)),
                direction_all=float(np.mean(np.sign(p)==np.sign(y))),
                n_strong=int(strong.sum()),direction_strong=float(np.mean(np.sign(p[strong])==np.sign(y[strong]))) if strong.any() else np.nan,
                increase_recall=float(np.mean(p[y>=.1]>0)) if (y>=.1).any() else np.nan,
                decrease_recall=float(np.mean(p[y<=-.1]<0)) if (y<=-.1).any() else np.nan,
                top20_abs_overlap=len(set(np.argsort(np.abs(y))[-n:])&set(np.argsort(np.abs(p))[-n:]))/n,
                true_mean=float(y.mean()),pred_mean=float(p.mean()))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    records={r['sample_id']:r for r in read_jsonl(READY)}
    family=pd.read_csv(FAMILY)
    fam_index=family.set_index('sample_id')
    rows=[]; inventory=[]; source_files={READY,PAIR,FAMILY}
    checkpoints=set()
    for pair in read_jsonl(PAIR):
        a_id,h_id=pair['apo']['sample_id'],pair['holo']['sample_id']
        f_a=records.get(a_id,{}).get('metadata',{}).get('fold')
        f_h=records.get(h_id,{}).get('metadata',{}).get('fold')
        inv=dict(pair_id=pair['pair_id'],pdb_id=pair['metadata']['pdb_id'],apo_sample=a_id,holo_sample=h_id,
                 apo_fold=f_a,holo_fold=f_h,independent_repeats_per_state=1,within_trajectory_windows=5,
                 dataset='DynaRepo_test53',family_components='',family_train_record_count=np.nan)
        if a_id not in records or h_id not in records:
            inv.update(status='excluded_missing_corrected_ready_record',reason='Cannot use legacy target/prediction as corrected Full OOF')
            inventory.append(inv);continue
        if f_a!=f_h:
            inv.update(status='excluded_different_folds',reason='Differences between two different models are not a paired prediction')
            inventory.append(inv);continue
        components=set(fam_index.loc[[a_id,h_id],'component_id'])
        leaks=family[family.component_id.isin(components)&family.original_fold.ne(f_a)]
        inv.update(family_components=';'.join(sorted(components)),family_train_record_count=len(leaks),family_clean=len(leaks)==0)
        pa,ph=EXPORT/(a_id+'.npz'),EXPORT/(h_id+'.npz')
        if not pa.exists() or not ph.exists():
            inv.update(status='excluded_missing_full_oof',reason='Missing corrected cached prediction')
            inventory.append(inv);continue
        a,h=np.load(pa),np.load(ph)
        # Legacy pair indices refer to the old parser's residue list.  Restored
        # chemical types can change offsets; remap the CURRENT full sequences.
        ia0=np.flatnonzero(np.isin(a['chain'],pair['apo']['target_chains']))
        ih0=np.flatnonzero(np.isin(h['chain'],pair['holo']['target_chains']))
        sa,sh=''.join(a['amino_acid'][ia0]),''.join(h['amino_acid'][ih0])
        if sa==sh:
            ia,ih=ia0,ih0
            identity=coverage=1.0
        else:
            match=global_sequence_alignment(sa,sh)
            ia,ih=ia0[match.indices_a],ih0[match.indices_b]
            identity,coverage=match.identity,match.coverage
        if identity<.9 or coverage<.8:
            raise ValueError('Corrected sequence mapping below thresholds: '+pair['pair_id'])
        inv.update(legacy_mapped_residues=len(pair['apo_indices']),corrected_mapped_residues=len(ia),
                   corrected_mapping_identity=identity,corrected_mapping_coverage=coverage,
                   corrected_target_chain_sequences_identical=sa==sh)
        assert len(ia)==len(ih)
        assert np.array_equal(a['amino_acid'][ia],h['amino_acid'][ih]),pair['pair_id']
        inv.update(status='evaluated',reason='same heldout fold; corrected Full absolute state subtraction; no response head')
        inventory.append(inv)
        source_files.update([pa,ph])
        checkpoints.add(CORE/f'trained_checkpoints/full/fold_{f_a}/epoch_040.pt')
        for j,metric in enumerate(['rmsf','chi_entropy']):
            ok=a['valid'][ia,j]&h['valid'][ih,j]
            # Exclude structurally uninformative AGP from side-chain analysis.
            if j==1: ok &= ~np.isin(a['amino_acid'][ia],['A','G','P'])
            for x,z in zip(ia[ok],ih[ok]):
                rows.append(dict(pair_id=pair['pair_id'],pdb_id=inv['pdb_id'],fold=f_a,metric=metric,
                                 apo_sample=a_id,holo_sample=h_id,apo_index=int(x),holo_index=int(z),
                                 chain=str(h['chain'][z]),author_residue_id=str(h['author_residue_id'][z]),
                                 amino_acid=str(h['amino_acid'][z]),family_clean=inv['family_clean'],
                                 family_train_record_count=len(leaks),md_apo=float(a['target'][x,j]),
                                 md_holo=float(h['target'][z,j]),pred_apo=float(a['mean'][x,j]),pred_holo=float(h['mean'][z,j]),
                                 md_delta=float(h['target'][z,j]-a['target'][x,j]),
                                 model_delta=float(h['mean'][z,j]-a['mean'][x,j]),
                                 label_block_sd_delta=float(np.sqrt(a['label_variance'][x,j]+h['label_variance'][z,j]))))
    inventory=pd.DataFrame(inventory); inventory.to_csv(OUT/'paired_inventory.csv',index=False)
    data=pd.DataFrame(rows); data.to_csv(OUT/'test53_corrected_full_oof_residue_deltas.csv.gz',index=False,compression='gzip')
    per=[]
    for key,g in data.groupby(['pair_id','pdb_id','fold','family_clean','metric']):
        per.append(dict(zip(['pair_id','pdb_id','fold','family_clean','metric'],key),**metrics(g)))
    per=pd.DataFrame(per);per.to_csv(OUT/'test53_per_pair_metrics.csv',index=False)
    # First average chain sides within PDB (4ETP occurs twice), then PDB equally.
    numeric=[x for x in per.columns if x not in ['pair_id','pdb_id','fold','family_clean','metric']]
    system=per.groupby(['pdb_id','fold','family_clean','metric'])[numeric].mean().reset_index()
    system.to_csv(OUT/'test53_per_system_metrics.csv',index=False)
    rng=np.random.default_rng(20261025);summary=[]
    for cohort,sub in [('all_same_fold_oof',system),('family_clean_under_original_oof',system[system.family_clean])]:
        for metric,g in sub.groupby('metric'):
            result=dict(cohort=cohort,metric=metric,systems=len(g))
            for name in numeric:
                values=g[name].dropna().to_numpy(float)
                result[name]=float(values.mean()) if len(values) else np.nan
                if len(values):
                    boot=rng.choice(values,(5000,len(values)),replace=True).mean(axis=1)
                    result[name+'_ci_low'],result[name+'_ci_high']=np.quantile(boot,[.025,.975]).tolist()
            result['systems_beating_zero']=int((g.mae_gain_over_zero>0).sum())
            summary.append(result)
    pd.DataFrame(summary).to_csv(OUT/'test53_macro_summary.csv',index=False)
    # Corrected KIX predictions already exported for the existing four states.
    kix_source=ROOT/'paper/figures_v59/03_KIX与泛素局部动态验证/source_archive/data/kix_fold1_all_residues.csv'
    k=pd.read_csv(kix_source);krows=[]
    for target,reference in [('myb','apo'),('mll','apo'),('both','apo'),('both','myb'),('both','mll')]:
        left=k[k.state.eq(reference)];right=k[k.state.eq(target)]
        merged=right.merge(left,on=['metric','author_residue_number','amino_acid','model_fold'],suffixes=('_holo','_apo'))
        merged['md_delta']=merged.md_value_holo-merged.md_value_apo
        merged['model_delta']=merged.model_value_holo-merged.model_value_apo
        merged['contrast']=target+'-'+reference
        krows.append(merged)
    kd=pd.concat(krows,ignore_index=True)
    kd=kd[np.isfinite(kd.md_delta)&np.isfinite(kd.model_delta)].copy()
    kd.to_csv(OUT/'kix_corrected_fold1_residue_deltas.csv',index=False)
    km=[]
    for (contrast,metric),g in kd.groupby(['contrast','metric']):
        km.append(dict(contrast=contrast,metric=metric,model_fold=1,independent_repeats_per_state=1,
                       independent_family=False,**metrics(g)))
    pd.DataFrame(km).to_csv(OUT/'kix_corrected_fold1_delta_metrics.csv',index=False)
    source_files.add(kix_source)
    provenance=dict(protocol='corrected_full_absolute_state_subtraction_epoch40_seed7',
                    definition='delta = holo - apo; no response head or fitted calibration',
                    residue_mapping='Rebuild on CURRENT corrected cache amino_acid/chain arrays; select legacy target chain IDs, require complete chain sequence identity or exact Needleman-Wunsch correspondences (identity>=.90, coverage>=.80). Never apply old parser indices to corrected residues.',
                    sidechain_mask='both state valid & amino acid not A/G/P',
                    strong_thresholds={'rmsf_angstrom':0.1,'chi_entropy_nats':0.1},
                    independent_unit='PDB; average the two 4ETP sides before macro/bootstrap',
                    bootstrap='5000 PDB-resampled macro means; seed20261025',
                    family_clean_definition='No ready training record in either state all-chain conservative family component lies outside the original held-out fold. This is a selected conservative clean subset of original OOF, not the retrained family-OOF run.',
                    family_graph_rule=json.loads((ROOT/'experiments/_sources/family_folds/training_freeze.json').read_text())['rule'],
                    files=[dict(path=str(f),sha256=sha(f)) for f in sorted(source_files)],
                    checkpoints=[dict(path=str(f),sha256=sha(f)) for f in sorted(checkpoints)])
    (OUT/'paired_provenance.json').write_text(json.dumps(provenance,indent=2,ensure_ascii=False),encoding='utf-8')
    print(pd.DataFrame(summary)[['cohort','metric','systems','mae_model','mae_zero','spearman','direction_strong','systems_beating_zero']].to_string(index=False))
    print(pd.DataFrame(km)[['contrast','metric','mae_model','mae_zero','spearman']].to_string(index=False))

if __name__=='__main__':main()
