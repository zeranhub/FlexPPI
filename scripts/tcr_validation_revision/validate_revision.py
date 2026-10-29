"""Independent numerical, protocol, and figure-export checks for the revision."""
from pathlib import Path
import hashlib,json,re,subprocess,sys,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/tcr_prediction_reference_revision_20260915'
BASE=ROOT/'artifacts/tcr_validation_revision/1g4_reference'
sys.path.insert(0,str(ROOT/'src'))
from flexppi.data.labels import _chainwise_rmsf,_kabsch_fit,_dihedral_degrees
from MDAnalysis.lib.distances import calc_dihedrals

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    # Check geometric convention without affecting the real data or selections.
    rng=np.random.default_rng(20260915);points=rng.normal(size=(100,4,3))
    a=_dihedral_degrees(points);b=np.degrees(calc_dihedrals(points[:,0],points[:,1],points[:,2],points[:,3]))
    np.testing.assert_allclose(np.abs((a-b+180)%360-180),0,atol=1e-4)
    bb=rng.normal(size=(10,24,3)).astype(np.float32);chain=np.r_[np.zeros(12),np.ones(12)]
    rotation=np.linalg.qr(rng.normal(size=(3,3)))[0]
    np.testing.assert_allclose(_chainwise_rmsf(bb,chain),_chainwise_rmsf(bb@rotation+[3,-5,8],chain),atol=1e-5)
    original=pd.read_csv(ROOT/'output/tcr_real_states_20260915/b423_existing_holo_md_comparison.csv',dtype={'author_id':str})
    d=pd.read_csv(OUT/'validation_paired_residues.csv',dtype={'author_id':str});v=d[d.evaluation_mask]
    both=d.merge(original,on=['chain','author_id','metric'],suffixes=('','_source'),validate='one_to_one')
    np.testing.assert_allclose(both.prediction,both.predicted_deposited_fold3_source,rtol=0,atol=1e-12)
    np.testing.assert_allclose(both.md_reference,both.md_reference_source,rtol=0,atol=1e-12)
    assert v.groupby('metric').size().to_dict()=={'hchi':340,'rmsf':428}
    m=pd.read_csv(OUT/'validation_metrics.csv')
    for metric,g in v.groupby('metric'):
        row=m[(m.scope=='TCR alpha beta')&(m.metric==metric)&(m.predictor=='FlexPPI')].iloc[0]
        np.testing.assert_allclose(abs(g.prediction-g.md_reference).mean(),row.mae)
        np.testing.assert_allclose(spearmanr(g.prediction,g.md_reference).statistic,row.spearman)
    a=json.loads((OUT/'1g4_validation_audit.json').read_text())
    assert len(a['runs'])==20 and set((r['state'],r['run']) for r in a['runs'])=={(s,i) for s in ['apo','holo'] for i in range(1,11)}
    for r in a['runs']:
        assert r['n_frames']==641 and r['snapshot_time_ps']==10000 and r['time_start_ps']==20000 and r['time_end_ps']==100000
        assert sha(r['snapshot'])==r['snapshot_sha256']
        assert sha(r['source_pdb'])==r['pdb_sha256'] and sha(r['source_xtc'])==r['xtc_sha256']
    runs=pd.read_csv(OUT/'1g4_run_paired_residues.csv',dtype={'author_id':str});valid=runs[runs.evaluation_mask]
    assert len(runs)==20*444*2 and len(valid)==20*(444+348)
    assert valid.groupby(['state','run','metric']).size().value_counts().to_dict()=={444:20,348:20}
    assert not valid.duplicated(['state','run','chain','author_id','metric']).any()
    delta=pd.read_csv(OUT/'1g4_delta_paired_residues.csv')
    np.testing.assert_allclose(delta.md_delta,delta.md_mean_holo-delta.md_mean_apo)
    np.testing.assert_allclose(delta.prediction_delta,delta.prediction_mean_holo-delta.prediction_mean_apo)
    split=ROOT/'output/atlas_hchi_retraining_v1_linux/payload/data/splits/state_rmsf_hchi_v2_atlas_hchi_v1'
    split_members={}
    for role in ['train','test']:
        records=[json.loads(x) for x in (split/f'fold_0_{role}.jsonl').read_text().splitlines()]
        split_members[role]={r['sample_id'] for r in records}
        pdbs={str(p).upper() for r in records for p in r.get('metadata',{}).get('pdb_ids',[])}
        if role=='train':assert not {'2BNQ','2BNR'}&pdbs
    assert not split_members['train']&split_members['test']
    config_audit={}
    for fold in [0,3]:
        path=ROOT/f'artifacts/performance_first_v2_atlas_hchi_v1/fold_{fold}/resolved_config.json'
        config=json.loads(path.read_text());assert not config['model']['use_condition_features']
        config_audit[str(fold)]=dict(path=str(path),sha256=sha(path),condition_features_enabled=False)
    # Per-run performance: never imply ensemble-mean accuracy is single-input accuracy.
    priors=pd.read_csv(ROOT/'artifacts/evaluation_atlas_hchi_v1/analysis/training_fold_aa_priors.csv')
    priors=priors[priors.fold==0].copy();priors.metric=priors.metric.replace({'chi_entropy':'hchi'})
    valid=valid.merge(priors[['metric','amino_acid','training_mean']],on=['metric','amino_acid'],validate='many_to_one')
    individual=[]
    for (state,run,metric),g in valid.groupby(['state','run','metric']):
        individual.append(dict(state=state,run=run,metric=metric,n=len(g),mae=float(abs(g.prediction-g.reference).mean()),
            spearman=float(spearmanr(g.prediction,g.reference).statistic),aa_only_mae=float(abs(g.training_mean-g.reference).mean())))
    pd.DataFrame(individual).to_csv(OUT/'1g4_individual_run_metrics.csv',index=False)
    # Resample whole trajectories, preserving each input/reference pair and all
    # correlated residues. Apo and holo sets are sampled independently.
    boot=[]
    for metric in ['rmsf','hchi']:
        values={}
        for state in ['apo','holo']:
            g=valid[(valid.state==state)&(valid.metric==metric)]
            values[state]={col:g.pivot(index='run',columns=['chain','author_id'],values=col).sort_index(axis=1).to_numpy() for col in ['prediction','reference']}
        results=[]
        for _ in range(2000):
            ia=rng.integers(0,10,size=10);ih=rng.integers(0,10,size=10)
            md=values['holo']['reference'][ih].mean(0)-values['apo']['reference'][ia].mean(0)
            pred=values['holo']['prediction'][ih].mean(0)-values['apo']['prediction'][ia].mean(0)
            results.append(abs(md).mean()-abs(pred-md).mean())
        g=delta[delta.metric==metric]
        boot.append(dict(metric=metric,mae_gain_over_zero=float(abs(g.md_delta).mean()-abs(g.prediction_delta-g.md_delta).mean()),
            percentile95_low=float(np.quantile(results,.025)),percentile95_high=float(np.quantile(results,.975)),n_bootstrap=2000,seed=20260915,
            note='Positive gain favors FlexPPI. Whole-run paired input/reference bootstrap, independent apo/holo sampling; descriptive sampling interval, not a residue-level p-value.'))
    pd.DataFrame(boot).to_csv(OUT/'1g4_delta_mae_gain_bootstrap.csv',index=False)
    exports=[]
    for f in json.loads((OUT/'figure_export_audit.json').read_text()):
        assert not f['outside_text'],f
        name=f['name']
        for ext in ['png','pdf','svg']:assert (OUT/(name+'.'+ext)).stat().st_size>5000
        svg=ET.parse(OUT/(name+'.svg'));nt=len(svg.findall('.//{http://www.w3.org/2000/svg}text'));assert nt>=15
        info=subprocess.check_output(['pdfinfo',str(OUT/(name+'.pdf'))],text=True)
        w,h=map(float,re.search(r'Page size:\s+([0-9.]+) x ([0-9.]+)',info).groups())
        assert abs(w*25.4/72-183)<.1 and h*25.4/72<247
        fonts=subprocess.check_output(['pdffonts',str(OUT/(name+'.pdf'))],text=True)
        assert 'TrueType' in fonts and 'Type 3' not in fonts
        im=Image.open(OUT/(name+'.png'));assert min(im.info.get('dpi',(0,0)))>=299
        exports.append(dict(name=name,svg_text=nt,width_mm=w*25.4/72,height_mm=h*25.4/72,
                            hashes={ext:sha(OUT/(name+'.'+ext)) for ext in ['png','pdf','svg']}))
    report=dict(status='passed',n_figures=len(exports),b423_counts={'rmsf':428,'hchi':340},ig4_runs=20,ig4_counts={'rmsf':444,'hchi':348},
        dihedral_convention_verified=True,rmsf_rigid_invariance=True,no_new_MD_or_training=True,ig4_bootstrap=boot,exports=exports,config_audit=config_audit,
        scope_warnings=['B4 absolute holo only, existing combination-level holdout, not untouched external validation.',
        '1G4 has no CD3 and uses fixed MD-derived snapshots; does not validate an MD-free PDB-only workflow.',
        '1G4 ensemble-mean figures average 10 predictions per state; per-run performance is separately supplied.',
        'All sampled proteins remain related to training systems; no unseen-family claim.'])
    (OUT/'deliverable_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='exports'},indent=2),flush=True)

if __name__=='__main__':main()
