"""Direct 1G4 validation with predeclared MD-derived static inputs.

Knapp 2019 publicly released 10 free-TCR and 10 TCR-pMHC trajectories.
Use every run; the snapshot nearest 10 ns is the sole input, and all stored
frames from 20 through 100 ns form the reference. No frame is selected by fit
to predictions. This is not CD3 validation or an MD-free deployment benchmark.
"""
from pathlib import Path
import argparse, hashlib, json, sys, tarfile, warnings
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors, calc_dihedrals

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'scripts'), str(ROOT/'scripts/tcr_real_states')]
from flexppi.data.labels import _chainwise_rmsf, _entropy
from flexppi.constants import THREE_TO_ONE, CHI_ATOMS
from flexppi.data.structures import parse_pdb
BASE = ROOT/'artifacts/tcr_validation_revision/1g4_reference'
OUT = ROOT/'output/tcr_prediction_reference_revision_20260915'
MD5 = 'b460344a3df16c93913ac0b239f5d295'
SNAPSHOT_PS = 10000.0
REFERENCE_START_PS = 20000.0

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def extract():
    archive=BASE/'1G4.tar.gz'
    assert hashlib.md5(archive.read_bytes()).hexdigest()==MD5
    dest=(BASE/'extracted').resolve();dest.mkdir(parents=True,exist_ok=True)
    with tarfile.open(archive) as t:
        for member in t.getmembers():
            target=(dest/member.name).resolve()
            assert dest in target.parents and (member.isdir() or member.isfile())
            if member.isdir():target.mkdir(parents=True,exist_ok=True)
            elif not target.exists():
                target.parent.mkdir(parents=True,exist_ok=True)
                with t.extractfile(member) as src, target.open('wb') as out:
                    import shutil
                    shutil.copyfileobj(src,out)
    return dest

def unwrap_chain_ca(ca, box):
    corrected=ca.astype(float).copy()
    if box is not None:
        corrected[1:]=corrected[0]+np.cumsum(minimize_vectors(np.diff(corrected,axis=0),box),axis=0)
    return corrected

def write_snapshot(u, path, anchor_chain):
    """Write deposited MD snapshot backbone with periodic images made whole.
    Atom/residue identities are unchanged; chain image choice is a box translation.
    """
    box=u.dimensions
    chains={str(s.segid):s.residues for s in u.segments}
    ca_by={c:unwrap_chain_ca(np.array([r.atoms.select_atoms('name CA').positions[0] for r in rs]),box) for c,rs in chains.items()}
    anchor=ca_by[anchor_chain].mean(0)
    lines=['REMARK 900 FIXED 10 NS MD SNAPSHOT; NOT AN EXPERIMENTAL APO STRUCTURE',
           'REMARK 900 BACKBONE ONLY; PERIODIC IMAGES MADE WHOLE; NO MINIMIZATION']
    serial=1
    for ch,rs in chains.items():
        ca=ca_by[ch];center=ca.mean(0)
        shift=minimize_vectors(np.array([center-anchor]),box)[0]-(center-anchor) if box is not None else np.zeros(3)
        for r,new_ca in zip(rs,ca):
            actual_ca=r.atoms.select_atoms('name CA').positions[0]
            names=set(r.atoms.names)
            for name in ['N','CA','C','O']:
                source_name='O1' if name=='O' and 'O' not in names and 'O1' in names else name
                if name=='O' and source_name not in names and 'OC1' in names:source_name='OC1'
                assert source_name in names,(ch,r.resid,name)
                pos=r.atoms.select_atoms('name '+source_name).positions[0]
                delta=pos-actual_ca
                if box is not None:delta=minimize_vectors(np.array([delta]),box)[0]
                xyz=new_ca+delta+shift
                lines.append(f'ATOM  {serial:5d} {name:^4s} {r.resname:>3s} {ch:1s}{int(r.resid):4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {name[0]:>2s}')
                serial+=1
        lines.append('TER')
    lines.append('END');path.write_text('\n'.join(lines)+'\n',encoding='ascii')

def reference_one(pdb):
    run=int(pdb.parent.name.split('run')[-1]);state='holo' if 'TCRpMHC' in pdb.parent.name else 'apo'
    target_chains=['D','E'] if state=='holo' else ['A','B']
    xtc=next(pdb.parent.glob('*.xtc'))
    cache=BASE/f'{state}_run{run:02d}_reference.npz'
    snapshot=BASE/f'{state}_run{run:02d}_input10ns.pdb'
    if cache.exists() and snapshot.exists():
        with np.load(cache) as z:
            old=json.loads(str(z['metadata']))
            assert old['snapshot_time_ps']==SNAPSHOT_PS and old['time_start_ps']==REFERENCE_START_PS
            assert sha(snapshot)==old['snapshot_sha256']
        return cache
    u=mda.Universe(str(pdb),str(xtc))
    residues=[r for s in u.segments if s.segid in target_chains for r in s.residues]
    assert len(residues)==444
    aa=np.array([THREE_TO_ONE[r.resname] for r in residues])
    chain=np.array(['alpha' if r.segid==target_chains[0] else 'beta' for r in residues])
    ids=np.array([str(r.resid) for r in residues]);ci=(chain=='beta').astype(int)
    caids=np.array([r.atoms.select_atoms('name CA').indices[0] for r in residues])
    quad=[];rchi={}
    for i,r in enumerate(residues):
        lookup={str(a.name):a.index for a in r.atoms};local=[]
        for names in CHI_ATOMS.get(r.resname,[]):
            if all(n in lookup for n in names):
                local.append(len(quad));quad.append([lookup[n] for n in names])
        if local and len(local)==len(CHI_ATOMS[r.resname]):rchi[i]=local
    quad=np.array(quad);frames=[];codes=[];times=[];max_ca_step=0.;snapshot_time=None
    all_times=np.array([ts.time for ts in u.trajectory])
    assert all_times[-1]>=99000 and all_times[0]<=SNAPSHOT_PS
    snap_idx=int(np.argmin(np.abs(all_times-SNAPSHOT_PS)));u.trajectory[snap_idx]
    write_snapshot(u,snapshot,target_chains[0]);snapshot_time=float(u.trajectory.ts.time)
    for ts in u.trajectory:
        if ts.time < REFERENCE_START_PS-1e-3:continue
        xyz=u.atoms.positions[caids].copy()
        for c in [0,1]:
            ix=np.flatnonzero(ci==c);xyz[ix]=unwrap_chain_ca(xyz[ix],u.dimensions)
            max_ca_step=max(max_ca_step,float(np.linalg.norm(np.diff(xyz[ix],axis=0),axis=1).max()))
        assert max_ca_step<5.,'Unresolved chain break/PBC: '+str(max_ca_step)
        frames.append(xyz);times.append(float(ts.time))
        q=u.atoms.positions[quad]
        angles=np.degrees(calc_dihedrals(q[:,0],q[:,1],q[:,2],q[:,3],box=u.dimensions))
        assert np.isfinite(angles).all()
        codes.append(np.floor(((angles+180.)%360.)/120.).astype(np.int8))
    arr=np.stack(frames);rmsf=_chainwise_rmsf(arr,ci)
    cc=np.stack(codes);entropy=np.full(len(residues),np.nan)
    for i,angles in rchi.items():
        joint=sum(cc[:,j].astype(np.int32)*(3**p) for p,j in enumerate(angles));entropy[i]=_entropy(joint)
    metadata=dict(state=state,run=run,n_frames=len(times),time_start_ps=times[0],time_end_ps=times[-1],
                  snapshot_time_ps=snapshot_time,pdb_sha256=sha(pdb),xtc_sha256=sha(xtc),snapshot_sha256=sha(snapshot),
                  max_adjacent_CA_A=max_ca_step,source_pdb=str(pdb),source_xtc=str(xtc),snapshot=str(snapshot),
                  target_chains=target_chains,ca_rmsf='Per-chain all-CA fit to first reference frame, then fluctuation around fitted temporal mean; project label implementation',
                  hchi='Joint rotamer Shannon entropy, 120-degree chi bins, natural log; all retained frames; A/G/P excluded only at evaluation')
    np.savez_compressed(cache,rmsf=rmsf,hchi=entropy,aa=aa,chain=chain,author=ids,metadata=json.dumps(metadata))
    print(state,run,len(times),'reference frames; snapshot',snapshot_time,'ps',flush=True)
    return cache

def infer(caches):
    import torch
    from flexppi.inference import _structure_example
    from flexppi.data.dataset import collate_states
    from export_atlas_hchi_v1_oof_predictions import load_fold
    from chunked_inference import forward_chunked
    split=ROOT/'output/atlas_hchi_retraining_v1_linux/payload/data/splits/state_rmsf_hchi_v2_atlas_hchi_v1'
    train=[json.loads(l) for l in (split/'fold_0_train.jsonl').read_text().splitlines()]
    test=[json.loads(l) for l in (split/'fold_0_test.jsonl').read_text().splitlines()]
    assert {'itrace_2bnq','itrace_2bnr'}.issubset({r['sample_id'] for r in test})
    assert not {'itrace_2bnq','itrace_2bnr'} & {r['sample_id'] for r in train}
    torch.set_num_threads(4)
    model,config,path,checkpoint_sha=load_fold(ROOT,0,torch.device('cpu'),True)
    assert model.condition_mlp is None and not config['model']['use_condition_features']
    rows=[];meta=[]
    for cache in caches:
        z=np.load(cache);m=json.loads(str(z['metadata']));snap=m['snapshot'];st=parse_pdb(snap)
        expected_aa=''.join(z['aa']);indices=[i for i,c in enumerate(st.chain_ids) if c in m['target_chains']]
        assert ''.join(st.sequence[i] for i in indices)==expected_aa
        pc=BASE/f"{m['state']}_run{m['run']:02d}_prediction_fold0.npz"
        if pc.exists():
            old=np.load(pc);assert str(old['checkpoint_sha'])==checkpoint_sha and str(old['input_sha'])==sha(snap)
            prediction=old['prediction']
        else:
            example=_structure_example(snap,m['target_chains'],310.,80.,0.)
            with torch.inference_mode():result=forward_chunked(model,collate_states([example]),48)
            prediction=result['mean'][0].numpy()[indices]
            assert result['mask'][0].numpy()[indices].all()
            np.savez_compressed(pc,prediction=prediction,checkpoint_sha=checkpoint_sha,input_sha=sha(snap))
        assert prediction.shape==(444,2) and np.isfinite(prediction).all()
        for i,(ch,author,aa) in enumerate(zip(z['chain'],z['author'],z['aa'])):
            for j,metric in enumerate(['rmsf','hchi']):
                reference=float(z[metric][i]);valid=np.isfinite(reference) and not (metric=='hchi' and aa in 'AGP')
                rows.append(dict(state=m['state'],run=m['run'],chain=ch,author_id=author,amino_acid=aa,metric=metric,reference=reference,
                    prediction=float(prediction[i,j]),evaluation_mask=valid,position_in_chain=i+1 if ch=='alpha' else i-202))
        m['checkpoint_sha256']=checkpoint_sha;m['fold']=0;meta.append(m)
        print('predicted',m['state'],m['run'],flush=True)
    frame=pd.DataFrame(rows);frame.to_csv(OUT/'1g4_run_paired_residues.csv',index=False)
    valid=frame.loc[frame.evaluation_mask]
    agg=valid.groupby(['state','chain','author_id','amino_acid','position_in_chain','metric'],sort=False).agg(
        md_mean=('reference','mean'),md_sd=('reference','std'),prediction_mean=('prediction','mean'),prediction_sd=('prediction','std'),n_runs=('run','size')).reset_index()
    assert (agg.n_runs==10).all();agg.to_csv(OUT/'1g4_state_summary.csv',index=False)
    keys=['chain','author_id','amino_acid','position_in_chain','metric']
    delta=agg[agg.state=='holo'].merge(agg[agg.state=='apo'],on=keys,suffixes=('_holo','_apo'),validate='one_to_one')
    delta['md_delta']=delta.md_mean_holo-delta.md_mean_apo
    delta['prediction_delta']=delta.prediction_mean_holo-delta.prediction_mean_apo
    # Unpaired 10-run state ensembles, not run-number paired biological repeats.
    delta['md_delta_se']=np.sqrt(delta.md_sd_holo**2/10+delta.md_sd_apo**2/10)
    delta['prediction_delta_se']=np.sqrt(delta.prediction_sd_holo**2/10+delta.prediction_sd_apo**2/10)
    delta.to_csv(OUT/'1g4_delta_paired_residues.csv',index=False)
    from scipy.stats import spearmanr
    summaries=[]
    for metric,g in delta.groupby('metric'):
        err=g.prediction_delta-g.md_delta
        summaries.append(dict(task='apo_to_holo_delta',metric=metric,n=len(g),mae=float(abs(err).mean()),rmse=float(np.sqrt((err**2).mean())),
            spearman=float(spearmanr(g.md_delta,g.prediction_delta).statistic),zero_change_mae=float(abs(g.md_delta).mean()),
            sign_agreement=float((np.sign(g.md_delta)==np.sign(g.prediction_delta)).mean())))
        for state,q in agg[agg.metric==metric].groupby('state'):
            e=q.prediction_mean-q.md_mean
            summaries.append(dict(task=state+'_absolute',metric=metric,n=len(q),mae=float(abs(e).mean()),rmse=float(np.sqrt((e**2).mean())),
                spearman=float(spearmanr(q.md_mean,q.prediction_mean).statistic)))
    pd.DataFrame(summaries).to_csv(OUT/'1g4_validation_metrics.csv',index=False)
    audit=dict(status='passed',source_doi='10.1371/journal.pcbi.1007338',archive_doi='10.6084/m9.figshare.8067743.v1',archive_md5=MD5,
        inclusion='All 20 published runs: 10 free TCR + 10 TCR-pMHC; no CD3 or membrane',input='Each run snapshot nearest 10 ns, fixed before comparison',
        reference='All stored 20–100 ns frames per run; matched RMSF implementation and joint Hchi; no new MD',
        scope='MD-derived conformation validation; not experimental-PDB-only deployment, not unseen-family validation',
        split='Frozen fold_0 excludes ITrace 2BNQ/2BNR samples; homologous training systems remain',
        experimental_apo_note='2BNU alpha has a sequence mismatch versus 2BNR/MD; not used as an exact matched apo input; no mutation/reversion performed',
        runs=meta,condition_features_enabled=False,
        runner_sha256=sha(__file__),chunked_inference_sha256=sha(ROOT/'scripts/tcr_real_states/chunked_inference.py'),
        statistics='Means and sample SD over 10 MD runs for each state; SE of difference from independent state variances; no residue-level significance claims')
    (OUT/'1g4_validation_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)

def main():
    BASE.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    warnings.filterwarnings('ignore',module='MDAnalysis')
    dest=extract();pdbs=sorted(dest.rglob('*.pdb'));assert len(pdbs)==20
    caches=[reference_one(p) for p in pdbs];infer(caches)

if __name__=='__main__':main()
