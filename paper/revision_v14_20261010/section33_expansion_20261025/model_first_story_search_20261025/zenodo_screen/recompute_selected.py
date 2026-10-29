"""Read original XTC sequentially; canonical names; new outputs only."""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
import sys,json,hashlib,argparse,importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.formats.libmdaxdr import XTCFile
CORE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
OUT=Path(__file__).resolve().parent
DATA=Path(r'E:\Project\FlexPPI\artifacts\zenodo_6638504\data')
sys.path.insert(0,str(CORE/'training_payload/src'))
from flexppi.constants import CHI_ATOMS,THREE_TO_ONE
from flexppi.data.structures import parse_pdb
from flexppi.data.labels import _dihedral_degrees,_entropy
spec=importlib.util.spec_from_file_location('correct',CORE/'prepare_corrected_labels.py')
correct=importlib.util.module_from_spec(spec);spec.loader.exec_module(correct)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(pids):
    records=[json.loads(x) for x in (DATA/'zenodo_6638504_states.jsonl').read_text().splitlines()]
    out=OUT/'raw_recomputed';out.mkdir(exist_ok=True)
    for r in records:
        if r['metadata']['pdb_id'] not in pids:continue
        sid=r['sample_id'];sub=out/sid;sub.mkdir(exist_ok=True)
        pdb=sub/'canonical_topology.pdb';changes=correct.canonical_pdb(r['structure_path'],pdb)
        s=parse_pdb(str(pdb));u=mda.Universe(str(pdb));n=s.num_residues
        residues=[q for q in u.residues if q.resname in THREE_TO_ONE and 'CA' in q.atoms.names]
        assert [(str(q.atoms.chainIDs[0]).strip() or '_',str(q.resid)+str(q.icode).strip(),q.resname) for q in residues]==list(zip(s.chain_ids,s.author_residue_ids,s.residue_names))
        quads=[];cols={};ca=[]
        for i,q in enumerate(residues):
            atoms={a.name:a.index for a in q.atoms};ca.append(atoms['CA']);defs=CHI_ATOMS.get(q.resname,[])
            if defs and all(all(a in atoms for a in quad) for quad in defs):
                cols[i]=list(range(len(quads),len(quads)+len(defs)));quads.extend([[atoms[a] for a in quad] for quad in defs])
        quads=np.array(quads);angles=[];coords=[];times=[];bmin=np.inf;bmax=0.
        xtc=Path(r['trajectory_paths'][0])
        with XTCFile(str(xtc),'r') as stream:
            for frame in stream:
                x=np.asarray(frame.x,dtype=np.float32)*10
                assert len(x)==len(u.atoms)
                q=x[quads];bonds=np.linalg.norm(np.diff(q,axis=1),axis=-1)
                bmin=min(bmin,float(bonds.min()));bmax=max(bmax,float(bonds.max()))
                if bmin<.95 or bmax>2.25:raise ValueError((sid,'Broken covalent geometry',bmin,bmax))
                angles.append(_dihedral_degrees(q));coords.append(x[ca]);times.append(float(frame.time))
        times=np.array(times);angles=np.array(angles,dtype=np.float32);coords=np.array(coords,dtype=np.float32)
        dt=np.diff(times);assert np.allclose(dt,dt[0],rtol=1e-4,atol=.02)
        wf=int(round(100000/dt[0]));nw=len(times)//wf;used=nw*wf
        assert nw>=3
        blocks=np.full((nw,n,2),np.nan)
        for b in range(nw):
            seg=slice(b*wf,(b+1)*wf);blocks[b,:,0]=correct.chain_fit(coords[seg],s.chain_ids)
            for i,cs in cols.items():
                bins=np.clip(np.floor((angles[seg][:,cs]+180)/120).astype(int),0,2)
                codes=(bins*3**np.arange(len(cs))).sum(1);blocks[b,i,1]=_entropy(codes)
        mean=np.nanmean(blocks,axis=0);mask=np.isfinite(mean)
        np.savez_compressed(sub/'labels_and_windows.npz',mean=mean,mask=mask,blocks=blocks,chain=np.array(s.chain_ids),author=np.array(s.author_residue_ids),aa=np.array(list(s.sequence)),time_ps=times[:used])
        np.savez_compressed(sub/'raw_chi_angles.npz',angles=angles[:used],residue=np.array(list(cols)),chi_start=np.array([c[0] for c in cols.values()]),chi_count=np.array([len(c) for c in cols.values()]))
        old=np.load(r['label_path']);assert old['values'].shape[0]==n
        discrep={m:float(np.max(np.abs(mean[mask[:,j],j]-old['values'][mask[:,j],j]))) for j,m in enumerate(['rmsf','hchi'])}
        audit=dict(status='complete',sample_id=sid,chemical_changes=changes,n_residues=n,frames=len(times),dt_ps=float(dt[0]),duration_ns=float((times[-1]-times[0])/1000),
            windows=nw,window_frames=wf,frames_used=used,chi_bond_min_A=bmin,chi_bond_max_A=bmax,legacy_label_max_difference=discrep,
            source=[dict(path=str(p),sha256=sha(p)) for p in [Path(r['structure_path']),xtc]],script_sha256=sha(__file__),
            interpretation='100ns contiguous blocks within one trajectory, not independent repeats; per-chain all-CA fitting; no silent PBC repair.')
        (sub/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
        print(json.dumps(audit),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('pids',nargs='+');args=ap.parse_args();main(args.pids)
