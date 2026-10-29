from pathlib import Path
import sys,os,json,hashlib,types,warnings
from collections import Counter
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='1'
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
SRC=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/src/flexppi')
for name,subdir in [('flexppi',''),('flexppi.data','data')]:
    mod=types.ModuleType(name);mod.__path__=[str(SRC/subdir)];sys.modules[name]=mod
from flexppi.constants import CHI_ATOMS,THREE_TO_ONE
from flexppi.data.labels import _dihedral_degrees,_entropy
import numpy as np
import pandas as pd
import MDAnalysis as mda
warnings.filterwarnings('ignore',message='DCDReader currently makes independent')

FRAME_INDICES=np.arange(7,110)
ALIASES={'HID':'HIS','HIE':'HIS','HIP':'HIS','HSD':'HIS','HSE':'HIS','HSP':'HIS','CYX':'CYS','CYM':'CYS','ASH':'ASP','GLH':'GLU','LYN':'LYS'}
STAGE=ROOT/'canonical_topology';STAGE.mkdir(exist_ok=True)
rows=[];audit=[]

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def canonical(source,dest):
    output=[];counts=Counter();order_source=[];order_dest=[]
    for line in source.read_text().splitlines(True):
        if line.startswith(('ATOM  ','HETATM')):
            original=line
            res=line[17:20]
            if res in ALIASES:
                counts[res+'->'+ALIASES[res]]+=1
                line=line[:17]+ALIASES[res]+line[20:]
            if res=='ILE' and line[12:16].strip()=='CD':
                counts['ILE_CD->CD1']+=1
                line=line[:12]+' CD1'+line[16:]
            # All source coordinates, serials, chain IDs, residue IDs, record order stay exact.
            assert original[:12]==line[:12] and original[20:]==line[20:]
            order_source.append((original[:12],original[20:54]))
            order_dest.append((line[:12],line[20:54]))
        output.append(line)
    assert order_source==order_dest
    dest.write_text(''.join(output),encoding='utf-8')
    return dict(counts),len(order_source)

def rmsf_chain(coords):
    centered=coords.astype(float)-coords.mean(axis=1,keepdims=True)
    ref=coords[0].astype(float)-coords[0].mean(axis=0)
    cov=np.einsum('tni,nj->tij',centered,ref,optimize=True)
    u,_,vt=np.linalg.svd(cov);rot=u@vt
    neg=np.linalg.det(rot)<0;u[neg,:,-1]*=-1;rot[neg]=u[neg]@vt[neg]
    aligned=np.einsum('tni,tij->tnj',centered,rot,optimize=True)
    return np.sqrt(np.mean(np.sum((aligned-aligned.mean(0))**2,axis=-1),axis=0))

for condition in ['WT','S32M']:
    folder=Path('E:/data/CatchBond')/condition/condition
    rawtop=folder/'system_renum.pdb'
    top=STAGE/f'{condition}_all_atoms_canonical.pdb'
    aliases,topatoms=canonical(rawtop,top)
    top_sha=sha(rawtop);prmtop_sha=sha(folder/'system.prmtop')
    for rep in range(1,11):
        path=folder/f'{condition}_reimaged_skip50_r{rep}.dcd'
        traj_sha=sha(path);stat=path.stat()
        u=mda.Universe(str(top),str(path))
        assert len(u.atoms)==topatoms
        residues=[r for r in u.residues if r.resname in THREE_TO_ONE and 'CA' in r.atoms.names]
        assert len(residues)==408
        chain=np.array([r.atoms.chainIDs[0] for r in residues])
        resid=np.array([str(r.resid)+(str(r.icode).strip() if hasattr(r,'icode') else '') for r in residues])
        aa=np.array([THREE_TO_ONE[r.resname] for r in residues])
        ca=np.array([r.atoms[r.atoms.names=='CA'].indices[0] for r in residues])
        quads=[];which={};chi_count=[]
        for i,r in enumerate(residues):
            definitions=CHI_ATOMS.get(r.resname,[]);look={a.name:a.index for a in r.atoms}
            if definitions and all(all(name in look for name in d) for d in definitions):
                which[i]=list(range(len(quads),len(quads)+len(definitions)))
                quads.extend([[look[name] for name in d] for d in definitions])
            chi_count.append(len(definitions))
        quads=np.array(quads,dtype=int)
        assert all(i in which for i,r in enumerate(residues) if r.resname=='PRO')
        assert len(u.trajectory)>int(FRAME_INDICES[-1])
        coords=[];angles=[];times=[]
        for frame in FRAME_INDICES:
            ts=u.trajectory[int(frame)]
            coords.append(ts.positions[ca].copy())
            angles.append(_dihedral_degrees(ts.positions[quads]))
            times.append(float(ts.time))
        totalframes=len(u.trajectory);dt=float(u.trajectory.dt)
        u.trajectory.close()
        coords=np.asarray(coords);angles=np.asarray(angles);times=np.asarray(times)
        bins=np.clip(np.floor((angles+180)/120).astype(np.int8),0,2)
        hchi=np.full(408,np.nan);values=np.full(408,np.nan)
        for ch in np.unique(chain):
            selected=np.flatnonzero(chain==ch)
            values[selected]=rmsf_chain(coords[:,selected])
        for i,columns in which.items():
            codes=sum(bins[:,col].astype(np.int32)*3**j for j,col in enumerate(columns))
            hchi[i]=_entropy(codes)
        assert np.all(np.isfinite(values))
        assert not np.any(np.isfinite(hchi[np.isin(aa,['A','G'])]))
        cache=ROOT/f'{condition}_rep{rep:02d}_all408.npz'
        np.savez_compressed(cache,chain=chain,author_residue_id=resid,amino_acid=aa,
            rmsf_A=values,Hchi_nats=hchi,frame_indices_zero_based=FRAME_INDICES,
            dcd_times_ps=times,readme_nominal_time_ns=FRAME_INDICES*10,
            ca_coordinates_A=coords,chi_angles_degrees=angles,chi_columns_json=json.dumps(which),
            chi_atom_indices_original_topology=quads,CA_atom_indices_original_topology=ca,
            source_trajectory_sha256=traj_sha)
        for i in range(408):
            rows.append({'condition':condition,'replica':rep,'chain':chain[i],'author_residue_id':resid[i],
                'aa':aa[i],'rmsf_A':values[i],'Hchi_nats':hchi[i],'n_chi':chi_count[i],
                'frames':len(FRAME_INDICES),'index_start_inclusive':7,'index_stop_exclusive':110})
        assert path.stat().st_size==stat.st_size and path.stat().st_mtime_ns==stat.st_mtime_ns
        audit.append({'condition':condition,'replica':rep,'original_topology':str(rawtop),'topology_sha256':top_sha,
            'prmtop':str(folder/'system.prmtop'),'prmtop_sha256':prmtop_sha,'canonical_topology':str(top),
            'canonical_topology_sha256':sha(top),'alias_counts':aliases,'topology_atoms':topatoms,
            'trajectory':str(path),'trajectory_sha256':traj_sha,'source_bytes':stat.st_size,'saved_frames':totalframes,
            'analysis_frames':len(FRAME_INDICES),'DCD_dt_ps':dt,'first_selected_DCD_time_ps':times[0],
            'last_selected_DCD_time_ps':times[-1],'protein_residues':408,'Hchi_valid_residues':int(np.isfinite(hchi).sum()),
            'output':str(cache),'source_stat_unchanged_after_read':True})
        print(condition,rep,'all408 RMSF/Hchi complete','frames',len(coords),'dt_ps',dt,flush=True)

df=pd.DataFrame(rows)
df.to_csv(ROOT/'all408_per_replica_labels.csv',index=False)
summary=df.groupby(['condition','chain','author_residue_id','aa'],sort=False).agg(
    rmsf_A_mean=('rmsf_A','mean'),rmsf_A_SD=('rmsf_A','std'),rmsf_A_n=('rmsf_A','count'),
    Hchi_nats_mean=('Hchi_nats','mean'),Hchi_nats_SD=('Hchi_nats','std'),Hchi_nats_n=('Hchi_nats','count')).reset_index()
summary['rmsf_A_SEM']=summary.rmsf_A_SD/np.sqrt(summary.rmsf_A_n)
summary['Hchi_nats_SEM']=summary.Hchi_nats_SD/np.sqrt(summary.Hchi_nats_n.replace(0,np.nan))
summary.to_csv(ROOT/'all408_ten_replica_mean_SD_SEM.csv',index=False)
focus=summary[(summary.chain=='C')|((summary.chain=='D')&summary.author_residue_id.isin([str(i) for i in range(25,36)]+['92','93']))|((summary.chain=='E')&summary.author_residue_id.isin(['96','97','98','99','100']))]
focus.to_csv(ROOT/'peptide_CDR1alpha_N92_N93_beta98_focus.csv',index=False)
(ROOT/'MD_recomputation_audit.json').write_text(json.dumps({'status':'passed','new_MD':False,'new_training':False,
    'protocol':'Each of 20 independent DCD trajectories uses zero-based frames7:110, 103frames; each protein chain independently C-alpha Kabsch fit to first selected frame. Joint chi3-bin entropy computed per residue from all existing chi; Pro2chi, A/G missing.',
    'units':{'RMSF':'Angstrom','Hchi':'nats','chi':'degrees'},
    'readme_sampling':'10ns per stored frame; exclude first7frames; time indices stored separately from DCD time metadata.',
    'readme_parameter_note':'Source README discloses incorrect Lennard-Jones radius for water oxygen; authors state relative comparisons and principal findings are unaffected. This is an external case study, not a benchmark standard.',
    'sources':audit},ensure_ascii=False,indent=2),encoding='utf-8')
print(focus.to_string(index=False),flush=True)
