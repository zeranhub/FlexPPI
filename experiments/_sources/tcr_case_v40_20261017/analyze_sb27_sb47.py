from pathlib import Path
import sys,os,json,hashlib,shutil,types
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
ROOT=Path(__file__).resolve().parent
RUNTIME=ROOT.parent/'label_sensitivity/runtime'
sys.path.insert(0,str(RUNTIME))
CORRECTED=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
for name,subdir in [('flexppi',''),('flexppi.data','data')]:
    module=types.ModuleType(name);module.__path__=[str(CORRECTED/'training_payload/src/flexppi'/subdir)];sys.modules[name]=module
from flexppi.constants import CHI_ATOMS,THREE_TO_ONE
from flexppi.data.labels import _dihedral_degrees,_entropy
import numpy as np
import pandas as pd
import MDAnalysis as mda
from scipy.spatial import cKDTree
OUT=ROOT/'sb27_sb47';OUT.mkdir(exist_ok=True)
STAGE=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010/tcr_v40_staged');STAGE.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=[json.loads(l) for l in (ROOT/'family_predictions/prediction_manifest.jsonl').read_text().splitlines()]
provenance=[];replica_results=[];timeseries=[];contact_rows=[];q7angles=[];q7states=[];outputs=[]
aliases={'HID':'HIS','HIE':'HIS','HIP':'HIS','HSD':'HIS','HSE':'HIS','HSP':'HIS'}

def canonical(source,dest):
    result=[]
    for line in Path(source).read_text().splitlines(True):
        if line.startswith(('ATOM  ','HETATM')):
            res=line[17:20]
            if res in aliases:line=line[:17]+aliases[res]+line[20:]
            if res=='ILE' and line[12:16].strip()=='CD':line=line[:12]+' CD1'+line[16:]
        result.append(line)
    dest.write_text(''.join(result),encoding='utf8')

def rmsf_chain(coords):
    centered=coords.astype(float)-coords.mean(axis=1,keepdims=True)
    ref=coords[0].astype(float)-coords[0].mean(axis=0)
    cov=np.einsum('tni,nj->tij',centered,ref,optimize=True)
    u,_,vt=np.linalg.svd(cov);rot=u@vt
    neg=np.linalg.det(rot)<0;u[neg,:,-1]*=-1;rot[neg]=u[neg]@vt[neg]
    aligned=np.einsum('tni,tij->tnj',centered,rot,optimize=True)
    return np.sqrt(np.mean(np.sum((aligned-aligned.mean(0))**2,axis=-1),axis=0))

def fine_metrics(ang):
    ent=[];width=[]
    for origin in [0,1.25,2.5,3.75]:
        hist=np.histogram((ang+180-origin)%360,bins=np.arange(0,365,5))[0]
        p=hist[hist>0]/hist.sum();ent.append(-np.sum(p*np.log(p)))
        width.append(np.sum(hist>=hist.max()/np.e)*5)
    return float(np.mean(ent)),float(np.mean(width))

for pdb,label in [('2ak4','SB27'),('4jry','SB47')]:
    prediction_path=ROOT/f'family_predictions/full/itrace_{pdb}.npz'
    pred=np.load(prediction_path,allow_pickle=False)
    mi=next(x for x in manifest if x['sample_id']=='itrace_'+pdb and x['variant']=='full')
    model_pdb=CORRECTED/f'training_payload/dataset/structures/itrace_{pdb}.pdb'
    assert sha(model_pdb)==mi['input_sha256']
    shutil.copy2(model_pdb,OUT/f'{pdb}_model_input.pdb')
    cids=np.flatnonzero(pred['chain']=='C');assert len(cids)==13
    expected=''.join(pred['amino_acid'][cids]);assert expected=='LPEPLPQGQLTAY'
    for rep in (1,2,3):
        rawtop=Path(f'E:/data/ITrace/trajectories/{pdb}_run{rep}_topology.pdb')
        rawxtc=Path(f'E:/data/ITrace/trajectories/{pdb}_run{rep}_md_processed.xtc')
        assert rawtop.exists() and rawxtc.exists()
        top=STAGE/f'{pdb}_{rep}_canonical.pdb';xtc=STAGE/f'{pdb}_{rep}.xtc'
        canonical(rawtop,top)
        if not xtc.exists() or xtc.stat().st_size!=rawxtc.stat().st_size:shutil.copy2(rawxtc,xtc)
        assert sha(xtc)==sha(rawxtc)
        provenance.append({'pdb':pdb,'replica':rep,'topology':str(rawtop),'topology_sha256':sha(rawtop),'trajectory':str(rawxtc),'trajectory_sha256':sha(rawxtc),'canonical_topology_sha256':sha(top),'prediction_sha256':sha(prediction_path),'model_input_sha256':sha(model_pdb),'checkpoint_sha256':mi['checkpoint_sha256']})
        u=mda.Universe(str(top),str(xtc))
        pep=[r for r in u.residues if r.atoms.chainIDs[0]=='C' and 'CA' in r.atoms.names]
        assert ''.join(THREE_TO_ONE[r.resname] for r in pep)==expected
        ca=np.array([r.atoms[r.atoms.names=='CA'].indices[0] for r in pep])
        heavy=lambda atom:not atom.name.lstrip('0123456789').startswith(('H','D'))
        ph=[np.array([a.index for a in r.atoms if heavy(a)]) for r in pep]
        flat_ph=np.concatenate(ph);starts=np.r_[0,np.cumsum([len(x) for x in ph])[:-1]]
        tcr_indices=np.array([a.index for a in u.atoms if a.chainID in ('D','E') and heavy(a)])
        q7_indices=np.array([a.index for a in pep[6].atoms if heavy(a) and a.name not in ('N','CA','C','O','OXT')])
        quads=[];which={}
        for i,r in enumerate(pep):
            look={a.name:a.index for a in r.atoms}
            definitions=CHI_ATOMS.get(r.resname,[])
            if definitions and all(all(n in look for n in d) for d in definitions):
                which[i]=list(range(len(quads),len(quads)+len(definitions)))
                quads.extend([[look[n] for n in d] for d in definitions])
        q7_atom_names=np.array([u.atoms[i].name for i in q7_indices])
        quads=np.array(quads,dtype=int);coords=[];angles=[];contact=[];side_dist=[];side_atom_dist=[];times=[]
        for ts in u.trajectory:
            xyz=ts.positions
            coords.append(xyz[ca].copy());angles.append(_dihedral_degrees(xyz[quads]));times.append(float(ts.time))
            tree=cKDTree(xyz[tcr_indices]);dist=tree.query(xyz[flat_ph],k=1)[0]
            contact.append(np.minimum.reduceat(dist,starts))
            per_atom=tree.query(xyz[q7_indices],k=1)[0]
            side_dist.append(float(per_atom.min()));side_atom_dist.append(per_atom)
        u.trajectory.close()
        times=np.asarray(times);coords=np.asarray(coords);angles=np.asarray(angles);contact=np.asarray(contact);side_dist=np.asarray(side_dist);side_atom_dist=np.asarray(side_atom_dist)
        assert len(times)>=1001 and abs(times[0])<1e-6 and abs(times[-1]-200000)<1e-3
        values=rmsf_chain(coords)
        statebins=np.clip(np.floor((angles+180)/120).astype(np.int8),0,2)
        H=np.full(13,np.nan)
        for i,columns in which.items():
            codes=sum(statebins[:,col].astype(np.int32)*3**j for j,col in enumerate(columns));H[i]=_entropy(codes)
        qang=angles[:,which[6]];qcode=sum(statebins[:,col].astype(np.int32)*3**j for j,col in enumerate(which[6]))
        q7angles.append(qang);q7states.append(qcode)
        np.savez_compressed(OUT/f'{pdb}_rep{rep}_peptide.npz',time_ps=times,ca_coordinates=coords,chi_angles_degrees=angles,chi_columns_json=json.dumps(which),minimum_TCR_heavy_distance_A=contact,q7_minimum_sidechain_TCR_distance_A=side_dist,q7_sidechain_atom_names=q7_atom_names,q7_sidechain_atom_TCR_distance_A=side_atom_dist,q7_chi_degrees=qang,q7_joint_state=qcode)
        hist=np.bincount(qcode,minlength=27)
        for code,n in enumerate(hist):outputs.append({'state':code,'pdb':pdb,'complex':label,'replica':rep,'count':int(n),'probability':n/len(qcode)})
        for i in range(13):
            replica_results.append({'pdb':pdb,'complex':label,'replica':rep,'position':i+1,'aa':expected[i],'rmsf_A':values[i],'Hchi_nats':H[i]})
            for cut in [4.0,4.5,5.0]:contact_rows.append({'pdb':pdb,'complex':label,'replica':rep,'position':i+1,'aa':expected[i],'cutoff_A':cut,'contact_fraction':float(np.mean(contact[:,i]<=cut))})
        ent,width=zip(*(fine_metrics(qang[:,j]) for j in range(3)))
        terminal=side_atom_dist[:,np.isin(q7_atom_names,['OE1','NE2'])].min(axis=1)
        timeseries.append({'pdb':pdb,'complex':label,'replica':rep,'frames':len(times),'duration_ns':(times[-1]-times[0])/1000,'Q7_Hchi_nats':H[6],'Q7_rmsf_A':values[6],'Q7_side_TCR_contact_fraction_4p5A':float(np.mean(side_dist<=4.5)),'Q7_terminal_TCR_contact_fraction_4p5A':float(np.mean(terminal<=4.5)),'Q7_side_TCR_distance_median_A':float(np.median(side_dist)),'Q7_terminal_TCR_distance_median_A':float(np.median(terminal)),'Q7_fine_entropy_mean_nats_per_chi':float(np.mean(ent)),'Q7_lowF_width_mean_deg_per_chi':float(np.mean(width))})
        print(pdb,rep,'Q7 H',round(H[6],6),'RMSF',round(values[6],6),'side TCR contact',round(np.mean(side_dist<=4.5),5),'lowFwidth',np.mean(width),flush=True)
    local=pd.DataFrame(replica_results);local=local[local.pdb==pdb]
    for i,j in enumerate(cids):
        r=local[local.position==i+1]
        np.testing.assert_allclose(r.rmsf_A.mean(),pred['target'][j,0],rtol=0,atol=3e-6)
        np.testing.assert_allclose(r.rmsf_A.var(ddof=1),pred['label_variance'][j,0],rtol=0,atol=3e-6)
        if pred['valid'][j,1]:
            np.testing.assert_allclose(r.Hchi_nats.mean(),pred['target'][j,1],rtol=0,atol=3e-6)
            np.testing.assert_allclose(r.Hchi_nats.var(ddof=1),pred['label_variance'][j,1],rtol=0,atol=3e-6)
pd.DataFrame(replica_results).to_csv(OUT/'all13_peptide_per_replica_labels.csv',index=False)
pd.DataFrame(contact_rows).to_csv(OUT/'all13_peptide_TCR_contact_fraction.csv',index=False)
pd.DataFrame(timeseries).to_csv(OUT/'Q7_per_replica_physical_summary.csv',index=False)
pd.DataFrame(outputs).to_csv(OUT/'Q7_all27_joint_state_occupancies.csv',index=False)
(OUT/'raw_trajectory_audit.json').write_text(json.dumps({'status':'passed','all13_peptide_original_reference_means_and_variances_recomputed':True,'absolute_tolerance':3e-6,'raw_topology_and_trajectory_files_unchanged':True,'frame_window':'all available1001–1003frames 0–200ns, nominal200ps; all source saved frames retained','alignment':'peptide C-alpha independently Kabsch-fitted to first frame for RMSF; contacts and chi are alignment invariant','contact_definition':'At least one heavy atom <=cutoff; Q7 sidechain excludes N,CA,C,O,OXT; terminal contact uses OE1 or NE2','source_files':provenance,'new_MD':False,'new_training':False},ensure_ascii=False,indent=2),encoding='utf8')
print(pd.DataFrame(timeseries).to_string(index=False),flush=True)
