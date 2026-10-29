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
SOURCE=ROOT.parent/'tcr_case_v40_20261017'
OUT=ROOT/'jm22_f5l_validation';OUT.mkdir(exist_ok=True)
STAGE=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010/tcr_v41_jm22_staged');STAGE.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=[json.loads(l) for l in (SOURCE/'family_predictions/prediction_manifest.jsonl').read_text().splitlines()]
provenance=[];replica_results=[];contact_rows=[];physical_rows=[];state_rows=[];predictions={};frames_audit=[]
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

for pdb,label,expected in [('1oga','JM22-WT-M1','GILGFVFTL'),('5hhm','JM22-F5L-M1','GILGLVFTL')]:
    prediction_path=SOURCE/f'family_predictions/full/itrace_{pdb}.npz'
    with np.load(prediction_path,allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
    mi=next(x for x in manifest if x['sample_id']=='itrace_'+pdb and x['variant']=='full')
    model_pdb=CORRECTED/f'training_payload/dataset/structures/itrace_{pdb}.pdb'
    assert sha(model_pdb)==mi['input_sha256']
    shutil.copy2(model_pdb,OUT/f'{pdb}_model_input.pdb')
    cids=np.flatnonzero(pred['chain']=='C');assert len(cids)==9
    assert ''.join(pred['amino_acid'][cids])==expected
    predictions[pdb]={k:pred[k][cids] for k in ['mean','target','label_variance','valid','amino_acid']}
    for rep in (1,2,3):
        rawtop=Path(f'E:/data/ITrace/trajectories/{pdb}_run{rep}_topology.pdb')
        rawxtc=Path(f'E:/data/ITrace/trajectories/{pdb}_run{rep}_md_processed.xtc')
        assert rawtop.exists() and rawxtc.exists()
        original_top_sha,original_xtc_sha=sha(rawtop),sha(rawxtc)
        top=STAGE/f'{pdb}_{rep}_canonical.pdb';xtc=STAGE/f'{pdb}_{rep}.xtc'
        canonical(rawtop,top)
        if not xtc.exists() or xtc.stat().st_size!=rawxtc.stat().st_size:shutil.copy2(rawxtc,xtc)
        assert sha(xtc)==original_xtc_sha
        provenance.append({'pdb':pdb,'replica':rep,'topology':str(rawtop),'topology_sha256':original_top_sha,'trajectory':str(rawxtc),'trajectory_sha256':original_xtc_sha,'canonical_topology_sha256':sha(top),'prediction_sha256':sha(prediction_path),'model_input_sha256':sha(model_pdb),'checkpoint_sha256':mi['checkpoint_sha256']})
        u=mda.Universe(str(top),str(xtc))
        pep=[r for r in u.residues if r.atoms.chainIDs[0]=='C' and 'CA' in r.atoms.names]
        assert ''.join(THREE_TO_ONE[r.resname] for r in pep)==expected
        ca=np.array([r.atoms[r.atoms.names=='CA'].indices[0] for r in pep])
        heavy=lambda atom:not atom.name.lstrip('0123456789').startswith(('H','D'))
        ph=[np.array([a.index for a in r.atoms if heavy(a)]) for r in pep]
        flat_ph=np.concatenate(ph);starts=np.r_[0,np.cumsum([len(x) for x in ph])[:-1]]
        flat_index={int(atom):i for i,atom in enumerate(flat_ph)}
        side_cols=[[flat_index[a.index] for a in r.atoms if heavy(a) and a.name not in ('N','CA','C','O','OXT')] for r in pep]
        tcr_indices=np.array([a.index for a in u.atoms if a.chainID in ('D','E') and heavy(a)])
        quads=[];which={}
        for i,r in enumerate(pep):
            look={a.name:a.index for a in r.atoms};definitions=CHI_ATOMS.get(r.resname,[])
            if definitions and all(all(n in look for n in d) for d in definitions):
                which[i]=list(range(len(quads),len(quads)+len(definitions)))
                quads.extend([[look[n] for n in d] for d in definitions])
        quads=np.array(quads,dtype=int);coords=[];angles=[];contact=[];side_contact=[];times=[]
        for ts in u.trajectory:
            xyz=ts.positions
            coords.append(xyz[ca].copy());angles.append(_dihedral_degrees(xyz[quads]));times.append(float(ts.time))
            dist=cKDTree(xyz[tcr_indices]).query(xyz[flat_ph],k=1)[0]
            contact.append(np.minimum.reduceat(dist,starts))
            side_contact.append([float(dist[c].min()) if len(c) else np.nan for c in side_cols])
        u.trajectory.close()
        times=np.asarray(times);coords=np.asarray(coords);angles=np.asarray(angles);contact=np.asarray(contact);side_contact=np.asarray(side_contact)
        assert len(times)>=1001 and abs(times[0])<1e-6 and abs(times[-1]-200000)<1e-3
        frames_audit.append({'pdb':pdb,'replica':rep,'frames':len(times),'first_ps':float(times[0]),'last_ps':float(times[-1])})
        values=rmsf_chain(coords)
        statebins=np.clip(np.floor((angles+180)/120).astype(np.int8),0,2)
        H=np.full(9,np.nan);fine_mean=np.full(9,np.nan);width_mean=np.full(9,np.nan)
        for i,columns in which.items():
            codes=sum(statebins[:,col].astype(np.int32)*3**j for j,col in enumerate(columns));H[i]=_entropy(codes)
            hist=np.bincount(codes,minlength=3**len(columns))
            for code,count in enumerate(hist):state_rows.append({'pdb':pdb,'complex':label,'replica':rep,'position':i+1,'aa':expected[i],'nchi':len(columns),'state':code,'count':int(count),'probability':count/len(codes)})
            fine=[];widths=[]
            for j,col in enumerate(columns):
                ent,width=fine_metrics(angles[:,col]);fine.append(ent);widths.append(width)
                physical_rows.append({'pdb':pdb,'complex':label,'replica':rep,'position':i+1,'aa':expected[i],'chi':j+1,'fine_entropy_nats':ent,'lowF_width_degrees':width,'frames':len(times)})
            fine_mean[i]=np.mean(fine);width_mean[i]=np.mean(widths)
        np.savez_compressed(OUT/f'{pdb}_rep{rep}_peptide.npz',time_ps=times,ca_coordinates=coords,chi_angles_degrees=angles,chi_columns_json=json.dumps(which),minimum_TCR_heavy_distance_A=contact,minimum_sidechain_TCR_heavy_distance_A=side_contact)
        for i in range(9):
            replica_results.append({'pdb':pdb,'complex':label,'replica':rep,'position':i+1,'aa':expected[i],'nchi':len(which.get(i,[])),'rmsf_A':values[i],'Hchi_nats':H[i],'fine_entropy_mean_nats_per_chi':fine_mean[i],'lowF_width_mean_degrees_per_chi':width_mean[i],'frames':len(times)})
            for cut in [4.0,4.5,5.0]:
                contact_rows.append({'pdb':pdb,'complex':label,'replica':rep,'position':i+1,'aa':expected[i],'cutoff_A':cut,'any_heavy_contact_fraction':float(np.mean(contact[:,i]<=cut)),'sidechain_contact_fraction':float(np.mean(side_contact[:,i]<=cut)) if len(side_cols[i]) else np.nan,'minimum_TCR_distance_median_A':float(np.median(contact[:,i])),'minimum_sidechain_TCR_distance_median_A':float(np.median(side_contact[:,i])) if len(side_cols[i]) else np.nan})
        assert sha(rawtop)==original_top_sha and sha(rawxtc)==original_xtc_sha
        print(pdb,rep,'RMSF',np.round(values,4),'Hchi',np.round(H,4),flush=True)
    local=pd.DataFrame(replica_results);local=local[local.pdb==pdb]
    for i,j in enumerate(cids):
        r=local[local.position==i+1]
        np.testing.assert_allclose(r.rmsf_A.mean(),pred['target'][j,0],rtol=0,atol=3e-6)
        np.testing.assert_allclose(r.rmsf_A.var(ddof=1),pred['label_variance'][j,0],rtol=0,atol=3e-6)
        if pred['valid'][j,1]:
            np.testing.assert_allclose(r.Hchi_nats.mean(),pred['target'][j,1],rtol=0,atol=3e-6)
            np.testing.assert_allclose(r.Hchi_nats.var(ddof=1),pred['label_variance'][j,1],rtol=0,atol=3e-6)

reps=pd.DataFrame(replica_results);reps.to_csv(OUT/'all9_peptide_per_replica_labels.csv',index=False)
pd.DataFrame(contact_rows).to_csv(OUT/'all9_peptide_TCR_contacts_per_replica.csv',index=False)
pd.DataFrame(physical_rows).to_csv(OUT/'all9_per_chi_fine_entropy_lowF_width.csv',index=False)
pd.DataFrame(state_rows).to_csv(OUT/'all9_joint_rotamer_state_occupancies.csv',index=False)
delta=[]
for i in range(9):
    row={'position':i+1,'aa_a':predictions['1oga']['amino_acid'][i],'aa_b':predictions['5hhm']['amino_acid'][i],'same_aa':predictions['1oga']['amino_acid'][i]==predictions['5hhm']['amino_acid'][i]}
    for metric,j in [('rmsf',0),('hchi',1)]:
        for pdb,side in [('1oga','a'),('5hhm','b')]:
            valid=bool(predictions[pdb]['valid'][i,j]);row[metric+'_valid_'+side]=valid
            row[metric+'_MD_'+side]=float(predictions[pdb]['target'][i,j]) if valid else np.nan
            row[metric+'_prediction_'+side]=float(predictions[pdb]['mean'][i,j]) if valid else np.nan
            row[metric+'_MD_SD_'+side]=float(np.sqrt(predictions[pdb]['label_variance'][i,j])) if valid else np.nan
        row[metric+'_MD_delta']=row[metric+'_MD_b']-row[metric+'_MD_a'];row[metric+'_prediction_delta']=row[metric+'_prediction_b']-row[metric+'_prediction_a']
        row[metric+'_delta_SEM']=np.sqrt((row[metric+'_MD_SD_a']**2+row[metric+'_MD_SD_b']**2)/3)
    for metric in ['fine_entropy_mean_nats_per_chi','lowF_width_mean_degrees_per_chi']:
        for pdb,side in [('1oga','a'),('5hhm','b')]:row[metric+'_'+side]=float(reps[(reps.pdb==pdb)&(reps.position==i+1)][metric].mean())
        row[metric+'_delta']=row[metric+'_b']-row[metric+'_a']
    delta.append(row)
df=pd.DataFrame(delta);df.to_csv(OUT/'all9_model_MD_physical_changes.csv',index=False)
df[df.same_aa].to_csv(OUT/'common8_amino_acid_changes.csv',index=False)
df[~df.same_aa].to_csv(OUT/'F5L_substitution_separate.csv',index=False)
summary=[]
for scope,mask in [('all9',np.ones(9,dtype=bool)),('common8',df.same_aa.to_numpy()),('F5L',~df.same_aa.to_numpy())]:
    for metric in ['rmsf','hchi']:
        f=df[mask&df[metric+'_valid_a']&df[metric+'_valid_b']];x=f[metric+'_MD_delta'].to_numpy();y=f[metric+'_prediction_delta'].to_numpy();sem=f[metric+'_delta_SEM'].to_numpy();robust=abs(x)>2*sem
        rho=float(np.corrcoef(pd.Series(x).rank(),pd.Series(y).rank())[0,1]) if len(f)>=4 else np.nan
        summary.append({'scope':scope,'metric':metric,'n':len(f),'rho_delta':rho,'MAE_delta':float(np.mean(abs(x-y))),'zero_change_MAE':float(np.mean(abs(x))),'relative_MAE_gain':float(1-np.mean(abs(x-y))/np.mean(abs(x))),'robust_n':int(robust.sum()),'robust_direction_fraction':float(np.mean(x[robust]*y[robust]>0)) if robust.any() else np.nan})
pd.DataFrame(summary).to_csv(OUT/'model_MD_change_summary.csv',index=False)
(OUT/'raw_trajectory_audit.json').write_text(json.dumps({'status':'passed','all9_peptide_original_reference_means_and_variances_recomputed':True,'absolute_tolerance':3e-6,'raw_topology_and_trajectory_files_unchanged':True,'frame_window':'all available frames 0–200ns; no subsampling','frames':frames_audit,'alignment':'peptide C-alpha independently Kabsch-fitted to first frame for RMSF; contacts and chi are alignment invariant','contact_definition':'Any peptide heavy atom, or sidechain heavy atom excluding N/CA/C/O/OXT, within stated cutoff of any TCR heavy atom; glycine sidechain contact undefined','fine_angle_definition':'Each chi 5 degree periodic bins, entropy and F*<=1 occupancy width averaged over four origins 0/1.25/2.5/3.75 degrees','comparison':'5HHM F5L minus1OGA WT. Substitution site5 analyzed separately. Common8 means identical amino acid on8 positions (only6 valid Hchi).','source_files':provenance,'new_MD':False,'new_training':False},ensure_ascii=False,indent=2),encoding='utf8')
print(pd.DataFrame(summary).to_string(index=False),flush=True)
