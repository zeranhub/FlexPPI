from pathlib import Path
import sys,os,json,hashlib,shutil,types
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'label_sensitivity/runtime'))
CORRECTED=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
for name,subdir in [('flexppi',''),('flexppi.data','data')]:
    module=types.ModuleType(name);module.__path__=[str(CORRECTED/'training_payload/src/flexppi'/subdir)];sys.modules[name]=module
from flexppi.constants import CHI_ATOMS,THREE_TO_ONE
from flexppi.data.labels import _dihedral_degrees,_entropy
import numpy as np
import pandas as pd
import MDAnalysis as mda
SOURCE=ROOT.parent/'tcr_case_v40_20261017'
STAGE=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010/tcr_v41_allprotein_staged');STAGE.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=[json.loads(l) for l in (SOURCE/'family_predictions/prediction_manifest.jsonl').read_text().splitlines()]
roles=json.loads(Path(r'E:\Project\FlexPPI\output\itrace_benchmark\raw\itrace_chain_roles.json').read_text())
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

def chain_fit(coordinates,chain_ids):
    # Exact double-precision implementation used in prepare_corrected_labels.py,
    # which generated the ITrace targets (the generic labels module uses float32).
    values=np.full(coordinates.shape[1],np.nan,dtype=np.float64);chain_ids=np.asarray(chain_ids)
    for chain in dict.fromkeys(chain_ids.tolist()):
        ids=np.flatnonzero(chain_ids==chain)
        if len(ids)<3:continue
        mobile=coordinates[:,ids].astype(np.float64);centered=mobile-mobile.mean(axis=1,keepdims=True);ref=mobile[0]-mobile[0].mean(axis=0)
        covariance=np.einsum('tni,nj->tij',centered,ref,optimize=True);u,_,vt=np.linalg.svd(covariance);rotation=u@vt;neg=np.linalg.det(rotation)<0;u[neg,:,-1]*=-1;rotation[neg]=u[neg]@vt[neg]
        aligned=np.einsum('tni,tij->tnj',centered,rotation,optimize=True);values[ids]=np.sqrt(np.mean(np.sum((aligned-aligned.mean(0))**2,axis=-1),axis=0))
    return values

def complex_labels(pdb,out):
    prediction_path=SOURCE/f'family_predictions/full/itrace_{pdb}.npz'
    with np.load(prediction_path,allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
    mi=next(x for x in manifest if x['sample_id']=='itrace_'+pdb and x['variant']=='full')
    model_pdb=CORRECTED/f'training_payload/dataset/structures/itrace_{pdb}.pdb'
    assert sha(model_pdb)==mi['input_sha256']
    shutil.copy2(model_pdb,out/f'{pdb}_model_input.pdb')
    n=len(pred['mean']);expected=''.join(pred['amino_acid']);chain_index=np.array([list(dict.fromkeys(pred['chain'])).index(c) for c in pred['chain']])
    chain_positions=[];seen={}
    for c in pred['chain']:seen[c]=seen.get(c,0)+1;chain_positions.append(seen[c])
    results=[];prov=[];frame_audit=[];sample_values=[]
    for rep in (1,2,3):
        rawtop=Path(f'E:/data/ITrace/trajectories/{pdb}_run{rep}_topology.pdb');rawxtc=Path(f'E:/data/ITrace/trajectories/{pdb}_run{rep}_md_processed.xtc')
        assert rawtop.exists() and rawxtc.exists()
        top_sha,xtc_sha=sha(rawtop),sha(rawxtc)
        top=STAGE/f'{pdb}_{rep}_canonical.pdb';xtc=STAGE/f'{pdb}_{rep}.xtc';canonical(rawtop,top)
        if not xtc.exists() or xtc.stat().st_size!=rawxtc.stat().st_size:shutil.copy2(rawxtc,xtc)
        assert sha(xtc)==xtc_sha
        u=mda.Universe(str(top),str(xtc))
        residues=[r for r in u.residues if r.resname in THREE_TO_ONE and 'CA' in r.atoms.names]
        assert len(residues)==n,(pdb,len(residues),n)
        assert ''.join(THREE_TO_ONE[r.resname] for r in residues)==expected,pdb
        np.testing.assert_array_equal(np.array([r.atoms.chainIDs[0] for r in residues]),pred['chain'])
        ca=np.array([r.atoms[r.atoms.names=='CA'].indices[0] for r in residues]);quads=[];which={}
        for i,r in enumerate(residues):
            look={a.name:a.index for a in r.atoms};definitions=CHI_ATOMS.get(r.resname,[])
            if definitions and all(all(name in look for name in d) for d in definitions):
                which[i]=list(range(len(quads),len(quads)+len(definitions)))
                quads.extend([[look[name] for name in d] for d in definitions])
        quads=np.array(quads,dtype=int);coords=[];angles=[];times=[]
        for ts in u.trajectory:
            xyz=ts.positions;coords.append(xyz[ca].copy());angles.append(_dihedral_degrees(xyz[quads]));times.append(float(ts.time))
        u.trajectory.close()
        coords=np.asarray(coords);angles=np.asarray(angles);times=np.asarray(times)
        assert len(times)>=1001 and abs(times[0])<1e-6 and abs(times[-1]-200000)<1e-3
        rmsf=chain_fit(coords,chain_index)
        statebins=np.clip(np.floor((angles+180)/120).astype(np.int8),0,2);H=np.full(n,np.nan)
        for i,columns in which.items():H[i]=_entropy(sum(statebins[:,col].astype(np.int32)*3**j for j,col in enumerate(columns)))
        values=np.stack([rmsf,H],axis=1);sample_values.append(values)
        for i,r in enumerate(residues):
            results.append({'pdb':pdb,'replica':rep,'index':i,'chain':str(pred['chain'][i]),'role':roles[pdb][str(pred['chain'][i])],'author_residue_id':str(pred['author_residue_id'][i]),'position_in_chain':chain_positions[i],'raw_MD_resid':int(r.resid),'aa':str(pred['amino_acid'][i]),'nchi':len(which.get(i,[])),'rmsf_A':float(rmsf[i]),'Hchi_nats':float(H[i]),'frames':len(times)})
        assert sha(rawtop)==top_sha and sha(rawxtc)==xtc_sha
        prov.append({'pdb':pdb,'replica':rep,'topology':str(rawtop),'topology_sha256':top_sha,'trajectory':str(rawxtc),'trajectory_sha256':xtc_sha,'canonical_topology_sha256':sha(top),'prediction_sha256':sha(prediction_path),'model_input_sha256':sha(model_pdb),'checkpoint_sha256':mi['checkpoint_sha256']})
        frame_audit.append({'pdb':pdb,'replica':rep,'frames':len(times),'first_ps':float(times[0]),'last_ps':float(times[-1]),'residues':n,'chi_angles':len(quads)})
        print(pdb,rep,'residues',n,'frames',len(times),'computed',flush=True)
    stacked=np.stack(sample_values);means=np.full((n,2),np.nan);variance=np.zeros((n,2));valid=np.isfinite(stacked).any(axis=0)
    for i in range(n):
        for j in range(2):
            x=stacked[:,i,j];x=x[np.isfinite(x)]
            if len(x):means[i,j]=x.mean();variance[i,j]=x.var(ddof=1) if len(x)>1 else 0
    np.testing.assert_array_equal(valid,pred['valid'][:,:2])
    metrics_audit=[]
    for j,metric in enumerate(['rmsf','hchi']):
        mask=valid[:,j]
        np.testing.assert_allclose(means[mask,j],pred['target'][mask,j],rtol=0,atol=3e-6)
        np.testing.assert_allclose(variance[mask,j],pred['label_variance'][mask,j],rtol=0,atol=3e-6)
        metrics_audit.append({'pdb':pdb,'metric':metric,'valid_residues':int(mask.sum()),'maximum_mean_absolute_error':float(np.max(abs(means[mask,j]-pred['target'][mask,j]))),'maximum_variance_absolute_error':float(np.max(abs(variance[mask,j]-pred['label_variance'][mask,j])))})
    np.savez_compressed(out/f'{pdb}_all_residue_replica_labels.npz',values=stacked,mean=means,sample_variance=variance,valid=valid,chain=pred['chain'],author_residue_id=pred['author_residue_id'],amino_acid=pred['amino_acid'],metric_names=np.array(['rmsf','hchi']))
    f=pd.DataFrame(results);f.to_csv(out/f'{pdb}_all_residue_per_replica_labels.csv',index=False)
    return f,pred,means,variance,prov,frame_audit,metrics_audit

def rho(x,y):return float(pd.Series(x).rank().corr(pd.Series(y).rank())) if len(x)>=4 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan

def pair_analysis(a,b,case):
    out=ROOT/case;out.mkdir(exist_ok=True)
    rows=[];provenance=[];frame_audit=[];metric_audit=[];data={}
    for pdb in [a,b]:
        f,p,m,v,prov,frames,audit=complex_labels(pdb,out);rows.append(f);provenance.extend(prov);frame_audit.extend(frames);metric_audit.extend(audit);data[pdb]=(p,m,v)
    pd.concat(rows,ignore_index=True).to_csv(out/'all_protein_per_residue_per_replica_labels.csv',index=False)
    frame=pd.read_csv(out/'original_pair_all_aligned_positions.csv',dtype={'id_a':str,'id_b':str,'chain_a':str,'chain_b':str})
    for side,pdb in [('a',a),('b',b)]:
        p,m,v=data[pdb];lookup={(str(c),str(i)):j for j,(c,i) in enumerate(zip(p['chain'],p['author_residue_id']))}
        indices=np.array([lookup[(str(c),str(i))] for c,i in zip(frame['chain_'+side],frame['id_'+side])]);np.testing.assert_array_equal(frame['aa_'+side],p['amino_acid'][indices])
        for j,metric in enumerate(['rmsf','hchi']):
            frame[metric+'_md_'+side]=m[indices,j];frame[metric+'_sd_'+side]=np.sqrt(v[indices,j]);frame[metric+'_full_'+side]=p['mean'][indices,j]
            for rep in (1,2,3):
                f=rows[0 if side=='a' else 1];frame[f'{metric}_rep{rep}_{side}']=f[f.replica==rep].set_index('index').loc[indices,'rmsf_A' if j==0 else 'Hchi_nats'].to_numpy()
    for metric in ['rmsf','hchi']:
        frame[metric+'_md_delta']=frame[metric+'_md_b']-frame[metric+'_md_a'];frame[metric+'_full_delta']=frame[metric+'_full_b']-frame[metric+'_full_a'];frame[metric+'_delta_sem']=np.sqrt((frame[metric+'_sd_a']**2+frame[metric+'_sd_b']**2)/3)
    frame.to_csv(out/'recomputed_pair_all_aligned_positions.csv',index=False)
    same=frame.same_aa;iface=frame.interface;mhc=frame.role=='mhc';pep=frame.role=='peptide';tcr=frame.role.isin(['tcr_alpha','tcr_beta']);cdr=frame.region_a.str.startswith('CDR')&frame.region_b.str.startswith('CDR');cdr3=frame.region_a.str.contains('CDR3')&frame.region_b.str.contains('CDR3')
    scopes={'all_common_AA':same,'all_interface_common_AA':same&iface,'MHC_interface_natural':mhc&iface,'MHC_interface_common_AA':same&mhc&iface,'peptide_all_natural':pep,'peptide_common_AA':same&pep,'MHC_plus_peptide_interface_natural':iface&(mhc|pep),'MHC_plus_peptide_interface_common_AA':same&iface&(mhc|pep),'MHC_groove_natural':mhc&(frame.position_a<=180)&(frame.position_b<=180),'MHC_groove_common_AA':same&mhc&(frame.position_a<=180)&(frame.position_b<=180),'TCR_interface_common_AA':same&iface&tcr,'TCR_variable_common_AA':same&tcr&(frame.position_a<=125)&(frame.position_b<=125),'TCR_CDR_common_AA':same&cdr,'CDR3_common_AA':same&cdr3}
    for name in ['CDR1α','CDR1β','CDR2α','CDR2β','CDR3α','CDR3β','α1 helix','α2 helix']:
        scopes[name+'_common_AA']=same&(frame.region_a==name)&(frame.region_b==name)
    stats=[]
    for scope,mask in scopes.items():
        for metric in ['rmsf','hchi']:
            f=frame[mask&frame[metric+'_valid']];x=f[metric+'_md_delta'].to_numpy();y=f[metric+'_full_delta'].to_numpy();sem=f[metric+'_delta_sem'].to_numpy();robust=abs(x)>2*sem
            if not len(f):continue
            stats.append({'a':a,'b':b,'scope':scope,'metric':metric,'n':len(f),'mutated_amino_acid_n':int((~f.same_aa).sum()),'rho_delta':rho(x,y),'MAE_delta':float(np.mean(abs(x-y))),'zero_change_MAE':float(np.mean(abs(x))),'relative_MAE_gain':float(1-np.mean(abs(x-y))/np.mean(abs(x))) if np.mean(abs(x))>0 else np.nan,'MD_mean_delta':float(np.mean(x)),'prediction_mean_delta':float(np.mean(y)),'MD_delta_SD':float(np.std(x)),'prediction_delta_SD':float(np.std(y)),'direction_fraction':float(np.mean(x*y>0)),'robust_n':int(robust.sum()),'robust_direction_fraction':float(np.mean(x[robust]*y[robust]>0)) if robust.any() else np.nan,'robust_relative_MAE_gain':float(1-np.mean(abs(x[robust]-y[robust]))/np.mean(abs(x[robust]))) if robust.any() else np.nan,'absolute_MAE_a':float(np.mean(abs(f[metric+'_md_a']-f[metric+'_full_a']))),'absolute_MAE_b':float(np.mean(abs(f[metric+'_md_b']-f[metric+'_full_b'])))})
    stat=pd.DataFrame(stats);stat.to_csv(out/'predefined_region_change_metrics.csv',index=False)
    contributions=[]
    for metric in ['rmsf','hchi']:
        f=frame[same&iface&frame[metric+'_valid']].copy();x=f[metric+'_md_delta'].to_numpy();y=f[metric+'_full_delta'].to_numpy();rx=pd.Series(x).rank().to_numpy();ry=pd.Series(y).rank().to_numpy();cross=(rx-rx.mean())*(ry-ry.mean())
        for group,mask in [('MHC',f.role=='mhc'),('peptide',f.role=='peptide'),('TCR',f.role.isin(['tcr_alpha','tcr_beta']))]:
            z=f[mask];xx=z[metric+'_md_delta'].to_numpy();yy=z[metric+'_full_delta'].to_numpy();rawmask=mask.to_numpy()
            contributions.append({'metric':metric,'component':group,'n':len(z),'fraction_of_interface_positions':len(z)/len(f),'sum_absolute_MD_change':float(np.sum(abs(xx))),'sum_absolute_prediction_error':float(np.sum(abs(xx-yy))),'sum_MAE_improvement_over_zero':float(np.sum(abs(xx)-abs(xx-yy))),'fraction_of_total_absolute_MD_change':float(np.sum(abs(xx))/np.sum(abs(x))),'fraction_of_global_rank_covariance_numerator':float(cross[rawmask].sum()/cross.sum()) if abs(cross.sum())>1e-12 else np.nan})
    pd.DataFrame(contributions).to_csv(out/'interface_component_decomposition.csv',index=False)
    audit={'status':'passed','pair':[a,b],'all_protein_residue_original_reference_means_and_sample_variances_recomputed':True,'absolute_tolerance':3e-6,'raw_topology_and_trajectory_files_unchanged':True,'frame_window':'all saved frames0–200ns, no subsampling','frames':frame_audit,'metric_cache_checks':metric_audit,'alignment':'Each chain CA independently Kabsch fitted to its first frame using double precision chain_fit identical to prepare_corrected_labels.py that generated ITrace targets','Hchi':'Original fixed3 equal120-degree bins per chi with joint rotamer entropy; invalid/no-chi positions remain missing','regions':'Original structural8A interface union and CDR annotations; all fields from original pair screen. Natural MHC/peptide positions separately analyzed from across-TCR identical-AA matches.','source_files':provenance,'new_MD':False,'new_training':False}
    (out/'raw_trajectory_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
    print(stat[stat.scope.isin(['all_interface_common_AA','MHC_interface_common_AA','peptide_common_AA','MHC_plus_peptide_interface_common_AA','TCR_interface_common_AA','CDR3_common_AA'])][['scope','metric','n','rho_delta','relative_MAE_gain','robust_n','robust_direction_fraction']].to_string(index=False),flush=True)

if __name__=='__main__':
    cases={'pik3ca':('7l1d','7rrg','pik3ca_validation'),'mel5':('3hg1','4jff','mel5_validation')}
    for name in (sys.argv[1:] or list(cases)):pair_analysis(*cases[name])
