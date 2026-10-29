from pathlib import Path
import sys, os, json, hashlib, importlib.util
for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):
    os.environ[key]='1'
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[2]
SOURCE=EXP/'tcr_expanded_screen_v41_20261017'
STAGE=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\tcr_v41_allprotein_staged')
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
spec=importlib.util.spec_from_file_location('old_recompute',SOURCE/'recompute_allprotein_pair_labels.py')
original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array, calc_bonds, calc_angles

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as file:
        for part in iter(lambda:file.read(1024*1024),b''):h.update(part)
    return h.hexdigest()

def align_chain(coords):
    x=np.asarray(coords,dtype=np.float64)
    centered=x-x.mean(1,keepdims=True)
    ref=x[0]-x[0].mean(0)
    covariance=np.einsum('tni,nj->tij',centered,ref,optimize=True)
    u,_,vt=np.linalg.svd(covariance);rotation=u@vt
    neg=np.linalg.det(rotation)<0;u[neg,:,-1]*=-1;rotation[neg]=u[neg]@vt[neg]
    return np.einsum('tni,tij->tnj',centered,rotation,optimize=True)

def fluct(x):
    return np.sqrt(np.mean(np.sum((x-x.mean(0))**2,axis=-1),axis=0))

def H(state):
    return float(original._entropy(np.asarray(state)))

audit0=json.loads((SOURCE/'mel5_validation/raw_trajectory_audit.json').read_text())
mapping=pd.read_csv(SOURCE/'mel5_validation/corrected_pair_all_aligned_positions.csv',dtype={'id_a':str,'id_b':str})
qrow=mapping[(mapping.role=='mhc')&(mapping.id_a=='72')].iloc[0]
rows=mapping[(mapping.role=='tcr_beta')&mapping.position_a.between(50,54)].copy()
assert len(rows)==5 and rows.aa_a.tolist()==list('SVGIG') and rows.aa_b.tolist()==list('WGPFG')
perframe=[];perrep=[];windows=[];states_out=[];provenance=[];identity=[]
for pdb,side in [('3hg1','a'),('4jff','b')]:
    with np.load(EXP/f'tcr_case_v40_20261017/family_predictions/full/itrace_{pdb}.npz',allow_pickle=False) as z:
        pred={key:z[key].copy() for key in z.files}
    for rep in (1,2,3):
        record=next(x for x in audit0['source_files'] if x['pdb']==pdb and x['replica']==rep)
        top=STAGE/f'{pdb}_{rep}_canonical.pdb';xtc=STAGE/f'{pdb}_{rep}.xtc'
        assert top.exists() and xtc.exists()
        assert sha(top)==record['canonical_topology_sha256']
        assert sha(xtc)==record['trajectory_sha256']
        rawtop=Path(record['topology']);rawxtc=Path(record['trajectory'])
        assert sha(rawtop)==record['topology_sha256'] and sha(rawxtc)==record['trajectory_sha256']
        u=mda.Universe(str(top),str(xtc))
        residues=[r for r in u.residues if r.resname in original.THREE_TO_ONE and 'CA' in r.atoms.names]
        assert ''.join(original.THREE_TO_ONE[r.resname] for r in residues)==''.join(pred['amino_acid'])
        np.testing.assert_array_equal([r.atoms.chainIDs[0] for r in residues],pred['chain'])
        lookup={(str(c),str(author)):idx for idx,(c,author) in enumerate(zip(pred['chain'],pred['author_residue_id']))}
        qi=lookup[(str(qrow['chain_'+side]),str(qrow['id_'+side]))];q=residues[qi]
        assert q.resname=='GLN' and q.atoms.chainIDs[0]=='A' and str(q.resid)=='72'
        qnames={a.name:a.index for a in q.atoms};quads=np.array([[qnames[name] for name in definition] for definition in original.CHI_ATOMS['GLN']])
        qsc=np.array([qnames[name] for name in ('CB','CG','CD','OE1','NE2')])
        beta_indices=np.flatnonzero(pred['chain']=='E');beta_res=[residues[idx] for idx in beta_indices]
        ca=np.array([r.atoms[r.atoms.names=='CA'].indices[0] for r in beta_res])
        target=[]
        for _,row in rows.iterrows():
            pi=lookup[(str(row['chain_'+side]),str(row['id_'+side]))]
            r=residues[pi];names={a.name:a.index for a in r.atoms}
            heavy=np.array([a.index for a in r.atoms if not a.name.startswith('H')])
            assert original.THREE_TO_ONE[r.resname]==row['aa_'+side]
            target.append({'residue':r,'names':names,'heavy':heavy,'pi':pi,'beta_index':int(np.where(beta_indices==pi)[0][0]),'common_position':int(row.position_a),'input_position':int(row['position_'+side]),'input_id':str(row['id_'+side]),'aa':str(row['aa_'+side]),'model_rmsf':float(row['rmsf_full_'+side]),'reference_rmsf':float(row[f'rmsf_rep{rep}_{side}'])})
        identity.append({'pdb':pdb,'replica':rep,'Q72_prediction_index':qi,'Q72_raw_resid':int(q.resid),'beta_positions':[{'paper_or_common_position':t['common_position'],'input_id':t['input_id'],'input_position':t['input_position'],'raw_topology_resid':int(t['residue'].resid),'aa':t['aa'],'prediction_index':t['pi']} for t in target]})
        allca=[];angles=[];time=[];contact=[];min_distance=[];specific_distance=[];hbonds=[]
        for ts in u.trajectory:
            xyz=ts.positions;box=ts.dimensions
            allca.append(xyz[ca].copy());angles.append(original._dihedral_degrees(xyz[quads]));time.append(float(ts.time)/1000)
            minima=[float(distance_array(xyz[qsc],xyz[t['heavy']],box=box).min()) for t in target]
            min_distance.append(minima);contact.append([[d<=4.0,d<=4.5] for d in minima])
            # Specific sidechain-backbone partners are specified before evaluating occupancies.
            phe=target[3];gly=target[4]
            qoe_pheN=float(calc_bonds(xyz[qnames['OE1']],xyz[phe['names']['N']],box=box))
            qne_glyO=float(calc_bonds(xyz[qnames['NE2']],xyz[gly['names']['O']],box=box))
            phe_dha=float(np.degrees(calc_angles(xyz[phe['names']['N']],xyz[phe['names']['HN']],xyz[qnames['OE1']],box=box)))
            q_dha=max(float(np.degrees(calc_angles(xyz[qnames['NE2']],xyz[qnames[hydrogen]],xyz[gly['names']['O']],box=box))) for hydrogen in ('HE21','HE22'))
            specific_distance.append([qoe_pheN,qne_glyO]);hbonds.append([qoe_pheN<=3.5 and phe_dha>=150,qne_glyO<=3.5 and q_dha>=150])
        u.trajectory.close()
        time=np.array(time);angles=np.array(angles);aligned=align_chain(allca)
        bins=np.clip(np.floor((angles+180)/120).astype(np.int8),0,2)
        state=sum(bins[:,j].astype(np.int32)*3**j for j in range(3))
        contact=np.array(contact,dtype=bool);min_distance=np.array(min_distance);specific_distance=np.array(specific_distance);hbonds=np.array(hbonds,dtype=bool)
        rmsf=fluct(aligned);h=H(state)
        np.testing.assert_allclose(h,float(qrow[f'hchi_rep{rep}_{side}']),atol=1e-12,rtol=0)
        for t in target:np.testing.assert_allclose(rmsf[t['beta_index']],t['reference_rmsf'],atol=1e-12,rtol=0)
        counts=np.bincount(state,minlength=27)
        for s,count in enumerate(counts):
            states_out.append({'pdb':pdb,'replica':rep,'joint_state':s,'chi1_bin':s%3,'chi2_bin':(s//3)%3,'chi3_bin':s//9,'count':int(count),'fraction':float(count/len(state)),'Q72_NE2_GlyO_contact_le4_A_fraction_in_state':float((specific_distance[state==s,1]<=4).mean()) if count else None,'Q72_NE2_GlyO_hbond_fraction_in_state':float(hbonds[state==s,1].mean()) if count else None})
        for f,tns in enumerate(time):
            item={'pdb':pdb,'replica':rep,'frame':f,'time_ns':float(tns),'Q72_chi1_deg':float(angles[f,0]),'Q72_chi2_deg':float(angles[f,1]),'Q72_chi3_deg':float(angles[f,2]),'Q72_joint_state':int(state[f]),'Q72_OE1_beta_I53_or_F54_N_A':float(specific_distance[f,0]),'Q72_NE2_beta_G54_or_G55_O_A':float(specific_distance[f,1]),'Q72_OE1_betaN_hbond':bool(hbonds[f,0]),'Q72_NE2_betaGlyO_hbond':bool(hbonds[f,1])}
            for j,t in enumerate(target):item[f'Q72_sidechain_beta_common{t["common_position"]}_minheavy_A']=float(min_distance[f,j])
            perframe.append(item)
        for j,t in enumerate(target):
            perrep.append({'pdb':pdb,'replica':rep,'frames':len(time),'first_ns':float(time[0]),'last_ns':float(time[-1]),'Q72_Hchi_nats':h,'beta_common_position':t['common_position'],'beta_input_position':t['input_position'],'beta_input_id':t['input_id'],'beta_aa':t['aa'],'beta_CA_RMSF_A':float(rmsf[t['beta_index']]),'model_beta_RMSF_A':t['model_rmsf'],'min_heavy_contact_le4_A_fraction':float(contact[:,j,0].mean()),'min_heavy_contact_le4p5_A_fraction':float(contact[:,j,1].mean()),'median_min_heavy_distance_A':float(np.median(min_distance[:,j])),'Q72_NE2_beta_GlyO_le3p5_A_fraction':float((specific_distance[:,1]<=3.5).mean()),'Q72_NE2_beta_GlyO_le4_A_fraction':float((specific_distance[:,1]<=4).mean()),'Q72_NE2_beta_GlyO_hbond_fraction':float(hbonds[:,1].mean()),'Q72_OE1_beta_Ile_or_PheN_le3p5_A_fraction':float((specific_distance[:,0]<=3.5).mean()),'Q72_OE1_beta_Ile_or_PheN_hbond_fraction':float(hbonds[:,0].mean())})
        for start in range(0,200,20):
            mask=(time>=start)&(time<start+20 if start<180 else time<=200)
            wr=fluct(aligned[mask]);wh=H(state[mask])
            for j,t in enumerate(target):
                windows.append({'pdb':pdb,'replica':rep,'start_ns':start,'end_ns':start+20,'frames':int(mask.sum()),'Q72_Hchi_nats':wh,'beta_common_position':t['common_position'],'beta_input_position':t['input_position'],'beta_aa':t['aa'],'beta_CA_RMSF_A':float(wr[t['beta_index']]),'min_heavy_contact_le4_A_fraction':float(contact[mask,j,0].mean()),'min_heavy_contact_le4p5_A_fraction':float(contact[mask,j,1].mean()),'Q72_NE2_beta_GlyO_le3p5_A_fraction':float((specific_distance[mask,1]<=3.5).mean()),'Q72_NE2_beta_GlyO_hbond_fraction':float(hbonds[mask,1].mean())})
        provenance.append({**record,'staged_topology':str(top),'staged_trajectory':str(xtc),'staged_topology_sha256':sha(top),'staged_trajectory_sha256':sha(xtc),'raw_files_unchanged':sha(rawtop)==record['topology_sha256'] and sha(rawxtc)==record['trajectory_sha256']})
        print(pdb,rep,'Hchi',round(h,4),'Gly RMSF',round(float(rmsf[target[4]['beta_index']]),4),'contact4',round(float(contact[:,4,0].mean()),4),'Gly Hbond',round(float(hbonds[:,1].mean()),4),flush=True)

pd.DataFrame(perframe).to_csv(ROOT/'per_frame_Q72_beta_contacts_and_states.csv',index=False)
repdf=pd.DataFrame(perrep);repdf.to_csv(ROOT/'per_replica_Q72_beta_contact_and_RMSF.csv',index=False)
wd=pd.DataFrame(windows);wd.to_csv(ROOT/'nonoverlapping_20ns_window_descriptors.csv',index=False)
pd.DataFrame(states_out).to_csv(ROOT/'Q72_joint_chi_state_occupancies.csv',index=False)
corr=[]
for (pdb,rep,pos),g in wd.groupby(['pdb','replica','beta_common_position']):
    for x,y in [('Q72_Hchi_nats','beta_CA_RMSF_A'),('min_heavy_contact_le4_A_fraction','beta_CA_RMSF_A'),('Q72_Hchi_nats','min_heavy_contact_le4_A_fraction')]:
        corr.append({'pdb':pdb,'replica':int(rep),'beta_common_position':int(pos),'windows':len(g),'x':x,'y':y,'Spearman_descriptive':float(g[x].rank().corr(g[y].rank())) if g[x].nunique()>1 and g[y].nunique()>1 else None,'Pearson_descriptive':float(g[x].corr(g[y])) if g[x].nunique()>1 and g[y].nunique()>1 else None})
pd.DataFrame(corr).to_csv(ROOT/'within_replica_20ns_descriptive_correlations.csv',index=False)
summary=repdf.groupby(['pdb','beta_common_position','beta_input_position','beta_aa']).agg(Q72_Hchi_mean=('Q72_Hchi_nats','mean'),Q72_Hchi_sd=('Q72_Hchi_nats','std'),beta_CA_RMSF_mean_A=('beta_CA_RMSF_A','mean'),beta_CA_RMSF_sd_A=('beta_CA_RMSF_A','std'),contact4_mean=('min_heavy_contact_le4_A_fraction','mean'),contact4_sd=('min_heavy_contact_le4_A_fraction','std'),contact4p5_mean=('min_heavy_contact_le4p5_A_fraction','mean'),Gly_O_Hbond_mean=('Q72_NE2_beta_GlyO_hbond_fraction','mean'),Gly_O_contact3p5_mean=('Q72_NE2_beta_GlyO_le3p5_A_fraction','mean'),Ile_or_Phe_N_Hbond_mean=('Q72_OE1_beta_Ile_or_PheN_hbond_fraction','mean')).reset_index()
summary.to_csv(ROOT/'three_replica_summary.csv',index=False)
audit={'status':'passed','new_MD':False,'new_inference':False,'source_window':'All saved frames from 0 to 200 ns, all 3 independent replicas for both states; no frame removal/subsampling. Window partition [0,20),...,[180,200] ns.','chain_alignment':'Each entire beta chain CA fitted independently to its first frame in double precision; exact original labels reproduced before any window descriptor. Window RMSF uses each window mean after the same entire-chain fitting.','Hchi':'Fixed 3 bins per chi, [-180,-60), [-60,60), [60,180]; joint state entropy natural logarithm. Q72 has 3 chis; original per-replica labels exactly reproduced. Finite 20ns windows are descriptive and can underestimate slowly exchanging state diversity.','contact':'Minimum distance from Q72 sidechain heavy atoms CB/CG/CD/OE1/NE2 to all heavy atoms of each target beta residue; <=4.0 A and <=4.5 A reported, minimum-image convention using trajectory box. This is geometric contact occupancy.','hydrogen_bond_geometry':'Q72 NE2 (either HE21/HE22) to beta Gly54/Gly55 backbone O; beta Ile53/Phe54 backbone N-HN to Q72 OE1; D-A <=3.5 A and D-H-A >=150 degrees. No energetic or causal assertion.','replication':'Three trajectories per state; no individual saved frame or window treated as an independent experimental replicate. Correlations are within-replica descriptive quantities only; no p-values/causality.','model_and_sequence_mapping':'Corrected all-position structural sequence alignment, WT beta50-54 SVGIG to engineered input beta51-55 WGPFG. Common/paper beta54 Gly is engineered model input beta55 Gly. Q72 is unchanged HLA chain A position72.','MDAnalysis_version':mda.__version__,'identity':identity,'source_files':provenance,'script_sha256':sha(__file__)}
(ROOT/'analysis_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(summary.to_string(index=False),flush=True)
