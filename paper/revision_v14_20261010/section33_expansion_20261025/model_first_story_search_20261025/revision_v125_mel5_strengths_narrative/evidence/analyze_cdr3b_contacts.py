from pathlib import Path
import sys,os,json,hashlib,time
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
EXP=Path(r'E:\Project\FlexPPI\experiments\_sources')
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
ROOT=Path(__file__).resolve().parent
VAL=EXP/'tcr_expanded_screen_v41_20261017/mel5_validation'
MAP=pd.read_csv(VAL/'corrected_pair_all_aligned_positions.csv',dtype={'id_a':str,'id_b':str})
REF=pd.read_csv(VAL/'plot_sources/mel5_all6_CDR_per_position_source.csv',float_precision='round_trip')
AUDIT=json.loads((VAL/'raw_trajectory_audit.json').read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
THREE={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def parse(path):
 d={};names={}
 for l in Path(path).read_text().splitlines():
  if not l.startswith(('ATOM  ','HETATM')) or l[17:20].strip() not in THREE:continue
  if l[16] not in [' ','A']:continue
  name=l[12:16].strip()
  if name.startswith('H') or l[76:78].strip() in ['H','D']:continue
  key=(l[21],l[22:27].strip())
  d.setdefault(key,[]).append([float(l[30:38]),float(l[38:46]),float(l[46:54])]);names[key]=THREE[l[17:20].strip()]
 return {k:np.array(x) for k,x in d.items()},names

S={};rows=MAP.to_dict('records');qlist=[i for i,r in enumerate(rows) if r['role']=='tcr_beta' and r['region']=='CDR3β'];assert len(qlist)==13
static={};candidate=set()
for state,side in [('3hg1','a'),('4jff','b')]:
 P,N=parse(VAL/(state+'_model_input.pdb'));S[state]=(P,N)
 for qi in qlist:
  qr=rows[qi];key=(qr['chain_'+side],qr['id_'+side]);assert N[key]==qr['aa_'+side]
  for vi,vr in enumerate(rows):
   k=(vr['chain_'+side],vr['id_'+side]);assert N[k]==vr['aa_'+side]
   if key[0]==k[0] and abs(int(qr['position_'+side])-int(vr['position_'+side]))<=2:continue
   md=float(np.sqrt(((P[key][:,None,:]-P[k][None,:,:])**2).sum(-1)).min());static[(state,qi,vi)]=md
   if md<=7.:candidate.add(vi)
vs=sorted(candidate);print('13 CDR3beta positions; candidate residues within7A of either start:',len(vs),flush=True)
np.savez(ROOT/'static_candidate_index.npz',query_indices=qlist,neighbor_indices=vs)
META={};pairs=[];repstats=[];identity=[];allR=[]
for state,side in [('3hg1','a'),('4jff','b')]:
 for rep in [1,2,3]:
  start=time.time();r=next(z for z in AUDIT['source_files'] if z['pdb']==state and z['replica']==rep)
  top=Path(r['topology']);traj=Path(r['trajectory']);assert sha(top)==r['topology_sha256'] and sha(traj)==r['trajectory_sha256']
  # Protein heavy-atom distance calculations do not require atom-name canonicalisation.
  u=mda.Universe(str(top),str(traj));protein=[z for z in u.residues if z.resname in THREE or z.resname in ['HID','HIE','HIP','HSD','HSE','HSP']]
  bykey={(str(z.atoms.chainIDs[0]),str(z.resid)):z for z in protein}
  ordered=[]
  for rr in rows:
   key=(rr['chain_'+side],rr['id_'+side]);rrr=bykey[key];aa=THREE.get(rrr.resname,'H');assert aa==rr['aa_'+side],(state,key,aa,rr['aa_'+side]);ordered.append(rrr)
  # Check peptide identities explicitly, including P9 Thr and P10 Val.
  pep=[z for rr,z in zip(rows,ordered) if rr['role']=='peptide'];assert ''.join(THREE[z.resname] for z in pep)=='ELAGIGILTV'
  def heavy(res):return np.array([a.index for a in res.atoms if not a.name.startswith('H') and a.element not in ['H','D']],dtype=int)
  qa=[heavy(ordered[i]) for i in qlist];va=[heavy(ordered[i]) for i in vs]
  qi=np.concatenate(qa);vi=np.concatenate(va);qstart=np.r_[0,np.cumsum([len(x) for x in qa])[:-1]];vstart=np.r_[0,np.cumsum([len(x) for x in va])[:-1]]
  beta=[z for z in protein if z.atoms.chainIDs[0]=='E'];beta_ca=np.array([z.atoms[z.atoms.names=='CA'].indices[0] for z in beta]);betaindex={str(z.resid):i for i,z in enumerate(beta)}
  target_ca_indices=[betaindex[rows[i]['id_'+side]] for i in qlist]
  ca_frames=[];distances=[];times=[]
  for ts in u.trajectory:
   mat=distance_array(ts.positions[qi],ts.positions[vi],box=ts.dimensions)
   # Two reductions turn atom-pair distances into residue-pair nearest distances.
   d=np.minimum.reduceat(np.minimum.reduceat(mat,qstart,axis=0),vstart,axis=1)
   distances.append(d.astype(np.float32));times.append(float(ts.time)/1000);ca_frames.append(ts.positions[beta_ca].copy())
  u.trajectory.close();ds=np.array(distances);times=np.array(times);assert len(times)>=1001 and times[0]==0 and abs(times[-1]-200)<1e-6
  x=np.asarray(ca_frames,dtype=np.float64);cen=x-x.mean(1,keepdims=True);ref=x[0]-x[0].mean(0);cov=np.einsum('tni,nj->tij',cen,ref,optimize=True);U,_,Vt=np.linalg.svd(cov);rot=U@Vt;neg=np.linalg.det(rot)<0;U[neg,:,-1]*=-1;rot[neg]=U[neg]@Vt[neg];aligned=np.einsum('tni,tij->tnj',cen,rot,optimize=True)
  rms=np.sqrt(((aligned-aligned.mean(0))**2).sum(-1).mean(0))[target_ca_indices]
  for j,idx in enumerate(qlist):
   rr=rows[idx];beforepos=int(rr['position_a']);reference=REF[(REF.region=='CDR3β') & (REF.position_a==beforepos)].iloc[0]
   np.testing.assert_allclose(rms[j],reference[f'rmsf_rep{rep}_{side}'],atol=1e-12,rtol=0)
   repstats.append({'pdb':state,'replica':rep,'common_position':beforepos,'aa':rr['aa_'+side],'RMSF_A':rms[j],'model_RMSF_A':rr['rmsf_full_'+side]})
   for k,vidx in enumerate(vs):
    nr=rows[vidx]
    if rr['chain_'+side]==nr['chain_'+side] and abs(int(rr['position_'+side])-int(nr['position_'+side]))<=2:continue
    pairs.append({'pdb':state,'replica':rep,'CDR3b_common_position':beforepos,'CDR3b_aa':rr['aa_'+side],'neighbor_role':nr['role'],'neighbor_region':nr['region'],'neighbor_chain':nr['chain_a'],'neighbor_common_position':nr['position_a'],'neighbor_aa':nr['aa_'+side],'nframes':len(times),'contact_le4_fraction':float((ds[:,j,k]<=4).mean()),'contact_le4p5_fraction':float((ds[:,j,k]<=4.5).mean()),'median_distance_A':float(np.median(ds[:,j,k])),'static_distance_A':static.get((state,idx,vidx),np.nan)})
  # Local loop-aligned RMSF separates loop-internal fluctuations from whole-loop displacement.
  loop=aligned[:,target_ca_indices];lc=loop-loop.mean(1,keepdims=True);lr=loop[0]-loop[0].mean(0);U,_,Vt=np.linalg.svd(np.einsum('tni,nj->tij',lc,lr,optimize=True));rot=U@Vt;neg=np.linalg.det(rot)<0;U[neg,:,-1]*=-1;rot[neg]=U[neg]@Vt[neg];la=np.einsum('tni,tij->tnj',lc,rot,optimize=True);lrms=np.sqrt(((la-la.mean(0))**2).sum(-1).mean(0))
  for j,idx in enumerate(qlist):allR.append({'pdb':state,'replica':rep,'common_position':int(rows[idx]['position_a']),'loop_aligned_RMSF_A':float(lrms[j]),'whole_beta_chain_aligned_RMSF_A':float(rms[j])})
  np.savez_compressed(ROOT/f'{state}_rep{rep}_all_candidate_distances.npz',distances_A=ds,time_ns=times,CDR3b_common_positions=[rows[i]['position_a'] for i in qlist],neighbor_common_positions=[rows[i]['position_a'] for i in vs],neighbor_chains=[rows[i]['chain_a'] for i in vs],neighbor_AA=[rows[i]['aa_'+side] for i in vs])
  identity.append({**r,'frames':len(times),'all_original_CDR3b_RMSF_labels_exactly_reproduced':True,'protein_residues':len(protein),'peptide_sequence':'ELAGIGILTV','candidate_neighbor_residues':len(vs)})
  print(state,rep,'frames',len(times),'betaRMSFmean',round(float(rms.mean()),5),'localLoopRMSFmean',round(float(lrms.mean()),5),'seconds',round(time.time()-start,1),flush=True)
pd.DataFrame(pairs).to_csv(ROOT/'all_CDR3b_nonlocal_contacts_per_replica.csv',index=False)
pd.DataFrame(repstats).to_csv(ROOT/'all_CDR3b_RMSF_per_replica.csv',index=False)
pd.DataFrame(allR).to_csv(ROOT/'all_CDR3b_loop_aligned_RMSF_per_replica.csv',index=False)
out=pd.DataFrame(pairs);key=['pdb','CDR3b_common_position','CDR3b_aa','neighbor_role','neighbor_region','neighbor_chain','neighbor_common_position','neighbor_aa'];summ=out.groupby(key,dropna=False).agg(contact4_mean=('contact_le4_fraction','mean'),contact4_SD=('contact_le4_fraction','std'),contact4p5_mean=('contact_le4p5_fraction','mean'),median_distance_replica_mean_A=('median_distance_A','mean'),static_distance_A=('static_distance_A','first')).reset_index();summ.to_csv(ROOT/'all_CDR3b_nonlocal_contacts_three_replica_summary.csv',index=False)
audit={'status':'passed','new_MD':False,'new_model_inference':False,'all_saved_frames':True,'frame_window_ns':[0,200],'replicas_per_state':3,'contact_definition':'minimum heavy-atom distance, minimum image convention, <=4A and <=4.5A; occupancy is not residence lifetime.','candidate_scope':'All mapped protein residues within7A of any CDR3beta residue in either bound starting structure. Same-chain sequence-neighbors separated by <=2 excluded. All13 loop positions retained.','peptide_sequence_verified':'ELAGIGILTV, P9 Thr P10 Val','RMSF_validation':'All13 positions in all6 trajectories reproduce archived chain-aligned double-precision RMSF to1e-12. Local-loop fitted RMSF is an additional descriptive metric, not model target.','source_files':identity,'mapping_file':str(VAL/'corrected_pair_all_aligned_positions.csv'),'mapping_sha256':sha(VAL/'corrected_pair_all_aligned_positions.csv'),'script_sha256':sha(__file__)}
(ROOT/'analysis_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print('COMPLETE',len(pairs),'residue-pair-replica observations',flush=True)
