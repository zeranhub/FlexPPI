from pathlib import Path
import sys,os,json,hashlib
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
EXP=Path(r'E:\Project\FlexPPI\experiments\_sources');sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
import numpy as np,pandas as pd,MDAnalysis as mda
from MDAnalysis.lib.distances import calc_bonds,calc_angles
R=Path(__file__).resolve().parent;V=EXP/'tcr_expanded_screen_v41_20261017/mel5_validation'
M=pd.read_csv(V/'corrected_pair_all_aligned_positions.csv',dtype={'id_a':str,'id_b':str});AUD=json.loads((V/'raw_trajectory_audit.json').read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
resrow=lambda chain,pos:M[(M.chain_a==chain)&(M.position_a==pos)].iloc[0]
rows=[];frames=[];checks=[]
for state,side in [('3hg1','a'),('4jff','b')]:
 for rep in [1,2,3]:
  a=next(x for x in AUD['source_files'] if x['pdb']==state and x['replica']==rep);assert sha(a['topology'])==a['topology_sha256'] and sha(a['trajectory'])==a['trajectory_sha256']
  u=mda.Universe(a['topology'],a['trajectory'])
  def atom(chain,pos,name):
   rr=resrow(chain,pos);ats=u.select_atoms('chainID '+str(rr['chain_'+side])+' and resid '+str(rr['id_'+side])+' and name '+name)
   assert len(ats)==1,(state,chain,pos,name);return ats.indices[0]
  pairs=[('E95_backboneNH_to_G101_O',atom('E',95,'N'),atom('E',95,'HN'),atom('E',101,'O'))]
  if state=='4jff':pairs.append(('W103_indoleNH_to_G99_O',atom('E',103,'NE1'),atom('E',103,'HE1'),atom('E',99,'O')))
  found={k:[] for k,_,_,_ in pairs};dist={k:[] for k,_,_,_ in pairs};ang={k:[] for k,_,_,_ in pairs}
  for ts in u.trajectory:
   for k,d,h,ac in pairs:
    distance=float(calc_bonds(ts.positions[d],ts.positions[ac],box=ts.dimensions));angle=float(np.degrees(calc_angles(ts.positions[d],ts.positions[h],ts.positions[ac],box=ts.dimensions)))
    hit=distance<=3.5 and angle>=150
    found[k].append(hit);dist[k].append(distance);ang[k].append(angle)
    frames.append({'pdb':state,'replica':rep,'frame':ts.frame,'time_ns':ts.time/1000,'interaction':k,'donor_acceptor_A':distance,'DHA_angle_deg':angle,'hydrogen_bond_geometry':hit})
  n=len(u.trajectory);u.trajectory.close()
  for k,_,_,_ in pairs:rows.append({'pdb':state,'replica':rep,'interaction':k,'frames':n,'Hbond_DA3p5_angle150_fraction':float(np.mean(found[k])),'DA_le3p5_fraction':float(np.mean(np.array(dist[k])<=3.5)),'median_DA_A':float(np.median(dist[k])),'median_DHA_deg':float(np.median(ang[k]))})
  if state=='4jff':
   w=np.asarray(found['W103_indoleNH_to_G99_O']);e=np.asarray(found['E95_backboneNH_to_G101_O']);rows.append({'pdb':state,'replica':rep,'interaction':'Both_intraloop_Hbond_geometry_same_frame','frames':n,'Hbond_DA3p5_angle150_fraction':float(np.mean(w&e))})
  print(state,rep,{k:round(float(np.mean(v)),4) for k,v in found.items()},flush=True)
pd.DataFrame(rows).to_csv(R/'intraloop_Hbond_geometry_per_replica.csv',index=False);pd.DataFrame(frames).to_csv(R/'intraloop_Hbond_geometry_all_frames.csv',index=False)
a={'status':'passed','new_MD':False,'new_model_inference':False,'frame_window':'All saved0-200ns frames of all3 trajectories for both states','definitions':'D-A<=3.5 A and D-H-A>=150 degrees, minimum-image convention. W103 uses NE1-HE1 donor and G99 O acceptor; E95 uses N-HN donor and G101 O acceptor.','WT_W103_comparison':'WT residue103 is Leu and has no indole NE1-HE1 donor; absence is chemistry, not an unmeasured zero Hbond occupancy. All-heavy G99-L103 comparison reported separately.','candidate_selection':'Exploratory follow-up to exhaustive13position/88neighbor contact screen; hypotheses generated from same trajectories, no independent holdout.','numbering':'All common/paper positions mapped through corrected pair table; mature input indices+1.','script_sha256':sha(__file__),'mapping_sha256':sha(V/'corrected_pair_all_aligned_positions.csv')}
(R/'intraloop_Hbond_analysis_audit.json').write_text(json.dumps(a,ensure_ascii=False,indent=2),encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False),flush=True)
