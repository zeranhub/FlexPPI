from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'case_analysis'
OUT.mkdir(exist_ok=True)
BASE=Path('E:/Project/FlexPPI/experiments/_case_validation_sources')
ANGLE=BASE/'ubiquitin/inputs/A003R_Q41_Q62_angles_all_5001_frames.csv'
PROVENANCE=BASE/'data/A003R_Q41_Q62_provenance_and_QA.json'
SUMMARY=BASE/'data/ubiquitin_Q41_Q62_validation_summary.csv'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
prov=json.loads(PROVENANCE.read_text(encoding='utf8'))
cache=Path(prov['authoritative_numeric_cache'])
assert sha(cache)==prov['cache_sha256']
plan={'question':'How do the existing same-type side-chain case descriptors relate to sampled angular distributions?',
      'case':'DynaRepo A003R / 2QHO ubiquitin Q41 and Q62; already chosen in v20',
      'input_sha256':{str(p):sha(p) for p in [ANGLE,PROVENANCE,SUMMARY,cache]},
      'frames':'Original 0..4999 frames, five non-overlapping 1000-frame label windows; interval 100ps',
      'bin_offsets_degrees':[0,30,60,90],'offset_definition':'wrap(theta-phase) to [-180,180), then apply original 120-degree bins',
      'chi_histogram_width_degrees':5,
      'projection':'F_star=-log(p_bin/p_max) over visited equal-width bins; empty bins remain missing',
      'temperature':'dimensionless thermal units; no assumption of known temperature or equilibrium',
      'target':'Frozen corrected joint Hchi unchanged; chi1 marginals reported separately'}
(OUT/'ubiquitin_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf8')
def bins(angles,phase=0):
    if phase==0:return np.clip(np.floor((angles+180)/120).astype(int),0,2)
    return np.floor(((angles+180-phase)%360)/120).astype(int)
def entropy(values,n):
    p=np.bincount(values,minlength=n).astype(float)/len(values)
    return float(-np.sum(p[p>0]*np.log(p[p>0])))
def jsd(a,b):
    m=(a+b)/2
    return float(.5*sum(np.sum(x[x>0]*np.log(x[x>0]/m[x>0])) for x in [a,b]))
angles=pd.read_csv(ANGLE)
values=pd.read_csv(SUMMARY).set_index('residue')
cache_arrays=np.load(cache,allow_pickle=False)
hist_rows=[];stat_rows=[];state_rows=[];transition_rows=[]
for residue in ['Q41','Q62']:
    d=angles[angles.residue==residue].sort_values('frame_zero_based')
    assert len(d)==5001 and np.array_equal(d.frame_zero_based,np.arange(5001))
    d=d[d.frame_zero_based<5000].copy()
    a=d[['chi1_degrees','chi2_degrees','chi3_degrees']].to_numpy()
    assert np.isfinite(a).all()
    assert np.allclose(np.diff(d.time_ns),.1)
    base=bins(a);codes=(base*np.array([1,3,9])).sum(axis=1)
    assert np.array_equal(codes,d.joint_code.to_numpy())
    index=40 if residue=='Q41' else 61
    label_h=float(np.mean([entropy(codes[k:k+1000],27) for k in range(0,5000,1000)]))
    assert abs(label_h-float(cache_arrays['target'][index,1]))<1e-6
    assert abs(float(cache_arrays['mean'][index,1])-values.loc[residue,'model_hchi'])<1e-6
    for phase in plan['bin_offsets_degrees']:
        b=bins(a,phase);joint=(b*np.array([1,3,9])).sum(axis=1)
        h=[entropy(joint[k:k+1000],27) for k in range(0,5000,1000)]
        h1=[entropy(b[k:k+1000,0],3) for k in range(0,5000,1000)]
        stat_rows.append({'residue':residue,'phase':phase,'joint_window_mean_h':np.mean(h),
                         'joint_pooled_h':entropy(joint,27),'chi1_window_mean_h':np.mean(h1),
                         'chi1_pooled_h':entropy(b[:,0],3),'joint_window_sd':np.std(h,ddof=1),
                         'md_rmsf_A':float(cache_arrays['target'][index,0]),
                         'pred_rmsf_A':float(cache_arrays['mean'][index,0]),
                         'frozen_joint_h':label_h,'pred_joint_h':float(cache_arrays['mean'][index,1])})
        occupancy=np.bincount(b[:,0],minlength=3)/len(b)
        nchange=int((np.diff(b[:,0])!=0).sum())
        transition_rows.append({'residue':residue,'phase':phase,'saved_frame_changes':nchange,
                                'saved_frame_change_fraction':nchange/(len(b)-1),'step_ps':100,
                                'first_last_half_occupancy_jsd':jsd(np.bincount(b[:2500,0],minlength=3)/2500,np.bincount(b[2500:,0],minlength=3)/2500)})
        for state,p in enumerate(occupancy):
            name=['g_minus','g_plus','trans'][state] if phase==60 else str(state)
            state_rows.append({'residue':residue,'phase':phase,'state':state,'name':name,'probability':p})
    for window in [0,1,2,3,4,5]:
        aq=a if window==0 else a[(window-1)*1000:window*1000]
        for angle_index in range(3):
            counts,edges=np.histogram(aq[:,angle_index],bins=np.arange(-180,180.0001,5))
            p=counts/counts.sum()
            f=np.full(len(counts),np.nan);f[p>0]=-np.log(p[p>0]/p.max())
            for j in range(len(counts)):
                hist_rows.append({'residue':residue,'window':window,'chi':angle_index+1,'angle_lower':edges[j],
                                  'angle_upper':edges[j+1],'angle_center':(edges[j]+edges[j+1])/2,
                                  'frames':len(aq),'count':int(counts[j]),'probability':p[j],
                                  'density_per_degree':p[j]/5,'F_star':f[j]})
    d.to_csv(OUT/(residue+'_angles_used.csv'),index=False)
pd.DataFrame(stat_rows).to_csv(OUT/'ubiquitin_bin_sensitivity.csv',index=False)
pd.DataFrame(state_rows).to_csv(OUT/'ubiquitin_chi1_state_occupancies.csv',index=False)
pd.DataFrame(hist_rows).to_csv(OUT/'ubiquitin_chi1_histograms_and_Fstar.csv',index=False)
pd.DataFrame(transition_rows).to_csv(OUT/'ubiquitin_saved_frame_changes.csv',index=False)
report={'status':'complete','plan_sha256':sha(OUT/'ubiquitin_plan.json'),
        'code_sha256':sha(__file__),'frames_per_residue':5000,'residues':2,
        'all_frames_used':True,'frozen_joint_label_and_prediction_reproduced':True,
        'visited_bins_only_F_star':True,'actual_new_trajectories':0,
        'summary':stat_rows,'saved_frame_changes':transition_rows}
(OUT/'ubiquitin_analysis_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(pd.DataFrame(stat_rows).to_string(index=False))
print(pd.DataFrame(transition_rows).to_string(index=False))
