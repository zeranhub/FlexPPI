"""Independent CPU verification of saved trajectory analysis outputs."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent


def digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):
            value.update(block)
    return value.hexdigest()


def entropy(codes):
    counts=np.bincount(codes)
    counts=counts[counts>0].astype(float)
    p=counts/counts.sum()
    return float(-(p*np.log(p)).sum())


def bins(angles,offset):
    if offset==0:
        transformed=angles
        return np.minimum(2,np.maximum(0,np.floor((transformed+np.float32(180))/np.float32(120)))).astype(int)
    transformed=(angles.astype(float)-offset+180)%360-180
    return np.floor((transformed+180)/120).astype(int)


def main():
    plan=json.loads((ROOT/'analysis_plan.json').read_text(encoding='utf-8'))
    complete=json.loads((ROOT/'completion_audit.json').read_text(encoding='utf-8'))
    assert complete['status']=='complete'
    assert digest(ROOT/'analysis_plan.json')==complete['plan_sha256']
    assert digest(ROOT/'run_angle_distribution_analysis.py')==plan['analysis_code_sha256']==complete['analysis_code_sha256']
    for item in plan['dependency_hashes'].values():
        assert digest(item['path'])==item['sha256']
    for item in complete['output_tables'].values():
        assert digest(item['path'])==item['sha256']
    original=json.loads(Path(plan['inherited_plan_path']).read_text(encoding='utf-8'))
    assert digest(plan['inherited_plan_path'])==plan['inherited_plan_sha256']
    assert [r['sample_id'] for r in original['records']]==[r['sample_id'] for r in plan['cohort']]
    tables=pd.read_csv(ROOT/'per_residue_trajectory.csv.gz',keep_default_na=False)
    repeat=pd.read_csv(ROOT/'repeat_jsd_per_residue.csv.gz')
    masks=pd.read_csv(ROOT/'valid_mask_audit.csv')
    assert len(tables)==87711*4
    assert len(tables.groupby(['sample_id','replica_index']))==203
    assert set(tables.shift_degrees)=={0,30,60,90}
    assert len(tables.groupby(['sample_id','shift_degrees']))==83*4
    assert (masks.old_only==0).all() and (masks.new_only==0).all()
    original_valid=masks[(masks.variant=='uniform_rep1')&(masks.shift_degrees==0)].groupby('metric')['common'].sum().to_dict()
    assert original_valid=={'chi_entropy':33493,'rmsf':38599}
    pro=tables.residue_name=='PRO'
    assert (tables.loc[pro,'chi1_conventional_eligible']==0).all()
    assert tables.loc[pro,'chi1_conventional_exclusion_reason'].str.contains('ring-constrained').all()
    for column in ['conventional_chi1_p_gminus','conventional_chi1_p_gplus','conventional_chi1_p_trans']:
        assert (tables.loc[pro,column]=='').all()
        assert (tables.loc[tables.shift_degrees!=60,column]=='').all()
    for column in ['chi1_p0','chi1_p1','chi1_p2','h_joint','h_chi1','half_joint_jsd','half_chi1_jsd',
                   'ca_total_variance_A2','ca_within_joint_variance_A2','ca_between_joint_variance_A2',
                   'ca_within_chi1_variance_A2','ca_between_chi1_variance_A2','full_minus_window_h_weighted']:
        tables[column]=pd.to_numeric(tables[column],errors='coerce')
    np.testing.assert_allclose(tables[['chi1_p0','chi1_p1','chi1_p2']].sum(axis=1),1,rtol=0,atol=1e-12)
    assert (tables.h_joint>=-1e-12).all() and (tables.h_joint<=tables.chi_count*np.log(3)+1e-12).all()
    assert (tables.h_chi1>=-1e-12).all() and (tables.h_chi1<=np.log(3)+1e-12).all()
    assert (tables.full_minus_window_h_weighted>=-1e-10).all()
    for frame,columns in [(tables,['half_joint_jsd','half_chi1_jsd']),(repeat,['joint_repeat_jsd','chi1_repeat_jsd'])]:
        for column in columns:
            assert (frame[column]>=0).all() and (frame[column]<=np.log(2)+1e-12).all()
    joint_error=np.abs(tables.ca_total_variance_A2-tables.ca_within_joint_variance_A2-tables.ca_between_joint_variance_A2)
    chi1_error=np.abs(tables.ca_total_variance_A2-tables.ca_within_chi1_variance_A2-tables.ca_between_chi1_variance_A2)
    assert (joint_error<=1e-8*np.maximum(tables.ca_total_variance_A2,1)).all()
    assert (chi1_error<=1e-8*np.maximum(tables.ca_total_variance_A2,1)).all()
    # Conventional chi1 center and trans-seam tests independently calculated.
    probe=np.asarray([-180,-179,179,180,-60,60,-120,0,120],dtype=np.float32)
    assert bins(probe,60).tolist()==[2,2,2,2,0,1,0,1,2]
    assert bins(np.asarray([-180,-60,60,180],dtype=np.float32),0).tolist()==[0,1,2,2]
    cache_reports=[]
    indexed=tables.set_index(['sample_id','replica_index','residue_index','shift_degrees'])
    for row in original['records']:
        for item in row['trajectories']:
            path=ROOT/'cache'/(row['sample_id']+'_rep'+str(item['replica_index'])+'.npz')
            assert digest(path)==complete['ca_and_chi_cache_sha256'][path.name]
            with np.load(path,allow_pickle=False) as cache:
                assert cache['aligned_ca_angstrom'].shape==(501,len(cache['chi_count']),3)
                assert cache['chi_degrees'].shape==(501,len(cache['chi_count']),4)
                np.testing.assert_array_equal(cache['selected_original_indices'],item['selected_indices'])
                np.testing.assert_allclose(cache['times_ps'],[item['times_ps'][i] for i in item['selected_indices']],rtol=0,atol=.05)
                assert cache['analysis_plan_sha256'].item()==complete['plan_sha256']
                # First/median/last eligible residue in every trajectory; independent
                # dihedral-bin entropy, saved-label and CA variance reproduction.
                eligible=np.flatnonzero(cache['chi_valid'])
                for ri in np.unique(eligible[[0,len(eligible)//2,-1]]):
                    q=int(cache['chi_count'][ri])
                    angles=cache['chi_degrees'][:,ri,:q]
                    for offset in (0,30,60,90):
                        state=bins(angles,offset)
                        code=(state*(3**np.arange(q))).sum(axis=1)
                        h=entropy(code)
                        actual=indexed.loc[(row['sample_id'],item['replica_index'],ri,offset),'h_joint']
                        assert abs(h-actual)<=1e-12
                        if offset==0:
                            assert abs(h-cache['uniform_hchi'][ri])<=1e-12
                    ca=cache['aligned_ca_angstrom'][:,ri]
                    rmsf=np.sqrt(np.sum((ca-ca.mean(axis=0))**2,axis=1).mean())
                    assert abs(rmsf-cache['uniform_rmsf'][ri])<=1e-12
                cache_reports.append({'sample_id':row['sample_id'],'replica_index':item['replica_index'],
                                      'cache_sha256':complete['ca_and_chi_cache_sha256'][path.name],'spot_residues':len(np.unique(eligible[[0,len(eligible)//2,-1]]))})
    result={'status':'passed','analysis_plan_sha256':complete['plan_sha256'],'verifier_code_sha256':digest(__file__),
            'records':83,'trajectories':203,'all_output_hashes_exact':True,'all_cache_hashes_exact':True,
            'rows':len(tables),'original_uniform_chi_valid':33493,'original_uniform_rmsf_valid':38599,
            'all_paired_masks_exact':True,'Pro_conventional_chi1_exclusions_explicit':True,
            'chi1_centers_and_circular_trans_seam_passed':True,
            'joint_total_variance_identity_max_abs_A2':float(joint_error.max()),
            'chi1_total_variance_identity_max_abs_A2':float(chi1_error.max()),'cache_spot_checks':cache_reports}
    (ROOT/'independent_verification_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='cache_spot_checks'},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
