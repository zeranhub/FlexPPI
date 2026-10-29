"""Independent reconstruction from saved NPZ and complete source CSV tables.

Does not import extraction code. Uses Counter/math for histogram and metrics,
then checks all saved slots, including physically absent NaN padding and Pro.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,math
from collections import Counter
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


plan=json.loads((HERE/'analysis_plan.json').read_text(encoding='utf8'))
audit=json.loads((HERE/'extraction_audit.json').read_text(encoding='utf8'))
assert audit['status']=='complete' and audit['plan_sha256']==sha(HERE/'analysis_plan.json')
for path,digest in plan['source_hash_files'].items():assert sha(path)==digest
for name,item in audit['outputs'].items():assert sha(HERE/name)==item['sha256']
z=np.load(HERE/'ubiquitin_all66_chi_501frames.npz',allow_pickle=False)
prediction=np.load(audit['prediction_source'],allow_pickle=False)
indices=np.flatnonzero((prediction['chain']=='A')&prediction['valid'][:,1])
np.testing.assert_array_equal(z['residue_index'],indices)
assert len(indices)==66
for name in ('chain','author_residue_id','amino_acid'):
    np.testing.assert_array_equal(z[name],prediction[name][indices])
np.testing.assert_allclose(z['predicted_Hchi_nats'],prediction['mean'][indices,1],rtol=0,atol=0)
np.testing.assert_allclose(z['frozen_MD_joint_Hchi_nats'],prediction['target'][indices,1],rtol=0,atol=0)
np.testing.assert_array_equal(z['source_frame_indices'],np.arange(0,1001,2))
np.testing.assert_allclose(z['times_ps'],np.arange(501)*200,rtol=0,atol=.05)
defined=z['chi_defined']
angles=z['chi_degrees']
assert angles.shape==(501,66,4) and defined.shape==(66,4)
assert defined.sum()==151 and (~defined).sum()==113
assert np.isfinite(angles[:,defined]).all() and np.isnan(angles[:,~defined]).all()
np.testing.assert_array_equal(defined.sum(axis=1),z['chi_count'])
for name in ('counts_origin0','probability_origin0','Fstar_origin0'):
    assert np.isnan(z[name][~defined]).all()
for name in ('width_per_origin_degrees','fine_entropy_per_origin_nats'):
    assert np.isnan(z[name][~defined]).all()

hist=pd.read_csv(HERE/'all66_histograms_origin0.csv')
metrics=pd.read_csv(HERE/'all66_per_angle_metrics.csv')
long=pd.read_csv(HERE/'all66_angles_long_501frames.csv')
summary=pd.read_csv(HERE/'all66_residue_summary.csv')
assert len(hist)==19008 and len(metrics)==264 and len(long)==132264 and len(summary)==66
checked_valid_hist_bins=0
max_Fstar_error=0.
max_metric_error=0.
checked_angle_values=0
for si,site in enumerate(z['residue']):
    sr=summary.iloc[si]
    assert sr.residue==site and bool(sr.all_defined_chi_valid)
    for qi in range(4):
        mr=metrics.iloc[si*4+qi]
        assert mr.residue==site and mr.chi=='chi'+str(qi+1) and bool(mr.chi_defined)==defined[si,qi]
        hr=hist.iloc[(si*4+qi)*72:(si*4+qi+1)*72]
        ar=long.iloc[(si*4+qi)*501:(si*4+qi+1)*501]
        assert (hr.residue==site).all() and (hr.chi=='chi'+str(qi+1)).all()
        assert (ar.residue==site).all() and (ar.chi=='chi'+str(qi+1)).all()
        np.testing.assert_array_equal(ar.frame_zero_based.to_numpy(),z['source_frame_indices'])
        np.testing.assert_allclose(ar.time_ps.to_numpy(),z['times_ps'],atol=0,rtol=0)
        np.testing.assert_allclose(ar.angle_degrees.to_numpy(),angles[:,si,qi].astype(float),rtol=0,atol=3e-14,equal_nan=True)
        checked_angle_values+=len(ar)
        if not defined[si,qi]:
            assert ar.angle_degrees.isna().all()
            assert hr[['count','probability','density_per_degree','F_star','low_Fstar_le1']].isna().all().all()
            valuecols=[c for c in metrics if c.startswith(('width_','fine_entropy_'))]
            assert mr[valuecols].isna().all()
            continue
        widths,entropies=[],[]
        for oi,origin in enumerate(z['bin_origins_degrees']):
            counter=Counter(math.floor((float(x)+180-float(origin))%360/5) for x in angles[:,si,qi])
            assert sum(counter.values())==501 and all(0<=b<72 for b in counter)
            maximum=max(counter.values())
            width=5*sum(math.log(maximum/n)<=1 for n in counter.values())
            entropy=math.fsum(-(n/501)*math.log(n/501) for n in counter.values())
            widths.append(width)
            entropies.append(entropy)
            assert z['width_per_origin_degrees'][si,qi,oi]==width
            error=abs(z['fine_entropy_per_origin_nats'][si,qi,oi]-entropy)
            max_metric_error=max(max_metric_error,error)
            assert error<1e-12
            if oi==0:
                for b in range(72):
                    count=counter.get(b,0)
                    row=hr.iloc[b]
                    assert int(row['count'])==count==z['counts_origin0'][si,qi,b]
                    probability=count/501
                    np.testing.assert_allclose([row.probability,z['probability_origin0'][si,qi,b]],probability,atol=1e-15,rtol=0)
                    assert abs(row.density_per_degree-probability/5)<1e-15
                    if count:
                        fstar=math.log(maximum/count)
                        np.testing.assert_allclose([row.F_star,z['Fstar_origin0'][si,qi,b]],fstar,atol=1e-14,rtol=0)
                        max_Fstar_error=max(max_Fstar_error,abs(row.F_star-fstar))
                        assert bool(row.low_Fstar_le1)==(fstar<=1)
                    else:
                        assert np.isnan(row.F_star) and np.isnan(z['Fstar_origin0'][si,qi,b])
                        assert row.probability==0 and row['count']==0 and not bool(row.low_Fstar_le1)
                    checked_valid_hist_bins+=1
        for column,expected in [('width_origin0_degrees',widths[0]),('width_four_origins_mean_degrees',math.fsum(widths)/4),
            ('fine_entropy_origin0_nats',entropies[0]),('fine_entropy_four_origins_mean_nats',math.fsum(entropies)/4)]:
            err=abs(mr[column]-expected)
            max_metric_error=max(max_metric_error,err)
            assert err<1e-12
    for column,source in [('mean_per_chi_width_origin0_degrees',z['width_per_origin_degrees'][si,:,0]),
        ('mean_per_chi_width_four_origins_degrees',z['width_per_origin_degrees'][si]),
        ('mean_per_chi_fine_entropy_origin0_nats',z['fine_entropy_per_origin_nats'][si,:,0]),
        ('mean_per_chi_fine_entropy_four_origins_nats',z['fine_entropy_per_origin_nats'][si])]:
        np.testing.assert_allclose(sr[column],np.nanmean(source),rtol=0,atol=1e-12)

pro=summary[summary.amino_acid=='P']
assert pro.residue.tolist()==['P19','P37','P38']
assert pro.frozen_MD_joint_Hchi_nats.eq(0).all()
assert pro.chi_count.eq(2).all() and pro.all_defined_chi_valid.all()
previous=HERE.parents[1]/'free_energy_curves_v27_20261014/case_curves'
old=pd.read_csv(previous/'angles_501frames_all6GLN.csv')
max_old_error=0.
for site,group in old.groupby('residue',sort=False):
    si=z['residue'].tolist().index(site)
    observed=angles[:,si,:3].astype(float)
    expected=group[[f'chi{i}_degrees' for i in (1,2,3)]].to_numpy(float)
    error=abs((observed-expected+180)%360-180)
    assert error.max()<1e-5
    max_old_error=max(max_old_error,float(error.max()))
historic=pd.read_csv(previous/'histograms_501frames.csv')
joined=historic.merge(hist,on=['residue','chi','angle_center_degrees'],suffixes=('_v27','_all66'),validate='one_to_one')
assert len(joined)==432
for column in ('count','probability','density_per_degree','F_star'):
    np.testing.assert_allclose(joined[column+'_v27'],joined[column+'_all66'],rtol=0,atol=1e-14,equal_nan=True)

result=dict(status='passed',verification_code_sha256=sha(__file__),source_files_sha256_match=True,
    extraction_audit_sha256=sha(HERE/'extraction_audit.json'),valid_Hchi_sites=66,defined_angle_positions=151,
    per_chi_positions=defined.sum(axis=0).tolist(),absent_angle_slots=113,missing_never_zero_filled=True,
    angle_values_including_missing_checked=checked_angle_values,valid_histogram_bins_checked=checked_valid_hist_bins,
    width_entropy_origin_combinations_checked=151*4,maximum_Fstar_numeric_error=max_Fstar_error,maximum_metric_numeric_error=max_metric_error,
    old_six_GLN_angle_values_checked=501*6*3,old_six_GLN_maximum_circular_error_degrees=max_old_error,
    old_two_site_histogram_rows_checked=432,all_three_zero_original_reference_Pro_preserved=True,
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in HERE.iterdir() if p.is_file() and p.name!='independent_verification.json'})
(HERE/'independent_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in result.items() if k!='files'},ensure_ascii=False,indent=2))
