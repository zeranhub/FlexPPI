"""Retrospective Q40 diagnostic from existing archived six-GLN angles.

Writes only this trajectory directory. Original 120-degree labels remain
unchanged; fine-bin occupancy measures are reported separately. No new MD,
training, inference, curve selection or inferential p-values.
"""
from __future__ import annotations
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,math
from collections import Counter
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
V27=BASE/'free_energy_curves_v27_20261014/case_curves'
V30=BASE/'rich_angle_figure_v30_20261014/experiment'
PREDICTION=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\neural_exports\evaluation\full\A003R.npz')
ORIGINS=(0.,1.25,2.5,3.75)


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):result.update(block)
    return result.hexdigest()


def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def write_json(path,value):
    Path(path).write_text(json.dumps(clean(value),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')


def entropy_counts(counts):
    counts=np.asarray(counts,dtype=float)
    p=counts[counts>0]/counts.sum()
    return float(-np.dot(p,np.log(p)))


def covering90(counts):
    counts=np.asarray(counts,dtype=int)
    ranked=np.sort(counts)[::-1]
    return int(np.searchsorted(np.cumsum(ranked),.9*counts.sum(),side='left')+1)


def fine_metrics(values):
    widths=np.zeros((3,4),float)
    entropies=np.zeros((3,4),float)
    for qi in range(3):
        for oi,origin in enumerate(ORIGINS):
            bins=np.floor(((values[:,qi].astype(float)+180-origin)%360)/5).astype(int)
            counts=np.bincount(bins,minlength=72)
            assert counts.sum()==len(values)
            widths[qi,oi]=5*((counts>0)&(counts>=counts.max()*np.exp(-1.))).sum()
            entropies[qi,oi]=entropy_counts(counts)
    return widths,entropies


def plan():
    HERE.mkdir(parents=True,exist_ok=True)
    old=json.loads((V27/'extraction_audit.json').read_text(encoding='utf8'))
    new=json.loads((V30/'extraction_audit.json').read_text(encoding='utf8'))
    assert old['status']==new['status']=='complete'
    sources={str(PREDICTION):'b8a764739993e3bce18dedb3cff493c50a1a52613a08baac16f1803c4edabb5a'}
    for name in ('angles_5000frames_all6GLN.csv','angles_501frames_all6GLN.csv','joint_entropy_original_five_windows.csv','selection_inventory_with_MD_outcomes.csv'):
        sources[str(V27/name)]=old['outputs'][name]['sha256']
    sources[str(V30/'all66_residue_summary.csv')]=new['outputs']['all66_residue_summary.csv']['sha256']
    for path,digest in sources.items():assert sha(path)==digest,path
    p=dict(protocol='Retrospective_six_GLN_original_joint_entropy_decomposition_v1',created_utc=datetime.now(timezone.utc).isoformat(),
        analysis_code_sha256=sha(__file__),source_hash_files=sources,
        motivation='User asks why Q40 joint Hchi prediction error is larger; descriptive diagnosis after the error was observed; no independent or prespecified hypothesis-testing claim',
        positions=['Q2','Q31','Q40','Q41','Q49','Q62'],
        original_joint_label='Three equal 120degree sectors per chi, [-180,-60),[-60,60),[60,180]; original float32 arithmetic and +180 final-sector rule; all three GLN chis form27 joint states',
        time_scopes='Original five nonoverlapping 1000frame windows0-99.9,100-199.9,...400-499.9ns. Report window arithmetic means/SD(ddof1)/range; also pooled5000 and common501(0-100ns,200ps) as separately labelled descriptive scopes',
        decomposition='Hchi1+Hchi2+Hchi3-Hjoint = total correlation, nonnegative descriptive dependence among coarse states; not a causal mechanical coupling energy and not equal to pairwise mutual information sum',
        additional_measures='Each-angle original-sector probabilities, marginal Shannon entropies;27jointstate counts/probabilities; occupied states/minimum number of largest states covering>=90%; pairwise coarse-state mutual information',
        fine_bins='Separate 5degree marginal Shannon entropy and Fstar<=1 total occupied-bin angular width, four origins0,1.25,2.5,3.75 averaged; exact v27/v30 definition; never relabel as original120degree joint entropy',
        inference='No pvalues or confidence intervals from five temporally consecutive windows; describe variation and error direction only; no new MD/training/model inference',
        source_preservation='All existing data/model files read-only; all source SHA256/size/mtime verified before and after; write only trajectory subdirectory')
    dest=HERE/'diagnostic_plan.json'
    if dest.exists():
        prior=json.loads(dest.read_text(encoding='utf8'))
        assert prior['analysis_code_sha256']==sha(__file__) and prior['source_hash_files']==sources
        p=prior
    else:write_json(dest,p)
    return p


def run(p):
    before={path:(Path(path).stat().st_size,Path(path).stat().st_mtime_ns) for path in p['source_hash_files']}
    source=pd.read_csv(V27/'angles_5000frames_all6GLN.csv')
    small=pd.read_csv(V27/'angles_501frames_all6GLN.csv')
    original=pd.read_csv(V27/'joint_entropy_original_five_windows.csv')
    oldsummary=pd.read_csv(V27/'selection_inventory_with_MD_outcomes.csv').set_index('residue')
    newsummary=pd.read_csv(V30/'all66_residue_summary.csv').set_index('residue')
    prediction=np.load(PREDICTION,allow_pickle=False)
    result_rows,marginal_rows,joint_rows,pair_rows,fine_rows,summary_rows=[],[],[],[],[],[]
    max_original_error=0.
    for site in p['positions']:
        raw=source[source.residue==site].sort_values('frame_zero_based')
        assert len(raw)==5000
        np.testing.assert_array_equal(raw.frame_zero_based.to_numpy(),np.arange(5000))
        values=raw[[f'chi{qi}_degrees' for qi in (1,2,3)]].to_numpy(np.float32)
        short=small[small.residue==site].sort_values('frame_zero_based')
        np.testing.assert_allclose(values[::2][:501],short[[f'chi{qi}_degrees' for qi in (1,2,3)]].to_numpy(),atol=0,rtol=0)
        pi=int(oldsummary.loc[site,'residue_index'])
        assert prediction['chain'][pi]=='A' and prediction['amino_acid'][pi]=='Q' and prediction['valid'][pi,1]
        pred=float(prediction['mean'][pi,1])
        scopes=[(f'window{wi+1}',wi+1,values[wi*1000:(wi+1)*1000]) for wi in range(5)]
        scopes += [('pooled5000',0,values),('common501',-1,values[::2][:501])]
        site_results=[]
        for scope,wi,angles in scopes:
            sector=np.clip(np.floor((angles+np.float32(180))/np.float32(120)).astype(int),0,2)
            code=(sector*(3**np.arange(3))).sum(axis=1)
            joint=np.bincount(code,minlength=27)
            hjoint=entropy_counts(joint)
            hsingle=[]
            marginals=[]
            for qi in range(3):
                counts=np.bincount(sector[:,qi],minlength=3)
                h=entropy_counts(counts)
                hsingle.append(h)
                marginals.append(counts)
                for si in range(3):
                    marginal_rows.append(dict(residue=site,scope=scope,window=wi,n_frames=len(angles),chi=f'chi{qi+1}',sector=si,
                        angle_lower_degrees=-180+120*si,angle_upper_degrees=-60+120*si,
                        count=int(counts[si]),probability=float(counts[si]/len(angles)),marginal_entropy_nats=h,
                        marginal_occupied_states=int((counts>0).sum()),marginal_states_cover90=covering90(counts)))
            totalcorr=float(sum(hsingle)-hjoint)
            assert totalcorr>=-1e-12
            for j in range(27):
                bins=(j%3,(j//3)%3,(j//9)%3)
                joint_rows.append(dict(residue=site,scope=scope,window=wi,n_frames=len(angles),joint_code=j,
                    chi1_sector=bins[0],chi2_sector=bins[1],chi3_sector=bins[2],count=int(joint[j]),probability=float(joint[j]/len(angles))))
            pairmis=[]
            for a,b in ((0,1),(0,2),(1,2)):
                pair=np.bincount(sector[:,a]+3*sector[:,b],minlength=9)
                mutual=hsingle[a]+hsingle[b]-entropy_counts(pair)
                assert mutual>=-1e-12
                pairmis.append(mutual)
                pair_rows.append(dict(residue=site,scope=scope,window=wi,chi_a=f'chi{a+1}',chi_b=f'chi{b+1}',pair_mutual_information_nats=mutual))
            observed_independence_entropy=sum(hsingle)
            maximum=3*math.log(3)
            result=dict(residue=site,scope=scope,window=wi,n_frames=len(angles),predicted_Hchi_nats=pred,
                joint_Hchi_nats=hjoint,Hchi1_marginal_nats=hsingle[0],Hchi2_marginal_nats=hsingle[1],Hchi3_marginal_nats=hsingle[2],
                marginal_entropy_sum_nats=observed_independence_entropy,total_correlation_nats=totalcorr,
                theoretical_max_27states_nats=maximum,independent_marginal_deficit_from_max_nats=maximum-observed_independence_entropy,
                total_deficit_from_27state_max_nats=maximum-hjoint,
                joint_occupied_states=int((joint>0).sum()),joint_states_cover90=covering90(joint),
                largest_joint_state_probability=float(joint.max()/len(angles)),
                predicted_minus_MD_nats=pred-hjoint,prediction_error_direction='over' if pred>hjoint else 'under' if pred<hjoint else 'equal')
            result_rows.append(result)
            if wi>0:
                expected=float(original[(original.residue==site)&(original.window==wi)].Hchi_nats.iloc[0])
                delta=abs(hjoint-expected)
                assert delta<=1e-12
                max_original_error=max(max_original_error,delta)
                site_results.append(result)
            widths,fineH=fine_metrics(angles)
            for qi in range(3):
                fine_rows.append(dict(residue=site,scope=scope,window=wi,n_frames=len(angles),chi=f'chi{qi+1}',
                    fine_entropy_origin0_nats=float(fineH[qi,0]),fine_entropy_four_origins_mean_nats=float(fineH[qi].mean()),
                    width_origin0_degrees=float(widths[qi,0]),width_four_origins_mean_degrees=float(widths[qi].mean())))
            if scope=='pooled5000':np.testing.assert_allclose(widths.mean(),oldsummary.loc[site,'MD_width_5000frames_deg'],atol=1e-12,rtol=0)
            if scope=='common501':
                np.testing.assert_allclose(widths.mean(),oldsummary.loc[site,'MD_width_501frames_deg'],atol=1e-12,rtol=0)
                np.testing.assert_allclose(fineH.mean(),newsummary.loc[site,'mean_per_chi_fine_entropy_four_origins_nats'],atol=1e-12,rtol=0)
        table=pd.DataFrame(site_results)
        means={name:float(table[name].mean()) for name in ('joint_Hchi_nats','Hchi1_marginal_nats','Hchi2_marginal_nats','Hchi3_marginal_nats',
            'marginal_entropy_sum_nats','total_correlation_nats','independent_marginal_deficit_from_max_nats','total_deficit_from_27state_max_nats',
            'joint_occupied_states','joint_states_cover90','largest_joint_state_probability')}
        mean=means['joint_Hchi_nats']
        np.testing.assert_allclose(mean,prediction['target'][pi,1],atol=1e-6,rtol=0)
        np.testing.assert_allclose(table.joint_Hchi_nats.var(ddof=1),prediction['label_variance'][pi,1],atol=1e-6,rtol=0)
        summary_rows.append(dict(residue=site,residue_index=pi,predicted_Hchi_nats=pred,frozen_MD_joint_Hchi_nats=float(prediction['target'][pi,1]),
            **{'mean_window_'+name:value for name,value in means.items()},
            MD_window_sample_SD_nats=float(table.joint_Hchi_nats.std(ddof=1)),MD_window_min_nats=float(table.joint_Hchi_nats.min()),
            MD_window_max_nats=float(table.joint_Hchi_nats.max()),MD_window_range_nats=float(table.joint_Hchi_nats.max()-table.joint_Hchi_nats.min()),
            predicted_minus_MD_mean_nats=pred-mean,prediction_absolute_error_nats=abs(pred-mean),
            windows_overestimated=int((table.predicted_minus_MD_nats>0).sum()),windows_underestimated=int((table.predicted_minus_MD_nats<0).sum()),
            error_same_direction_in_all_windows=bool((table.predicted_minus_MD_nats>0).all() or (table.predicted_minus_MD_nats<0).all()),
            prediction_within_observed_window_min_max=bool(table.joint_Hchi_nats.min()<=pred<=table.joint_Hchi_nats.max()),
            deficit_fraction_from_marginal_concentration=means['independent_marginal_deficit_from_max_nats']/means['total_deficit_from_27state_max_nats'],
            deficit_fraction_from_coarse_state_dependence=means['total_correlation_nats']/means['total_deficit_from_27state_max_nats']))
    metrics=pd.DataFrame(fine_rows)
    finecompare=[]
    for site,rows in metrics.groupby('residue',sort=False):
        long=rows[rows.scope=='pooled5000']
        short=rows[rows.scope=='common501']
        finecompare.append(dict(residue=site,MD_fine_entropy_501_four_origins_mean_nats=float(short.fine_entropy_four_origins_mean_nats.mean()),
            MD_fine_entropy_5000_four_origins_mean_nats=float(long.fine_entropy_four_origins_mean_nats.mean()),
            MD_width_501_four_origins_mean_degrees=float(short.width_four_origins_mean_degrees.mean()),
            MD_width_5000_four_origins_mean_degrees=float(long.width_four_origins_mean_degrees.mean())))
    output={}
    for name,rows in [('original120_entropy_by_window_and_scope.csv',result_rows),('original120_per_chi_state_probabilities.csv',marginal_rows),
        ('original120_all_joint_state_probabilities.csv',joint_rows),('original120_pairwise_mutual_information.csv',pair_rows),
        ('six_GLN_prediction_window_summary.csv',summary_rows),('fine5degree_per_angle_by_window_and_scope.csv',fine_rows),
        ('fine5degree_501_vs5000_summary.csv',finecompare)]:
        path=HERE/name
        pd.DataFrame(rows).to_csv(path,index=False)
        output[name]=dict(sha256=sha(path),rows=len(rows),bytes=path.stat().st_size)
    for path,digest in p['source_hash_files'].items():
        assert sha(path)==digest
        assert (Path(path).stat().st_size,Path(path).stat().st_mtime_ns)==before[path]
    audit=dict(status='complete',analysis_plan_sha256=sha(HERE/'diagnostic_plan.json'),analysis_code_sha256=sha(__file__),
        source_SHA256_size_mtime_unchanged=True,all6GLN_all5_original_windows_verified=True,maximum_original_window_Hchi_error=max_original_error,
        all6_original_joint_mean_and_sample_variance_match_frozen_cache_atol=1e-6,
        all6_common501_and_pooled5000_widths_match_v27_atol=1e-12,all6_common501_fineentropy_match_v30_atol=1e-12,
        outputs=output,summary=summary_rows,fine_comparison=finecompare)
    write_json(HERE/'diagnostic_audit.json',audit)
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print(pd.DataFrame(finecompare).to_string(index=False))


if __name__=='__main__':run(plan())
