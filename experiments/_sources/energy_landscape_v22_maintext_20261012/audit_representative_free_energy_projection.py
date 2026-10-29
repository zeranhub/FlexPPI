"""Read-only v21 input audit; all new products written only to v22."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
V21=OUT.parent/'energy_landscape_v21_20261012'
CASE=V21/'case_analysis'
TRAJ=V21/'trajectory_analysis'
SOURCE_CASES={'2py5_A':{'Y254':253,'Y521':520},'A0022':{'I78':77,'I21':20},
              'itrace_6p64':{'T149':721,'T25':408}}


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):
            result.update(block)
    return result.hexdigest()


def sectors(angles):
    return np.floor((((angles.astype(float)-60)+180)%360)/120).astype(int)


def connected_circular_regions(centers,counts,fstar,threshold=1):
    selected=np.isfinite(fstar)&(fstar<=threshold)
    visited=set()
    regions=[]
    for index in np.flatnonzero(selected):
        if index in visited:
            continue
        stack=[int(index)]
        component=[]
        while stack:
            current=stack.pop()
            if current in visited or not selected[current]:
                continue
            visited.add(current)
            component.append(current)
            stack.extend([(current-1)%len(centers),(current+1)%len(centers)])
        component=sorted(component)
        wrap=0 in component and len(centers)-1 in component
        if wrap:
            positive=[i for i in component if centers[i]>0]
            negative=[i for i in component if centers[i]<0]
            lower=float(centers[min(positive)]-2.5)
            upper=float(centers[max(negative)]+2.5)
        else:
            lower=float(centers[min(component)]-2.5)
            upper=float(centers[max(component)]+2.5)
        regions.append(dict(angle_lower_degrees=lower,angle_upper_degrees=upper,circular_seam_wrapped=wrap,
            observed_bins=len(component),observed_frames=int(counts[component].sum()),
            occupancy_probability=float(counts[component].sum()/counts.sum()),
            bin_centers_degrees=[float(centers[i]) for i in component]))
    return regions


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    files=[CASE/'representative_cases_selection.json',CASE/'representative_cases_audit.json',CASE/'representative_chi1_Fstar.csv',
           CASE/'representative_chi1_timeseries.csv',CASE/'representative_selected_residues.csv',TRAJ/'analysis_plan.json',
           TRAJ/'run_angle_distribution_analysis.py',TRAJ/'completion_audit.json',TRAJ/'per_residue_trajectory.csv.gz',
           TRAJ/'distribution_record_macro.csv',TRAJ/'distribution_source_summary.csv']
    before={str(p):sha(p) for p in files}
    metadata=json.loads((CASE/'representative_cases_audit.json').read_text(encoding='utf-8'))
    frozen=json.loads((TRAJ/'completion_audit.json').read_text(encoding='utf-8'))
    assert metadata['status']==frozen['status']=='complete'
    assert sha(CASE/'representative_cases_selection.json')==metadata['selection_sha256']
    histogram=pd.read_csv(CASE/'representative_chi1_Fstar.csv')
    times=pd.read_csv(CASE/'representative_chi1_timeseries.csv')
    selected=pd.read_csv(CASE/'representative_selected_residues.csv',dtype={'author_residue_id':str,'chain':str})
    assert len(selected)==6 and len(histogram)==14*72 and len(times)==14*501
    edges=np.arange(-180,185,5)
    centers=(edges[:-1]+edges[1:])/2
    cache_paths={Path(x['path']).stem:x for x in metadata['cache_sources']}
    verifications=[]
    residue_summary=[]
    components=[]
    for (sample,replica,residue),series in times.groupby(['sample_id','replica','residue']):
        assert sample in SOURCE_CASES and residue in SOURCE_CASES[sample]
        ri=SOURCE_CASES[sample][residue]
        assert len(series)==501 and series.input_index.eq(ri).all()
        series=series.sort_values('time_ns')
        np.testing.assert_allclose(series.time_ns.to_numpy(),np.arange(501)/5,rtol=0,atol=1e-12)
        entry=cache_paths[sample+'_rep'+str(replica)]
        assert sha(entry['path'])==entry['sha256']
        with np.load(entry['path'],allow_pickle=False) as cache:
            assert cache['chi_valid'][ri] and cache['ca_valid'][ri]
            assert cache['residue_name'][ri]!='PRO'
            assert cache['amino_acid'][ri]==residue[0] and cache['author_residue_id'][ri]==residue[1:]
            angle=cache['chi_degrees'][:,ri,0]
            np.testing.assert_allclose(series.chi1_degrees.to_numpy(),angle,rtol=0,atol=1e-12)
            chain=str(cache['chain'][ri])
        part=histogram[(histogram.sample_id==sample)&(histogram.replica==replica)&(histogram.residue==residue)].sort_values('chi1_degrees')
        assert len(part)==72 and part.input_index.eq(ri).all()
        np.testing.assert_array_equal(part.chi1_degrees.to_numpy(),centers)
        counts,_=np.histogram(angle,bins=edges)
        assert int(counts.sum())==501
        np.testing.assert_array_equal(part['count'].to_numpy(),counts)
        p=counts/501
        np.testing.assert_allclose(part.probability.to_numpy(),p,rtol=0,atol=1e-15)
        f=np.full(72,np.nan)
        f[counts>0]=-np.log(p[counts>0]/p.max())
        np.testing.assert_allclose(part.F_star.to_numpy(),f,rtol=0,atol=1e-14,equal_nan=True)
        assert part.loc[part['count']==0,'F_star'].isna().all()
        assert np.nanmin(f)==0
        verifications.append(dict(sample_id=sample,replica=int(replica),residue=residue,frames=501,
            cache_sha256=entry['sha256'],times_match_cache=True,all72histogram_counts_exact=True,
            probabilities_match_counts=True,all_visited_Fstar_match_minus_log_p_over_pmax=True,empty_bins_missing=True))
        if replica!=1:
            continue
        coding=sectors(angle)
        occ=np.bincount(coding,minlength=3)
        peaks=np.flatnonzero(counts==counts.max())
        low_regions=connected_circular_regions(centers,counts,f)
        old=selected[(selected.sample_id==sample)&(selected.input_index==ri)].iloc[0]
        row=dict(sample_id=sample,dataset=old.dataset,replica=1,residue=residue,chain=chain,input_index=ri,
            frames=501,window_start_ns=0,window_end_ns=100,interval_ps=200,
            chi1_gminus_frames=int(occ[0]),chi1_gplus_frames=int(occ[1]),chi1_trans_frames=int(occ[2]),
            chi1_gminus_probability=float(occ[0]/501),chi1_gplus_probability=float(occ[1]/501),chi1_trans_probability=float(occ[2]/501),
            maximum_5degree_bin_frames=int(counts.max()),maximum_5degree_bin_probability=float(p.max()),
            maximum_bin_centers_degrees=';'.join(format(centers[i],'g') for i in peaks),
            observed_5degree_bins=int((counts>0).sum()),low_Fstar_le1_circular_regions=len(low_regions),
            low_Fstar_le1_occupancy_probability=float(sum(r['occupancy_probability'] for r in low_regions)),
            chi1_conventional_sector_entropy_nats=float(-np.sum((occ[occ>0]/501)*np.log(occ[occ>0]/501))),
            reference_RMSF_A=float(old.rmsf_reference),predicted_RMSF_A=float(old.rmsf_prediction),
            reference_joint_Hchi_nats=float(old.hchi_reference),predicted_joint_Hchi_nats=float(old.hchi_prediction))
        for sector_index,label in enumerate(['gminus','gplus','trans']):
            indexes=np.flatnonzero(sectors(centers)==sector_index)
            if counts[indexes].max()>0:
                top=indexes[counts[indexes]==counts[indexes].max()]
                row[label+'_maximum_bin_centers_degrees']=';'.join(format(centers[i],'g') for i in top)
                row[label+'_maximum_bin_probability']=float(p[top[0]])
            else:
                row[label+'_maximum_bin_centers_degrees']=''
                row[label+'_maximum_bin_probability']=0.0
        residue_summary.append(row)
        for number,region in enumerate(low_regions,start=1):
            components.append(dict(sample_id=sample,residue=residue,chain=chain,replica=1,region=number,Fstar_threshold=1,
                **{k:v for k,v in region.items() if k!='bin_centers_degrees'},
                bin_centers_degrees=';'.join(format(c,'g') for c in region['bin_centers_degrees'])))
    assert len(verifications)==14 and len(residue_summary)==6
    pd.DataFrame(residue_summary).to_csv(OUT/'representative_rep1_chi1_occupancy_and_Fstar.csv',index=False,float_format='%.17g')
    pd.DataFrame(components).to_csv(OUT/'representative_rep1_low_Fstar_regions.csv',index=False,float_format='%.17g')
    ca_file=TRAJ/'per_residue_trajectory.csv.gz'
    assert sha(ca_file)==frozen['output_tables']['per_residue_trajectory']['sha256']
    fields=['sample_id','dataset','replica_index','residue_index','shift_degrees','chi_count','chi1_conventional_eligible',
            'ca_total_variance_A2','ca_within_joint_variance_A2','ca_between_joint_variance_A2','ca_between_joint_fraction']
    ca=pd.read_csv(ca_file,usecols=fields,low_memory=False)
    primary=ca[(ca.replica_index==1)&ca.shift_degrees.isin([0,60])].copy()
    assert len(primary)==33493*2
    mask0=primary[primary.shift_degrees==0].sort_values(['sample_id','residue_index'])
    mask60=primary[primary.shift_degrees==60].sort_values(['sample_id','residue_index'])
    np.testing.assert_array_equal(mask0[['sample_id','residue_index']].to_numpy(),mask60[['sample_id','residue_index']].to_numpy())
    np.testing.assert_allclose(mask0.ca_total_variance_A2,mask60.ca_total_variance_A2,rtol=0,atol=1e-12)
    assert primary.ca_between_joint_fraction.notna().all()
    identity_error=np.abs(primary.ca_total_variance_A2-primary.ca_within_joint_variance_A2-primary.ca_between_joint_variance_A2)
    assert (identity_error<=1e-8*np.maximum(primary.ca_total_variance_A2,1)).all()
    np.testing.assert_allclose(primary.ca_between_joint_fraction,primary.ca_between_joint_variance_A2/primary.ca_total_variance_A2,rtol=0,atol=1e-12)
    ca_record=primary.groupby(['sample_id','dataset','shift_degrees']).agg(
        valid_residues=('residue_index','size'),between_fraction_record_mean=('ca_between_joint_fraction','mean'),
        ca_total_variance_residue_sum=('ca_total_variance_A2','sum'),ca_between_variance_residue_sum=('ca_between_joint_variance_A2','sum')).reset_index()
    assert ca_record.groupby('shift_degrees').sample_id.nunique().to_dict()=={0:83,60:83}
    ca_record.to_csv(OUT/'CA_joint_state_decomposition_per_record.csv',index=False,float_format='%.17g')
    summaries=[]
    for shift,rows in primary.groupby('shift_degrees'):
        for source in ['all','atlas','dynarepo','itrace']:
            part=rows if source=='all' else rows[rows.dataset==source]
            rec=ca_record[(ca_record.shift_degrees==shift)&(ca_record.dataset.isin(part.dataset))]
            summaries.append(dict(dataset=source,shift_degrees=int(shift),records=int(rec.sample_id.nunique()),
                effective_residues=len(part),minimum_residues_per_record=int(rec.valid_residues.min()),
                maximum_residues_per_record=int(rec.valid_residues.max()),
                mean_of_record_residue_mean_between_fraction=float(rec.between_fraction_record_mean.mean()),
                record_sample_SD_between_fraction=float(rec.between_fraction_record_mean.std(ddof=1)),
                pooled_residue_mean_between_fraction=float(part.ca_between_joint_fraction.mean()),
                ratio_of_pooled_between_to_total_variance=float(part.ca_between_joint_variance_A2.sum()/part.ca_total_variance_A2.sum())))
    pd.DataFrame(summaries).to_csv(OUT/'CA_joint_state_decomposition_summary_and_denominators.csv',index=False,float_format='%.17g')
    for file,actual in before.items():
        assert sha(file)==actual,'v21sourcechanged during read-only audit'
    for file in [CASE/'representative_chi1_Fstar.csv',CASE/'representative_chi1_timeseries.csv']:
        assert sha(file)==before[str(file)]
    interpretation={
        '2py5_A':'2py5_A的两个Tyr位点分别集中在g−和g+占据区。Y254的χ1最大5°概率箱位于−57.5°，低F*区域横跨约−85°至−50°；Y521的主峰位于+67.5°，对应低F*区约+60°至+80°。两者主要差别包括占据区中心及局部宽度；原联合Hχ还同时包含其他χ角的信息。',
        'A0022':'A0022的I78主要处于g−占据区（95.2%），并存在少量g+占据（4.8%）；I21在本窗口内均处于g−。其χ1概率峰分别位于−57.5°和−72.5°，I78的低F*区约−70°至−40°，I21约−90°至−65°，呈现同一主要占据区内不同的局部位置和展宽。',
        'itrace_6p64':'ITrace 6p64的E链T149同时分布于g+和trans占据区（47.9%与51.9%），分别在约+47.5°和−177.5°呈现概率峰与两个分离的低F*区域；D链T25则集中在trans（100%），主峰位于−177.5°。圆周图上±180°两侧连成同一trans区域。T149较小的原参考RMSF与两个主要侧链占据区并存，T25较大的原参考RMSF与单一主要侧链占据区并存，具体体现两项读数所强调的不同运动特征。',
        'joint_CA_decomposition':'以每残基的联合χ离散状态作为条件，将链内对齐Cα的总均方波动分解为状态内波动和状态均值之间的波动。主分析83条记录、33493个同时具有有效联合χ与Cα坐标的位置中，先逐位计算状态均值间波动占总波动的比例，再对每条记录内残基取平均，最后对83条记录等权平均；原分区及整体平移60°等宽分区的结果分别为6.24%和7.77%。这两个比例描述与侧链状态相关的Cα位置差异，分区仅提供条件标签，未识别全蛋白自由能盆地。'}
    result=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),
        v21_inputs_sha256=before,source_files_unchanged_after_audit=True,
        source_case_selection='Three source records were frozen by original OOF two-target MAE rank; residue pairs frozen by same-AA reference label rank contrast; this audit does not select new cases',
        complete_histogram_time_cache_checks=verifications,
        Fstar=dict(definition='-ln(p_5degree_bin/p_max) for nonzero counts',dimensionless=True,zero_count_bins='missing, not interpolated',
                   histogram_bins=72,bin_width_degrees=5,frames_per_series=501,frame_interval_ps=200,first_replica_window_ns=[0,100],
                   peak='Highest observed5degree-bin probability; not a fitted continuous density mode',
                   low_region_description='Contiguous circular occupied bins with F*<=1, i.e.p>=pmax/e; descriptive summary threshold only',
                   physical_scope='Empirical relative occupancy projection along chi1; a physical free-energy reading requires adequate stationary equilibrium sampling and that temperature/coordinate definition; no activation barrier or true transition rate inferred'),
        residue_summary=residue_summary,low_Fstar_regions=components,short_maintext_chinese=interpretation,
        CA_decomposition=dict(total_variance='mean_t ||x_t-mean_t x||² after per-chainCA Kabsch to first sampled frame',
            within='sum_s p_s mean_{t in s} ||x_t-mean_s x||²',between='sum_s p_s ||mean_s x-overall_mean x||²',
            grouping='All defined legacy joint chi-state codes (includesPro); same angles/all4chis shifted together for60degrees; descriptive categorical conditioning',
            maximum_identity_error_A2=float(identity_error.max()),same_effective_residue_positions_for_0_and_60=True,
            denominator='33493primary-replica residues with valid jointχ andCA; atlas5660,dynarepo6384,itrace21449; all83records,52–735eligiblepositions/record; not the38599RMSF-only-validpositions',
            averaging='Per-residue between/total ratio -> arithmetic residue mean within record ->83record equal mean',
            source_summaries=summaries,
            interpretation='A small state-conditioned fraction ofCAmean-square variation, not a fractional free-energy barrier or causal sidechain-to-backbone contribution'),
        predictions='Original frozen model predicts scalarRMSF/jointHchi, not chi1 F* curve; these6positions do not establish complete recovery of all empirical landscape differences',
        new_csv_outputs={})
    for name in ['representative_rep1_chi1_occupancy_and_Fstar.csv','representative_rep1_low_Fstar_regions.csv',
                 'CA_joint_state_decomposition_per_record.csv','CA_joint_state_decomposition_summary_and_denominators.csv']:
        result['new_csv_outputs'][name]=sha(OUT/name)
    (OUT/'representative_Fstar_and_CA_maintext_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (OUT/'representative_Fstar_maintext_explanation_zh.md').write_text('\n\n'.join(interpretation.values())+'\n\nF*按5°等宽角箱的经验占据概率定义为−ln(p/pmax)，是χ1坐标的一维相对占据投影；低F*区指本窗口较常访问的角区间。峰位置对应概率最大的5°角箱，低F*范围仅按F*≤1的观测角箱汇总。原联合Hχ由所有定义χ角的联合状态计算，因此与χ1一维投影相互补充。\n',encoding='utf-8')
    print(pd.DataFrame(residue_summary).to_string(index=False))
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__=='__main__':
    main()
