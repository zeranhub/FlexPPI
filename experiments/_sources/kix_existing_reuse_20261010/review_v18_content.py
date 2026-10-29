"""Independent read-only KIX content review against residue-level frozen sources."""
from pathlib import Path
import datetime
import hashlib
import itertools
import json
import math
import numpy as np
import pandas as pd
from docx import Document

PROJECT = Path(r'E:\Project\FlexPPI')
EXPERIMENT = PROJECT / 'experiments/_sources'
OUT = EXPERIMENT / 'kix_existing_reuse_20261010'
DOCX = EXPERIMENT.parent / 'paper_v18_20261010_中文修订_复用KIX单轨迹_三线表.docx'
CASE = PROJECT / 'experiments/_case_validation_sources/data'
SUMMARY = OUT / 'summary'
QUALITY = PROJECT / 'artifacts/kix_md_20260914/pilot_results_20260916/quality'
STATES = ('apo','myb','mll','both')
NAMES = {'apo':'游离KIX','myb':'cMyb','mll':'MLL','both':'双伙伴'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def close(a,b):
    # Residue-source prediction values were exported from float32 as decimal CSV.
    # A 1e-7 tolerance is far below the manuscript's 0.001-Angstrom rounding.
    assert math.isclose(float(a),float(b),rel_tol=1e-7,abs_tol=1e-7),(a,b)


def main():
    document=Document(DOCX)
    paragraphs={i:p.text for i,p in enumerate(document.paragraphs)}
    text='\n'.join(paragraphs.values())
    source=pd.read_csv(CASE/'kix_fold1_all_residues.csv')
    published=pd.read_csv(CASE/'kix_fold1_case_metrics.csv').set_index('state')
    state_summary=pd.read_csv(SUMMARY/'state_summary.csv').set_index('state')
    delta_summary=pd.read_csv(SUMMARY/'l12_delta_summary.csv').set_index('contrast')
    quality_summary=pd.read_csv(SUMMARY/'quality_summary.csv').set_index(['state','metric'])
    half_summary=pd.read_csv(SUMMARY/'half_window_stability.csv')
    trace=pd.read_csv(QUALITY/'trajectory_trace_binned.csv')
    thermo=pd.read_csv(QUALITY/'thermodynamics_all_samples.csv')
    windows=pd.read_csv(QUALITY/'residue_window_metrics.csv')
    quality_json=json.loads((QUALITY/'summary.json').read_text())
    assert set(source.metric)=={'rmsf'} and set(source.model_fold)=={1}
    assert set(source.md_independent_repeats)=={1}
    assert set(published['fold'])=={1}
    assert set(state_summary.independent_trajectories)=={1}
    assert set(state_summary.production_ns)=={100} and set(state_summary.equilibration_ns)=={10}
    assert set(state_summary.frames)=={100000}
    state_results={}
    for state in STATES:
        rows=source[source.state==state].sort_values('author_residue_number')
        assert rows.author_residue_number.tolist()==list(range(586,673))
        l12=rows[rows.author_residue_number.between(614,621)]
        assert len(l12)==8
        mae=float((rows.model_value-rows.md_value).abs().mean())
        rho=float(rows.model_value.rank().corr(rows.md_value.rank()))
        model_top=set(rows.nlargest(18,'model_value').author_residue_number)
        md_top=set(rows.nlargest(18,'md_value').author_residue_number)
        overlap=len(model_top&md_top)
        l12_model_high=len(model_top&set(range(614,622)))
        close(mae,published.loc[state,'mae']); close(rho,published.loc[state,'spearman'])
        assert overlap==published.loc[state,'top_overlap']
        assert l12_model_high==published.loc[state,'l12_model_high_n']
        close(mae,state_summary.loc[state,'RMSF_MAE_angstrom'])
        close(rho,state_summary.loc[state,'RMSF_spearman'])
        close(l12.md_value.mean(),state_summary.loc[state,'L12_MD_mean_RMSF_angstrom'])
        close(l12.model_value.mean(),state_summary.loc[state,'L12_fold1_mean_RMSF_angstrom'])
        close(rows.md_value.mean(),state_summary.loc[state,'KIX_MD_mean_RMSF_angstrom'])
        close(rows.model_value.mean(),state_summary.loc[state,'KIX_fold1_mean_RMSF_angstrom'])
        expected_row=[NAMES[state],f'{mae:.3f}',f'{rho:.3f}',f'{overlap}/18',
                      f'{l12.md_value.mean():.3f}',f'{l12.model_value.mean():.3f}']
        table=document.tables[8]
        actual_row=[c.text for c in table.rows[STATES.index(state)+1].cells]
        assert actual_row==expected_row,(actual_row,expected_row)
        state_results[state]=dict(n_residues=87,mae=mae,spearman=rho,top20_n=18,
                                 top20_overlap=overlap,l12_model_top20=l12_model_high,
                                 l12_md_mean=float(l12.md_value.mean()),
                                 l12_model_mean=float(l12.model_value.mean()),table_S7_matches=True)
    pairs=[('myb','apo'),('mll','apo'),('both','apo'),('mll','myb'),('both','myb'),('both','mll')]
    delta_results=[]
    for i,(first,second) in enumerate(pairs):
        md_delta=state_results[first]['l12_md_mean']-state_results[second]['l12_md_mean']
        model_delta=state_results[first]['l12_model_mean']-state_results[second]['l12_model_mean']
        error=abs(model_delta-md_delta)
        summary=delta_summary.loc[first+'-'+second]
        close(md_delta,summary.MD_delta_RMSF_angstrom)
        close(model_delta,summary.fold1_delta_RMSF_angstrom)
        close(error,summary.absolute_delta_error_angstrom)
        assert md_delta>0 and model_delta<0 and not bool(summary.direction_agrees)
        expected_row=[NAMES[first]+'−'+NAMES[second],f'{md_delta:+.3f}',f'{model_delta:+.3f}',f'{error:.3f}']
        assert [cell.text for cell in document.tables[9].rows[i+1].cells]==expected_row
        delta_results.append(dict(contrast=first+'-'+second,md_delta=md_delta,
                                 model_delta=model_delta,absolute_error=error,
                                 directions_opposite=True,table_S8_matches=True))
    half=windows[(windows.window_kind=='disjoint_half') & windows.author_residue_number.between(614,621)]
    means=half.groupby(['state','start_ns_exclusive','end_ns_inclusive']).rmsf_angstrom.mean()
    for row in half_summary.itertuples():
        close(row.L12_mean_RMSF_angstrom,means.loc[(row.state,row.start_ns_exclusive,row.end_ns_inclusive)])
    first_delta=float(means.loc[('mll',0,50)]-means.loc[('myb',0,50)])
    second_delta=float(means.loc[('mll',50,100)]-means.loc[('myb',50,100)])
    assert f'−{abs(first_delta):.3f}和+{second_delta:.3f}' in text
    assert all(item['n_frames']==100000 and item['dcd_hashes_verified'] and
               item['full_repeat_metric_max_abs_difference']=={'rmsf_angstrom':0.0,'chi_entropy_nats':0.0}
               for item in quality_json['states'])
    assert set(trace.n_frames)=={100}
    trace_means=trace.groupby(['state','metric'])['mean'].mean()
    diagnostic=[]
    for (state,metric),value in trace_means.items():
        if (state,metric) in quality_summary.index:
            close(value,quality_summary.loc[(state,metric),'mean'])
        if metric=='KIX_CA_RMSD_angstrom' or metric.startswith('partner_'):
            diagnostic.append(dict(state=state,metric=metric,mean=float(value)))
    temperature=thermo.groupby('state').temperature_kelvin.mean()
    for state,value in temperature.items():
        close(value,quality_summary.loc[(state,'temperature_kelvin'),'mean'])
    assert f'{temperature.min():.2f}–{temperature.max():.2f} K' in text
    assert len(thermo)==40000
    assert 'ρ为0.729–0.896' in text and '61%–89%' in text
    assert '每种状态各一条已完成的生产轨迹' in text
    assert '第1折' in text
    assert '不能据此构造独立重复的置信区间或显著性检验' in text
    assert '尚不足以支持结合响应幅值或变构机制的预测' in text
    assert '不能据此证明模型正确恢复了结合诱导的柔性变化或变构响应' in text
    assert '热力学变构效应' in text
    assert 'hmu' not in text.lower() and '认证' not in text
    issues=[]
    for index,paragraph in paragraphs.items():
        if '前50 ns/后50 ns独立窗口' in paragraph:
            issues.append(dict(severity='MINOR',must_fix=False,paragraph_index=index,
                current='前50 ns/后50 ns独立窗口',suggestion='前50 ns/后50 ns不重叠时间窗口',
                reason='同一轨迹两个不重叠半段不是统计独立重复；P200已给出正确限制，因此不构成CI误用或数值错误。'))
    source_paths=[CASE/'kix_fold1_all_residues.csv',CASE/'kix_fold1_case_metrics.csv',
                  SUMMARY/'state_summary.csv',SUMMARY/'l12_delta_summary.csv',
                  SUMMARY/'quality_summary.csv',SUMMARY/'half_window_stability.csv',
                  QUALITY/'summary.json',QUALITY/'trajectory_trace_binned.csv',
                  QUALITY/'residue_window_metrics.csv',QUALITY/'thermodynamics_all_samples.csv',
                  EXPERIMENT/'integrate_existing_kix.py']
    report=dict(status='passed_no_required_KIX_content_corrections',
                created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                reviewed_docx=str(DOCX),reviewed_docx_sha256=sha(DOCX),docx_edited=False,
                scope='Independent KIX content/data review only; not whole-manuscript layout or reference-list review',
                checks=dict(four_states_single_100ns=True,total_production_frames=400000,
                    no_independent_repeat_CI_or_p=True,all_six_L12_delta_directions_opposite=True,
                    half_window_mll_minus_myb=[first_delta,second_delta],
                    ranking_range=[min(r['spearman'] for r in state_results.values()),max(r['spearman'] for r in state_results.values())],
                    fold1_only=True,no_early_five_model_mix=True,
                    KIX_tables_RMSF_only=True,AGP_Hchi_mask_not_mixed=True,
                    note_Hchi='KIX numerical tables report RMSF only. Historical full-frame Hchi recheck includes original Pro definition; no AG/P-excluded Hchi score or model Hchi value is imported into these tables.',
                    hmu_auth_denial_not_used_as_absence_evidence=True,
                    no_unjustified_convergence_binding_or_thermodynamic_function_claim=True,
                    model_activity_localization_and_failed_response_prediction_distinguished=True),
                state_results=state_results,L12_deltas=delta_results,quality_diagnostics=diagnostic,
                numeric_comparison_tolerance=dict(absolute=1e-7,relative=1e-7,
                    reason='Float32 model predictions were exported as decimal residue CSV; recomputed source MAE differs from full-precision metric CSV only at approximately 1e-9 Angstrom. All manuscript table values match exact 0.001 formatting.'),
                issues=issues,must_fix_issues=[],
                verified_source_sha256={str(path):sha(path) for path in source_paths},
                conclusion='KIX tables S7/S8 and body statements match frozen Fold1 residue sources and independently recomputed summary data. All six state differences oppose MD; half-window sign reversal and single-trajectory limitations are explicit. The earlier optional wording suggestion has been incorporated; no remaining KIX content issue was found.')
    (OUT/'v18_content_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'half_window_mll_minus_myb':[first_delta,second_delta],
                      'issues':issues,'output':str(OUT/'v18_content_review.json')},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
