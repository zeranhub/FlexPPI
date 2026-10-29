"""Transparent derived summaries; the frozen analysis outputs are unchanged."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for data in iter(lambda:handle.read(4*1024*1024),b''):
            result.update(data)
    return result.hexdigest()


def main():
    audit=json.loads((ROOT/'completion_audit.json').read_text(encoding='utf-8'))
    assert audit['status']=='complete'
    input_path=ROOT/'per_residue_trajectory.csv.gz'
    assert sha(input_path)==audit['output_tables']['per_residue_trajectory']['sha256']
    frame=pd.read_csv(input_path,dtype={'sample_id':str,'chain':str,'author_residue_id':str,'amino_acid':str,'residue_name':str},low_memory=False)
    nonpro=frame[(frame.chi1_conventional_eligible==1)].copy()
    metrics=['h_chi1','chi1_p0','chi1_p1','chi1_p2','apparent_chi1_changes','chi1_boundary_10deg_fraction',
             'half_chi1_jsd','ca_between_chi1_fraction','conventional_chi1_p_gminus','conventional_chi1_p_gplus','conventional_chi1_p_trans']
    keys=['sample_id','dataset','fold','replica_index','shift_degrees']
    record=nonpro.groupby(keys)[metrics].mean().reset_index()
    count=nonpro.groupby(keys).size().reset_index(name='nonpro_chi1_residues')
    record=record.merge(count,on=keys,validate='one_to_one')
    record.to_csv(ROOT/'nonpro_chi1_per_record_trajectory.csv',index=False)
    primary=record[record.replica_index==1].copy().assign(replica_mode='replica1').drop(columns='replica_index')
    available=record.groupby(['sample_id','dataset','fold','shift_degrees'])[metrics+['nonpro_chi1_residues']].mean().reset_index().assign(replica_mode='available_replica_mean')
    macro=pd.concat([primary,available],ignore_index=True)
    macro.to_csv(ROOT/'nonpro_chi1_record_macro.csv',index=False)
    combined=pd.concat([macro,macro.copy().assign(dataset='all')])
    summary=combined.groupby(['replica_mode','dataset','shift_degrees'])[metrics+['nonpro_chi1_residues']].agg(['mean','std','count'])
    summary.columns=[a+'_'+b for a,b in summary.columns]
    summary.reset_index().to_csv(ROOT/'nonpro_chi1_source_summary.csv',index=False)
    excluded=frame[(frame.shift_degrees==60)&(frame.chi1_conventional_eligible==0)][
        ['sample_id','dataset','fold','replica_index','residue_index','chain','author_residue_id','residue_name','chi1_conventional_exclusion_reason']]
    excluded.to_csv(ROOT/'conventional_chi1_exclusions.csv',index=False)
    repeat=pd.read_csv(ROOT/'repeat_jsd_record_macro.csv')
    combined=pd.concat([repeat,repeat.copy().assign(dataset='all')])
    summary_repeat=combined.groupby(['dataset','shift_degrees'])[['joint_repeat_jsd','chi1_repeat_jsd']].agg(['mean','std','count'])
    summary_repeat.columns=[a+'_'+b for a,b in summary_repeat.columns]
    summary_repeat.reset_index().to_csv(ROOT/'repeat_jsd_source_summary_flat.csv',index=False)
    # Means and errors for manuscript-ready long table. Units retained explicitly.
    distribution=pd.read_csv(ROOT/'distribution_source_summary.csv')
    prediction=pd.read_csv(ROOT/'prediction_effect_source_summary.csv')
    selected=distribution[(distribution.replica_mode=='replica1')]
    manuscript=[]
    for _,row in selected.iterrows():
        match=prediction[(prediction.variant=='uniform_rep1')&(prediction.dataset==row.dataset)&
                         (prediction.metric=='chi_entropy')&(prediction.shift_degrees==row.shift_degrees)].iloc[0]
        manuscript.append(dict(dataset=row.dataset,records=int(row.h_joint_count),shift_degrees=int(row.shift_degrees),
            joint_H_mean_nats=row.h_joint_mean,joint_H_record_sd_nats=row.h_joint_std,
            chi1_equalwidth_H_mean_nats=row.h_chi1_mean,chi1_boundary10_fraction=row.chi1_boundary_10deg_fraction_mean,
            apparent_chi1_changes_per501_mean=row.apparent_chi1_changes_mean,
            half_joint_JSD_mean_nats=row.half_joint_jsd_mean,
            full_H_minus_weighted_window_H_mean_nats=row.full_minus_window_h_weighted_mean,
            joint_CA_between_fraction_mean=row.ca_between_joint_fraction_mean,
            fixed_OOF_MAE_original_label_nats=match.old_mae_mean,
            fixed_OOF_MAE_uniform_offset_label_nats=match.new_mae_mean))
    pd.DataFrame(manuscript).to_csv(ROOT/'manuscript_summary_table.csv',index=False)
    schema={
        'cache':{'pattern':'cache/{sample_id}_rep{replica_index}.npz','chi_degrees':'float32 [501,N,4]; degrees, undefined/missing chis NaN',
            'aligned_ca_angstrom':'float64 [501,N,3]; angstrom; source-chain Kabsch to first selected frame; each chain centered at origin',
            'chi_count':'int8 [N]; number of chemically defined chi angles',
            'chi_atom_available':'bool [N,4]; individual quadruplet availability',
            'chi_valid/ca_valid':'bool [N]; inherited target/frame validity and mapped geometry',
            'times_ps':'float64 [501]; 0..100000 ps at200 ps',
            'selected_original_indices':'int [501]; actual original frame indices',
            'source_residue_indices':'int [N]; source-to-canonical mapping, -1 if missing',
            'chain/author_residue_id/amino_acid/residue_name':'Unicode [N]; canonical identities; read CSV author IDs as string',
            'uniform_rmsf/uniform_hchi':'float64 [N]; exact source of reproduced uniform labels, Å/nats',
            'analysis_plan_sha256/raw_trajectory_sha256/topology_sha256':'Unicode scalar provenance'},
        'CSV':{
            'per_residue_trajectory.csv.gz':'350844 rows =87711 residue-replica observations x4 shifts; all joint/chi1 descriptors, paired frozen OOF errors, state-conditioned alignedCA variance; Pro explicitly flagged for conventional chi1',
            'repeat_jsd_per_residue.csv.gz':'325308 rows; all within-record replica pairs for ATLAS/ITrace, all4 shifts; DynaRepo does not have replica JSD',
            'per_record_trajectory.csv':'812 rows =203 trajectories x4 shifts; arithmetic residue mean within each trajectory',
            'distribution_record_macro.csv':'664 rows =83records x4shifts x2replica modes; replica1 primary and available-replica descriptor mean supplementary',
            'distribution_source_summary.csv':'record mean/sampleSD/count by source/all, shift and replica mode',
            'prediction_effect_per_record.csv':'1328rows =83records x4shifts x2replica modes x2targets; target rmsf0/chi1 indices explicit; unchanged frozen OOF mean re-scored',
            'valid_mask_audit.csv':'paired residue-count intersections and old/new-only counts for each scoring row',
            'nonpro_chi1_record_macro.csv':'Pro excluded from conventional chi1 interpretation; all4shifts still equal-width descriptions; physically named g-/g+/trans probabilities only populated atshift60',
            'conventional_chi1_exclusions.csv':'One row per Pro residue-replica; explicit exclusion reason',
            'repeat_jsd_source_summary_flat.csv':'flat single-row-header version of repeated JSD summaries',
            'manuscript_summary_table.csv':'Primary replica1 record macro by source/all and shift; nats and apparent counts explicitly named'}
    }
    (ROOT/'output_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    created=['nonpro_chi1_per_record_trajectory.csv','nonpro_chi1_record_macro.csv','nonpro_chi1_source_summary.csv',
             'conventional_chi1_exclusions.csv','repeat_jsd_source_summary_flat.csv','manuscript_summary_table.csv','output_schema.json']
    derivative={'status':'complete','frozen_plan_sha256':audit['plan_sha256'],'frozen_per_residue_input_sha256':sha(input_path),
        'code_sha256':sha(__file__),'derivation':'Only aggregation and explicit previously planned Pro exclusion; no new record selection, no labels/model edits',
        'outputs':{f:sha(ROOT/f) for f in created}}
    (ROOT/'delivery_summary_audit.json').write_text(json.dumps(derivative,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(pd.DataFrame(manuscript).to_string(index=False))


if __name__=='__main__':
    main()
