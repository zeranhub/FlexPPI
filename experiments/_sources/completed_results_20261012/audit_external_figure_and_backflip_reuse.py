"""Independent read-only old/new external labels, official BF2 and fig4 audit."""
import ast
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
TASK=ROOT.parent
OLD=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
OLD_FIG=OLD.parents[1]/'figures/corrected_external_results_preparation/completed/20261003T163520Z_03a3b308'
NEW=TASK/'family_folds/external_2430_evaluation'
FIG=ROOT/'figures'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text(encoding='utf8'))
def csv(path):return pd.read_csv(path,keep_default_na=False,float_precision='round_trip',dtype={'chain':str,'author_residue_id':str})


class CapturedAxes:
    def __init__(self):self.lines=[];self.scatters=[];self.limits=None;self.scale='linear';self.tick_labels=[]
    def plot(self,x,y,**kwargs):self.lines.append((np.asarray(x),np.asarray(y,dtype=float)))
    def scatter(self,x,y,**kwargs):self.scatters.append((np.asarray(x),np.asarray(y,dtype=float),kwargs))
    def set_xticks(self,ticks,labels):self.tick_labels=labels
    def set_xlim(self,*args):pass
    def set_ylim(self,*args):self.limits=args
    def set_yscale(self,scale):self.scale=scale
    def set_yticks(self,*args):pass
    def set_ylabel(self,*args):pass


class CapturedFigure:
    def subplots_adjust(self,**kwargs):pass
    def text(self,*args,**kwargs):pass


def main():
    new_audit=read(NEW/'evaluation_audit.json')
    new=csv(NEW/'per_record_metrics.csv')
    old_all=csv(OLD_FIG/'all_32_per_record_metrics.csv')
    old_bf=old_all[(old_all.variant=='backflip2_official')&(old_all.metric=='rmsf')]
    bf=old_bf[old_bf.protocol=='strict_single_seed7'].sort_values('sample_id').reset_index(drop=True)
    current_rmsf=new[(new.variant=='full')&(new.metric=='rmsf')].sort_values('sample_id').reset_index(drop=True)
    assert len(bf)==bf.sample_id.nunique()==len(current_rmsf)==32
    assert bf.sample_id.tolist()==current_rmsf.sample_id.tolist()
    assert np.array_equal(bf['count'],current_rmsf['count']) and int(bf['count'].sum())==18748
    copies=[]
    for protocol,rows in old_bf.groupby('protocol'):
        rows=rows.sort_values('sample_id').reset_index(drop=True)
        assert rows.sample_id.tolist()==bf.sample_id.tolist()
        assert np.array_equal(rows[['count','mae','rmse','bias','spearman','pearson']].to_numpy(),bf[['count','mae','rmse','bias','spearman','pearson']].to_numpy())
        copies.append(dict(protocol=protocol,records=32,identical_fixed_official_prediction=True))
    old_records={row['sample_id']:row for row in [json.loads(line) for line in (OLD/'external_thermopcd/ready_test.jsonl').read_text().splitlines()]}
    new_inputs={row['sample_id']:row for row in new_audit['inputs']}
    assert set(old_records)==set(new_inputs)==set(bf.sample_id)
    residue_table=csv(OLD/'backflip2_corrected_scoring/external_per_residue.csv.gz')
    assert len(residue_table)==18748 and residue_table.sample_id.nunique()==32
    with np.load(OLD/'backflip2_corrected_scoring/external_aligned.npz',allow_pickle=False) as old_cache:
        old_prediction=old_cache['prediction'].copy();old_target=old_cache['target'].copy();old_valid=old_cache['valid'].copy()
    assert np.array_equal(old_target,residue_table.target_rmsf.to_numpy(np.float32))
    assert np.array_equal(old_valid,residue_table.valid_rmsf.to_numpy(bool))
    assert np.array_equal(old_prediction,residue_table.backflip2_official.to_numpy(float))
    manifest=csv(OLD/'backflip2_corrected_inputs/input_manifest.csv')
    manifest=manifest[manifest.role=='external']
    mapping=csv(OLD/'backflip2_corrected_inputs/residue_mapping.csv')
    mapping=mapping[mapping.role=='external']
    official_audit=read(OLD/'backflip2_corrected_scoring/scoring_audit.json')
    output_hashes={entry['input_id']:entry['sha256'] for entry in official_audit['output_hashes']}
    official_prediction={}
    for row in manifest.itertuples():
        input_path=OLD/'backflip2_corrected_inputs'/row.pdb_file
        prediction_path=OLD/'backflip2_corrected_completed/predictions'/(row.input_id+'.npz')
        assert sha(input_path)==row.input_sha256
        assert sha(prediction_path)==output_hashes[row.input_id]
        with np.load(prediction_path,allow_pickle=False) as saved:
            prediction=saved['rmsf_angstrom'].copy();cov=saved['covariance_nm2']
            assert prediction.shape==(row.n_residues,)
            assert np.allclose(prediction,10.0*np.sqrt(np.maximum(np.trace(cov,axis1=1,axis2=2),0.0)),atol=2e-6,rtol=2e-6)
        official_prediction[row.input_id]=prediction
    record_checks=[]
    for sample in bf.sample_id:
        old_record=old_records[sample];item=new_inputs[sample]
        # Manifest strings contain doubled Windows separators; Path accepts
        # those existing paths. Raw label hashes prove unchanged source labels.
        assert sha(old_record['label_path'])==old_record['metadata']['label_sha256']==item['label_sha256']
        assert sha(old_record['structure_path'])==old_record['metadata']['input_sha256']==item['structure_sha256']
        old_neural_path=OLD/'neural_exports/strict_external_evaluation/strict_full/fold_0'/(sample+'.npz')
        with np.load(old_neural_path,allow_pickle=False) as old_neural, np.load(NEW/'full'/(sample+'.npz'),allow_pickle=False) as new_neural:
            for key in ['amino_acid','chain','author_residue_id']:assert np.array_equal(old_neural[key],new_neural[key])
            assert np.array_equal(old_neural['target'][:,:2],new_neural['target'],equal_nan=True)
            assert np.array_equal(old_neural['valid'][:,:2],new_neural['valid'])
            rows=residue_table[residue_table.sample_id==sample].sort_values('residue_index')
            assert np.array_equal(rows.residue_index,np.arange(item['residues']))
            assert rows.chain.tolist()==new_neural['chain'].astype(str).tolist()
            assert rows.author_residue_id.tolist()==new_neural['author_residue_id'].astype(str).tolist()
            assert rows.amino_acid.tolist()==new_neural['amino_acid'].astype(str).tolist()
            assert np.array_equal(rows.target_rmsf.to_numpy(np.float32),new_neural['target'][:,0])
            assert np.array_equal(rows.valid_rmsf.to_numpy(bool),new_neural['valid'][:,0])
        predictions=np.full(len(rows),np.nan);visits=np.zeros(len(rows),int)
        for input_id,chain_rows in mapping[mapping.sample_id==sample].groupby('input_id'):
            ids=chain_rows.record_residue_index.to_numpy(int);indices=chain_rows.input_index.to_numpy(int)
            predictions[ids]=official_prediction[input_id][indices];visits[ids]+=1
        assert (visits==1).all() and np.array_equal(predictions,rows.backflip2_official.to_numpy(float))
        valid=rows.valid_rmsf.to_numpy(bool);target=rows.target_rmsf.to_numpy(np.float32)[valid];prediction=predictions[valid]
        mae=float(np.mean(abs(target-prediction)))
        rho=float(pd.Series(target).rank().corr(pd.Series(prediction).rank()))
        old_metric=bf[bf.sample_id==sample].iloc[0]
        assert abs(mae-old_metric.mae)<1e-12 and abs(rho-old_metric.spearman)<1e-12
        record_checks.append(dict(sample_id=sample,rmsf_count=int(valid.sum()),mae=mae,spearman=rho,
            raw_label_sha256=item['label_sha256'],raw_structure_sha256=item['structure_sha256'],
            labels_masks_residue_identity_old_new_exact=True,official_chain_predictions_reassembled_exact=True))
    assert abs(bf.mae.mean()-.7485120335786064)<1e-12 and abs(bf.spearman.mean()-.5930374471464412)<1e-12
    # Execute only the actual fig4 AST with recording plot operations. File
    # writers/exporters are intercepted, so no script, figure or source CSV is
    # changed and no Matplotlib rendering or inference is requested.
    script_path=FIG/'redraw_external_2430.py';script=script_path.read_text(encoding='utf8')
    node=next(node for node in ast.parse(script).body if isinstance(node,ast.FunctionDef) and node.name=='fig4')
    axes=[CapturedAxes(),CapturedAxes()];captured={}
    def capture_write(frame,name):captured['exported_points']=frame.copy();captured['points_filename']=name
    def capture_export(fig,stem,caption,panels,stats,sources):captured.update(stem=stem,caption=caption,panels=panels,stats=pd.DataFrame(stats),sources=sources)
    metrics=['rmsf','chi_entropy'];methods=['full','intrachain','static_sequence_complex','aa_only']
    environment=dict(pd=pd,np=np,M=methods,ML=dict(zip(methods,['Full','Intra','Static','AA'])),MC=dict(zip(methods,['b','o','g','grey'])),
        TASKS=metrics,UN={'rmsf':'Å','chi_entropy':'nats'},TN={'rmsf':'RMSF','chi_entropy':'Hχ'},BLACK='black',
        write=capture_write,grid=lambda *a,**k:(CapturedFigure(),axes),numeric_ticks=lambda *a,**k:None,
        label=lambda *a,**k:None,export=capture_export,EXT=NEW)
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(script_path),'exec'),environment)
    source_points=pd.read_csv(NEW/'per_record_metrics.csv')
    source_points['protocol']='strict_single_seed7';source_points['family_components']=1
    environment['fig4'](source_points)
    plotted=csv(FIG/'fig4_all_32_strict_single_model_metrics.csv')
    keys=['sample_id','variant','metric']
    source_sorted=new.sort_values(keys).reset_index(drop=True);plotted_sorted=plotted.sort_values(keys).reset_index(drop=True)
    assert len(plotted)==256 and not plotted.duplicated(keys).any()
    assert source_sorted[keys].equals(plotted_sorted[keys])
    assert np.array_equal(source_sorted['count'],plotted_sorted['count'])
    assert np.allclose(source_sorted[['mae','rmse','bias','spearman','pearson']],plotted_sorted[['mae','rmse','bias','spearman','pearson']],atol=1e-14,rtol=0)
    stats=csv(FIG/'fig4_strict_external_paired_statistics.csv')
    macro=csv(NEW/'record_macro_metrics.csv')
    panel_checks=[]
    for index,metric in enumerate(metrics):
        ax=axes[index];wide=source_points[source_points.metric==metric].pivot(index='sample_id',columns='variant',values='mae').reindex(columns=methods).sort_index()
        assert len(ax.lines)==32 and len(ax.scatters)==8
        for (_,values),(_,expected_row) in zip(ax.lines,wide.iterrows()):assert np.array_equal(values,expected_row.to_numpy())
        assert ax.limits[0]<wide.to_numpy().min() and ax.limits[1]>wide.to_numpy().max()
        means=[]
        for position,method in enumerate(methods):
            x,y,settings=ax.scatters[2*position];mx,my,mean_settings=ax.scatters[2*position+1]
            assert np.array_equal(x,np.full(32,position)) and np.array_equal(y,wide[method].to_numpy())
            assert mean_settings['marker']=='D' and int(mx)==position
            raw=new[(new.metric==metric)&(new.variant==method)];mean=float(raw.mae.mean())
            assert abs(float(my)-mean)<1e-14
            stored=stats[(stats.metric==metric)&(stats.variant==method)].iloc[0]
            published=macro[(macro.metric==metric)&(macro.variant==method)].iloc[0]
            assert abs(stored.mean_mae-mean)<1e-14 and abs(published.mae-mean)<1e-14
            assert stored.n_records==32 and stored.valid_residues==int(raw['count'].sum())
            assert abs(stored.median_mae-raw.mae.median())<1e-14 and abs(stored.min_mae-raw.mae.min())<1e-14 and abs(stored.max_mae-raw.mae.max())<1e-14
            assert ax.tick_labels[position].endswith(format(mean,'.3f'))
            means.append(dict(variant=method,records=32,valid_residues=int(raw['count'].sum()),mean_mae=mean,
                all_scatter_points_and_diamond_exact_to_source=True))
        panel_checks.append(dict(metric=metric,paired_lines=32,individual_points=128,mean_diamonds=4,
            yscale=ax.scale,limits_include_all_points=True,method_macros=means))
    manifest_path=FIG/'external_2430_figure_manifest.json';manifest=read(manifest_path)
    assert manifest['training_records']==2430 and manifest['excluded_records']==633 and manifest['new_training']==manifest['new_inference']==0
    figure_manifest=manifest['figures']['fig4_strict_external_paired']
    assert figure_manifest['input_sha256'][str(NEW/'per_record_metrics.csv')]==sha(NEW/'per_record_metrics.csv')
    assert manifest['audit_sha256']==sha(NEW/'evaluation_audit.json')
    for extension,item in figure_manifest['outputs'].items():assert sha(item['path'])==item['sha256']
    paths=[OLD_FIG/'all_32_per_record_metrics.csv',OLD/'backflip2_corrected_scoring/external_aligned.npz',
        OLD/'backflip2_corrected_scoring/external_per_residue.csv.gz',OLD/'backflip2_corrected_scoring/scoring_audit.json',
        NEW/'per_record_metrics.csv',NEW/'record_macro_metrics.csv',NEW/'evaluation_audit.json',script_path,
        FIG/'fig4_all_32_strict_single_model_metrics.csv',FIG/'fig4_strict_external_paired_statistics.csv',manifest_path]
    report=dict(status='passed_external_figure_and_backflip_official_reuse',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        evidence_sha256={str(path):sha(path) for path in paths},backflip_reuse=dict(allowed=True,records=32,rmsf_valid=18748,
            all_record_ids_counts_identical=True,all_raw_input_label_hashes_identical=True,all_rmsf_masks_targets_and_residue_identity_exact=True,
            also_two_target_old_new_masks_labels_exact=True,official_external_chain_outputs=95,all95_input_output_hashes_verified=True,
            all95_covariance_to_angstrom_conversion_verified=True,all_official_predictions_reassembled_exact=True,
            record_macro_mae=float(bf.mae.mean()),record_macro_spearman=float(bf.spearman.mean()),
            same_official_row_repeated_across_old_comparison_protocols=copies,record_checks=record_checks,
            released_checkpoint_sha256=official_audit['official_inference_metadata']['checkpoint_sha256'],
            official_source_commit=official_audit['official_inference_metadata']['source_commit'],
            source_note='Fixed released official independent-chain weights, not retrained on2430. Train purge changes FlexPPI/Static models, not this fixed comparator; source weights and external labels/predictions are unchanged.',
            limitation='Released comparator pretraining overlap and MD protocol differences remain disclosed; do not describe BF2 as trained with2430 records.'),
        external_figure=dict(training_records=2430,purged_records=633,all32_records_all4_methods_all2_metrics=True,
            plotted_point_rows=256,actual_fig4_function_executed_with_capture_only=True,no_figure_or_CSV_rewritten=True,
            no_sampling_no_outlier_removal_no_clipping=True,panels=panel_checks,record_macro_means_match_table_and_stats=True,
            figure_manifest_input_and_image_hashes_verified=True,visual_qa_scope='Root has visually inspected exported image; this audit independently checks data operations and existing artifact hashes.',
            descriptive_one_external_family=True,no_ci_or_pvalues=True),new_training=False,new_inference=False,completion_certificates_not_modified=True)
    output=ROOT/'audit/external_figure_and_backflip_reuse.json';output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':report['status'],'path':str(output),'BF2_MAE':report['backflip_reuse']['record_macro_mae'],
        'BF2_Spearman':report['backflip_reuse']['record_macro_spearman'],'external_records':32,'RMSF_positions':18748,'figure_panels':panel_checks},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
