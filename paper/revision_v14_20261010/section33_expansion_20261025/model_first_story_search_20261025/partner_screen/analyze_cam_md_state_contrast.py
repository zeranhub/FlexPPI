"""CaM NHE1 / IQCG structural-state contrast, all common residues.

The PDB states differ in calcium occupancy as well as partner; never a causal
partner-only intervention. Fixed mature numbering and lobe ranges precede
regional summaries. Model fold0 is in-sample for both states.
"""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import joblib
from threadpoolctl import threadpool_limits
from screen_exact_sequence_contexts import score
OUT=Path(__file__).resolve().parent;R=OUT.parent.parent
CORE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
REF='ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK'
SIDS=['A004C','A007U']

def main():
    contrast=pd.read_csv(OUT/'fixed_fold0_candidate_all_residue_deltas.csv.gz')
    contrast=contrast[contrast.pair_id.eq('A004C__A007U')].copy()
    inv=[]
    for sid in SIDS:
        z=np.load(CORE/'neural_exports/evaluation/full'/(sid+'.npz'))
        target='B' if sid=='A004C' else 'A';idx=np.flatnonzero(z['chain']==target)
        assert ''.join(z['amino_acid'][idx[:145]])==REF[3:]
        for j,i in enumerate(idx[:145]):
            inv.append(dict(sample_id=sid,chain=target,residue_index=int(i),local_author=z['author_residue_id'][i],mature_cam_number=j+4,aa=z['amino_acid'][i]))
    inv=pd.DataFrame(inv);inv.to_csv(OUT/'cam_md_mature_number_mapping.csv',index=False)
    lookup=inv.set_index(['sample_id','residue_index'])
    contrast['mature_cam_number']=[int(lookup.loc[(s,int(i)),'mature_cam_number']) for s,i in zip(contrast.sample_a,contrast.index_a)]
    assert all(contrast.amino_acid==[REF[x-1] for x in contrast.mature_cam_number])
    assert all(contrast.mature_cam_number==[int(lookup.loc[(s,int(i)),'mature_cam_number']) for s,i in zip(contrast.sample_b,contrast.index_b)])
    contrast['region']=np.select([contrast.mature_cam_number.le(75),contrast.mature_cam_number.le(81)],['N_lobe_4_75','linker_76_81'],default='C_lobe_82_148')
    # Apply the EXISTING fold0 static models to both states, no fitting.
    contract=json.loads((R/'feature_baseline/baseline_contract.json').read_text())
    feat=pd.read_parquet(R/'feature_baseline/features_all.parquet',filters=[('sample_id','in',SIDS)])
    base={};manifest=[]
    with threadpool_limits(limits=2):
        for model in ['aa_sasa','aa_sasa_contacts','static47_sasa']:
            cols=contract['feature_sets'][model]
            for metric in ['rmsf','hchi']:
                path=R/'feature_baseline/models'/f'{model}_fold0_{metric}.joblib'
                fitted=joblib.load(path);values=fitted.predict(feat[cols].to_numpy(dtype=np.float32))
                tab=pd.Series(values,index=pd.MultiIndex.from_frame(feat[['sample_id','residue_index']]))
                base[model,metric]=tab
                manifest.append(dict(model=model,metric=metric,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        for model in ['aa_sasa','aa_sasa_contacts','static47_sasa']:
            contrast[model+'0_delta']=[float(base[model,m].loc[(s2,int(i2))]-base[model,m].loc[(s1,int(i1))]) for m,s1,s2,i1,i2 in zip(contrast.metric,contrast.sample_a,contrast.sample_b,contrast.index_a,contrast.index_b)]
    # State values as well as differences, to make the effect concrete.
    for model in ['full','intrachain']:
        pred={sid:np.load(OUT/'training_member_predictions'/(model+'_'+sid+'.npz'))['mean'] for sid in SIDS}
        prefix='full0' if model=='full' else 'intra0'
        for state in ['a','b']:
            contrast[prefix+'_'+state]=[float(pred[sid][int(i),int(m=='hchi')]) for sid,i,m in zip(contrast['sample_'+state],contrast['index_'+state],contrast.metric)]
    contrast.to_csv(OUT/'cam_md_state_contrast_all_common_residues.csv',index=False)
    metrics=[]
    for m,g in contrast.groupby('metric'):
        for name,gg in [('all',g),('interface_union',g[g.interface]),('noninterface',g[~g.interface])]+list(g.groupby('region')):
            if len(gg)<3:continue
            for model in ['full0','intra0','aa_sasa0','aa_sasa_contacts0','static47_sasa0']:
                metrics.append(dict(metric=m,region=name,model=model,**score(gg.md_delta.to_numpy(),gg[model+'_delta'].to_numpy())))
    pd.DataFrame(metrics).to_csv(OUT/'cam_md_state_contrast_region_metrics.csv',index=False)
    regionmean=[]
    for (m,region),g in contrast.groupby(['metric','region']):
        row=dict(metric=m,region=region,n=len(g),md_a=g.md_a.mean(),md_b=g.md_b.mean(),md_delta=g.md_delta.mean())
        for model in ['full0','intra0']:
            row[model+'_a']=g[model+'_a'].mean();row[model+'_b']=g[model+'_b'].mean();row[model+'_delta']=g[model+'_delta'].mean()
        regionmean.append(row)
    pd.DataFrame(regionmean).to_csv(OUT/'cam_md_state_contrast_region_means.csv',index=False)
    # Independent trajectory recomputation, five chronological blocks per state.
    windows=[];blockrows=[]
    for metric,file,col in [('hchi','all_sidechains_per_window.csv','hchi'),('rmsf','all_residues_rmsf_per_window.csv','rmsf')]:
        tabs={sid:pd.read_csv(OUT/'verified_case_trajectories'/sid/file).set_index(['residue_index','block']) for sid in SIDS}
        g=contrast[contrast.metric.eq(metric)]
        for _,row in g.iterrows():
            for block in range(1,6):
                a=tabs[SIDS[0]].loc[int(row.index_a),block];b=tabs[SIDS[1]].loc[int(row.index_b),block]
                values=dict(metric=metric,mature_cam_number=int(row.mature_cam_number),aa=row.amino_acid,region=row.region,block=block,
                            md_a=float(a[col]),md_b=float(b[col]),md_delta=float(b[col]-a[col]),mean_md_delta=row.md_delta,full0_delta=row.full0_delta,intra0_delta=row.intra0_delta)
                if metric=='hchi':
                    for quantity in ['fine_entropy_mean_per_chi','width_mean_per_chi','shifted60_joint_entropy']:
                        values[quantity+'_delta']=float(b[quantity]-a[quantity])
                windows.append(values)
        w=pd.DataFrame(windows);w=w[w.metric.eq(metric)]
        for block,q in w.groupby('block'):
            for region,sub in [('all',q)]+list(q.groupby('region')):
                for model in ['full0','intra0']:
                    blockrows.append(dict(metric=metric,region=region,block=block,model=model,**score(sub.md_delta.to_numpy(),sub[model+'_delta'].to_numpy())))
    windows=pd.DataFrame(windows);windows.to_csv(OUT/'cam_md_state_contrast_all_residues_five_blocks.csv',index=False)
    pd.DataFrame(blockrows).to_csv(OUT/'cam_md_state_contrast_five_block_metrics.csv',index=False)
    robust=[]
    for (m,pos),g in windows.groupby(['metric','mature_cam_number']):
        row=dict(metric=m,mature_cam_number=pos,aa=g.aa.iloc[0],region=g.region.iloc[0],md_delta=g.mean_md_delta.iloc[0],
                 full0_delta=g.full0_delta.iloc[0],intra0_delta=g.intra0_delta.iloc[0],min_block_delta=g.md_delta.min(),max_block_delta=g.md_delta.max(),
                 blocks_same_sign_as_mean=int((np.sign(g.md_delta)==np.sign(g.mean_md_delta)).sum()))
        for q in ['fine_entropy_mean_per_chi','width_mean_per_chi','shifted60_joint_entropy']:
            if m=='hchi':row[q+'_delta']=g[q+'_delta'].mean()
        robust.append(row)
    robust=pd.DataFrame(robust);robust.to_csv(OUT/'cam_md_state_contrast_residue_stability.csv',index=False)
    # Complete chemically defined methionine panel, not selected successful sites.
    contrast[contrast.amino_acid.eq('M')].to_csv(OUT/'cam_md_all_nine_methionines.csv',index=False)
    audit=dict(status='complete',contrast='IQCG / apo-CaM 4LZX minus NHE1 / Ca-bound CaM 2YGG',
               caveat='Both partner and calcium occupancy differ; not partner-only causal inference. RMSF fitted per complete CaM chain, not independently per lobe.',
               common_mature_numbers='4 through 148; NHE1 CaM expression-tail Ala149 excluded',
               fixed_regions={'N_lobe':[4,75],'linker':[76,81],'C_lobe':[82,148]},
               baseline_checkpoints=manifest,model_scope='Frozen corrected fold0 epoch40, both states in fold0 training set; no heldout-generalization claim.')
    (OUT/'cam_md_state_contrast_analysis_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print('\nREGION MEANS\n',pd.DataFrame(regionmean).to_string(index=False))
    print('\nMETRICS\n',pd.DataFrame(metrics).query("region=='all'").to_string(index=False))
    print('\nMETHIONINES\n',contrast[contrast.amino_acid.eq('M')][['metric','mature_cam_number','md_a','md_b','md_delta','full0_a','full0_b','full0_delta','intra0_delta','contact_a','contact_b']].to_string(index=False))
    print('\nSTABILITY\n',robust.groupby('metric').apply(lambda g:pd.Series(dict(n=len(g),strong=int((g.md_delta.abs()>=.1).sum()),strong_4of5=int(((g.md_delta.abs()>=.1)&(g.blocks_same_sign_as_mean>=4)).sum()),strong_5of5=int(((g.md_delta.abs()>=.1)&(g.blocks_same_sign_as_mean==5)).sum()))),include_groups=False).to_string())
if __name__=='__main__':main()
