"""Read existing joint-model predictions and previously recomputed labels only."""
from pathlib import Path
import csv, hashlib, json, math
import numpy as np

ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent
OLD=EXP/'tcr_expanded_screen_v41_20261017'
PRED=EXP/'tcr_case_v40_20261017/family_predictions/full'
OUT=ROOT/'data_audit';OUT.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def yes(x):return str(x).lower()=='true'
def num(x):return float(x) if x not in ('',None) else None
def rank(x):
    order=np.argsort(x,kind='stable');ans=np.empty(len(x),float);i=0
    while i<len(x):
        j=i+1
        while j<len(x) and x[order[j]]==x[order[i]]:j+=1
        ans[order[i:j]]=(i+j-1)/2;i=j
    return ans
def rho(a,b):
    if len(a)<4 or np.ptp(a)==0 or np.ptp(b)==0:return None
    return float(np.corrcoef(rank(np.array(a)),rank(np.array(b)))[0,1])
sources=[];allrows=[];checks=[];target_errors=[]
npzs={}
def getnpz(pdb):
    if pdb not in npzs:
        p=PRED/f'itrace_{pdb}.npz'
        npzs[pdb]=dict(np.load(p,allow_pickle=False))
        sources.append({'path':str(p),'sha256':sha(p)})
    return npzs[pdb]
def canonid(x):
    try:return str(int(x))
    except ValueError:return x
def match(pdb,chain,resid,aa,metric,value,target_value):
    z=getnpz(pdb)
    inds=[i for i,(ch,rs) in enumerate(zip(z['chain'],z['author_residue_id'])) if str(ch)==chain and canonid(str(rs))==canonid(resid)]
    assert len(inds)==1,(pdb,chain,resid,len(inds))
    i=inds[0];assert str(z['amino_acid'][i])==aa
    k=0 if metric=='rmsf' else 1
    err=abs(float(z['mean'][i,k])-value)
    assert err<3e-6,(pdb,chain,resid,metric,err)
    target_err=abs(float(z['target'][i,k])-target_value)
    assert target_err<3e-6,(pdb,chain,resid,metric,'target',target_err)
    target_errors.append(target_err)
    return err
datasets=[
 ('MEL5_parent_to_alpha24beta17','3hg1','4jff',OLD/'mel5_validation/corrected_pair_all_aligned_positions.csv'),
 ('PIK3CA_two_TCRs','7l1d','7rrg',OLD/'pik3ca_validation/recomputed_pair_all_aligned_positions.csv'),
 ('SB27_to_SB47','2ak4','4jry',EXP/'tcr_case_v40_20261017/figure6/all13_model_and_MD_changes.csv')]
for case,pa,pb,path in datasets:
    rr=read(path);sources.append({'path':str(path),'sha256':sha(path)})
    errors=[]
    for row in rr:
        for metric in ['rmsf','hchi']:
            if not yes(row[metric+'_valid']):continue
            for state,pdb in [('a',pa),('b',pb)]:
                errors.append(match(pdb,row['chain_'+state],row['id_'+state],row['aa_'+state],metric,float(row[f'{metric}_full_{state}']),float(row[f'{metric}_md_{state}'])))
        row['case']=case
        allrows.append(row)
    checks.append({'case':case,'aligned_rows':len(rr),'checked_prediction_values':len(errors),'max_absolute_NPZ_to_CSV_error':max(errors),'passed':True})

jpath=OLD/'jm22_f5l_validation/all9_model_MD_physical_changes.csv'
sources.append({'path':str(jpath),'sha256':sha(jpath)})
errors=[]
for r in read(jpath):
    row={'case':'JM22_F5L','a':'1oga','b':'5hhm','role':'peptide','chain_a':'C','chain_b':'C',
      'id_a':r['position'],'id_b':r['position'],'position_a':r['position'],'position_b':r['position'],
      'aa_a':r['aa_a'],'aa_b':r['aa_b'],'same_aa':r['same_aa'],'region':'peptide','region_a':'peptide','region_b':'peptide','interface':'True'}
    for metric in ['rmsf','hchi']:
        row[metric+'_valid']=str(yes(r[metric+'_valid_a']) and yes(r[metric+'_valid_b']))
        for state,pdb in [('a','1oga'),('b','5hhm')]:
            row[f'{metric}_full_{state}']=r[f'{metric}_prediction_{state}']
            row[f'{metric}_md_{state}']=r[f'{metric}_MD_{state}']
            if yes(row[metric+'_valid']):
                errors.append(match(pdb,'C',r['position'],r['aa_'+state],metric,float(r[f'{metric}_prediction_{state}']),float(r[f'{metric}_MD_{state}'])))
        row[f'{metric}_full_delta']=r[f'{metric}_prediction_delta']
        row[f'{metric}_md_delta']=r[f'{metric}_MD_delta']
    allrows.append(row)
checks.append({'case':'JM22_F5L','aligned_rows':9,'checked_prediction_values':len(errors),'max_absolute_NPZ_to_CSV_error':max(errors),'passed':True})

def region(row):return row.get('region',row.get('region_a',''))
def scopes(case):
    if case=='MEL5_parent_to_alpha24beta17':
        return [('all_six_CDRs',lambda r:'CDR' in region(r)),
          ('both_CDR3',lambda r:'CDR3' in region(r)),
          ('CDR2beta',lambda r:region(r)=='CDR2β'),
          ('CDR3beta',lambda r:region(r)=='CDR3β'),
          ('same_AA_both_CDR3',lambda r:'CDR3' in region(r) and yes(r['same_aa'])),
          ('same_AA_pMHC_interface',lambda r:r['role'] in ['mhc','peptide'] and yes(r['interface']) and yes(r['same_aa']))]
    if case=='PIK3CA_two_TCRs':
        return [('same_AA_pMHC_interface',lambda r:r['role'] in ['mhc','peptide'] and yes(r['interface']) and yes(r['same_aa'])),
          ('same_peptide_all',lambda r:r['role']=='peptide')]
    return [('all_peptide',lambda r:r['role']=='peptide'),('same_AA_peptide',lambda r:r['role']=='peptide' and yes(r['same_aa']))]
summary=[]
for case in [d[0] for d in datasets]+['JM22_F5L']:
    rows=[r for r in allrows if r['case']==case]
    for scope,select in scopes(case):
        subset=[r for r in rows if select(r)]
        for metric,unit in [('rmsf','Å'),('hchi','nats')]:
            vv=[r for r in subset if yes(r[metric+'_valid'])]
            pred_a=np.array([float(r[f'{metric}_full_a']) for r in vv]);pred_b=np.array([float(r[f'{metric}_full_b']) for r in vv])
            md_a=np.array([float(r[f'{metric}_md_a']) for r in vv]);md_b=np.array([float(r[f'{metric}_md_b']) for r in vv])
            dp=pred_b-pred_a;dm=md_b-md_a
            summary.append({'case':case,'scope':scope,'metric':metric,'unit':unit,'n':len(vv),
              'n_changed_AA':sum(not yes(r['same_aa']) for r in vv),
              'predicted_mean_before':float(pred_a.mean()),'predicted_mean_after':float(pred_b.mean()),
              'predicted_mean_change':float(dp.mean()),'n_predicted_decrease':int((dp<0).sum()),'n_predicted_increase':int((dp>0).sum()),
              'existing_MD_mean_before':float(md_a.mean()),'existing_MD_mean_after':float(md_b.mean()),
              'change_spearman_with_existing_MD':rho(dp,dm),'change_MAE_with_existing_MD':float(abs(dp-dm).mean()),
              'zero_change_MAE_with_existing_MD':float(abs(dm).mean()),
              'improvement_over_zero_change':float(1-abs(dp-dm).mean()/abs(dm).mean()) if abs(dm).mean()>0 else None,
              'interpretation':'Descriptive exploratory comparison of existing predictions and independent MD labels; not an affinity calculation or new generalization test.'})

normalized=[]
for r in allrows:
    for metric,unit in [('rmsf','Å'),('hchi','nats')]:
        if not yes(r[metric+'_valid']):continue
        normalized.append({'case':r['case'],'pdb_before':r['a'],'pdb_after':r['b'],'role':r['role'],
          'chain_before':r['chain_a'],'chain_after':r['chain_b'],'residue_before':r['id_a'],'residue_after':r['id_b'],
          'aa_before':r['aa_a'],'aa_after':r['aa_b'],'same_AA':yes(r['same_aa']),'region':region(r),
          'interface':yes(r['interface']),'metric':metric,'unit':unit,
          'predicted_before':float(r[f'{metric}_full_a']),'predicted_after':float(r[f'{metric}_full_b']),
          'existing_MD_before':float(r[f'{metric}_md_a']),'existing_MD_after':float(r[f'{metric}_md_b'])})

def writecsv(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
writecsv(OUT/'candidate_complete_region_summary.csv',summary)
writecsv(OUT/'all_aligned_absolute_values.csv',normalized)
context_path=OLD/'mel5_context/regional_context_comparison.csv'
sources.append({'path':str(context_path),'sha256':sha(context_path)})
context=read(context_path)
context_rows_path=OLD/'mel5_context/all_natural_pMHC_context_comparison.csv'
sources.append({'path':str(context_rows_path),'sha256':sha(context_rows_path)})
context_rows=read(context_rows_path)
context_arrays={}
for pdb in ['3hg1','4jff']:
    for kind in ['full','pMHC_only']:
        path=OLD/f'mel5_context/{pdb}_{kind}_context_prediction.npz'
        sources.append({'path':str(path),'sha256':sha(path)})
        zz=dict(np.load(path,allow_pickle=False))
        context_arrays[(pdb,kind)]={ (str(ch),canonid(str(rs))):(str(aa),mu) for ch,rs,aa,mu in zip(zz['chain'],zz['author_residue_id'],zz['amino_acid'],zz['mean'])}
context_errors=[]
for rr in context_rows:
    for side,pdb in [('a','3hg1'),('b','4jff')]:
        key=(rr['chain_'+side],canonid(rr['id_'+side]))
        for kind in ['full','pMHC_only']:
            aa,mu=context_arrays[(pdb,kind)][key]
            assert aa==rr['aa_'+side]
            for k,metric in enumerate(['rmsf','hchi']):
                context_errors.append(abs(float(mu[k])-float(rr[f'{metric}_{kind}_{side}'])))
assert max(context_errors)<1e-10
for metric in ['hchi','rmsf']:
    vv=[r for r in context_rows if yes(r['interface']) and yes(r[metric+'_valid']) and yes(r['same_aa'])]
    md=[float(r[f'{metric}_md_b'])-float(r[f'{metric}_md_a']) for r in vv]
    for kind in ['full','pMHC_only']:
        pred=[float(r[f'{metric}_{kind}_b'])-float(r[f'{metric}_{kind}_a']) for r in vv]
        original=next(r for r in context if r['scope']=='natural_pMHC_interface' and r['metric']==metric and r['input_kind']==kind)
        assert len(vv)==int(original['n'])
        assert abs(rho(pred,md)-float(original['rho']))<1e-12
        assert abs(float(np.mean(abs(np.array(pred)-md)))-float(original['MAE']))<1e-12
checks.append({'case':'MEL5_existing_TCR_input_removal_control','checked_raw_context_values':len(context_errors),'max_absolute_NPZ_to_CSV_error':max(context_errors),'both_readout_interface_summary_recomputed':True,'passed':True})
report={'status':'passed','date':'2026-10-19','model_checkpoint_sha256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',
 'training_ITrace_records':0,'training_records':2430,'new_inference':0,'new_training':0,'new_MD':0,
 'mapping':'MEL5 globally aligned corrected_ mapping only; PIK3CA same-sequence pMHC primary; SB27 and JM22 natural peptide positions.',
 'checks':checks,'sources':sources,'summary':summary,'max_existing_MD_mean_to_NPZ_target_error':max(target_errors),
 'MEL5_context_existing_control':[r for r in context if r['scope']=='natural_pMHC_interface' and r['metric']=='hchi'],
 'recommendation':'MEL5 parent vs alpha24beta17 is the primary complementary two-readout candidate: TCR loop RMSF plus same-sequence pMHC interface Hchi with an existing TCR-input-removal control. SB27/SB47 is the most coherent same-peptide two-readout alternative. JM22 is a clean single-point contrast with small model Hchi response. PIK3CA is a secondary broad-interface contrast.',
 'boundaries':['RMSF and Hchi need not decrease together or on identical residues.','Model input removal measures model context dependence, not biological unbinding or causal affinity.','Do not use changed amino-acid Hchi values as same-chemistry environment validation.','Do not turn bound-state side-chain entropy into binding entropy or affinity.','Keep complete predefined regions, including weaker positions; exploratory case selection is not independent model validation.']}
(OUT/'candidate_numerical_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':report['status'],'checks':checks,'summary':[r for r in summary if (r['scope'] in ['both_CDR3','CDR2beta','CDR3beta','same_AA_pMHC_interface','same_AA_peptide'])]},ensure_ascii=False,indent=2))
