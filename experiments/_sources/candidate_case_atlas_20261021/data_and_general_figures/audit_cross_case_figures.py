"""Independent numeric/coverage audit of figures produced by different authors."""
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parent
ATLAS=ROOT.parent;EXP=ATLAS.parent
sys.path.insert(0,str(EXP/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
META=pd.read_csv(EXP/'itrace_case_search_only_20261017/all244_ITrace_metadata_and_sequences.csv').set_index('pdb')
results=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cache(p):return np.load(EXP/'tcr_case_v40_20261017/family_predictions/full'/f'itrace_{p}.npz',allow_pickle=True)
def finiteassert(x,y,where,tol=3e-6):
    x=np.asarray(x,float);y=np.asarray(y,float)
    assert np.allclose(x,y,atol=tol,rtol=0,equal_nan=True),(where,np.nanmax(abs(x-y)))
    return float(np.nanmax(abs(x-y))) if np.isfinite(x-y).any() else 0.

# The transferred MEL5 figure retains canonical parent numbering and all six CDRs.
out=ATLAS/'figures/C01'
if out.exists():
    f=pd.read_csv(out/'source_plotted_CDR_RMSF.csv')
    assert len(f)==45 and len(f.region.unique())==6
    na=cache('3hg1');nb=cache('4jff');checks=[]
    for row in f.itertuples():
        role='tcr_alpha' if 'α' in row.region else 'tcr_beta'
        ca=str(META.loc['3hg1',role+'_chain']);cb=str(META.loc['4jff',role+'_chain'])
        ia=np.where(na['chain']==ca)[0][int(row.common_TCR_position)-1]
        ib=np.where(nb['chain']==cb)[0][int(row.derivative_input_position)-1]
        checks.append(finiteassert([row.MEL5_model_RMSF_A,row.alpha24beta17_model_RMSF_A],[na['mean'][ia,0],nb['mean'][ib,0]],'C01 CDR',1e-13))
    h=pd.read_csv(out/'source_plotted_interface_Hchi.csv',dtype={'position_a':str,'chain_a':str})
    assert len(h)==23
    for row in h.itertuples():
        role_col='mhc_chain' if row.role=='mhc' else 'peptide_chain'
        ca=str(META.loc['3hg1',role_col]);cb=str(META.loc['4jff',role_col])
        ia=np.where((na['chain']==ca)&(na['author_residue_id'].astype(str)==row.position_a))[0]
        ib=np.where((nb['chain']==cb)&(nb['author_residue_id'].astype(str)==row.position_a))[0]
        assert len(ia)==len(ib)==1
        checks.append(finiteassert([row.hchi_full_a,row.hchi_full_b],[na['mean'][ia[0],1],nb['mean'][ib[0],1]],'C01 Hchi',1e-13))
    a=json.loads((out/'data_audit.json').read_text(encoding='utf8'))
    for key,value in a['source_files'].items():assert sha(value['path'])==value['SHA256'],('C01 source hash',key)
    results.append(dict(case_id='C01',status='passed',RMSF_sites=45,Hchi_sites=23,raw_model_values_checked=136,max_difference=max(checks),PNG=str(out/'fig4_mel5_joint.png'),PNG_SHA256=sha(out/'fig4_mel5_joint.png')))

three=ATLAS/'figure_design/three_case_figures_manifest.json'
if three.exists():
    for entry in json.loads(three.read_text(encoding='utf8')):
        cid=entry['case_id'];out=Path(entry['directory'])
        stem={'C02':('2ak4','4jry'),'C03':('7ow5','7ow6'),'C04':('1oga','5hhm')}[cid]
        f=pd.read_csv(out/'source_complete_peptide.csv',float_precision='round_trip')
        pp=pd.read_csv(out/'source_absolute_RMSF.csv',float_precision='round_trip')
        hh=pd.read_csv(out/'source_absolute_Hchi.csv',float_precision='round_trip')
        assert len(f)==len(pp)=={'C02':13,'C03':10,'C04':9}[cid]
        assert len(hh)=={'C02':11,'C03':5,'C04':7}[cid]
        errors=[];n=0
        for suffix,p in zip(['a','b'],stem):
            z=cache(p);chain=str(META.loc[p,'peptide_chain']);idx=np.where(z['chain']==chain)[0]
            assert len(idx)==len(f)
            assert np.array_equal(z['amino_acid'][idx],f['aa_'+suffix].to_numpy(str))
            for metric,j in [('rmsf',0),('hchi',1)]:
                col=metric+('_prediction_' if cid=='C04' else '_full_')+suffix
                actual=f[col].to_numpy(float);valid=np.isfinite(actual)
                errors.append(finiteassert(actual[valid],z['mean'][idx,j][valid],(cid,p,metric),1e-13));n+=valid.sum()
            plotted=pp['first' if suffix=='a' else 'second']
            errors.append(finiteassert(plotted,z['mean'][idx,0],(cid,'RMSF artist source'),1e-13))
        assert entry['PNG_SHA256']==sha(entry['PNG'])
        delivery=json.loads((out/'delivery_manifest.json').read_text(encoding='utf8'))
        for name,h in delivery['files'].items():assert sha(out/name)==h,(cid,'file hash',name)
        results.append(dict(case_id=cid,status='passed',RMSF_sites=len(pp),Hchi_sites=len(hh),raw_model_values_checked=int(n),max_difference=max(errors),PNG=entry['PNG'],PNG_SHA256=sha(entry['PNG'])))

# Independent coverage arithmetic for the generalized discovery-stage figures.
general=json.loads((ROOT/'general_figure_manifest.json').read_text(encoding='utf8'))
catalog=json.loads((ATLAS/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
bycase={r['case_id']:r for r in catalog['cases']}
for entry in general['figures']:
    cid=entry['case_id'];case=bycase[cid];a,b=case['contrasts'][0]
    f=pd.read_csv(ROOT/cid/'data'/f'{a}_{b}_all_chain_alignment.csv')
    p=pd.read_csv(ROOT/cid/'figure/panel_b_complete_aligned_peptide.csv')
    h=pd.read_csv(ROOT/cid/'figure/panel_c_complete_comparable_pMHC_interface.csv')
    assert len(p)==f.role.eq('peptide').sum()
    expected=f.role.isin(['mhc','peptide'])&f.interface_union&f.paired&f.hchi_valid&f.same_aa
    assert len(h)==expected.sum()
    for metric in ['rmsf','hchi']:
        g=pd.read_csv(ROOT/cid/'figure'/f'panel_d_{metric}_all_contrasts.csv')
        expected_n=0
        for aa,bb in case['contrasts']:
            whole=pd.read_csv(ROOT/cid/'data'/f'{aa}_{bb}_all_chain_alignment.csv')
            mask=(whole.role.eq('peptide')&whole.paired&whole.rmsf_valid) if metric=='rmsf' else (whole.role.isin(['mhc','peptide'])&whole.interface_union&whole.paired&whole.hchi_valid&whole.same_aa)
            expected_n+=int(mask.sum())
        assert len(g)==expected_n,(cid,metric,len(g),expected_n)
    results.append(dict(case_id=cid,status='passed',RMSF_sites=len(p),Hchi_sites=len(h),all_defined_contrasts_covered_in_panel_d=True,PNG=entry['png'],PNG_SHA256=sha(entry['png'])))

result=dict(status='passed_for_completed_figures',audited_case_count=len(results),cases=results,
            remaining_case_ids=[f'C{i:02}' for i in range(5,18)],
            scope='Numerical and complete-region coverage audit; visual QA is separately recorded. No manuscript modification.')
(ROOT/'cross_case_figure_numeric_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(dict(status=result['status'],audited_case_count=len(results),remaining_case_ids=result['remaining_case_ids']),ensure_ascii=False,indent=2))
