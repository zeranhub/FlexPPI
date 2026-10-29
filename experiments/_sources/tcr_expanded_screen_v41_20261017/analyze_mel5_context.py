from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012'/'plot_runtime'))
import numpy as np
import pandas as pd
read=lambda p:pd.read_csv(p,float_precision='round_trip')
source=read(ROOT/'mel5_validation'/'corrected_pair_all_aligned_positions.csv')
pred=read(ROOT/'mel5_context'/'mel5_context_predictions.csv')
out=ROOT/'mel5_context'
audit={'interface_definition':'C-alpha to C-alpha distance <=8 Å, union of either-state interface; same definition as evaluate_itrace_benchmark.py:interface_mask.','state_interface_flag_mismatches':{}}
def atom_ca(path):
    data=[];seen=set()
    for l in path.read_text().splitlines():
        if l.startswith('ATOM') and l[12:16].strip()=='CA':
            key=(l[21],l[22:27].strip())
            if key in seen:continue
            seen.add(key)
            data.append((key,np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])))
    return data
for state,pdb in [('a','3hg1'),('b','4jff')]:
    data=atom_ca(ROOT/'mel5_validation'/f'{pdb}_model_input.pdb')
    tcr=np.array([v for (ch,_),v in data if ch in 'DE'])
    pmhc=np.array([v for (ch,_),v in data if ch in 'ABC'])
    flags={}
    for (ch,rid),v in data:
        other=tcr if ch in 'ABC' else pmhc
        flags[(ch,rid)]=float(np.sqrt(((other-v)**2).sum(axis=1)).min())<=8
    mismatch=[]
    for _,r in source.iterrows():
        key=(str(r[f'chain_{state}']),str(r[f'id_{state}']))
        if bool(flags[key])!=bool(r[f'interface_{state}']):mismatch.append([*key,bool(flags[key]),bool(r[f'interface_{state}'])])
    audit['state_interface_flag_mismatches'][pdb]=mismatch
    assert not mismatch,(pdb,mismatch)

natural=source[source.role.isin(['mhc','peptide'])&source.same_aa].copy()
for kind in ['full','pMHC_only']:
    for state,pdb in [('a','3hg1'),('b','4jff')]:
        tab=pred[pred.pdb.eq(pdb)&pred.input_kind.eq(kind)].set_index(['chain','residue_id'])
        for metric,col in [('rmsf','rmsf_pred'),('hchi','Hchi_pred')]:
            natural[f'{metric}_{kind}_{state}']=[float(tab.loc[(str(ch),int(rid)),col]) for ch,rid in zip(natural[f'chain_{state}'],natural[f'id_{state}'])]
    for metric in ['rmsf','hchi']:
        natural[f'{metric}_{kind}_delta']=natural[f'{metric}_{kind}_b']-natural[f'{metric}_{kind}_a']
rows=[]
for scope,frame in [('natural_pMHC_all',natural),('natural_pMHC_interface',natural[natural.interface]),('MHC_interface',natural[natural.interface&natural.role.eq('mhc')]),('peptide_all',natural[natural.role.eq('peptide')])]:
    for metric in ['rmsf','hchi']:
        f=frame[frame[f'{metric}_valid']].copy()
        y=f[f'{metric}_md_delta'].to_numpy()
        for kind in ['full','pMHC_only']:
            p=f[f'{metric}_{kind}_delta'].to_numpy()
            rows.append(dict(scope=scope,metric=metric,input_kind=kind,n=len(f),rho=float(np.corrcoef(pd.Series(y).rank(),pd.Series(p).rank())[0,1]),MAE=float(np.mean(abs(y-p))),zero_change_MAE=float(np.mean(abs(y))),MAE_gain=float(1-np.mean(abs(y-p))/np.mean(abs(y))),mean_abs_model_delta=float(np.mean(abs(p))),mean_abs_MD_delta=float(np.mean(abs(y)))))
summary=pd.DataFrame(rows)
natural.to_csv(out/'all_natural_pMHC_context_comparison.csv',index=False)
summary.to_csv(out/'regional_context_comparison.csv',index=False)
q=natural[natural.role.eq('mhc')&natural.id_a.eq(68)]
assert len(q)==1
audit['K68']={c:float(q.iloc[0][c]) for c in ['hchi_md_a','hchi_md_b','hchi_md_delta','hchi_full_a','hchi_full_b','hchi_full_delta','hchi_pMHC_only_a','hchi_pMHC_only_b','hchi_pMHC_only_delta']}
audit['status']='passed'
audit['interpretation']='Static input removal tests input dependence of this joint model, not causal biological entropy or affinity change; pMHC coordinates differ between structures and are retained in both inputs.'
(out/'context_comparison_and_interface_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
print(summary[summary.scope.isin(['natural_pMHC_interface','MHC_interface'])].to_string(index=False))
print(json.dumps(audit,indent=2))
