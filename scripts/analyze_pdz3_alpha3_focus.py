"""Residue-level alpha3 analysis: no fitted response or fabricated dynamics."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from Bio.PDB import PDBParser

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/pdz3_alpha3_focus_20260914'
OUT=ROOT/'output/flexppi_pdz3_alpha3_focus_20260914'
METRICS=['rmsf','chi_entropy']
REGIONS={'carboxylate_loop':(318,323),'beta2_beta3_loop':(330,336),'middle_region':(341,356),'core_without_cut':(306,390)}

def dump(path,obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    specs=json.loads((BASE/'states.json').read_text())
    sm={s['id']:s for s in specs}
    d=pd.read_csv(BASE/'residue_deltas.csv')
    a=pd.read_csv(BASE/'absolute_predictions.csv')
    inter=pd.read_csv(BASE/'mutation_alpha3_interactions.csv')
    summary=d.groupby(['background','context','kind','mutant','author_number','aa'],sort=False).agg(
        rmsf_mean=('delta_rmsf','mean'),rmsf_sd=('delta_rmsf','std'),
        hchi_mean=('delta_chi_entropy','mean'),hchi_sd=('delta_chi_entropy','std'),
        rmsf_n_pos=('delta_rmsf',lambda x:(x>1e-5).sum()),hchi_n_pos=('delta_chi_entropy',lambda x:(x>1e-5).sum()),
        rmsf_n_neg=('delta_rmsf',lambda x:(x< -1e-5).sum()),hchi_n_neg=('delta_chi_entropy',lambda x:(x< -1e-5).sum()),
        n_models=('fold','nunique')).reset_index()
    summary['hchi_valid']=~summary.aa.isin(['A','G','P'])
    summary['boundary_near']=summary.author_number>=391
    summary.to_csv(BASE/'residue_summary.csv',index=False)
    parser=PDBParser(QUIET=True)
    coords={};geometry=[];contacts=[]
    for bg in ['apo1bfe','strip5heb','strip1be9']:
        st=parser.get_structure(bg,sm[f'{bg}_full_WT']['path'])[0]['A']
        rs=list(st.get_residues());coords[bg]={r.id[1]:{at.name:at.coord for at in r if at.element not in ('H','D')} for r in rs}
        amap=dict(zip(a[a.state==f'{bg}_full_WT'].author_number,a[a.state==f'{bg}_full_WT'].aa))
        alpha=[(p,k,v) for p,ats in coords[bg].items() if 396<=p<=402 for k,v in ats.items()]
        for p,ats in coords[bg].items():
            if p>=396:continue
            keys=list(ats);mat=cdist(np.stack(list(ats.values())),np.stack([x[2] for x in alpha]));ii,jj=np.unravel_index(mat.argmin(),mat.shape)
            geometry.append(dict(background=bg,author_number=p,aa=amap[p],alpha3_distance=float(mat[ii,jj]),
                atom=keys[ii],alpha3_residue=alpha[jj][0],alpha3_atom=alpha[jj][1],
                sidechain_complete=not(bg=='apo1bfe' and p==328)))
        for p,q in [(337,400),(337,397),(328,397),(328,400),(334,399),(355,401),(358,396),(337,334),(337,328),(337,327)]:
            kp,kq=list(coords[bg][p]),list(coords[bg][q]);mat=cdist(np.stack(list(coords[bg][p].values())),np.stack(list(coords[bg][q].values())))
            i,j=np.unravel_index(mat.argmin(),mat.shape)
            contacts.append(dict(background=bg,residue1=p,residue2=q,atom1=kp[i],atom2=kq[j],distance=float(mat[i,j]),
                sidechain_complete=not(bg=='apo1bfe' and 328 in (p,q))))
    geometry=pd.DataFrame(geometry);geometry.to_csv(BASE/'alpha3_residue_geometry.csv',index=False)
    pd.DataFrame(contacts).to_csv(BASE/'selected_static_contacts.csv',index=False)
    deletion=summary[(summary.kind=='deletion')&(summary.mutant=='WT')&(summary.context.isin(['d7','d8']))].merge(geometry,on=['background','author_number','aa'])
    deletion.to_csv(BASE/'deletion_residue_evidence.csv',index=False)
    # All residue-level region values stay in d; these summaries are only secondary.
    rrows=[]
    for name,(lo,hi) in REGIONS.items():
        sub=d[d.author_number.between(lo,hi)]
        for metric in METRICS:
            v=sub[sub.hchi_valid] if metric=='chi_entropy' else sub
            for keys,g in v.groupby(['background','context','kind','mutant','fold'],sort=False):
                rrows.append(dict(zip(['background','context','kind','mutant','fold'],keys),region=name,metric=metric,n=len(g),mean=g['delta_'+metric].mean(),mean_abs=g['delta_'+metric].abs().mean()))
    regional=pd.DataFrame(rrows);regional.to_csv(BASE/'regional_secondary_summary.csv',index=False)
    # Predeclared F337 saturation: affected OTHER core sites only, equal target mask.
    sat=d[(d.kind=='mutation')&(d.mutant.str.match(r'^F337[A-Z]$'))&(d.context.isin(['full','d8']))&d.author_number.between(312,390)]
    srows=[]
    for keys,g in sat.groupby(['background','context','mutant','fold']):
        for metric in METRICS:
            vals=g[g.hchi_valid]['delta_'+metric] if metric=='chi_entropy' else g['delta_'+metric]
            srows.append(dict(zip(['background','context','mutant','fold'],keys),metric=metric,n_targets=len(vals),mean_abs_response=vals.abs().mean(),response_norm=float(np.sqrt(np.mean(vals*vals)))))
    pd.DataFrame(srows).to_csv(BASE/'f337_saturation_summary.csv',index=False)
    # 2026 original experiment has its own energy metric; do NOT rescale it to Hchi.
    exp=pd.read_csv(BASE/'experiment_2026_data2.csv')
    exp3=pd.read_csv(BASE/'experiment_2026_data3.csv')
    exp3=exp3[exp3.Pos.notna()].copy()
    expsites=exp[(exp.assay=='Binding')&(exp.energy_term=='dddg_delC')].groupby('Pos').agg(
        exp_abs_energy=('total_mean_kcal_mol',lambda v:v.abs().median()),exp_signed_energy=('total_mean_kcal_mol','median'),exp_n=('id','size')).reset_index().rename(columns={'Pos':'author_number'})
    es=summary[(summary.kind=='deletion')&(summary.context=='d8')&(summary.mutant=='WT')&summary.hchi_valid&summary.author_number.between(312,377)].merge(expsites,on='author_number',validate='many_to_one')
    es.to_csv(BASE/'experimental_site_comparison.csv',index=False)
    stats=[]
    for bg,g in es.groupby('background'):
        stats.append(dict(background=bg,comparison='Absolute deletion Hchi vs median absolute experimental alpha3 binding-energy coupling, all eligible measured sites',n=len(g),rho=float(spearmanr(g.hchi_mean.abs(),g.exp_abs_energy).statistic)))
    for (bg,ctx),g in deletion[deletion.hchi_valid & ~deletion.boundary_near].groupby(['background','context']):
        stats.append(dict(background=bg,context=ctx,comparison='alpha3 distance vs absolute Hchi response; core<=390',n=len(g),rho=float(spearmanr(g.alpha3_distance,g.hchi_mean.abs()).statistic)))
    dump(BASE/'descriptive_statistics.json',stats)
    # Direct numerical audits independent from figure generation.
    checks={}
    checks['states']=len(specs);checks['all_single_chain_no_ligand']=all(s['chains']==['A'] for s in specs)
    checks['absolute_unique']=not a.duplicated(['state','fold','author_number']).any()
    checks['models_per_residue']=sorted(summary.n_models.unique().tolist())
    checks['delta_arithmetic_max_error']=max(float(np.max(np.abs(d['delta_'+m]-(d['after_'+m]-d['before_'+m])))) for m in METRICS)
    checks['mutated_site_rows']=int(sum(d[(d.kind=='mutation')&(d.mutant==s['mutant'])&(d.background==s['background'])&(d.context==s['context'])].author_number.isin([int(m[1:-1]) for m in s['mutations']]).sum() for s in specs if s['mutations']))
    maxdiff=0.;bad_backbones=[]
    for s in specs:
        st=parser.get_structure(s['id'],s['path'])[0]['A'];ref=coords[s['background']]
        for r in st:
            for at in r:
                if at.name in ('N','CA','C','O'):
                    maxdiff=max(maxdiff,float(np.max(np.abs(at.coord-ref[r.id[1]][at.name]))))
    checks['max_fixed_backbone_deviation_A']=maxdiff
    checks['verified_cache_count']=0
    cp={p['fold']:p for p in json.loads((BASE/'checkpoint_provenance.json').read_text())}
    for s in specs:
        assert hashlib.sha256(Path(s['path']).read_bytes()).hexdigest()==s['input_sha256']
        for fold in range(5):
            p=BASE/'prediction_cache'/f"{s['id']}_fold{fold}.npz"
            with np.load(p,allow_pickle=False) as z:
                prov=json.loads(str(z['provenance']))
                assert prov['input']==s['input_sha256']
                assert prov['checkpoint']==cp[fold]['sha256']
                assert prov['implementation']==cp[fold]['implementation_sha256']
                assert prov['protocol']=='pdz3_alpha3_no_ligand_v1'
                assert z['prediction'].shape==(s['n_residues'],2)
                checks['verified_cache_count']+=1
    expwide=exp.pivot(index=['assay','id'],columns='energy_term',values='total_mean_kcal_mol')
    checks['experiment_delC_additive_identity_max_error_kcal']=float((expwide.ddg_Npdz-expwide.ddg_NpdzC-expwide.dddg_delC).abs().max())
    checks['experiment_readme_threshold_discrepancy']='README says FDR<0.5; article Methods says0.05. Original supplied significance retained, no newly inferred experimental stars.'
    assert checks['all_single_chain_no_ligand'] and checks['absolute_unique']
    assert checks['models_per_residue']==[5] and maxdiff==0 and checks['mutated_site_rows']==0
    assert checks['delta_arithmetic_max_error']<1e-6
    dump(BASE/'numerical_audit.json',checks)
    print(json.dumps(checks,indent=2))
    print('Main mutant responses:')
    print(regional[(regional.region=='beta2_beta3_loop')&(regional.metric=='rmsf')&(regional.context=='full')&regional.mutant.isin(['V328A','V328I','F337A','Y397A','F400A','R399Q','E401Q'])].groupby(['background','mutant'])['mean'].agg(['mean','std']).round(4).to_string())

if __name__=='__main__':main()
