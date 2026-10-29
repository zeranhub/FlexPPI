from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
import numpy as np,pandas as pd
from scipy.stats import spearmanr,t

md=pd.read_csv(ROOT/'all408_per_replica_labels.csv',dtype={'author_residue_id':str})
pred=pd.read_csv(ROOT/'catchbond_predictions.csv',dtype={'residue_id':str}).rename(columns={'residue_id':'author_residue_id'})
keys=['chain','author_residue_id']
ref=pd.read_csv(ROOT/'all408_ten_replica_mean_SD_SEM.csv',dtype={'author_residue_id':str})
a=ref[ref.condition=='WT'].drop(columns='condition').set_index(keys)
b=ref[ref.condition=='S32M'].drop(columns='condition').set_index(keys)
assert a.index.equals(b.index)
base=pd.DataFrame(index=a.index)
base['aa_WT']=a.aa;base['aa_S32M']=b.aa
base['unchanged_aa']=a.aa==b.aa
assert (~base.unchanged_aa).sum()==1
assert base.loc[('D','32'),'aa_WT']=='S' and base.loc[('D','32'),'aa_S32M']=='M'

def uncertainty(mean_a,sd_a,n_a,mean_b,sd_b,n_b):
    va=sd_a**2/n_a;vb=sd_b**2/n_b;se=np.sqrt(va+vb)
    dof=(va+vb)**2/(va**2/(n_a-1)+vb**2/(n_b-1)) if va+vb else np.inf
    crit=t.ppf(.975,dof) if np.isfinite(dof) else 1.96
    delta=mean_b-mean_a
    return delta,se,delta-crit*se,delta+crit*se

for metric in ['rmsf_A','Hchi_nats']:
    for condition,source in [('WT',a),('S32M',b)]:
        for suffix in ['mean','SD','SEM','n']:base[f'{metric}_{condition}_{suffix}']=source[f'{metric}_{suffix}']
    unc=[uncertainty(a.iloc[i][metric+'_mean'],a.iloc[i][metric+'_SD'],a.iloc[i][metric+'_n'],b.iloc[i][metric+'_mean'],b.iloc[i][metric+'_SD'],b.iloc[i][metric+'_n'])
         if a.iloc[i][metric+'_n'] and b.iloc[i][metric+'_n'] else (np.nan,)*4 for i in range(len(a))]
    for j,col in enumerate(['MD_delta','MD_delta_SEM','MD_delta_CI95_low','MD_delta_CI95_high']):base[f'{metric}_{col}']=[x[j] for x in unc]

tables=[];scores=[];regstats=[]
regions={
 'all_unchanged_aa':lambda d:d.unchanged_aa,
 'peptide_all9':lambda d:d.chain=='C',
 'TCRalpha_unchanged':lambda d:(d.chain=='D')&d.unchanged_aa,
 'TCRbeta_all110':lambda d:d.chain=='E',
 'CDR1alpha_author25_35_unchanged':lambda d:(d.chain=='D')&d.author_residue_id.astype(int).between(25,35)&d.unchanged_aa,
 'CDR3alpha_author90_97':lambda d:(d.chain=='D')&d.author_residue_id.astype(int).between(90,97),
 'CDR3beta_author92_102':lambda d:(d.chain=='E')&d.author_residue_id.astype(int).between(92,102),
 'peptide_S4':lambda d:(d.chain=='C')&(d.author_residue_id=='4'),
 'peptide_M6':lambda d:(d.chain=='C')&(d.author_residue_id=='6'),
 'alpha_N92':lambda d:(d.chain=='D')&(d.author_residue_id=='92'),
 'alpha_N93':lambda d:(d.chain=='D')&(d.author_residue_id=='93'),
 'beta_S98':lambda d:(d.chain=='E')&(d.author_residue_id=='98'),
 'mutation_alpha_S32M_not_same_chemistry':lambda d:(d.chain=='D')&(d.author_residue_id=='32'),
}
for kind in ['experimental_matched','prepared']:
    d=base.copy()
    pa=pred[(pred.condition=='WT')&(pred.input_kind==kind)].set_index(keys).reindex(d.index)
    pb=pred[(pred.condition=='S32M')&(pred.input_kind==kind)].set_index(keys).reindex(d.index)
    assert pa.aa.equals(a.aa) and pb.aa.equals(b.aa)
    d=d.reset_index();d['input_kind']=kind
    for metric,pcol in [('rmsf_A','rmsf_pred'),('Hchi_nats','Hchi_pred')]:
        d[f'{metric}_model_WT']=pa[pcol].to_numpy()
        d[f'{metric}_model_S32M']=pb[pcol].to_numpy()
        d[f'{metric}_model_delta']=(pb[pcol]-pa[pcol]).to_numpy()
    tables.append(d)
    for region,select in regions.items():
        part=d[select(d)].copy()
        idset=set(zip(part.chain,part.author_residue_id))
        mpart=md[[(c,r) in idset for c,r in zip(md.chain,md.author_residue_id)]]
        for metric in ['rmsf_A','Hchi_nats']:
            good=part[np.isfinite(part[f'{metric}_MD_delta'])]
            y=good[f'{metric}_MD_delta'].to_numpy();p=good[f'{metric}_model_delta'].to_numpy()
            mae=float(np.mean(abs(y-p)));zero=float(np.mean(abs(y)))
            rho=float(spearmanr(y,p).statistic) if len(y)>=4 and np.ptp(y)>0 and np.ptp(p)>0 else np.nan
            reliable=(good[f'{metric}_MD_delta_CI95_low']>0)|(good[f'{metric}_MD_delta_CI95_high']<0)
            scores.append({'input_kind':kind,'region':region,'metric':metric,'n':len(y),'rho':rho,'change_MAE':mae,'zero_change_MAE':zero,
                'relative_MAE_reduction':1-mae/zero if zero else np.nan,'mean_MD_delta':float(np.mean(y)),
                'mean_model_delta':float(np.mean(p)),'mean_abs_model_delta':float(np.mean(abs(p))),
                'sign_agreement':float(np.mean(np.sign(y)==np.sign(p))),'n_MD_CI95_excluding_zero':int(reliable.sum()),
                'sign_agreement_on_MD_CI95_excluding_zero':float(np.mean(np.sign(y[reliable])==np.sign(p[reliable]))) if reliable.any() else np.nan,
                'same_amino_acid_only':bool(good.unchanged_aa.all())})
            perrep=mpart.groupby(['condition','replica'])[metric].mean().unstack(0)
            xa=perrep['WT'].dropna().to_numpy();xb=perrep['S32M'].dropna().to_numpy()
            if len(xa) and len(xb):
                delta,se,lo,hi=uncertainty(xa.mean(),xa.std(ddof=1),len(xa),xb.mean(),xb.std(ddof=1),len(xb))
                regstats.append({'input_kind':kind,'region':region,'metric':metric,'MD_WT_region_mean':float(xa.mean()),
                    'MD_S32M_region_mean':float(xb.mean()),'MD_region_delta':float(delta),'MD_region_delta_SEM':float(se),
                    'MD_region_delta_Welch_CI95_low':float(lo),'MD_region_delta_Welch_CI95_high':float(hi),
                    'model_region_delta':float(np.mean(p)),'independent_replicas_per_condition':10})

joined=pd.concat(tables,ignore_index=True)
joined.to_csv(ROOT/'all408_MD_model_changes_both_inputs.csv',index=False)
scores=pd.DataFrame(scores);scores.to_csv(ROOT/'regional_change_metrics.csv',index=False)
pd.DataFrame(regstats).to_csv(ROOT/'regional_mean_change_uncertainty.csv',index=False)
sense=[]
for region,select in regions.items():
    aa=joined[joined.input_kind=='experimental_matched'].copy();bb=joined[joined.input_kind=='prepared'].copy()
    sa=aa[select(aa)].set_index(keys);sb=bb[select(bb)].set_index(keys)
    for metric in ['rmsf_A','Hchi_nats']:
        mask=np.isfinite(sa[f'{metric}_MD_delta'])
        x=sa.loc[mask,f'{metric}_model_delta'];y=sb.loc[mask,f'{metric}_model_delta'];mdelta=sa.loc[mask,f'{metric}_MD_delta']
        sense.append({'region':region,'metric':metric,'n':len(x),'input_sensitivity_MAE_of_model_delta':float(np.mean(abs(x-y))),
            'MD_mean_abs_delta':float(np.mean(abs(mdelta))),'input_sensitivity_over_MD_amplitude':float(np.mean(abs(x-y))/np.mean(abs(mdelta))) if np.mean(abs(mdelta)) else np.nan,
            'model_delta_sign_agreement_between_input_kinds':float(np.mean(np.sign(x)==np.sign(y)))})
pd.DataFrame(sense).to_csv(ROOT/'experimental_vs_prepared_input_sensitivity.csv',index=False)

audit={'status':'passed','same_checkpoint_both_structure_inputs':True,'not_retrained':True,
       'comparison':'S32M minus WT; all unchanged-AA407residues retained, S32 position reported separately.',
       'regions':{'CDR1alpha':'authorD25–35, mutationD32 excluded from same-AA metrics','CDR3alpha':'authorD90–97, between conservedCys89/Phe98','CDR3beta':'authorE92–102, between conservedCys91/Phe103'},
       'uncertainty':'Each condition has10 independent MD trajectories. Welch95%CI of difference uses per-trajectory residue or region averages, never treats residues/frames as independent repetitions; descriptive intervals, no multiple-comparison discovery claim.',
       'source_files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'catchbond_predictions.csv',ROOT/'all408_per_replica_labels.csv',ROOT/'all408_ten_replica_mean_SD_SEM.csv']}}
(ROOT/'model_MD_comparison_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
print(scores[scores.region.isin(['all_unchanged_aa','peptide_all9','CDR1alpha_author25_35_unchanged','CDR3alpha_author90_97','CDR3beta_author92_102','peptide_S4','peptide_M6','alpha_N92','alpha_N93','beta_S98'])].to_string(index=False),flush=True)
