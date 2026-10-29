"""Descriptive, all-case analyses; no hypothesis selection hidden in plotting."""
from __future__ import annotations

import json
import re
import shlex
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/apo_holo_exploration_20260913'
ANALYSIS = OUT / 'analysis'
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'scripts')]
from flexppi.constants import THREE_TO_ONE
from flexppi.data.structures import parse_pdb
from flexppi.data.pairing import global_sequence_alignment

BM_LABEL = {15183:'CaMKK',15184:'CaMKI',15185:'eNOS',15186:'smMLCK',15187:'PDE',15188:'free',15191:'nNOS'}


def corr(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x)&np.isfinite(y)
    if ok.sum()<4 or np.std(x[ok])<1e-10 or np.std(y[ok])<1e-10:
        return np.nan
    return float(spearmanr(x[ok],y[ok]).statistic)


def partial_rank_correlation(x,y,distance,aa):
    """Exploratory association after linear adjustment of ranks for distance/AA.

    This is not a trained predictive baseline or an out-of-sample result.
    """
    x,y,distance=np.asarray(x,float),np.asarray(y,float),np.asarray(distance,float)
    aa=np.asarray(aa)
    ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(distance)
    x,y,distance,aa=x[ok],y[ok],distance[ok],aa[ok]
    if len(x)<10:
        return np.nan
    controls=np.column_stack([np.ones(len(x)),rankdata(distance),pd.get_dummies(aa,drop_first=True).to_numpy(float)])
    rx,ry=rankdata(x),rankdata(y)
    rx=rx-controls@np.linalg.lstsq(controls,rx,rcond=None)[0]
    ry=ry-controls@np.linalg.lstsq(controls,ry,rcond=None)[0]
    return float(np.corrcoef(rx,ry)[0,1]) if min(np.std(rx),np.std(ry))>1e-10 else np.nan


def loops(text, prefix):
    for match in re.finditer(r'(?m)^\s*loop_\s*$([\s\S]*?)^\s*stop_\s*$', text):
        lines = [line.strip() for line in match[1].splitlines() if line.strip()]
        tags = [line for line in lines if line.startswith('_')]
        if not tags or not tags[0].startswith(prefix):
            continue
        tokens = shlex.split('\n'.join(lines[len(tags):]))
        if len(tokens)%len(tags):
            raise ValueError('STAR loop token count is not a multiple of its tag count')
        for start in range(0,len(tokens),len(tags)):
            yield {key.split('.',1)[1]:value for key,value in zip(tags,tokens[start:start+len(tags)])}


def parse_bmrb():
    rows, entities, audit = [], {}, []
    for number, label in BM_LABEL.items():
        path = OUT/'raw'/f'bmr{number}_3.str'
        text = path.read_text()
        ca_id, ca_seq, peptide_seq = None, None, None
        for block in re.findall(r'(?m)^save_\S+\s*$([\s\S]*?)^save_\s*$',text):
            name = re.search(r'(?m)^\s*_Entity.Name\s+(.+)$',block)
            entity_id = re.search(r'(?m)^\s*_Entity.ID\s+(\S+)',block)
            sequence = re.search(r'_Entity.Polymer_seq_one_letter_code\s*\n;\s*\n([\s\S]*?)\n;',block)
            if name and entity_id and sequence:
                seq = re.sub(r'\s+', '', sequence[1])
                if 'calmodulin' in name[1].lower():
                    ca_id, ca_seq = entity_id[1], seq
                else:
                    peptide_seq = seq
        if ca_id is None:
            raise ValueError(f'Calmodulin entity not found in {number}')
        entities[label] = dict(calmodulin_sequence=ca_seq, peptide_sequence=peptide_seq, bmrb=number)
        for row in loops(text, '_Order_param.'):
            if row['Entity_ID'] != ca_id or row['Order_param_val'] in ('.','?'):
                continue
            position = int(row['Comp_index_ID'])
            aa = THREE_TO_ONE[row['Comp_ID']]
            if ca_seq[position-1] != aa:
                raise ValueError(f'BMRB sequence/observable mismatch {number}:{position}')
            rows.append(dict(bmrb=number, label=label, seq_id=position, amino_acid=aa, atom=row['Atom_ID'],
                observable='backbone' if row['Atom_ID']=='N' else 'methyl', s2=float(row['Order_param_val']),
                fit_error=float(row['Order_param_val_fit_err']) if row['Order_param_val_fit_err'] not in ('.','?') else np.nan))
        audit.append(dict(bmrb=number, label=label, entity_id=ca_id, ca_length=len(ca_seq), peptide_sequence=peptide_seq,
                          note='15184 title mislabeled; mapping uses entity CaMK1p and peptide sequence' if number==15184 else ''))
    df = pd.DataFrame(rows)
    df.to_csv(ANALYSIS/'calmodulin_nmr_order_parameters.csv',index=False)
    (ANALYSIS/'calmodulin_entities.json').write_text(json.dumps(audit,indent=2))
    return df, entities


def md_evaluate(df):
    rows = []
    for case, group in df[df.direction=='md_delta'].groupby('case_id'):
        for metric in ('rmsf','chi_entropy'):
            for predictor in ('delta','context','geometry','zero'):
                y = group['true_delta_'+metric].to_numpy()
                p = np.zeros(len(group)) if predictor=='zero' else group[predictor+'_'+metric].to_numpy()
                ok = np.isfinite(p)&np.isfinite(y)
                y,p = y[ok],p[ok]
                active = np.abs(y)>=.1
                rows.append(dict(case_id=case, group=group.group.iloc[0], metric=metric,predictor=predictor,n=len(y),
                    spearman=corr(p,y),mae=np.mean(np.abs(p-y)),rmse=np.sqrt(np.mean((p-y)**2)),
                    true_sd=np.std(y),predicted_sd=np.std(p),
                    active_count=int(active.sum()), active_sign_accuracy=np.mean(np.sign(p[active])==np.sign(y[active])) if active.any() else np.nan))
    result = pd.DataFrame(rows)
    result.to_csv(ANALYSIS/'md_delta_per_case.csv',index=False)
    result.groupby(['group','metric','predictor']).mean(numeric_only=True).reset_index().to_csv(ANALYSIS/'md_delta_per_complex.csv',index=False)
    return result


def aa_difference(g, value, inside, outside, match_sasa=False):
    differences = []
    for _, residue in g[inside].iterrows():
        control = g[outside & (g.amino_acid==residue.amino_acid)]
        if match_sasa and len(control):
            logdiff = np.abs(np.log(control.isolated_sasa/residue.isolated_sasa))
            control = control.loc[logdiff[logdiff<=np.log(2)].nsmallest(3).index]
        if len(control) and np.isfinite(residue[value]):
            differences.append(residue[value]-control[value].mean())
    return float(np.nanmean(differences)) if differences else np.nan, len(differences)


def matched_entropy_area(g,inside,outside):
    """Ratio-of-sums after same-AA and approximately matched exposure controls."""
    target_h,target_a,control_h,control_a=[],[],[],[]
    for _,residue in g[inside].iterrows():
        control=g[outside&(g.amino_acid==residue.amino_acid)]
        logdiff=np.abs(np.log(control.isolated_sasa/residue.isolated_sasa))
        control=control.loc[logdiff[logdiff<=np.log(2)].nsmallest(3).index]
        if len(control):
            target_h.append(residue.iso_chi_entropy);target_a.append(residue.isolated_sasa)
            control_h.append(control.iso_chi_entropy.mean());control_a.append(control.isolated_sasa.mean())
    if not target_h:
        return np.nan
    return 100*(sum(target_h)/sum(target_a)-sum(control_h)/sum(control_a))


def preorganization(df):
    rows = []
    for cid, group in df[df.direction=='preorganization'].groupby('case_id'):
        for surface_cut in (1.,5.,10.):
            for burial_cut in (0.5,1.,5.):
                for omit_p in (False,True):
                    g=group[group.iso_chi_entropy.notna() & (~group.amino_acid.isin(['P']) if omit_p else True)].copy()
                    surface=g.isolated_sasa>=surface_cut
                    inside=surface&(g.burial>=burial_cut)
                    outside=surface&(g.burial<burial_cut)
                    if inside.sum()<3 or outside.sum()<3:
                        continue
                    aa,n = aa_difference(g,'iso_chi_entropy',inside,outside)
                    aasa,nasa = aa_difference(g,'iso_chi_entropy',inside,outside,True)
                    rows.append(dict(case_id=cid,surface_cut=surface_cut,burial_cut=burial_cut,omit_proline=omit_p,
                        interface_n=int(inside.sum()),surface_control_n=int(outside.sum()),
                        raw_difference=g.loc[inside,'iso_chi_entropy'].mean()-g.loc[outside,'iso_chi_entropy'].mean(),
                        entropy_per_area_difference=100*(g.loc[inside,'iso_chi_entropy'].sum()/g.loc[inside,'isolated_sasa'].sum()-g.loc[outside,'iso_chi_entropy'].sum()/g.loc[outside,'isolated_sasa'].sum()),
                        aa_matched_difference=aa,aa_matched_n=n,aa_sasa_matched_difference=aasa,aa_sasa_matched_n=nasa,
                        aa_sasa_matched_entropy_per_area_difference=matched_entropy_area(g,inside,outside)))
    result=pd.DataFrame(rows)
    result.to_csv(ANALYSIS/'preorganization_all_controls.csv',index=False)
    return result[(result.surface_cut==5)&(result.burial_cut==1)&(~result.omit_proline)]


def redistribution(df):
    rows=[]
    for cid,g in df[df.direction=='redistribution'].groupby('case_id'):
        for metric in ('rmsf','chi_entropy'):
            for cutoff in (4.,5.,6.):
                inside=g.partner_distance<=cutoff
                distal=g.partner_distance>12.
                value='delta_'+metric
                aa,n=aa_difference(g,value,inside,distal)
                for region,mask in [('interface',inside),('distal',distal),('whole',np.ones(len(g),bool))]:
                    rows.append(dict(case_id=cid,group=g.group.iloc[0],metric=metric,cutoff=cutoff,region=region,
                         mean_delta=g.loc[mask,value].mean(),n=int(g.loc[mask,value].notna().sum()),
                         context_mean=g.loc[mask,'context_'+metric].mean(),geometry_mean=g.loc[mask,'geometry_'+metric].mean(),
                         aa_adjusted_interface_minus_distal=aa,aa_matched_n=n))
    result=pd.DataFrame(rows)
    result.to_csv(ANALYSIS/'redistribution_regions.csv',index=False)
    return result[result.cutoff==5]


def cam_comparison(df,nmr,entities):
    result,paired,checks=[],[],[]
    free=nmr[nmr.label=='free']
    free_seq=entities['free']['calmodulin_sequence']
    apo=parse_pdb(str(OUT/'raw/1CLL.pdb'))
    amap=global_sequence_alignment(free_seq,apo.chain_sequence('A'))
    author_to_free={apo.author_residue_ids[b]:a+1 for a,b in zip(amap.indices_a,amap.indices_b)}
    mapping_dir=OUT/'mapping'
    for cid,g in df[df.direction=='backbone_sidechain'].groupby('case_id'):
        label=cid.removeprefix('cam_')
        ent=entities[label]
        alignment=global_sequence_alignment(ent['calmodulin_sequence'],free_seq)
        to_free={a+1:b+1 for a,b in zip(alignment.indices_a,alignment.indices_b)}
        bound=nmr[nmr.label==label].copy()
        bound['seq_id']=bound.seq_id.map(to_free)
        merged=bound.merge(free,on=['seq_id','amino_acid','atom','observable'],suffixes=('_bound','_free'))
        provenance=json.loads((mapping_dir/(cid+'.json')).read_text())
        holo=parse_pdb(provenance['case']['holo']['path'])
        partners=[c for c in provenance['case']['holo']['chains'] if c not in provenance['case']['target_chains']]
        pseq=''.join(holo.chain_sequence(c) for c in partners)
        paln=global_sequence_alignment(ent['peptide_sequence'],pseq)
        pcov=min(len(paln.indices_a)/len(ent['peptide_sequence']),len(paln.indices_b)/len(pseq))
        checks.append(dict(case_id=cid,peptide_sequence_pdb=pseq,peptide_sequence_bmrb=ent['peptide_sequence'],
                           peptide_identity=paln.identity,peptide_bidirectional_coverage=pcov))
        predicted=g.copy()
        predicted['seq_id']=predicted.apo_author_id.astype(str).map(author_to_free)
        # Average equivalent methyl probes only AFTER matched-probe bound/free subtraction.
        merged['minus_delta_s2']=-(merged.s2_bound-merged.s2_free)
        matched=merged.groupby(['seq_id','amino_acid','observable']).agg(
            minus_delta_s2=('minus_delta_s2','mean'),probe_count=('atom','size')).reset_index()
        matched=matched.merge(predicted,on=['seq_id','amino_acid'])
        for observable,metric in [('backbone','rmsf'),('methyl','chi_entropy')]:
            m=matched[(matched.observable==observable)&matched['delta_'+metric].notna()].copy()
            m['predicted_change']=m['delta_'+metric]
            paired.append(m[['case_id','seq_id','apo_author_id','amino_acid','observable','minus_delta_s2','predicted_change','probe_count']])
            result.append(dict(case_id=cid,observable=observable,n=len(m),
                spearman=corr(m.predicted_change,m.minus_delta_s2),
                mean_predicted_change=m.predicted_change.mean(),mean_minus_delta_s2=m.minus_delta_s2.mean(),
                sd_minus_delta_s2=m.minus_delta_s2.std(),peptide_match_pass=paln.identity>=.95 and pcov>=.85,
                distance_baseline_spearman=corr(m.partner_distance,m.minus_delta_s2),
                context_only_spearman=corr(m['context_'+metric],m.minus_delta_s2),
                geometry_only_spearman=corr(m['geometry_'+metric],m.minus_delta_s2),
                partial_rank_aa_distance=partial_rank_correlation(m.predicted_change,m.minus_delta_s2,m.partner_distance,m.amino_acid)))
    pd.DataFrame(checks).to_csv(ANALYSIS/'calmodulin_structure_experiment_match.csv',index=False)
    pd.concat(paired).to_csv(ANALYSIS/'calmodulin_prediction_nmr_paired.csv',index=False)
    result=pd.DataFrame(result)
    result.to_csv(ANALYSIS/'calmodulin_nmr_comparison.csv',index=False)
    return result


def pdz_response(df):
    g=df[df.direction=='remote_response'].copy()
    responsive='L18 I20 V22 V26 V30 A39 V40 V61 V64 L66 A69 L78 T81 V85'.split()
    actual=set(g.amino_acid+g.author_id.astype(str))
    if not set(responsive).issubset(actual):
        raise ValueError('PDZ experimental position/sequence mapping mismatch')
    g['known_responsive']=(g.amino_acid+g.author_id.astype(str)).isin(responsive)
    g['abs_delta_chi_entropy']=g.delta_chi_entropy.abs()
    g['abs_delta_rmsf']=g.delta_rmsf.abs()
    g['interface_proximity']=-g.partner_distance
    # Compare within the methyl-bearing AA universe; unknown != negative.
    eligible=g.amino_acid.isin(list('ILVAMT'))
    rows=[]
    for value in ('abs_delta_chi_entropy','abs_delta_rmsf','interface_proximity'):
        inside=eligible&g.known_responsive&g[value].notna()
        outside=eligible&~g.known_responsive&g[value].notna()
        effect,n=aa_difference(g,value,inside,outside)
        fractions=[]
        for _,r in g[inside].iterrows():
            control=g[outside&(g.amino_acid==r.amino_acid)][value]
            if len(control):
                fractions.append(np.mean(control<r[value])+.5*np.mean(control==r[value]))
        rows.append(dict(predictor=value,aa_matched_difference=effect,matched_responsive_n=n,
            mean_within_aa_percentile=np.mean(fractions),
            caution='Unannotated comparison sites have unknown experimental status; this is NOT ROC/AUC'))
    g.to_csv(ANALYSIS/'pdz2_residue_response.csv',index=False)
    result=pd.DataFrame(rows)
    result.to_csv(ANALYSIS/'pdz2_response_enrichment.csv',index=False)
    return result


def plots(md,red,pre,cam):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    for metric,color in [('rmsf','#53698a'),('chi_entropy','#6fa89b')]:
        g=md[md.metric==metric].groupby(['group','predictor']).mae.mean().unstack()
        axes[0,0].scatter(g.zero,g.delta,label={'rmsf':'RMSF (Å)','chi_entropy':'Hχ (nats)'}[metric],color=color,s=35)
    lim=max(*axes[0,0].get_xlim(),*axes[0,0].get_ylim())
    axes[0,0].plot([0,lim],[0,lim],'--',color='.6',lw=1)
    axes[0,0].set(xlabel='Zero-change baseline MAE',ylabel='FlexPPI change MAE',title='a  Paired MD: below diagonal is better')
    axes[0,0].legend(frameon=False)
    for j,metric in enumerate(['rmsf','chi_entropy']):
        g=red[red.metric==metric].groupby(['group','region']).mean_delta.mean().unstack()
        d=(g.interface-g.distal).dropna()
        axes[0,1].scatter(np.repeat(j,len(d))+np.linspace(-.14,.14,len(d)),d,color=['#53698a','#6fa89b'][j],s=25)
    axes[0,1].axhline(0,color='.6',ls='--',lw=1)
    axes[0,1].set(xticks=[0,1],xticklabels=['RMSF (Å)','Hχ (nats)'],ylabel='Interface minus distal mean change',title='b  Redistribution: one point per complex')
    for j,column in enumerate(['raw_difference','aa_matched_difference','aa_sasa_matched_difference']):
        g=pre[column].dropna()
        axes[1,0].scatter(np.repeat(j,len(g))+np.linspace(-.16,.16,len(g)),g,s=22,color='#809888')
    axes[1,0].axhline(0,color='.6',ls='--',lw=1)
    axes[1,0].set(xticks=[0,1,2],xticklabels=['Raw','AA-matched','AA + exposure\nmatched'],ylabel='Interface minus other surface Hχ (nats)',title='c  Preorganization: controls can change conclusion')
    labels=cam.case_id.unique()
    for i,obs in enumerate(['backbone','methyl']):
        g=cam[cam.observable==obs].set_index('case_id').reindex(labels)
        valid=g.peptide_match_pass.to_numpy(bool)
        positions=np.arange(len(labels))+(i-.5)*.14
        axes[1,1].scatter(positions[valid],g.spearman[valid],label=obs,s=40,color=['#53698a','#ba8760'][i])
        axes[1,1].scatter(positions[~valid],g.spearman[~valid],marker='x',s=40,color='.6')
    methyl=cam[cam.observable=='methyl'].set_index('case_id').reindex(labels)
    passed=methyl.peptide_match_pass.to_numpy(bool)
    axes[1,1].scatter(np.arange(len(labels))[passed]+.22,methyl.distance_baseline_spearman[passed],
                      marker='D',s=25,facecolors='none',edgecolors='.35',label='methyl: distance only')
    axes[1,1].axhline(0,color='.6',ls='--',lw=1)
    axes[1,1].set(xticks=range(len(labels)),xticklabels=[x[4:] for x in labels],ylabel='Spearman: predicted change vs −ΔS²',ylim=(-1,1),title='d  CaM: crosses = peptide mismatch (excluded)')
    axes[1,1].legend(frameon=False)
    fig.suptitle('FlexPPI apo–holo exploratory screen | five-checkpoint ensemble\nCaM has near-identical training homologues: NOT independent validation',fontsize=12)
    fig.savefig(ANALYSIS/'exploratory_overview.png',dpi=180)
    fig.savefig(ANALYSIS/'exploratory_overview.svg')
    plt.close(fig)


def main():
    ANALYSIS.mkdir(exist_ok=True)
    raw=pd.read_csv(OUT/'residue_predictions.csv',dtype={'author_id':str,'apo_author_id':str})
    # Average models first. Folds are not biological replicates.
    keys=['case_id','direction','group','amino_acid','chain','author_id','apo_chain','apo_author_id','apo_index','holo_index','region']
    df=raw.groupby(keys,dropna=False).mean(numeric_only=True).reset_index().drop(columns='fold')
    df.to_csv(ANALYSIS/'ensemble_residue_predictions.csv',index=False)
    nmr,entities=parse_bmrb()
    md=md_evaluate(df)
    red=redistribution(df)
    pre=preorganization(df)
    cam=cam_comparison(df,nmr,entities)
    pdz=pdz_response(df)
    plots(md,red,pre,cam)
    print('MD macro summaries:\n',md.groupby(['metric','predictor'])[['spearman','mae','rmse']].mean().to_string())
    print('Preorganization negative counts:\n',(pre[['raw_difference','aa_matched_difference','aa_sasa_matched_difference']]<0).sum().to_string(), 'n=',len(pre))
    print('CaM:\n',cam.to_string(index=False))
    print('PDZ:\n',pdz.to_string(index=False))


if __name__=='__main__':
    main()
