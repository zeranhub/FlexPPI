"""Three separately scored, budget-matched exploratory screening uses.

No neural inference/training, no model-selection loop, no remote access. Pool
outcomes are held by the evaluation driver; policies receive observed labels
only. Full-pool feature normalization is explicitly label-free/transductive.
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import argparse,csv,json,hashlib,math
import numpy as np

HERE=Path(__file__).resolve().parent
CHECKPOINT='b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f'
SEEDS=list(range(20261022,20261042))
MODEL_BLOCKS={'sequence':['sequence'],'static':['static'],'baseline':['sequence','static'],
              'baseline_R':['sequence','static','R'],'baseline_H':['sequence','static','H'],
              'baseline_RH':['sequence','static','R','H']}
EXTRA_STRATEGIES=[]  # Optional independently bound pre-evaluation diagnostic controls.
SITE_STRATEGIES=['random_positions','static_contact_degree','static_relative_SASA','parent_R_percentile',
                 'parent_H_percentile','parent_joint_percentile','parent_joint_low','parent_partner_context_response']
AA=set('ACDEFGHIKLMNPQRSTVWY')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def dump_json(path,x):Path(path).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
def read_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h))
def dump_csv(path,rows):
    rows=list(rows);keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=keys);w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()})
def yes(v):
    if v is True or str(v).lower() in {'true','1','yes','passed'}:return True
    if v is False or str(v).lower() in {'false','0','no','failed'}:return False
    raise ValueError('Invalid boolean: '+repr(v))
def finite(v):
    try:r=float(v)
    except (TypeError,ValueError):return None
    return r if math.isfinite(r) else None
def stable_order(values,indices,ids):return sorted(indices,key=lambda i:(-float(values[i]),ids[i]))

def rank_ties(v):
    order=np.argsort(v,kind='stable');ranks=np.empty(len(v));start=0
    while start<len(v):
        end=start+1
        while end<len(v) and v[order[end]]==v[order[start]]:end+=1
        ranks[order[start:end]]=(start+end-1)/2+1;start=end
    return ranks

def diagnostic_rho(y,p):
    y,p=np.asarray(y),np.asarray(p)
    if len(y)<3 or len(set(y))<2 or len(set(p))<2:return None
    return float(np.corrcoef(rank_ties(y),rank_ties(p))[0,1])

def pairwise_sq(x):
    x=np.asarray(x,dtype=float)
    assert x.ndim==2 and np.isfinite(x).all()
    d=np.maximum(np.sum(x*x,axis=1)[:,None]+np.sum(x*x,axis=1)[None,:]-2*x@x.T,0.)
    np.fill_diagonal(d,0.)
    return d

def distance_kernel(d):
    n=len(d);positive=d[np.triu_indices(n,k=1)];positive=positive[positive>1e-12]
    bandwidth=float(np.median(positive)) if len(positive) else 1.
    k=np.exp(-d/bandwidth)
    np.fill_diagonal(k,1.)
    return k,dict(median_positive_distance=bandwidth,constant_block=not len(positive),label_values_used=False)

def construct_kernels(sequences,static,R,H,Rscale,Hscale):
    n=len(sequences);lengths={len(s) for s in sequences};assert len(lengths)==1
    assert all(set(s)<=AA for s in sequences)
    chars=np.array([list(s) for s in sequences])
    hd=np.mean(chars[:,None,:]!=chars[None,:,:],axis=2)
    kseq,sseq=distance_kernel(hd)
    static=np.asarray(static,dtype=float);assert static.shape==(n,2)
    scale=np.sqrt(np.mean(static**2,axis=0));scale[scale==0]=1.
    ks,ss=distance_kernel(pairwise_sq(static/scale));ss['unlabeled_full_pool_RMS_scales']=scale.tolist()
    blocks={'sequence':kseq,'static':ks};settings={'sequence':sseq,'static':ss}
    for name,x,scale in [('R',R,Rscale),('H',H,Hscale)]:
        x=np.asarray(x,dtype=float);scale=np.asarray(scale,dtype=float)
        assert x.ndim==2 and x.shape[0]==n and x.shape[1]>0 and scale.shape==(x.shape[1],)
        assert np.isfinite(x).all() and np.isfinite(scale).all() and np.all(scale>0)
        k,s=distance_kernel(pairwise_sq(x/scale));blocks[name]=k;s['profile_dimensions']=x.shape[1]
        settings[name]=s
    kernels={name:sum(blocks[b] for b in names)/len(names) for name,names in MODEL_BLOCKS.items()}
    return kernels,settings

def ridge_predict(kernel,observed,observed_y,test):
    """No pool outcomes are accepted; fixed mean-loss penalty lambda0.1."""
    observed=np.asarray(observed,dtype=int);test=np.asarray(test,dtype=int);y=np.asarray(observed_y,dtype=float)
    assert len(observed)==len(y)>0 and len(set(observed))==len(observed)
    assert not set(observed).intersection(test) and np.isfinite(y).all()
    intercept=float(np.mean(y));centered=y-intercept
    alpha=np.linalg.solve(kernel[np.ix_(observed,observed)]+len(observed)*.1*np.eye(len(observed)),centered)
    return intercept+kernel[np.ix_(test,observed)]@alpha

def novelty_one(kernel,observed,available,ids):
    """Kernel-induced coverage distance, not predictive uncertainty."""
    if not observed:return min(available,key=lambda i:ids[i])
    # Normalized kernels have diagonal1. Farthest min distance is smallest
    # maximum kernel similarity to any already measured/selected candidate.
    return min(available,key=lambda i:(float(np.max(kernel[i,observed])),ids[i]))

def greedy_diversity(kernel,initial,available,budget,ids):
    chosen=[];available=set(available);covered=list(initial)
    for _ in range(min(budget,len(available))):
        i=novelty_one(kernel,covered,sorted(available),ids);chosen.append(i);covered.append(i);available.remove(i)
    return chosen

def choose_batch(strategy,kernels,observed,observed_y,available,budget,ids,rng,feedback):
    """Selection API deliberately has no hidden y vector, hits or threshold."""
    available=sorted(available);budget=min(budget,len(available))
    if strategy=='random':return rng.permutation(available)[:budget].tolist()
    if strategy=='sequence_diversity':return greedy_diversity(kernels['sequence'],observed,available,budget,ids)
    kernel=kernels[strategy]
    pred=ridge_predict(kernel,observed,observed_y,available)
    score=np.full(len(ids),np.nan);score[available]=pred
    take=budget-1 if feedback and budget>1 else budget
    chosen=stable_order(score,available,ids)[:take]
    if feedback and budget>1:
        remaining=[i for i in available if i not in chosen]
        chosen.append(novelty_one(kernel,list(observed)+chosen,remaining,ids))
    return chosen

def outcome_record(panel,mode,strategy,seed,initial,chosen,y,hits,ids,budget=None):
    all_ids=list(initial)+list(chosen)
    first_new=None
    if mode=='mode3':
        for end in range(5,len(chosen)+1,5):
            if np.any(hits[np.asarray(chosen[:end],dtype=int)]):first_new=len(initial)+end;break
    else:
        for i,index in enumerate(chosen,1):
            if hits[index]:first_new=len(initial)+i;break
    return dict(panel_id=panel['panel_id'],family_id=panel['family_id'],task_object=panel['task_object'],endpoint=panel['endpoint'],
        mode=mode,strategy=strategy,seed=seed,initial_n=len(initial),extra_assays=len(chosen),total_assays=len(all_ids),
        initial_hits=int(np.sum(hits[initial])),new_hits=int(np.sum(hits[chosen])) if chosen else 0,
        cumulative_hits=int(np.sum(hits[all_ids])),best_extra=float(np.max(y[chosen])) if chosen else None,
        cumulative_best=float(np.max(y[all_ids])),first_new_hit_total_assay_count=first_new,
        threshold_found_in_initial=bool(np.any(hits[initial])),
        first_threshold_total_assay_count=len(initial) if np.any(hits[initial]) else first_new,
        threshold_not_found_by_current_budget=not bool(np.any(hits[all_ids])),
        first_new_threshold_right_censored=first_new is None,
        initial_candidate_ids=[ids[i] for i in initial],selected_candidate_ids=[ids[i] for i in chosen],
        hit_definition=panel['hit_definition'],hit_is_above_parent_improvement=panel.get('hit_is_above_parent_improvement',False),
        independent_replication='biological_family;20seeds are paired sampling variability only')

def make_hits(panel,label_rows,y):
    kind=panel['hit_kind']
    assert yes(panel['threshold_fixed_without_model_outputs'])
    if kind=='source_binary':
        values=[r['is_hit'] for r in label_rows];assert all(v in {'0','1','False','True','false','true'} for v in values)
        return np.array([yes(v) for v in values],dtype=bool),dict(kind=kind,cutoff=None)
    if kind=='fixed_threshold':
        threshold=float(panel['declared_hit_threshold'])
        return y>=threshold,dict(kind=kind,cutoff=threshold)
    if kind=='relative_top_decile':
        assert not panel.get('hit_is_above_parent_improvement',False)
        cutoff=float(np.quantile(y,.9))
        return y>=cutoff,dict(kind=kind,cutoff=cutoff,not_absolute_improvement=True,only_offline_scoring_knows_cutoff=True)
    raise ValueError('Undefined prespecified hit rule '+str(kind))

def load_panel(panel):
    lp=Path(panel['labels_csv']);fp=Path(panel['features_npz'])
    assert sha(lp)==panel['labels_sha256'] and sha(fp)==panel['features_sha256']
    labels=read_csv(lp);labels=[r for r in labels if r['panel_id']==panel['panel_id']]
    assert labels and len({r['candidate_id'] for r in labels})==len(labels)
    assert yes(panel['complete_candidate_pool']) and not yes(panel['positive_only_selection'])
    if any(not yes(r['label_qc_pass']) for r in labels):raise ValueError('SourceQC exclusions require a prebound common cohort, not runtime outcome selection')
    if any(r['censoring'] not in {'none','none_reported','uncensored'} for r in labels):raise ValueError('Censored quantitative pool: provide a separately prespecified source-binary endpoint or mark this panel not_evaluable; no boundary imputation')
    byid={r['candidate_id']:r for r in labels}
    with np.load(fp,allow_pickle=False) as z:
        assert str(z['checkpoint_sha256'].item())==CHECKPOINT
        ids=[str(x) for x in z['candidate_ids']];seq=[str(x) for x in z['sequences']]
        assert len(ids)==len(set(ids)) and set(ids)==set(byid)
        ordered=[byid[c] for c in ids]
        assert all(r['sequence']==s for r,s in zip(ordered,seq))
        y=np.array([float(r['value_oriented_higher_better']) for r in ordered]);assert np.isfinite(y).all()
        assert all(c!=r['parent_candidate_id'] for c,r in zip(ids,ordered)),'Parent not a selectable candidate'
        kernels,settings=construct_kernels(seq,z['static_values'],z['R_pattern'],z['H_pattern'],z['R_training_scale'],z['H_training_scale'])
        assert bool(z['unchanged_spectator_U_verified'].item())
        assert bool(z['all_library_variable_positions_excluded'].item())
        scores=np.array(z['zero_shot_scores'],dtype=float)
        names=[str(x) for x in z['zero_shot_score_names']]
        assert names==['HERMES_fixed','ProteinMPNN'] and scores.shape==(len(ids),2) and np.isfinite(scores).all()
        assert bool(z['zero_shot_scores_oriented_higher_better'].item())
        bindings={k:str(z[k].item()) for k in ['template_sha256','regions_sha256','training_reference_sha256']}
    hits,hitinfo=make_hits(panel,ordered,y)
    return ids,y,hits,kernels,settings,scores,bindings,hitinfo

def evaluate_candidate_panel(panel,output):
    ids,y,hits,kernels,settings,zero_scores,bindings,hitinfo=load_panel(panel);n=len(ids)
    records=[];selections=[];gates=[];trace=[];diagnostics=[]
    for mode,minn in [('mode2',60),('mode3',100)]:
        if not yes(panel.get(mode+'_eligible',False)) or n<minn:
            gates.append(dict(panel_id=panel['panel_id'],mode=mode,status='not_evaluable',reason='Curator data contract ineligible or candidate count below the prespecified gate',N=n));continue
        initial_n=min(40,max(12 if mode=='mode2' else 20,math.ceil(.1*n)))
        budget=min(30,n-initial_n) if mode=='mode2' else 40
        if n-initial_n<budget:
            gates.append(dict(panel_id=panel['panel_id'],mode=mode,status='not_evaluable',reason='Insufficient unseen pool for all8 complete5-candidate rounds',N=n));continue
        favorable_usable=np.sum(hits)>=5 and np.sum(~hits)>=5
        gates.append(dict(panel_id=panel['panel_id'],family_id=panel['family_id'],mode=mode,status='evaluable',N=n,initial_n=initial_n,extra_budget=budget,n_favorable=int(np.sum(hits)),n_nonfavorable=int(np.sum(~hits)),favorable_discovery_claim_eligible=bool(favorable_usable)))
        strategies=['random','sequence_diversity']+list(MODEL_BLOCKS) if mode=='mode2' else ['random','sequence_diversity','baseline','baseline_R','baseline_H','baseline_RH']
        if mode=='mode2':strategies+=['published_HERMES_fixed','published_ProteinMPNN']
        strategies+=EXTRA_STRATEGIES
        for seed in SEEDS:
            initial=np.sort(np.random.default_rng(seed).choice(n,size=initial_n,replace=False)).tolist()
            for strategy in strategies:
                observed=list(initial);chosen=[];available=set(range(n))-set(initial)
                rng=np.random.default_rng(seed+104729)
                if mode=='mode2':
                    if strategy.startswith('published_'):
                        col=0 if strategy=='published_HERMES_fixed' else 1
                        chosen=stable_order(zero_scores[:,col],available,ids)[:budget]
                    else:chosen=choose_batch(strategy,kernels,initial,y[initial],available,budget,ids,rng,False)
                    assert len(set(chosen))==budget and not set(chosen).intersection(initial)
                    for k in sorted({k for k in [5,10,20,30,budget] if k<=budget}):
                        records.append(outcome_record(panel,mode,strategy,seed,initial,chosen[:k],y,hits,ids))
                    hidden=sorted(available)
                    if strategy in MODEL_BLOCKS:
                        estimates=ridge_predict(kernels[strategy],initial,y[initial],hidden)
                        diagnostics.append(dict(panel_id=panel['panel_id'],mode=mode,strategy=strategy,seed=seed,
                            n_hidden=len(hidden),MAE=float(np.mean(np.abs(estimates-y[hidden]))),Spearman=diagnostic_rho(y[hidden],estimates),
                            metric_not_used_for_policy_or_model_selection=True))
                    elif strategy.startswith('published_'):
                        estimates=zero_scores[hidden,0 if strategy=='published_HERMES_fixed' else 1]
                        diagnostics.append(dict(panel_id=panel['panel_id'],mode=mode,strategy=strategy,seed=seed,
                            n_hidden=len(hidden),MAE=None,Spearman=diagnostic_rho(y[hidden],estimates),
                            raw_score_not_in_phenotype_units=True,metric_not_used_for_policy_or_model_selection=True))
                else:
                    records.append(outcome_record(panel,mode,strategy,seed,initial,[],y,hits,ids))
                    for roundid in range(1,9):
                        # Hidden outcomes are not passed to policy. The source
                        # array is indexed only at already measured candidates.
                        batch=choose_batch(strategy,kernels,observed,y[observed],available,5,ids,rng,True)
                        assert len(batch)==5 and len(set(batch))==5 and not set(batch).intersection(observed)
                        trace.append(dict(panel_id=panel['panel_id'],mode=mode,strategy=strategy,seed=seed,round=roundid,
                            observed_candidate_ids_before=[ids[i] for i in observed],selected_batch=[ids[i] for i in batch],
                            no_within_batch_feedback=True,observed_label_count=len(observed)))
                        chosen.extend(batch);observed.extend(batch);available-=set(batch)
                        records.append(outcome_record(panel,mode,strategy,seed,initial,chosen,y,hits,ids))
                selections.append(dict(panel_id=panel['panel_id'],family_id=panel['family_id'],mode=mode,strategy=strategy,seed=seed,
                    initial_candidate_ids=[ids[i] for i in initial],selected_candidate_ids=[ids[i] for i in chosen],
                    initial_n=initial_n,extra_budget=budget,total_budget=initial_n+budget))
    folder=output/panel['panel_id'];folder.mkdir(parents=True,exist_ok=True)
    dump_csv(folder/'budget_metrics.csv',records);dump_csv(folder/'selected_candidates.csv',selections);dump_csv(folder/'feedback_information_trace.csv',trace)
    dump_csv(folder/'one_shot_hidden_pool_diagnostics.csv',diagnostics)
    curves=defaultdict(list)
    for r in records:
        if r['mode']=='mode3':curves[r['strategy'],r['seed']].append(r)
    auc=[]
    for (strategy,seed),rr in curves.items():
        rr=sorted(rr,key=lambda r:r['extra_assays']);xx=np.array([r['extra_assays'] for r in rr]);yy=np.array([r['new_hits'] for r in rr])
        auc.append(dict(panel_id=panel['panel_id'],family_id=panel['family_id'],strategy=strategy,seed=seed,
            extra_assay_budget=40,cumulative_hit_curve_area=float(np.trapz(yy,xx)),
            only_batch_endpoints_observed=True,descriptive_seed_variability_not_biological_inference=True))
    dump_csv(folder/'feedback_cumulative_hit_area.csv',auc)
    dump_json(folder/'data_eligibility.json',gates)
    dump_json(folder/'kernel_binding.json',dict(kernel_settings=settings,source_bindings=bindings,hit_rule=hitinfo,
        labels_sha256=panel['labels_sha256'],features_sha256=panel['features_sha256'],all_pool_scaling_label_free=True))
    return records,gates

def evaluate_site_panel(panel,output):
    assert panel['task_object']=='TCR_mutation' or panel['task_object']=='antigen_peptide_mutation'
    labels=read_csv(panel['site_labels_csv']);scores=read_csv(panel['site_scores_csv'])
    assert sha(panel['site_labels_csv'])==panel['site_labels_sha256'] and sha(panel['site_scores_csv'])==panel['site_scores_sha256']
    labels=[r for r in labels if r['panel_id']==panel['panel_id']];scores={r['position_key']:r for r in scores if r['panel_id']==panel['panel_id']}
    positions=defaultdict(dict)
    for r in labels:
        if not yes(r['label_qc_pass']) or r['censoring'] not in {'none','none_reported','uncensored'}:continue
        pos=r['position_key'];aa=r['mutant_aa'];assert aa in AA and aa!=r['parent_aa']
        assert aa not in positions[pos];positions[pos][aa]=r
    complete={p:rr for p,rr in positions.items() if set(rr)==AA-{next(iter(rr.values()))['parent_aa']}}
    if len(complete)>=6:positions=complete;shared=None;substitution_n=19
    else:
        shared=set.intersection(*(set(rr) for rr in positions.values())) if positions else set()
        if len(positions)<6 or len(shared)<10:
            return [],[dict(panel_id=panel['panel_id'],mode='mode1',status='not_evaluable',reason='No equal-budget saturation coverage: need6 complete19-substitution positions or6 positions sharing>=10 substitutions.')]
        positions={p:{aa:r for aa,r in rr.items() if aa in shared} for p,rr in positions.items()};substitution_n=len(shared)
    assert all(p in scores for p in positions)
    posids=sorted(positions);all_rows=[r for p in posids for r in positions[p].values()]
    y=np.array([float(r['value_oriented_higher_better']) for r in all_rows]);hits,hitinfo=make_hits(panel,all_rows,y)
    outcome={}
    for r,value,hit in zip(all_rows,y,hits):outcome[r['candidate_id']]=(float(value),bool(hit))
    primary={}
    for p in posids:
        s=scores[p];R=float(s['parent_R_percentile']);assert 0<=R<=1
        H=.5 if s['parent_aa'] in {'G','A'} else float(s['parent_H_percentile']);assert 0<=H<=1
        joint=R if s['parent_aa'] in {'G','A'} else (R+H)/2
        primary[p]={'static_contact_degree':float(s['static_contact_degree']),'static_relative_SASA':float(s['static_relative_SASA']),
            'parent_R_percentile':R,'parent_H_percentile':H,'parent_joint_percentile':joint,'parent_joint_low':1-joint,
            'parent_partner_context_response':finite(s.get('parent_partner_context_response'))}
    records=[];gates=[]
    for k in [2,5]:
        if k==5 and len(posids)<10:continue
        for strategy in SITE_STRATEGIES:
            if strategy!='random_positions' and any(primary[p][strategy] is None for p in posids):
                gates.append(dict(panel_id=panel['panel_id'],mode='mode1',strategy=strategy,status='companion_not_evaluable',reason='Predeclared partner-removed data unavailable; no substitute score chosen.'));continue
            seeds=SEEDS if strategy=='random_positions' else [None]
            for seed in seeds:
                chosen=np.random.default_rng(seed).choice(posids,size=k,replace=False).tolist() if seed else sorted(posids,key=lambda p:(-primary[p][strategy],p))[:k]
                selected=[r for p in chosen for r in positions[p].values()];values=[outcome[r['candidate_id']][0] for r in selected];hh=[outcome[r['candidate_id']][1] for r in selected]
                records.append(dict(panel_id=panel['panel_id'],family_id=panel['family_id'],task_object=panel['task_object'],endpoint=panel['endpoint'],mode='mode1',strategy=strategy,seed=seed,
                    primary_strategy=strategy=='parent_joint_percentile',primary_budget=k==2,n_positions=k,substitutions_per_position=substitution_n,
                    assay_budget=k*substitution_n,new_hits=sum(hh),best_effect=max(values),mean_effect=sum(values)/len(values),
                    fraction_of_all_eligible_hits_recovered=sum(hh)/int(np.sum(hits)) if np.sum(hits) else None,
                    selected_positions=chosen,selected_candidate_ids=[r['candidate_id'] for r in selected],
                    hit_definition=panel['hit_definition'],not_KD_endpoint='enrichment' in panel['endpoint'].lower() or panel.get('endpoint_is_proxy',False)))
    folder=output/panel['panel_id'];folder.mkdir(parents=True,exist_ok=True)
    dump_csv(folder/'site_nomination_metrics.csv',records)
    dump_json(folder/'site_scope.json',dict(position_keys=posids,shared_substitutions=sorted(shared) if shared else 'all19_nonparent',substitution_n=substitution_n,hit_rule=hitinfo,
        scores_sha256=panel['site_scores_sha256'],labels_sha256=panel['site_labels_sha256'],no_model_score_used_for_eligibility=True))
    gates.append(dict(panel_id=panel['panel_id'],family_id=panel['family_id'],mode='mode1',status='evaluable',n_positions=len(posids),substitutions_per_position=substitution_n))
    return records,gates

def summarize_candidate_records(records,gates):
    index=defaultdict(dict)
    for r in records:index[r['panel_id'],r['mode'],r['extra_assays'],r['seed']][r['strategy']]=r
    paired=[]
    for (panel,mode,budget,seed),rows in sorted(index.items()):
        if budget==0 or 'baseline_RH' not in rows:continue
        joint=rows['baseline_RH']
        for comparator in ['baseline','baseline_R','baseline_H','random','sequence_diversity','published_HERMES_fixed','published_ProteinMPNN']:
            if comparator not in rows:continue
            other=rows[comparator]
            assert joint['initial_candidate_ids']==other['initial_candidate_ids']
            paired.append(dict(panel_id=panel,family_id=joint['family_id'],mode=mode,extra_assays=budget,seed=seed,comparator=comparator,
                joint_minus_comparator_new_hits=joint['new_hits']-other['new_hits'],joint_minus_comparator_best=joint['cumulative_best']-other['cumulative_best'],
                equal_initial_labels_and_total_budget=True))
    groups=defaultdict(list)
    for r in paired:groups[r['panel_id'],r['mode'],r['extra_assays'],r['comparator']].append(r)
    summary=[]
    for (p,m,b,c),rr in sorted(groups.items()):
        diff=[r['joint_minus_comparator_new_hits'] for r in rr];best=[r['joint_minus_comparator_best'] for r in rr]
        assert len(rr)==20
        summary.append(dict(panel_id=p,family_id=rr[0]['family_id'],mode=m,extra_assays=b,comparator=c,
            paired_seeds=20,mean_joint_minus_comparator_new_hits=float(np.mean(diff)),median_paired_new_hit_difference=float(np.median(diff)),
            seed_wins=sum(d>0 for d in diff),seed_ties=sum(d==0 for d in diff),seed_losses=sum(d<0 for d in diff),
            mean_joint_minus_comparator_cumulative_best=float(np.mean(best)),
            seed_percentile10=float(np.quantile(diff,.1)),seed_percentile90=float(np.quantile(diff,.9)),
            variability_is_sampling_seeds_not_independent_biology=True))
    return paired,summary

def success_checks(candidate_summary,site_records,gates):
    perpanel=[]
    for gate in gates:
        if gate['status']!='evaluable':continue
        panel=gate['panel_id'];mode=gate['mode']
        if mode in {'mode2','mode3'}:
            final=gate['extra_budget']
            rows={r['comparator']:r for r in candidate_summary if r['panel_id']==panel and r['mode']==mode and r['extra_assays']==final}
            required=['baseline','random','sequence_diversity']
            eligible=gate['favorable_discovery_claim_eligible']
            consistent=eligible and all(c in rows and rows[c]['mean_joint_minus_comparator_new_hits']>0 for c in required)
            complementary=consistent and all(c in rows and rows[c]['mean_joint_minus_comparator_new_hits']>0 for c in ['baseline_R','baseline_H'])
            perpanel.append(dict(panel_id=panel,family_id=gate['family_id'],mode=mode,primary_budget=final,
                favorable_discovery_assessment_eligible=eligible,joint_better_than_baseline_random_and_sequence_diversity=consistent,
                both_readouts_add_complementary_yield=complementary,
                no_posthoc_budget_or_winning_method_selection=True))
        else:
            rr=[r for r in site_records if r['panel_id']==panel and r['n_positions']==2]
            if not rr:continue
            joint=next(r for r in rr if r['strategy']=='parent_joint_percentile')
            random=[r['new_hits'] for r in rr if r['strategy']=='random_positions']
            static=[r['new_hits'] for r in rr if r['strategy'] in {'static_contact_degree','static_relative_SASA'}]
            assert len(random)==20 and len(static)==2
            perpanel.append(dict(panel_id=panel,family_id=gate['family_id'],mode=mode,primary_budget=joint['assay_budget'],
                primary_joint_high_better_than_both_static_and_mean_random=joint['new_hits']>max(static+[float(np.mean(random))]),
                joint_high_hits=joint['new_hits'],random_mean_hits=float(np.mean(random)),static_hits=static,
                task_object=joint['task_object'],endpoint=joint['endpoint'],
                no_companion_substitution_for_primary=True))
    mode_checks=[]
    for mode in ['mode1','mode2','mode3']:
        rr=[r for r in perpanel if r['mode']==mode]
        key='primary_joint_high_better_than_both_static_and_mean_random' if mode=='mode1' else 'joint_better_than_baseline_random_and_sequence_diversity'
        mode_checks.append(dict(mode=mode,evaluable_panels=len(rr),evaluable_biological_families=len({r['family_id'] for r in rr}),
            predefined_consistency_criterion_met_in_every_panel=bool(rr) and all(r[key] for r in rr),
            no_generalization_claim_from_seed_repetitions=True))
    return dict(per_panel_primary=perpanel,mode_consistency=mode_checks,
        three_modes_are_separate_hypotheses=True,no_validated_design_claim_without_new_experiments=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--contracts',required=True);p.add_argument('--output',required=True)
    args=p.parse_args();output=Path(args.output)
    assert not output.exists() or not any(output.iterdir()),'Use a new output directory; do not overwrite outcomes'
    output.mkdir(parents=True,exist_ok=True)
    protocol_path=HERE/'multimode_protocol.json';binding=read_json(HERE/'protocol_freeze_binding.json')
    assert sha(protocol_path)==binding['current_protocol_sha256']
    data=read_json(args.contracts);panels=data['panels'];assert len({r['panel_id'] for r in panels})==len(panels)
    records=[];sites=[];gates=[];failures=[]
    for panel in panels:
        assert panel['family_id'] and panel['endpoint'] and panel['task_object'] in {'TCR_mutation','antigen_peptide_mutation'}
        assert yes(panel['threshold_fixed_without_model_outputs'])
        if yes(panel.get('mode1_eligible',False)):
            rr,gg=evaluate_site_panel(panel,output);sites.extend(rr);gates.extend(gg)
        if yes(panel.get('mode2_eligible',False)) or yes(panel.get('mode3_eligible',False)):
            rr,gg=evaluate_candidate_panel(panel,output);records.extend(rr);gates.extend(gg)
        if not any(yes(panel.get(m+'_eligible',False)) for m in ['mode1','mode2','mode3']):
            gates.append(dict(panel_id=panel['panel_id'],family_id=panel['family_id'],mode='all',status='not_evaluable',reason=panel.get('ineligibility_reason','Source coverage or endpoint contract not suitable')))
    paired,summary=summarize_candidate_records(records,gates)
    dump_csv(output/'all_candidate_budget_metrics.csv',records);dump_csv(output/'all_site_nomination_metrics.csv',sites)
    dump_csv(output/'paired_seed_differences.csv',paired);dump_csv(output/'paired_seed_summary.csv',summary)
    dump_json(output/'all_data_gates.json',gates)
    dump_json(output/'predeclared_success_checks.json',success_checks(summary,sites,gates))
    dump_json(output/'execution_audit.json',dict(status='completed',timestamp_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(protocol_path),contracts_path=args.contracts,contracts_sha256=sha(args.contracts),
        evaluator_sha256=sha(__file__),checkpoint_sha256=CHECKPOINT,n_panels=len(panels),n_families=len({p['family_id'] for p in panels}),
        seed_ids=SEEDS,statistical_replication='biological family, not20randomseeds',
        no_neural_inference_or_training=True,no_hyperparameter_or_feature_selection=True,no_remote_access=True,
        all_eligible_and_ineligible_panels_preserved=True,three_roles_not_pooled=True,
        no_candidate_std_as_epistemic_uncertainty=True,whole_pool_bandwidths_and_feature_scales_label_free_transductive=True,
        exploration_policy='4highest current predicted mean+1farthest coverage; outcomes revealed only after each5-candidate batch',
        hypotheses_unvalidated_until_new_independent_results=True))
    print(json.dumps(dict(status='completed',panels=len(panels),families=len({p['family_id'] for p in panels}),candidate_metric_rows=len(records),site_metric_rows=len(sites),output=str(output))))

if __name__=='__main__':main()
