"""Find construct-tolerant biological comparisons; score held-out predictions."""
from pathlib import Path
from functools import lru_cache
from collections import Counter
import sys,json,itertools,hashlib
ROOT=Path(__file__).resolve().parent;EXP=ROOT.parent;OLD=EXP/'tcr_case_v40_20261017'
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
OUT=ROOT/'construct_tolerant_screen';OUT.mkdir(exist_ok=True)
inventory=json.loads((OLD/'local_pair_screen/sequence_inventory.json').read_text(encoding='utf8'))
al=PairwiseAligner();al.mode='global';al.match_score=2;al.mismatch_score=-1;al.open_gap_score=-6;al.extend_gap_score=-.5
al.target_end_gap_score=0;al.query_end_gap_score=0
@lru_cache(None)
def align(a,b):
    if len(a)==len(b):return tuple((i,i) for i in range(len(a)))
    x=al.align(a,b)[0]
    return tuple((int(i),int(j)) for i,j in zip(*x.indices) if i>=0 and j>=0)
def comparison(a,b,limit=None):
    pairs=align(a,b)
    if limit:pairs=[(i,j) for i,j in pairs if i<limit and j<limit]
    matches=sum(a[i]==b[j] for i,j in pairs);n=len(pairs)
    minimum=min(len(a),len(b),limit or 100000)
    return {'identity':matches/n if n else 0,'edits':n-matches,'coverage':n/minimum if minimum else 0,'n':n}
pairs=[]
for a,b in itertools.combinations(sorted(inventory),2):
    ar,br=inventory[a]['roles'],inventory[b]['roles']
    if not all(role in ar and role in br for role in ['mhc','peptide','tcr_alpha','tcr_beta']):continue
    seq=lambda role:(ar[role]['sequence'],br[role]['sequence'])
    m=comparison(*seq('mhc'),180);ta=comparison(*seq('tcr_alpha'),125);tb=comparison(*seq('tcr_beta'),125)
    pa,pb=seq('peptide');p=comparison(pa,pb);samepep=pa==pb
    mhc_equal=m['identity']>=.995 and m['coverage']>=.95
    tcr_equal=all(c['identity']>=.99 and c['edits']<=1 and c['coverage']>=.8 for c in [ta,tb])
    category=None
    if mhc_equal and tcr_equal and not samepep and p['coverage']>=.7 and p['identity']>=.4:category='antigen_variant_same_TCR_groove'
    elif mhc_equal and samepep and ta['coverage']>=.8 and tb['coverage']>=.8:
        if ta['edits']+tb['edits']==0:category='same_variable_domains_construct_control'
        elif min(ta['identity'],tb['identity'])>=.9 and ta['edits']+tb['edits']<=12:category='TCR_variable_variant_same_pMHC'
        else:category='different_TCR_same_peptide_groove'
    elif tcr_equal and samepep and m['coverage']>=.95 and 0<m['edits']<=8:category='MHC_groove_variant_same_TCR_peptide'
    if not category:continue
    rows={'a':a,'b':b,'category':category,'peptide_a':pa,'peptide_b':pb,'peptide_aligned_identity':p['identity'],'peptide_edits':p['edits'],'mhc_groove_identity':m['identity'],'mhc_groove_edits':m['edits'],'tcr_alpha_variable_identity':ta['identity'],'tcr_beta_variable_identity':tb['identity'],'tcr_alpha_variable_edits':ta['edits'],'tcr_beta_variable_edits':tb['edits'],'mhc_full_equal':ar['mhc']['sequence']==br['mhc']['sequence'],'tcr_alpha_full_equal':ar['tcr_alpha']['sequence']==br['tcr_alpha']['sequence'],'tcr_beta_full_equal':ar['tcr_beta']['sequence']==br['tcr_beta']['sequence']}
    for side,pdb in [('a',a),('b',b)]:
        meta=inventory[pdb]['metadata'];rows['title_'+side]=meta.get('title','');rows['doi_'+side]=meta.get('citation',{}).get('doi','')
    pairs.append(rows)
pf=pd.DataFrame(pairs);pf.to_csv(OUT/'expanded_candidate_pairs.csv',index=False)
oldpairs=pd.read_csv(OLD/'local_pair_screen/biological_candidate_pairs.csv')
oldset={tuple(sorted((r.a,r.b))) for r in oldpairs.itertuples()}
pf['new_vs_initial_screen']=[(r.a,r.b) not in oldset for r in pf.itertuples()]
pf.to_csv(OUT/'expanded_candidate_pairs.csv',index=False)
print('Pairs',len(pf),Counter(pf.category),'new',pf.new_vs_initial_screen.sum(),flush=True)
@lru_cache(None)
def load(pdb):
    z=np.load(OLD/f'family_predictions/full/itrace_{pdb}.npz')
    f=pd.DataFrame({'chain':z['chain'],'aa':z['amino_acid'],'author_id':z['author_residue_id'],'index':np.arange(len(z['mean']))})
    rolemap={v['chain']:k for k,v in inventory[pdb]['roles'].items()};f['role']=f.chain.map(rolemap)
    for metric,j in [('rmsf',0),('hchi',1)]:
        f[metric+'_md']=z['target'][:,j];f[metric+'_pred']=z['mean'][:,j];f[metric+'_var']=z['label_variance'][:,j];f[metric+'_valid']=z['valid'][:,j]
    return f
def rho(x,y):
    return float(np.corrcoef(pd.Series(x).rank(),pd.Series(y).rank())[0,1]) if len(x)>=4 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
summary=[];details=[]
for pi,pair in enumerate(pf.itertuples()):
    ar,br=load(pair.a),load(pair.b);rows=[]
    for role in ['mhc','b2m','peptide','tcr_alpha','tcr_beta']:
        a=ar[ar.role.eq(role)].reset_index(drop=True);b=br[br.role.eq(role)].reset_index(drop=True)
        if a.empty or b.empty:continue
        sa,sb=''.join(a.aa),''.join(b.aa)
        for i,j in align(sa,sb):
            x,y=a.iloc[i],b.iloc[j]
            row=dict(a=pair.a,b=pair.b,category=pair.category,role=role,position_a=i+1,position_b=j+1,id_a=x.author_id,id_b=y.author_id,aa_a=x.aa,aa_b=y.aa,same_aa=x.aa==y.aa)
            for metric in ['rmsf','hchi']:
                valid=x[metric+'_valid'] and y[metric+'_valid'];row[metric+'_valid']=bool(valid)
                for side,r in [('a',x),('b',y)]:
                    for m in ['md','pred','var']:row[f'{metric}_{m}_{side}']=float(r[f'{metric}_{m}']) if valid else np.nan
                row[metric+'_md_delta']=row[metric+'_md_b']-row[metric+'_md_a'];row[metric+'_pred_delta']=row[metric+'_pred_b']-row[metric+'_pred_a']
                row[metric+'_sem']=np.sqrt((row[metric+'_var_a']+row[metric+'_var_b'])/3)
            rows.append(row)
    frame=pd.DataFrame(rows);details.append(frame)
    scopes={'peptide_unchanged':frame.role.eq('peptide')&frame.same_aa,'peptide_all':frame.role.eq('peptide'),'TCR_variable_unchanged':frame.role.isin(['tcr_alpha','tcr_beta'])&frame.same_aa&(frame.position_a<=125)&(frame.position_b<=125),'MHC_groove_unchanged':frame.role.eq('mhc')&frame.same_aa&(frame.position_a<=180)&(frame.position_b<=180)}
    for scope,select in scopes.items():
        for metric in ['rmsf','hchi']:
            f=frame[select&frame[metric+'_valid']];x=f[metric+'_md_delta'].to_numpy();y=f[metric+'_pred_delta'].to_numpy();sem=f[metric+'_sem'].to_numpy()
            if not len(f):continue
            zero=np.mean(abs(x));mae=np.mean(abs(x-y));robust=abs(x)>2*sem
            summary.append(dict(a=pair.a,b=pair.b,category=pair.category,new_vs_initial_screen=pair.new_vs_initial_screen,scope=scope,metric=metric,n=len(x),rho_delta=rho(x,y),MAE_delta=mae,zero_change_MAE=zero,relative_MAE_gain=1-mae/zero if zero>1e-8 else np.nan,direction_fraction=float(np.mean(x*y>0)),robust_n=int(robust.sum()),robust_direction_fraction=float(np.mean(x[robust]*y[robust]>0)) if robust.any() else np.nan,MD_abs_delta=zero,pred_abs_delta=float(np.mean(abs(y))),absolute_MAE_a=float(np.mean(abs(f[metric+'_md_a']-f[metric+'_pred_a']))),absolute_MAE_b=float(np.mean(abs(f[metric+'_md_b']-f[metric+'_pred_b'])))))
    if (pi+1)%100==0:print('scored',pi+1,flush=True)
pd.DataFrame(summary).to_csv(OUT/'expanded_pair_metrics.csv',index=False)
pd.concat(details,ignore_index=True).to_csv(OUT/'expanded_pair_all_positions.csv.gz',index=False,compression='gzip')
audit={'status':'complete','n_structures':len(inventory),'n_pairs':len(pf),'new_pairs_vs_original':int(pf.new_vs_initial_screen.sum()),'categories':dict(Counter(pf.category)),'source_predictions':'Existing Full joint2output model; family fold0; all244 ITrace held out','sequence_rule':'MHC first180 positions and TCR first125 positions; exact same-length positional alignment, affine-gap global alignment for unequal lengths, terminal gaps unpenalized; both variable-domain coverages>=80%; complete identity fields retained for manual biological review','retained_substitution_positions':'Both peptide_all and unchanged-AA results; chi-entropy comparisons of changed chemistry must not be interpreted as local-environment-only evidence','new_inference_or_MD':0}
(OUT/'expanded_screen_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
print(json.dumps(audit,indent=2),flush=True)
