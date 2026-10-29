"""Prepare complete, gap-preserving chain alignments from existing NPZ outputs.

No inference, training or trajectory processing occurs in this script.
"""
from pathlib import Path
from functools import lru_cache
import json, hashlib, sys, math, warnings
warnings.filterwarnings('ignore',message="The attribute.*was renamed")
ROOT=Path(__file__).resolve().parent
ATLAS=ROOT.parent
EXP=ATLAS.parent
sys.path.insert(0,str(EXP/'completed_results_20261012/plot_runtime'))
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
_dumps=json.dumps
json.dumps=lambda *args,**kw:_dumps(*args,default=lambda x:x.item() if hasattr(x,'item') else str(x),**kw)

CAT=json.loads((ATLAS/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
ANNOT=Path(r'E:\Project\FlexPPI\output\itrace_residue_biology\raw\matched_residue_context_responses.csv')
USE=['pdb_id','residue_index','chain','chain_role','author_residue_id','amino_acid_1','region','is_tcr_pmhc_interface_8A']
ANN=pd.read_csv(ANNOT,usecols=USE,dtype={'author_residue_id':str},low_memory=False)
MAN=pd.read_csv('E:/data/ITrace/tables/manifest.csv')
ORDER=['mhc','beta2m','peptide','tcr_alpha','tcr_beta']
MEL5={'3hg1','4jff','4jfd','4jfe','4qok','6eqa','7q99','7q9a'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def make_aligner(peptide=False):
    x=PairwiseAligner();x.mode='global';x.match_score=2;x.mismatch_score=-1
    x.open_gap_score=-6;x.extend_gap_score=-.5
    if not peptide:x.target_end_gap_score=0;x.query_end_gap_score=0
    return x

@lru_cache(None)
def load(pdb):
    src=EXP/'tcr_case_v40_20261017/family_predictions/full'/f'itrace_{pdb}.npz'
    if pdb in {'9nmu','9nmx'}:
        state='WT' if pdb=='9nmu' else 'S32M'
        base=EXP/'tcr_expanded_screen_v41_20261017/catchbond_validation'
        src=base/f'catchbond_{state}_experimental_matched_prediction.npz'
        z=np.load(src,allow_pickle=True)
        f=pd.DataFrame({'index':np.arange(len(z['mean'])),'chain':z['chain'],'id':z['author_residue_id'].astype(str),'aa':z['amino_acid']})
        f['role']=f.chain.map({'A':'mhc','C':'peptide','D':'tcr_alpha','E':'tcr_beta'})
        f['region']='Unannotated'
        rid=pd.to_numeric(f.id,errors='coerce')
        for ch,lo,hi,name in [('D',25,35,'CDR1α'),('D',90,97,'CDR3α'),('E',92,102,'CDR3β')]:
            f.loc[f.chain.eq(ch)&rid.between(lo,hi),'region']=name
        f['interface']=False
        lab=pd.read_csv(base/'all408_ten_replica_mean_SD_SEM.csv',dtype={'author_residue_id':str})
        lab=lab[lab.condition.eq(state)].set_index(['chain','author_residue_id'])
        for metric,col,j in [('rmsf','rmsf_A',0),('hchi','Hchi_nats',1)]:
            f[metric+'_full']=z['mean'][:,j]
            f[metric+'_md']=[float(lab.loc[(r.chain,r.id),col+'_mean']) for r in f.itertuples()]
            f[metric+'_var']=[float(lab.loc[(r.chain,r.id),col+'_SD'])**2 for r in f.itertuples()]
            f[metric+'_nrep']=[int(lab.loc[(r.chain,r.id),col+'_n']) for r in f.itertuples()]
            f[metric+'_valid']=np.isfinite(f[metric+'_md'])&(f[metric+'_nrep']>0)
        f['annotation_source']='CatchBond prepared-chain mapping; author-defined CDR1α/CDR3 only, other regions unannotated'
        inputp=base/f'catchbond_{state}_experimental_matched.pdb'
    else:
        z=np.load(src,allow_pickle=True)
        af=ANN[ANN.pdb_id.eq(pdb)].sort_values('residue_index')
        assert len(af)==len(z['mean']),(pdb,len(af),len(z['mean']))
        assert np.array_equal(af.chain.to_numpy(str),z['chain'])
        assert np.array_equal(af.amino_acid_1.to_numpy(str),z['amino_acid'])
        assert np.array_equal(af.residue_index.to_numpy(),np.arange(len(af)))
        f=pd.DataFrame({'index':np.arange(len(af)),'chain':z['chain'],'id':z['author_residue_id'].astype(str),'aa':z['amino_acid'],
                        'role':af.chain_role.to_numpy(),'region':af.region.to_numpy(),'interface':af.is_tcr_pmhc_interface_8A.to_numpy(bool)})
        mdreps=MAN[MAN.pdb_id.eq(pdb)].traj_id.nunique()
        assert mdreps==3,(pdb,mdreps)
        for metric,j in [('rmsf',0),('hchi',1)]:
            f[metric+'_full']=z['mean'][:,j]
            f[metric+'_md']=z['target'][:,j]
            f[metric+'_var']=z['label_variance'][:,j]
            f[metric+'_valid']=z['valid'][:,j]&np.isfinite(z['target'][:,j])
            f[metric+'_nrep']=mdreps
        f['annotation_source']=pdb+' original full-chain annotation'
        inputp=Path(f'E:/data/ITrace/trajectories/{pdb}_run1_topology.pdb')
    # Explicit float64 export preserves the exact binary float32 value from NPZ;
    # pandas otherwise writes float32 columns using a shortened decimal display.
    for metric in ['rmsf','hchi']:
        for value in ['full','md','var']:
            f[metric+'_'+value]=f[metric+'_'+value].astype('float64')
    assert f.role.notna().all(),(pdb,'unmapped chain role')
    f['region_legacy']=f.region
    if pdb in MEL5 and pdb!='3hg1':
        w,_=load('3hg1')
        for role in ['tcr_alpha','tcr_beta']:
            wa=w[w.role.eq(role)].reset_index(drop=True);fa=f[f.role.eq(role)].reset_index(drop=True)
            hit=make_aligner().align(''.join(wa.aa),''.join(fa.aa))[0]
            for i,j in zip(*hit.indices):
                if i>=0 and j>=0:
                    ix=int(fa.iloc[j]['index']);f.loc[ix,'region']=wa.iloc[i].region
                    f.loc[ix,'annotation_source']='3HG1 WT CDR boundaries globally projected'
    return f,dict(pdb=pdb,prediction=str(src),prediction_sha256=sha(src),structure=str(inputp),structure_sha256=sha(inputp),
                  chain_roles={str(ch):str(g.role.iloc[0]) for ch,g in f.groupby('chain',sort=False)},
                  label_mean_source='NPZ target first two columns' if pdb not in {'9nmu','9nmx'} else 'ten-replica mean CSV',
                  label_variance_source='NPZ sample variance across 3 independent MD runs' if pdb not in {'9nmu','9nmx'} else 'SD squared across 10 independent MD runs')

def rho(x,y):
    return float(pd.Series(x).rank().corr(pd.Series(y).rank())) if len(x)>2 and np.ptp(x)>0 and np.ptp(y)>0 else None

def write_pair(case,a,b,out):
    fa,sa=load(a);fb,sb=load(b)
    rs=[];audit=[]
    roles=sorted(set(fa.role)|set(fb.role),key=lambda r:ORDER.index(r) if r in ORDER else 99)
    for role in roles:
        aa=fa[fa.role.eq(role)].reset_index(drop=True);bb=fb[fb.role.eq(role)].reset_index(drop=True)
        if aa.empty or bb.empty:raise ValueError((a,b,role,'missing role'))
        sx=''.join(aa.aa);sy=''.join(bb.aa)
        al=make_aligner(role=='peptide').align(sx,sy)[0]
        pairs=list(zip(*al.indices))
        assert sorted(i for i,j in pairs if i>=0)==list(range(len(aa)))
        assert sorted(j for i,j in pairs if j>=0)==list(range(len(bb)))
        jointly=[(i,j) for i,j in pairs if i>=0 and j>=0]
        var=[(i,j) for i,j in jointly if i<125 and j<125]
        ident=sum(sx[i]==sy[j] for i,j in jointly)/len(jointly)
        vid=sum(sx[i]==sy[j] for i,j in var)/len(var) if var else 0
        near=not role.startswith('tcr') or (vid>=.9 and len(var)/min(125,len(aa),len(bb))>=.85)
        audit.append(dict(role=role,length_a=len(aa),length_b=len(bb),aligned_columns=len(pairs),paired_positions=len(jointly),
                          identity=ident,variable_identity_first125=vid,near_related_for_CDR_comparison=near,
                          substitution_n=sum(sx[i]!=sy[j] for i,j in jointly),gap_a=sum(i<0 for i,j in pairs),gap_b=sum(j<0 for i,j in pairs),
                          score=float(al.score),sequence_a=sx,sequence_b=sy,alignment=str(al)))
        for col,(i,j) in enumerate(pairs):
            ar=aa.iloc[i] if i>=0 else None;br=bb.iloc[j] if j>=0 else None
            r=dict(case_id=case['case_id'],a=a,b=b,role=role,alignment_column=col+1,gap_a=i<0,gap_b=j<0,
                   paired=i>=0 and j>=0,same_aa=bool(i>=0 and j>=0 and ar.aa==br.aa),
                   region_a=str(ar.region) if ar is not None else '',region_b_legacy=str(br.region_legacy) if br is not None else '',
                   region_b_projected=str(ar.region) if ar is not None and br is not None else '',
                   interface_union=bool((ar is not None and ar.interface) or (br is not None and br.interface)),
                   near_related_for_CDR_comparison=near,CDR_annotation_source=str(ar.annotation_source) if ar is not None else '',
                   label='')
            for suffix,x,pos in [('a',ar,i),('b',br,j)]:
                for key in ['index','chain','id','aa']:
                    r[key+'_'+suffix]=x[key] if x is not None else ''
                r['position_'+suffix]=pos+1 if pos>=0 else ''
                for metric in ['rmsf','hchi']:
                    for v in ['full','md','var','nrep','valid']:
                        r[f'{metric}_{v}_{suffix}']=x[f'{metric}_{v}'] if x is not None else (False if v=='valid' else np.nan)
            r['label']=f"{r['aa_a']}{r['id_a']}" if r['same_aa'] else f"{r['aa_a'] or '-'}{r['id_a']}→{r['aa_b'] or '-'}{r['id_b']}"
            for metric in ['rmsf','hchi']:
                r[metric+'_valid']=bool(r[metric+'_valid_a'] and r[metric+'_valid_b'])
                r[metric+'_full_delta']=r[metric+'_full_b']-r[metric+'_full_a']
                r[metric+'_md_delta']=r[metric+'_md_b']-r[metric+'_md_a']
                r[metric+'_delta_sem']=math.sqrt(r[metric+'_var_a']/r[metric+'_nrep_a']+r[metric+'_var_b']/r[metric+'_nrep_b']) if r[metric+'_valid'] else np.nan
            rs.append(r)
    f=pd.DataFrame(rs)
    stem=f'{a}_{b}'
    f.to_csv(out/(stem+'_all_chain_alignment.csv'),index=False,encoding='utf-8-sig')
    # Convenience aliases follow the established manuscript source schema.
    f['region_b']=f.region_b_projected
    f[f.role.eq('peptide')].to_csv(out/(stem+'_complete_peptide.csv'),index=False,encoding='utf-8-sig')
    f[f.role.isin(['tcr_alpha','tcr_beta'])&f.region_a.str.startswith('CDR')].to_csv(out/(stem+'_all_CDR_projected.csv'),index=False,encoding='utf-8-sig')
    f[f.role.isin(['mhc','peptide'])&f.interface_union].to_csv(out/(stem+'_complete_pMHC_interface.csv'),index=False,encoding='utf-8-sig')
    metrics=[]
    for scope,sel in [('whole_peptide',f.role.eq('peptide')),('natural_pMHC_interface',f.role.isin(['mhc','peptide'])&f.interface_union),('all_CDR',f.role.isin(['tcr_alpha','tcr_beta'])&f.region_a.str.startswith('CDR'))]:
        for m in ['rmsf','hchi']:
            selected=sel&f.paired&f[m+'_valid']&((f.same_aa) if m=='hchi' else True)
            q=f[selected];x=q[m+'_md_delta'].to_numpy(float);y=q[m+'_full_delta'].to_numpy(float)
            zero=float(np.mean(abs(x))) if len(x) else None;mae=float(np.mean(abs(x-y))) if len(x) else None
            metrics.append(dict(scope=scope,metric=m,n=len(q),rho=rho(x,y),MAE=mae,zero_change_MAE=zero,relative_MAE_gain=1-mae/zero if zero else None,
                                near_related_for_CDR_comparison=bool(q.near_related_for_CDR_comparison.all()) if len(q) else None,
                                caution='Cross-unrelated CDR mapping is exploratory' if scope=='all_CDR' and len(q) and not q.near_related_for_CDR_comparison.all() else 'MD comparison is corroboration, not experimental truth'))
    aj=dict(case_id=case['case_id'],a=a,b=b,alignment_method='Complete-chain global PairwiseAligner, affine gap; all gaps and residues retained',
            scores=dict(match=2,mismatch=-1,gap_open=-6,gap_extend=-.5,protein_terminal_gap=0,peptide_terminal_gap='penalized'),
            source_a=sa,source_b=sb,chain_alignments=audit,total_rows=len(f),paired_rows=int(f.paired.sum()),
            all_residues_accounted_once=True,CDR_rule='a-state WT boundaries projected through complete-chain alignment; canonical 3HG1 for MEL5 family',
            metrics=metrics)
    (out/(stem+'_alignment_audit.json')).write_text(json.dumps(aj,ensure_ascii=False,indent=2),encoding='utf8')
    return aj

def main():
    manifest=[]
    for case in CAT['cases']:
        if int(case['case_id'][1:])<5:continue
        out=ROOT/case['case_id']/'data';out.mkdir(parents=True,exist_ok=True)
        for p in case['pdbs']:
            f,s=load(p);f.to_csv(out/(p+'_absolute_state_values.csv'),index=False,encoding='utf-8-sig')
            (out/(p+'_source_audit.json')).write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf8')
        pairs=[write_pair(case,a,b,out) for a,b in case['contrasts']]
        item=dict(case_id=case['case_id'],name=case['name'],tier=case['tier'],case_directory=str(out.parent),data_directory=str(out),
                  pair_stems=[f'{a}_{b}' for a,b in case['contrasts']],all_states=case['pdbs'],pairs=pairs)
        (out.parent/'data_audit.json').write_text(json.dumps(item,ensure_ascii=False,indent=2),encoding='utf8')
        manifest.append(item)
        print(case['case_id'],len(case['pdbs']),'states',len(pairs),'complete pairs',flush=True)
    (ROOT/'data_manifest.json').write_text(json.dumps(dict(status='passed',cases=manifest,annotation_source=str(ANNOT),annotation_sha256=sha(ANNOT),new_MD=0,new_training=0,new_inference=0),ensure_ascii=False,indent=2),encoding='utf8')

if __name__=='__main__':main()
