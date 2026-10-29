from pathlib import Path
from functools import lru_cache
from difflib import SequenceMatcher
from collections import Counter
import json,sys,hashlib,itertools
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
OUT=ROOT/('change_screen_focus' if '--focus' in sys.argv else 'change_screen');OUT.mkdir(exist_ok=True)
inventory=json.loads((ROOT/'local_pair_screen/sequence_inventory.json').read_text(encoding='utf8'))
roles=json.loads(Path('E:/Project/FlexPPI/output/itrace_benchmark/raw/itrace_chain_roles.json').read_text())

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

meta_cols=['pdb_id','condition','chain','chain_role','author_residue_id','amino_acid_1','region','is_tcr_pmhc_interface_8A','min_distance_to_peptide_angstrom','min_distance_to_tcr_angstrom']
annotation_path=Path('E:/Project/FlexPPI/output/itrace_residue_biology/raw/matched_residue_context_responses.csv')
annotations=[]
for chunk in pd.read_csv(annotation_path,usecols=meta_cols,dtype={'pdb_id':str,'chain':str,'author_residue_id':str},chunksize=100000):
    annotations.append(chunk[chunk.condition=='full'].drop(columns='condition'))
annotations=pd.concat(annotations,ignore_index=True)
assert not annotations.duplicated(['pdb_id','chain','author_residue_id']).any()

@lru_cache(None)
def load(pdb):
    with np.load(ROOT/f'family_predictions/full/itrace_{pdb}.npz',allow_pickle=False) as z:
        a={k:z[k].copy() for k in z.files}
    with np.load(ROOT/f'family_predictions/intrachain/itrace_{pdb}.npz',allow_pickle=False) as z:
        b={k:z[k].copy() for k in z.files}
    assert np.array_equal(a['chain'],b['chain']) and np.array_equal(a['author_residue_id'],b['author_residue_id'])
    assert np.array_equal(a['target'],b['target'],equal_nan=True)
    rows=pd.DataFrame({'pdb_id':pdb,'index':np.arange(len(a['mean'])),'chain':a['chain'],'author_residue_id':a['author_residue_id'],'aa':a['amino_acid']})
    rows['role']=rows.chain.map(roles[pdb]);assert rows.role.notna().all()
    for k,j in [('rmsf',0),('hchi',1)]:
        rows[k+'_md']=a['target'][:,j];rows[k+'_sd']=np.sqrt(a['label_variance'][:,j])
        rows[k+'_full']=a['mean'][:,j];rows[k+'_intra']=b['mean'][:,j];rows[k+'_valid']=a['valid'][:,j]
    rows=rows.merge(annotations,on=['pdb_id','chain','author_residue_id'],how='left',validate='one_to_one')
    assert (rows.aa==rows.amino_acid_1).all()
    rows['region']=rows.region.fillna('Unannotated')
    rows['is_tcr_pmhc_interface_8A']=rows.is_tcr_pmhc_interface_8A.fillna(False)
    return rows

@lru_cache(None)
def align(sa,sb):
    if len(sa)==len(sb):return tuple(zip(range(len(sa)),range(len(sb))))
    result=[]
    for tag,a0,a1,b0,b1 in SequenceMatcher(None,sa,sb,autojunk=False).get_opcodes():
        if tag=='equal' or tag=='replace' and a1-a0==b1-b0:result+=list(zip(range(a0,a1),range(b0,b1)))
    return tuple(result)

def rho(x,y):
    if len(x)<4 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
    return np.corrcoef(pd.Series(x).rank(method='average'),pd.Series(y).rank(method='average'))[0,1]

pairs=pd.read_csv(ROOT/'local_pair_screen/biological_candidate_pairs.csv')
# These independent literature-defined contrasts also qualify as useful checks when construct lengths differ.
extra=[('1ao7','1qrn'),('1ao7','1qse'),('1ao7','1qsf'),('1ao7','3h9s'),('1ao7','3pwp'),('2bnq','2bnr')]
known={tuple(sorted((r.a,r.b))) for r in pairs.itertuples()}
for a,b in extra:
    if tuple(sorted((a,b))) not in known:
        pairs=pd.concat([pairs,pd.DataFrame([{'a':a,'b':b,'category':'literature_defined_antigen_contrast'}])],ignore_index=True)
knownseries=[('7ow5','7ow6'),('8rlt','8rlu'),('8rlt','8rlv'),('8rlu','8rlv')]+extra
if '--focus' in sys.argv:
    knownset={tuple(sorted(p)) for p in knownseries}
    pairs=pairs[pairs.apply(lambda row:tuple(sorted((row.a,row.b))) in knownset,axis=1)]
details=[];summary=[]
for pi,p in enumerate(pairs.itertuples()):
    a,b=load(p.a),load(p.b)
    parts=[]
    for role in ['mhc','b2m','peptide','tcr_alpha','tcr_beta']:
        aa=a[a.role==role].reset_index(drop=True);bb=b[b.role==role].reset_index(drop=True)
        if not len(aa) or not len(bb):continue
        sa=''.join(aa.aa);sb=''.join(bb.aa)
        for i,j in align(sa,sb):
            ra,rb=aa.iloc[i],bb.iloc[j]
            same_aa=ra.aa==rb.aa
            row={'a':p.a,'b':p.b,'category':p.category,'role':role,'chain_a':ra.chain,'chain_b':rb.chain,'id_a':ra.author_residue_id,'id_b':rb.author_residue_id,'position_a':i+1,'position_b':j+1,'aa_a':ra.aa,'aa_b':rb.aa,'same_aa':same_aa,'region_a':ra.region,'region_b':rb.region,'interface':bool(ra.is_tcr_pmhc_interface_8A or rb.is_tcr_pmhc_interface_8A),'distance_peptide_a':ra.min_distance_to_peptide_angstrom,'distance_peptide_b':rb.min_distance_to_peptide_angstrom}
            for metric in ['rmsf','hchi']:
                valid=ra[metric+'_valid'] and rb[metric+'_valid']
                row[metric+'_valid']=valid
                for cond,r in [('a',ra),('b',rb)]:
                    for method in ['md','sd','full','intra']:row[f'{metric}_{method}_{cond}']=float(r[metric+'_'+method]) if valid else np.nan
                for method in ['md','full','intra']:row[f'{metric}_{method}_delta']=row[f'{metric}_{method}_b']-row[f'{metric}_{method}_a']
                row[f'{metric}_delta_sem']=np.sqrt((row[f'{metric}_sd_a']**2+row[f'{metric}_sd_b']**2)/3) if valid else np.nan
            parts.append(row)
    frame=pd.DataFrame(parts);details.append(frame)
    scopes={'all_unchanged':frame.same_aa,'peptide_unchanged':frame.same_aa&(frame.role=='peptide'),'tcr_interface_unchanged':frame.same_aa&frame.interface&frame.role.isin(['tcr_alpha','tcr_beta']),'cdr3_unchanged':frame.same_aa&frame.region_a.str.contains('CDR3')&frame.region_b.str.contains('CDR3'),'mhc_groove_unchanged':frame.same_aa&(frame.role=='mhc')&(frame.position_a<=180)&(frame.position_b<=180),'all_interface_unchanged':frame.same_aa&frame.interface}
    for scope,mask in scopes.items():
        for metric in ['rmsf','hchi']:
            f=frame[mask&frame[metric+'_valid']]
            x=f[metric+'_md_delta'].to_numpy();sem=f[metric+'_delta_sem'].to_numpy()
            if not len(f):continue
            for method in ['full','intra']:
                y=f[f'{metric}_{method}_delta'].to_numpy();robust=np.abs(x)>2*sem
                summary.append({'a':p.a,'b':p.b,'category':p.category,'scope':scope,'metric':metric,'method':method,'n':len(f),'rho_delta':rho(x,y),'mae_delta':np.mean(np.abs(x-y)),'md_mean_delta':np.mean(x),'model_mean_delta':np.mean(y),'direction_fraction':np.mean(x*y>0),'md_robust_n':int(robust.sum()),'direction_fraction_robust':float(np.mean(x[robust]*y[robust]>0)) if robust.any() else np.nan,'mean_abs_md_delta':np.mean(abs(x)),'mean_abs_model_delta':np.mean(abs(y))})
    if (p.a,p.b) in knownseries or (p.b,p.a) in knownseries:
        frame.to_csv(OUT/f'pair_{p.a}_{p.b}_all_positions.csv',index=False,encoding='utf-8-sig')
    if (pi+1)%50==0:print('pairs',pi+1,flush=True)
pd.DataFrame(summary).to_csv(OUT/'all_candidate_change_metrics.csv',index=False,encoding='utf-8-sig')
pd.concat(details,ignore_index=True).to_csv(OUT/'all_candidate_changes.csv.gz',index=False,compression='gzip')
primary=pd.DataFrame(summary)
print(primary[primary.apply(lambda r:(r.a,r.b) in knownseries or (r.b,r.a) in knownseries,axis=1)&primary.scope.isin(['peptide_unchanged','cdr3_unchanged'])].to_string(index=False),flush=True)
(OUT/'screen_audit.json').write_text(json.dumps({'status':'complete','pairs':len(pairs),'candidate_counts':dict(Counter(pairs.category)),'independent_biology_candidates':knownseries,'model':'Full and Intra family-disjoint fold0 seed7 epoch40; all244 ITrace excluded from training','source_annotation_scope':'Only chain role, exact identity, CDR and contact annotations; no legacy model values used','annotation_source_sha256':sha(annotation_path),'all_alignment_positions_saved':True,'selection_note':'Descriptive screening across all sequence-defined candidates. No inference or training in this screen. Differences on same-amino-acid positions are primary; substituted sites are retained in detailed tables but excluded from same-type summaries.'},ensure_ascii=False,indent=2),encoding='utf8')
