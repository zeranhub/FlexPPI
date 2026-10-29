"""Fixed-series ITrace screen. No training, no manuscript modifications.

Only the existing two-output fold-0 Full model is primary. Intrachain is a
predeclared control. All region definitions precede scoring and all outcomes
are retained. The cache's label variance is the variance across three MD runs.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import sys, json, hashlib

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
EXP = ROOT.parent
sys.path.insert(0, str(EXP / 'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

OLD = EXP / 'tcr_case_v40_20261017'
META = pd.read_csv(ROOT / 'all244_ITrace_metadata_and_sequences.csv').set_index('pdb')
ANNOTATION = Path(r'E:\Project\FlexPPI\output\itrace_residue_biology\raw\matched_residue_context_responses.csv')
anno = pd.read_csv(ANNOTATION, usecols=['pdb_id', 'residue_index', 'chain', 'chain_role', 'author_residue_id', 'amino_acid_1', 'region', 'is_tcr_pmhc_interface_8A'], dtype={'author_residue_id': str}, low_memory=False)
PLAN = json.loads((ROOT / 'analysis_plan.json').read_text(encoding='utf8'))
SERIES = dict(PLAN['series'])
SERIES.update(PLAN.get('additional_biology_verified_comparisons', {}))
SERIES['single_MHC_change_and_binding_thermodynamics'] = ['1lp9', '2j8u', '2jcc', '2uwe']
# Biology determines the reference state, not ordering by PDB identifier.
SERIES['MAGE_A10_engineering'] = ['7pdw', '7pbc', '7qpj']
ALIGNER = PairwiseAligner()
ALIGNER.mode = 'global'
ALIGNER.match_score = 2
ALIGNER.mismatch_score = -1
ALIGNER.open_gap_score = -6
ALIGNER.extend_gap_score = -.5
ALIGNER.target_end_gap_score = 0
ALIGNER.query_end_gap_score = 0
PEPTIDE_ALIGNER = PairwiseAligner()
PEPTIDE_ALIGNER.mode = 'global'
PEPTIDE_ALIGNER.match_score = 2
PEPTIDE_ALIGNER.mismatch_score = -1
PEPTIDE_ALIGNER.open_gap_score = -6
PEPTIDE_ALIGNER.extend_gap_score = -.5
ROLES = ['mhc', 'peptide', 'tcr_alpha', 'tcr_beta']
CDRS = ['CDR1α', 'CDR2α', 'CDR3α', 'CDR1β', 'CDR2β', 'CDR3β']
MODELS = ['full', 'old_strict', 'intrachain']
FRESH = {'full': ROOT/'family_full_no_ITrace', 'old_strict':ROOT/'old_strict_full_no_ITrace'}
MEL5_CORRESPONDENCE = {'3hg1','4jff','4jfd','4jfe','4qok','6eqa','7q99','7q9a'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

@lru_cache(None)
def alignment(sa, sb, peptide=False):
    # Always align complete chains, including equal-length chains.
    hit = (PEPTIDE_ALIGNER if peptide else ALIGNER).align(sa, sb)[0]
    return tuple((int(i), int(j)) for i, j in zip(*hit.indices) if i >= 0 and j >= 0)

@lru_cache(None)
def load(pdb):
    z = np.load(OLD / 'family_predictions/full' / f'itrace_{pdb}.npz')
    iz = np.load(OLD / 'family_predictions/intrachain' / f'itrace_{pdb}.npz')
    for key in ['chain', 'author_residue_id', 'amino_acid', 'valid']:
        assert np.array_equal(z[key], iz[key]), (pdb, key)
    for key in ['target', 'label_variance']:
        assert np.allclose(z[key], iz[key], atol=1e-8, equal_nan=True), (pdb, key)
    af = anno[anno.pdb_id.eq(pdb)].set_index('residue_index').sort_index()
    assert len(af) == len(z['mean']), (pdb, 'annotation length', len(af), len(z['mean']))
    assert np.array_equal(af.index, np.arange(len(af))), (pdb, 'annotation indices')
    assert np.array_equal(af.chain.to_numpy(dtype=str), z['chain']), (pdb, 'annotation chain')
    assert np.array_equal(af.amino_acid_1.to_numpy(dtype=str), z['amino_acid']), (pdb, 'annotation AA')
    f = pd.DataFrame({'chain': z['chain'], 'id': z['author_residue_id'], 'aa': z['amino_acid'], 'index': np.arange(len(z['mean'])), 'role': af.chain_role.to_numpy(), 'region': af.region.to_numpy(), 'interface': af.is_tcr_pmhc_interface_8A.to_numpy(dtype=bool)})
    fresh_means = {}
    for model, source in FRESH.items():
        fresh = np.load(source/f'itrace_{pdb}.npz')
        lookup = {(str(ch),str(rid),str(aa)):i for i,(ch,rid,aa) in enumerate(zip(fresh['chain'],fresh['author_residue_id'],fresh['amino_acid']))}
        assert len(lookup)==len(fresh['mean'])==len(f), (pdb, model, 'unique residue mapping')
        order = [lookup[(str(row.chain),str(row.id),str(row.aa))] for row in f.itertuples()]
        fresh_means[model] = fresh['mean'][order]
    assert np.max(abs(fresh_means['full']-z['mean']))<5e-5, (pdb,'fresh primary vs cache')
    for m, j in [('rmsf', 0), ('hchi', 1)]:
        f[m+'_md'] = z['target'][:, j]
        f[m+'_var'] = z['label_variance'][:, j]
        f[m+'_valid'] = z['valid'][:, j].astype(bool)
        for model,means in fresh_means.items():
            f[m+'_'+model] = means[:,j]
        f[m+'_intrachain'] = iz['mean'][:, j]
    f['legacy_region'] = f.region
    f['CDR_annotation_source'] = pdb+' legacy annotation'
    # The original 4JFF beta-chain CDR labels were off by one owing to an
    # N-terminal construct residue. Use the same biological WT boundaries for
    # every MEL5 state, including the engineered receptor and its peptide variants.
    if pdb in MEL5_CORRESPONDENCE and pdb!='3hg1':
        wt = load('3hg1')
        for role in ['tcr_alpha','tcr_beta']:
            w = wt[wt.role.eq(role)].reset_index(drop=True)
            r = f[f.role.eq(role)].reset_index(drop=True)
            for i,j in alignment(''.join(w.aa),''.join(r.aa)):
                idx = int(r.iloc[j]['index'])
                f.loc[idx,'region'] = w.iloc[i].region
                f.loc[idx,'CDR_annotation_source'] = '3hg1 WT boundaries globally aligned to '+pdb
    return f

def rho(x, y):
    return float(np.corrcoef(pd.Series(x).rank(), pd.Series(y).rank())[0, 1]) if len(x) >= 4 and np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan

summary, positions, alignment_audit, bounds = [], [], [], []
for series, pdbs in SERIES.items():
    for a, b in combinations(pdbs, 2):
        fa, fb = load(a), load(b)
        rows, tcr_stats = [], []
        for role in ROLES:
            x = fa[fa.role.eq(role)].reset_index(drop=True)
            y = fb[fb.role.eq(role)].reset_index(drop=True)
            if x.empty or y.empty:
                continue
            sa, sb = ''.join(x.aa), ''.join(y.aa)
            pairs = alignment(sa, sb, role == 'peptide')
            variable = [(i, j) for i, j in pairs if i < 125 and j < 125]
            identity = sum(sa[i] == sb[j] for i, j in variable) / len(variable) if variable else 0
            coverage = len(variable) / min(125, len(sa), len(sb))
            if role.startswith('tcr'):
                tcr_stats.append((identity, coverage))
            alignment_audit.append(dict(series=series, a=a, b=b, role=role, chain_a=x.chain.iloc[0], chain_b=y.chain.iloc[0], length_a=len(x), length_b=len(y), n_aligned=len(pairs), n_aa_changes=sum(sa[i] != sb[j] for i, j in pairs), variable_identity=identity, variable_coverage=coverage, alignment_method='full-chain global affine gap; peptide end gaps penalized, protein construct end gaps unpenalized', sequence_a=sa, sequence_b=sb))
            if role.startswith('tcr'):
                for region in CDRS:
                    rr = [(i, j) for i, j in pairs if x.iloc[i].region == region]
                    if not rr:
                        continue
                    bounds.append(dict(series=series, a=a, b=b, role=role, CDR=region, a_positions=';'.join(str(i+1) for i,j in rr), b_positions=';'.join(str(j+1) for i,j in rr), a_author_ids=';'.join(str(x.iloc[i].id) for i,j in rr), b_author_ids=';'.join(str(y.iloc[j].id) for i,j in rr), sequence_a=''.join(sa[i] for i,j in rr), sequence_b=''.join(sb[j] for i,j in rr), aligned_n=len(rr), a_region_total=int(x.region.eq(region).sum()), b_legacy_region_total=int(y.legacy_region.eq(region).sum()), annotation_source=x.iloc[rr[0][0]].CDR_annotation_source, rule='a-state annotation projected through complete-chain alignment; MEL5 series uses canonical 3HG1 WT boundaries; b-state legacy CDR labels not used for region selection'))
            for i,j in pairs:
                ar, br = x.iloc[i], y.iloc[j]
                rec = dict(series=series, a=a, b=b, role=role, chain_a=ar.chain, chain_b=br.chain, position_a=i+1, position_b=j+1, index_a=int(ar['index']), index_b=int(br['index']), id_a=ar.id, id_b=br.id, aa_a=ar.aa, aa_b=br.aa, same_aa=ar.aa == br.aa, region_a=ar.region, region_b_legacy=br.legacy_region, region_b=br.region, CDR_annotation_source=ar.CDR_annotation_source, interface_a=bool(ar.interface), interface_b=bool(br.interface), interface_union=bool(ar.interface or br.interface))
                for m in ['rmsf','hchi']:
                    valid = bool(ar[m+'_valid'] and br[m+'_valid'])
                    rec[m+'_valid'] = valid
                    for s, r in [('a', ar), ('b', br)]:
                        for value in ['md', 'var']+MODELS:
                            rec[f'{m}_{value}_{s}'] = float(r[f'{m}_{value}']) if valid else np.nan
                    rec[m+'_md_delta'] = rec[m+'_md_b'] - rec[m+'_md_a']
                    rec[m+'_delta_sem'] = float(np.sqrt((rec[m+'_var_a']+rec[m+'_var_b'])/3))
                    for model in MODELS:
                        rec[f'{m}_{model}_delta'] = rec[f'{m}_{model}_b'] - rec[f'{m}_{model}_a']
                rows.append(rec)
        f = pd.DataFrame(rows)
        near_related = len(tcr_stats)==2 and all(identity >= .9 and coverage >= .85 for identity,coverage in tcr_stats)
        f['TCR_near_related_alignment'] = near_related
        positions.append(f)
        tcr = f.role.isin(['tcr_alpha','tcr_beta'])
        scopes = {'all_CDR': tcr & f.region_a.isin(CDRS), 'CDR3': tcr & f.region_a.str.contains('CDR3', na=False), 'whole_peptide':f.role.eq('peptide'), 'natural_MHC_groove': f.role.eq('mhc') & (f.position_a<=180) & (f.position_b<=180), 'natural_pMHC_interface': f.role.isin(['mhc','peptide']) & f.interface_union}
        scopes.update({cdr:tcr & f.region_a.eq(cdr) for cdr in CDRS})
        for scope, selector in scopes.items():
            is_tcr_scope = scope == 'all_CDR' or scope.startswith('CDR')
            interpretable = not is_tcr_scope or near_related
            chemistries = ['same_AA'] if scope.startswith('natural_') else ['including_substitutions', 'same_AA']
            for chemistry in chemistries:
                selected = selector & (f.same_aa if chemistry=='same_AA' else True)
                for metric in ['rmsf','hchi']:
                    g = f[selected & f[metric+'_valid']]
                    for model in MODELS:
                        n = len(g)
                        rec = dict(series=series,a=a,b=b,direction='b minus a',scope=scope,chemistry=chemistry,metric=metric,model=model,n_aligned_selected=int(selected.sum()),n=n,n_substituted=int((~g.same_aa).sum()),TCR_near_related_alignment=near_related,region_interpretable=interpretable,interpretation='Cross-unrelated TCR loops not evaluated as complete homologous CDR regions' if not interpretable else ('Changed chemistry retained; Hchi is not pure environment response' if chemistry=='including_substitutions' and metric=='hchi' and (~g.same_aa).any() else 'Complete fixed region'))
                        x = g[metric+'_md_delta'].to_numpy()
                        y = g[f'{metric}_{model}_delta'].to_numpy()
                        sem = g[metric+'_delta_sem'].to_numpy()
                        if n and interpretable:
                            zero = float(np.mean(abs(x))); mae = float(np.mean(abs(x-y))); robust = abs(x)>2*sem
                            rec.update(rho_delta=rho(x,y),MAE_delta=mae,zero_change_MAE=zero,relative_MAE_gain=1-mae/zero if zero>1e-12 else np.nan,direction_fraction=float(np.mean(x*y>0)),MD_abs_delta=zero,model_abs_delta=float(np.mean(abs(y))),MD_mean_delta=float(np.mean(x)),model_mean_delta=float(np.mean(y)),mean_delta_SEM=float(np.mean(sem)),descriptor_above_2SEM_n=int(robust.sum()),descriptor_above_2SEM_direction=float(np.mean(x[robust]*y[robust]>0)) if robust.any() else np.nan,descriptor_above_2SEM_MAE=float(np.mean(abs(x[robust]-y[robust]))) if robust.any() else np.nan,descriptor_above_2SEM_zero_MAE=float(np.mean(abs(x[robust]))) if robust.any() else np.nan,absolute_MAE_a=float(np.mean(abs(g[metric+'_md_a']-g[f'{metric}_{model}_a']))),absolute_MAE_b=float(np.mean(abs(g[metric+'_md_b']-g[f'{metric}_{model}_b']))))
                        summary.append(rec)
        print(series,a,b,'aligned',len(f),'near_related_TCR',near_related,flush=True)

res = pd.DataFrame(summary)
res.to_csv(OUT/'all_prespecified_region_model_metrics.csv',index=False,encoding='utf-8-sig')
pd.concat(positions,ignore_index=True).to_csv(OUT/'all_aligned_position_changes.csv.gz',index=False,compression='gzip')
pd.DataFrame(alignment_audit).to_csv(OUT/'complete_chain_alignment_audit.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(bounds).to_csv(OUT/'all_CDR_boundaries_and_sequences.csv',index=False,encoding='utf-8-sig')
main = res[(res.model=='full') & res.region_interpretable & ((res.metric=='rmsf') | res.chemistry.eq('same_AA')) & ((res.chemistry=='including_substitutions') | res.scope.str.startswith('natural_') | res.metric.eq('hchi'))]
main.to_csv(OUT/'primary_full_complete_regions.csv',index=False,encoding='utf-8-sig')
audit = dict(status='complete',n_PDB=len({p for group in SERIES.values() for p in group}),n_pairs=res[['a','b']].drop_duplicates().shape[0],n_series=len(SERIES),series=SERIES,model_rule='Existing Full joint2output family fold0 primary; intrachain fixed control. Training exclusion independently audited by another agent.',checkpoint_sha256='ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',CDR_rule='All six legacy a-state CDR regions projected through full-chain global alignment to b. Complete aligned RMSF positions retain substitutions. Distinct same-AA and including-substitution results. No standalone CDR metrics for unrelated TCR (both variable chains require >=.90 identity and >=.85 coverage).',interface_rule='Union of existing annotations is_tcr_pmhc_interface_8A on either state; existing C-alpha 8-Angstrom interface annotations, not newly recomputed coordinates.',groove_rule='Full-chain alignment; positions <=180 in each MHC, identical AA only.',SEM_rule='sqrt((label_variance_a + label_variance_b)/3), three existing MD replicas per state. abs(delta)>2SEM is descriptive and not formal significance.',alignment=dict(mode='global',match=2,mismatch=-1,gap_open=-6,gap_extend=-.5,protein_end_gap=0,peptide_end_gap='same as internal; avoids clipping unrelated equal-length antigens'),new_training=0,new_MD=0,manuscript_edits=0,annotation_source=str(ANNOTATION),annotation_sha256=sha(ANNOTATION),source_cache_sha256={f'{model}/itrace_{p}.npz':sha(OLD/'family_predictions'/model/f'itrace_{p}.npz') for model in ['full','intrachain'] for p in sorted({p for group in SERIES.values() for p in group})})
(OUT/'rescreen_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
audit['fresh_predictions']={model:{p:sha(source/f'itrace_{p}.npz') for p in sorted({p for group in SERIES.values() for p in group})} for model,source in FRESH.items()}
audit['additional_existing_model']={'name':'old_strict','checkpoint_sha256':'b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f','fresh_inference':str(FRESH['old_strict'])}
audit['CDR_correction']='MEL5 WT/engineered/multiple-antigen series uses 3HG1 canonical WT six CDR boundaries before pairwise projection; fixes 4JFF beta N-terminal numbering offset.'
(OUT/'rescreen_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print('DONE',audit['n_pairs'],'pairs',len(res),'region/model/chemistry/metric rows',flush=True)
