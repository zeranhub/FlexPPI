"""Independent read-only verification of already saved predictions and tables."""
from pathlib import Path
from collections import Counter
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

OUT = Path(__file__).resolve().parent
BASE = Path(r'E:\Project\FlexPPI\experiments\_sources')
GROUP = BASE / 'local_model_benchmark_20261023/grouping_protein_complex_20261024'
COHORT = BASE / 'local_model_benchmark_20261023/cohort'
STRUCT = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures')
LABEL = STRUCT.parent / 'labels'
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
report = {}
issues = []

def check(name, condition, detail=None):
    report[name] = bool(condition)
    if not condition:
        issues.append({'check': name, 'detail': detail})

def rho(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 4 or np.ptp(a[ok]) == 0 or np.ptp(b[ok]) == 0:
        return np.nan
    return float(pd.Series(a[ok]).rank().corr(pd.Series(b[ok]).rank()))

def aa_residual(values, aa):
    values = np.asarray(values, float)
    result = np.full(len(values), np.nan)
    for typ in set(aa):
        sel = (aa == typ) & np.isfinite(values)
        if sel.sum() >= 4:
            result[sel] = values[sel] - values[sel].mean()
    return result

def calc(d):
    output = {}
    for metric in ('rmsf', 'chi_entropy'):
        selected = d[d['valid_' + metric]]
        y = selected['target_' + metric].to_numpy(float)
        f = selected['pred_' + metric].to_numpy(float)
        i = selected['intra_' + metric].to_numpy(float)
        aa = selected.amino_acid.to_numpy(str)
        output[metric + '_n'] = len(y)
        if len(y):
            full_mae = float(np.abs(y - f).mean())
            intra_mae = float(np.abs(y - i).mean())
            output.update({metric + '_full_mae': full_mae,
                           metric + '_intra_mae': intra_mae,
                           metric + '_full_rho': rho(y, f),
                           metric + '_intra_rho': rho(y, i),
                           metric + '_full_type_centered_rho': rho(aa_residual(y, aa), aa_residual(f, aa)),
                           metric + '_intra_type_centered_rho': rho(aa_residual(y, aa), aa_residual(i, aa)),
                           metric + '_mae_gain_fraction': 1 - full_mae / intra_mae if intra_mae else np.nan,
                           metric + '_fraction_sites_full_lower_error': float((np.abs(y-f) < np.abs(y-i)).mean())})
    return output

def ca_residues(path):
    residues = {}
    for line in Path(path).read_text().splitlines():
        if line[:6] not in ('ATOM  ', 'HETATM') or line[12:16].strip() != 'CA' or line[16:17] not in (' ', 'A'):
            continue
        name = line[17:20].strip()
        if name not in AA:
            continue
        key = (line[21:22].strip() or '_', line[22:26].strip()+line[26:27].strip())
        residues[key] = (AA[name], tuple(float(line[a:b]) for a,b in ((30,38),(38,46),(46,54))))
    return [(key[0], key[1], item[0], item[1]) for key,item in residues.items()]

def compare_table(table, mode, profiles):
    max_error = 0.0
    mismatches = []
    grouped = {sid: d for sid,d in profiles.groupby('sample_id', sort=False)}
    for _, row in table.iterrows():
        d = grouped[row.sample_id]
        if mode == 'chain':
            d = d[d.chain == row.chain]
        if mode == 'region':
            if row.region == 'interface_ca8': d = d[d.interface_ca8]
            elif row.region == 'outside_interface': d = d[~d.interface_ca8]
            elif row.region == 'remote_ca12': d = d[d.partner_ca_distance > 12]
            else: raise ValueError(row.region)
        for name, value in calc(d).items():
            actual = row[name]
            if np.isnan(value) and pd.isna(actual): continue
            error = abs(float(value)-float(actual))
            if not np.isfinite(error) or error > 1e-10:
                mismatches.append((row.sample_id, name, value, actual))
            else: max_error = max(max_error, error)
    return {'rows': len(table), 'max_absolute_error': max_error, 'mismatches': mismatches[:20], 'mismatch_count': len(mismatches)}

profiles = pd.read_csv(OUT/'all_non_tcr_profiles.csv.gz', dtype={'author_residue_id':str, 'chain':str})
summary = pd.read_csv(OUT/'all_616_record_screen.csv')
chain_table = pd.read_csv(OUT/'all_chain_metrics.csv')
regions = pd.read_csv(OUT/'all_region_metrics.csv')
members = pd.read_csv(GROUP/'membership_frozen.csv').fillna('')
records = pd.read_csv(COHORT/'records.csv').set_index('sample_id')
identities = pd.read_csv(GROUP/'chain_biological_identity_evidence.csv').fillna('')
folds = pd.read_csv(BASE/'family_folds/record_family_folds.csv').set_index('sample_id')
per_record = pd.read_csv(BASE/'family_folds/evaluation/per_record_metrics.csv')
wanted = set(members.loc[(members.benchmark_group=='protein_complex') & ~members.sample_id.str.contains(r'^\d[a-z0-9]{3}$', regex=True), 'sample_id'])
wanted &= set(per_record.loc[per_record.dataset != 'itrace', 'sample_id'])
check('requested_and_output_sample_ids_exactly_match_616', len(wanted)==616 and set(profiles.sample_id)==set(summary.sample_id)==wanted)
check('profile_sample_residue_keys_unique', not profiles.duplicated(['sample_id','residue_index']).any())
check('summary_sample_keys_unique', not summary.sample_id.duplicated().any())
check('no_itrace_rows', not (profiles.dataset=='itrace').any())

raw_parts=[]
for chunk in pd.read_csv(next((OUT/'archive_extract').rglob('oof_per_residue.csv.gz')), chunksize=200000,
                         dtype={'author_residue_id':str,'chain':str}):
    raw_parts.append(chunk[chunk.sample_id.isin(wanted)])
raw = pd.concat(raw_parts, ignore_index=True)
key = ['sample_id','residue_index']
full = raw[raw.variant=='full'].set_index(key).sort_index()
intra = raw[raw.variant=='intrachain'].set_index(key).sort_index()
check('raw_full_intra_residue_keys_equal_unique', full.index.equals(intra.index) and full.index.is_unique and intra.index.is_unique)
metadata_checks = {}
for column in ('dataset','fold','amino_acid','chain','author_residue_id','valid_rmsf','valid_chi_entropy','target_rmsf','target_chi_entropy'):
    a,b=full[column],intra[column]
    if pd.api.types.is_numeric_dtype(a):
        equal=np.allclose(a,b,rtol=0,atol=0,equal_nan=True)
    else: equal=a.fillna('<NA>').equals(b.fillna('<NA>'))
    metadata_checks[column] = bool(equal)
check('raw_full_intra_all_metadata_and_labels_identical', all(metadata_checks.values()), metadata_checks)
report['raw_metadata_field_checks'] = metadata_checks
p = profiles.set_index(key).sort_index()
check('profile_keys_equal_raw_full', p.index.equals(full.index))
source_checks = {}
for column in ('fold','chain','author_residue_id','amino_acid','valid_rmsf','valid_chi_entropy','target_rmsf','target_chi_entropy','pred_rmsf','pred_chi_entropy'):
    a,b=p[column],full[column]
    if pd.api.types.is_numeric_dtype(a): equal=np.allclose(a,b,rtol=0,atol=1e-12,equal_nan=True)
    else: equal=a.fillna('<NA>').equals(b.fillna('<NA>'))
    source_checks[column]=bool(equal)
for metric in ('rmsf','chi_entropy'):
    source_checks['intra_'+metric] = bool(np.allclose(p['intra_'+metric],intra['pred_'+metric],rtol=0,atol=1e-12,equal_nan=True))
check('saved_profiles_exactly_match_raw_source', all(source_checks.values()), source_checks)
check('all_profile_folds_match_frozen_record_folds', all(set(d.fold)=={int(folds.loc[sid,'fold'])} for sid,d in profiles.groupby('sample_id')))

mapping_failures=[]
geometry_failures=[]
identity_failures=[]
length_failures=[]
for sid,d in profiles.groupby('sample_id',sort=False):
    d=d.sort_values('residue_index')
    r=records.loc[sid]
    if len(d)!=r.n_residues or d.chain.nunique()!=r.n_chains:
        length_failures.append((sid,len(d),int(r.n_residues),d.chain.nunique(),int(r.n_chains)))
    residues=ca_residues(d.structure_path.iloc[0])
    observed=list(zip(d.chain,d.author_residue_id,d.amino_acid))
    if observed != [x[:3] for x in residues]:
        mapping_failures.append(sid)
        continue
    identity_map=identities[identities.sample_id==sid].drop_duplicates('input_chain',keep='first').set_index('input_chain').biological_identity_evidence.to_dict()
    if any(not identity_map.get(chain) for chain in d.chain.unique()):
        identity_failures.append((sid,'missing evidence'))
        continue
    if list(d.partner_identity) != [identity_map[c] for c in d.chain]:
        identity_failures.append((sid,'profile evidence mismatch'))
    xyz=np.array([x[3] for x in residues])
    ident=d.partner_identity.to_numpy(str)
    nearest=np.full(len(d),np.inf)
    partner_neighbors=np.zeros(len(d),int)
    for typ in set(ident):
        target=ident==typ
        if (~target).sum():
            tree=cKDTree(xyz[~target])
            nearest[target]=tree.query(xyz[target],k=1)[0]
            partner_neighbors[target]=[len(items) for items in tree.query_ball_point(xyz[target],r=8)]
    tree=cKDTree(xyz)
    all_neighbors=np.array([len(items)-1 for items in tree.query_ball_point(xyz,r=8)])
    if not np.allclose(nearest,d.partner_ca_distance,rtol=0,atol=1e-10) or not np.array_equal(nearest<=8,d.interface_ca8) or not np.array_equal(partner_neighbors,d.partner_ca_neighbors8) or not np.array_equal(all_neighbors,d.all_ca_neighbors8):
        geometry_failures.append(sid)
check('all_616_record_lengths_and_chain_counts_match_cohort', not length_failures, length_failures[:20])
check('all_616_PDB_chain_author_residue_and_sequence_mapping_exact', not mapping_failures, mapping_failures[:20])
check('all_616_biological_identity_evidence_complete_and_matches', not identity_failures, identity_failures[:20])
ambiguous_identity_groups=identities[identities.sample_id.isin(wanted)].groupby(['sample_id','input_chain']).biological_identity_evidence.nunique()
report['ambiguous_identity_evidence_groups']=[{'sample_id':sid,'input_chain':chain,'candidate_identities':int(n)} for (sid,chain),n in ambiguous_identity_groups.items() if n>1]
check('all_616_interface_and_density_geometry_reproduced_independently', not geometry_failures, geometry_failures[:20])

report['table_recomputations']={}
for name, table, mode in [('records',summary,'record'),('chains',chain_table,'chain'),('regions',regions,'region')]:
    result=compare_table(table,mode,profiles)
    report['table_recomputations'][name]=result
    check('recomputed_'+name+'_metrics_match',result['mismatch_count']==0,result['mismatches'])

case_details={}
for sid in ('A004F','A007U'):
    d=profiles[profiles.sample_id==sid].sort_values('residue_index')
    labels=np.load(LABEL/(sid+'.npz'),allow_pickle=False)
    names=list(labels['metric_names'])
    label_results={}
    for metric in ('rmsf','chi_entropy'):
        index=names.index(metric)
        label_results[metric+'_values']=bool(np.array_equal(labels['values'][:,index],d['target_'+metric].to_numpy(dtype=labels['values'].dtype),equal_nan=True))
        label_results[metric+'_mask']=bool(np.array_equal(labels['mask'][:,index],d['valid_'+metric]))
    centered_counts=[]
    subsets=[('record',d)]+[('chain_'+chain,cd) for chain,cd in d.groupby('chain')]+[
        ('interface_ca8',d[d.interface_ca8]),('outside_interface',d[~d.interface_ca8]),('remote_ca12',d[d.partner_ca_distance>12])]
    for name,sub in subsets:
        for metric in ('rmsf','chi_entropy'):
            valid=sub[sub['valid_'+metric]]
            aa=valid.amino_acid.to_numpy(str)
            residual=aa_residual(valid['target_'+metric].to_numpy(float),aa)
            counts=Counter(aa)
            centered_counts.append({'subset':name,'metric':metric,'unadjusted_valid_n':len(valid),
                'sameAA_effective_n':int(np.isfinite(residual).sum()),'retained_aa_types':sum(n>=4 for n in counts.values())})
    case_details[sid]={'pdb_id':d.pdb_id.iloc[0],'folds':sorted(set(d.fold)),
        'rows':len(d),'chains':d.groupby('chain').size().to_dict(),
        'identities':d.groupby('chain').partner_identity.first().to_dict(),
        'NPZ_direct_label_checks':label_results,'NPZ_comparison_rule':'CSV targets cast back to original NPZ float32; require exact numeric equality including NaN positions',
        'NPZ_csv_decimal_rounding_max_abs':{m:float(np.nanmax(np.abs(labels['values'][:,names.index(m)]-d['target_'+m].to_numpy(float)))) for m in ('rmsf','chi_entropy')},
        'sameAA_effective_counts':centered_counts}
    check(sid+'_direct_label_NPZ_matches',all(label_results.values()),label_results)
report['case_details']=case_details
report['profile_rows']=len(profiles)
report['records']=profiles.sample_id.nunique()
report['chain_rows']=len(chain_table)
report['region_rows']=len(regions)
report['issues']=issues
(OUT/'independent_numeric_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
