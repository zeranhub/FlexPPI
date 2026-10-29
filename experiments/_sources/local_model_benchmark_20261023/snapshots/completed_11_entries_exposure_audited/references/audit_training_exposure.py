"""Local ID-only exposure audit; never import torch, fit, predict or select records."""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import re
from collections import Counter, defaultdict

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DYNA = ROOT/'vendor/dynaprot-main'
LISTS = DYNA/'dynaprot/data/preprocessing/protein_lists'
BF1 = Path('E:/Project/FlexPPI_baselines')
ROLE = ('train', 'val', 'test')
SOURCES = set()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def text(path):
    SOURCES.add(Path(path))
    return Path(path).read_text(encoding='utf-8-sig')


def csv_read(path):
    SOURCES.add(Path(path))
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def csv_write(name, rows):
    with (HERE/name).open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def normalized_id(value):
    found = re.fullmatch(r'([0-9A-Za-z]{4})[._]([A-Za-z0-9]+)', str(value))
    assert found, value
    # PDB accession case is immaterial; author-chain case is deliberately preserved.
    return found[1].lower()+'_'+found[2]


def main():
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    records = csv_read(ROOT/'cohort/records.csv')
    chains = csv_read(ROOT/'cohort/chains.csv')
    record_by_id = {r['sample_id']:r for r in records}
    assert len(records) == len(record_by_id) == 3063 and len(chains) == 4973
    assert set(c['sample_id'] for c in chains) == set(record_by_id)
    split_sets = {}
    official_rows = []
    for family in ('atlas', 'flexpert'):
        split_sets[family] = {}
        for role in ROLE:
            path = LISTS/(family+'_chains_'+role+'.npy')
            SOURCES.add(path)
            array = np.load(path, allow_pickle=False)
            assert array.ndim == 1 and array.dtype.kind == 'U'
            ids = [normalized_id(x) for x in array.tolist()]
            assert len(ids) == len(set(ids))
            split_sets[family][role] = set(ids)
            official_rows.extend(dict(split_family=family, split=role, pdb_chain_id=x,
                source=str(path), source_sha256=sha(path)) for x in ids)
        assert all(not split_sets[family][a] & split_sets[family][b]
                   for a,b in [('train','val'),('train','test'),('val','test')])
    flexpert_path = ROOT/'vendor/Flexpert-main/data/atlas_splits.json'
    flexpert_original = json.loads(text(flexpert_path))
    cross_check = {}
    for role in ROLE:
        key = 'validation' if role == 'val' else role
        original = {normalized_id(x) for x in flexpert_original[key]}
        cross_check[role] = dict(dynaprot_npy_ids=len(split_sets['flexpert'][role]),
                                official_flexpert_json_ids=len(original), exact_set_equal=True)
        assert split_sets['flexpert'][role] == original
    all_path = LISTS/'atlas_chains.npy'
    SOURCES.add(all_path)
    all_atlas = {normalized_id(x) for x in np.load(all_path, allow_pickle=False).tolist()}
    atlas_union = set().union(*split_sets['atlas'].values())
    data_cfg = yaml.safe_load(text(ROOT/'weights/dynaprot/dynaprot-M/data_config.yaml'))
    assert Path(data_cfg['protein_chains_path']).name == 'atlas_chains.npy'
    text(ROOT/'weights/dynaprot/dynaprot-M/README.md')
    text(DYNA/'configs/README.md')
    text(DYNA/'configs/data/atlas_config_randomcropping.yaml')
    text(DYNA/'dynaprot/data/datasets.py')
    text(DYNA/'dynaprot/config.py')
    text(DYNA/'dynaprot/train.py')
    text(BF1/'BackFlip/README.md')
    text(ROOT/'cohort/README.md')
    text(ROOT/'cohort/annotation_architecture_candidates.csv')
    bindings = {}
    dyna_meta = json.loads(text(ROOT/'dynaprot/predictions/metadata.json'))
    bindings['DynaProt-M'] = dict(split_family='atlas',
        tier='checkpoint_config_bound_published_training_list',
        checkpoint_sha256=dyna_meta['checkpoint_sha256'],
        data_list_base=data_cfg['protein_chains_path'],
        rule='Released checkpoint folder data_config plus official config mapping; Dataset appends _train/_val/_test and resolves packaged list basename.',
        actual_weight_training_history_verified=False)
    for tag in ('backflip-1.0', 'backflip-1.0-seq'):
        config = yaml.safe_load(text(BF1/'model_cache/backflip'/tag/'config.yaml'))
        paths = {role:config['data']['dataset'][role+'_csv_path'] for role in ROLE}
        assert all(Path(path).name == 'flexpert_'+role+'.csv' for role,path in paths.items())
        meta = json.loads(text(ROOT/'backflip/results'/tag/'run_metadata.json'))
        bindings[tag] = dict(split_family='flexpert',
            tier='checkpoint_split_name_bound_official_candidate_ids_csv_not_available',
            checkpoint_sha256=meta['checkpoint_sha256'], configured_paths=paths,
            rule='Released config and README bind Flexpert split name; exact release CSV files are absent in searched local sources. DynaProt NPY agrees exactly with original Flexpert JSON, but does not byte-bind the absent checkpoint CSV.',
            actual_weight_training_history_verified=False)
    meta = json.loads(text(ROOT/'backflip/results/backflip-2.1/run_metadata.json'))
    config = yaml.safe_load(text(ROOT/'backflip/weights/backflip-2.1/config.yaml'))
    bindings['backflip-2.1'] = dict(split_family='', tier='unknown_no_bound_id_manifest',
        checkpoint_sha256=meta['checkpoint_sha256'],
        configured_paths={role:config['data']['dataset'][role+'_csv_path'] for role in ROLE},
        rule='Config names alphaflow_*.csv; files not located. Naming resemblance does not establish equality with DynaProt atlas lists.',
        actual_weight_training_history_verified=False)
    details = []
    for model, binding in bindings.items():
        family = binding['split_family']
        lists = split_sets.get(family, {})
        for c in chains:
            r = record_by_id[c['sample_id']]
            assert re.fullmatch('pdb:[A-Za-z0-9]{4}', r['biological_system'])
            pdb_id = r['biological_system'][4:].lower()
            if r['dataset'] == 'atlas':
                resolved = normalized_id(r['sample_id'])
                assert resolved[:4] == pdb_id and r['n_chains'] == '1'
                resolution = 'resolved_original_chain_from_ATLAS_sample_id'
                source_chain = resolved[5:]
            else:
                # Frozen MD labels are valid within-record keys, not a certified
                # author-chain mapping back to the deposited PDB. Do not promote.
                resolved, source_chain = '', ''
                resolution = 'unresolved_MD_chain_to_original_PDB_author_chain'
            label_candidate = pdb_id+'_'+c['chain'] if c['chain'] != '_' else ''
            exact_roles = [role for role,ids in lists.items() if resolved and resolved in ids]
            same_pdb_roles = [role for role,ids in lists.items() if any(x[:4] == pdb_id for x in ids)]
            label_roles = [role for role,ids in lists.items() if label_candidate and label_candidate in ids]
            if not family:
                match_status = 'unknown_checkpoint_training_ID_manifest'
            elif resolved:
                match_status = 'exact_PDB_chain_ID_match' if exact_roles else 'no_exact_ID_in_this_manifest'
            elif same_pdb_roles:
                match_status = 'same_PDB_only_original_chain_unresolved'
            else:
                match_status = 'PDB_absent_from_this_manifest_original_chain_unresolved'
            documented = ('yes' if 'train' in exact_roles else 'no_exact_train_ID_in_this_manifest') if (
                model == 'DynaProt-M' and resolved) else 'unresolved_or_not_checkpoint_bound'
            details.append(dict(model=model, sample_id=r['sample_id'], dataset=r['dataset'],
                group=r['main_group_candidate'], input_id=c['input_id'], cohort_chain=c['chain'],
                n_residues=c['n_residues'], pdb_id=pdb_id,
                original_author_chain=source_chain, original_PDB_chain_ID=resolved,
                identifier_resolution=resolution, split_family=family,
                checkpoint_binding=binding['tier'], exact_ID_split=';'.join(exact_roles),
                same_PDB_split=';'.join(same_pdb_roles),
                current_MD_label_candidate_ID=label_candidate,
                current_MD_label_candidate_split_NOT_exact=';'.join(label_roles),
                match_status=match_status, documented_training_list_exposure=documented))
    csv_write('training_exposure_chain_details.csv', details)
    record_details = []
    for (model,sample), items in group_by(details, lambda x:(x['model'],x['sample_id'])).items():
        r = record_by_id[sample]
        record_details.append(dict(model=model, sample_id=sample, dataset=r['dataset'],
            group=r['main_group_candidate'], pdb_id=items[0]['pdb_id'],
            cohort_chains=';'.join(x['cohort_chain'] for x in items),
            all_original_chain_IDs_resolved=all(x['original_PDB_chain_ID'] for x in items),
            original_PDB_chain_IDs=';'.join(x['original_PDB_chain_ID'] for x in items if x['original_PDB_chain_ID']),
            split_family=items[0]['split_family'], checkpoint_binding=items[0]['checkpoint_binding'],
            exact_train_chains=sum('train' in x['exact_ID_split'].split(';') for x in items),
            exact_val_chains=sum('val' in x['exact_ID_split'].split(';') for x in items),
            exact_test_chains=sum('test' in x['exact_ID_split'].split(';') for x in items),
            same_PDB_split=';'.join(sorted(set(y for x in items for y in x['same_PDB_split'].split(';') if y))),
            match_status=';'.join(sorted(set(x['match_status'] for x in items)))))
    csv_write('training_exposure_record_details.csv', record_details)
    csv_write('training_exposure_official_split_ids.csv', official_rows)
    summary = []
    for (model,dataset), items in group_by(record_details, lambda x:(x['model'],x['dataset'])).items():
        summary.append(dict(model=model, dataset=dataset, records=len(items),
            original_chain_ID_resolved_records=sum(bool(x['all_original_chain_IDs_resolved']) for x in items),
            exact_train_ID_records=sum(x['exact_train_chains'] > 0 for x in items),
            exact_val_ID_records=sum(x['exact_val_chains'] > 0 for x in items),
            exact_test_ID_records=sum(x['exact_test_chains'] > 0 for x in items),
            same_PDB_only_unresolved_records=sum(bool(x['same_PDB_split']) and not x['all_original_chain_IDs_resolved'] for x in items),
            checkpoint_binding=items[0]['checkpoint_binding']))
    csv_write('training_exposure_summary.csv', summary)
    own_audit = audit_own_folds(records)
    membership = []
    for r in record_details:
        sample = record_by_id[r['sample_id']]
        binding = bindings[r['model']]
        roles = [role for role in ROLE if r['exact_'+role+'_chains']>0]
        assert len(roles)<=1
        if not binding['split_family']:
            documented_split, match = 'unknown', 'unknown_training_list'
        elif roles:
            documented_split, match = roles[0], 'yes'
        elif r['all_original_chain_IDs_resolved']:
            documented_split, match = 'unlisted', 'no'
        elif r['same_PDB_split']:
            documented_split, match = 'unresolved', 'unresolved'
        else:
            documented_split, match = 'unlisted', 'no_PDB_in_list_original_chain_unresolved'
        membership.append(dict(method=r['model'],checkpoint_sha256=binding['checkpoint_sha256'],
            sample_id=r['sample_id'],dataset=r['dataset'],group=r['group'],
            documented_split=documented_split,binding_strength=binding['tier'],
            exact_identifier_match=match,identifier_resolution=r['match_status'],
            pdb_id=r['pdb_id'],original_PDB_chain_IDs=r['original_PDB_chain_IDs'],
            official_split_family=binding['split_family'],same_PDB_split=r['same_PDB_split'],
            original_fold=sample['original_fold'],family_fold=sample['family_fold'],
            own_same_family_in_train='',snapshot_utc=timestamp))
    with (HERE/'training_exposure_own_oof_details.csv').open(encoding='utf-8-sig',newline='') as handle:
        for r in csv.DictReader(handle):
            membership.append(dict(method='Full' if r['variant']=='full' else 'Intra',
                checkpoint_sha256=r['checkpoint_sha256'],sample_id=r['sample_id'],dataset=r['dataset'],group=r['group'],
                documented_split='record_biological_group_OOF',
                binding_strength='prediction_manifest_and_frozen_train_test_files_confirmed',
                exact_identifier_match='heldout_record_not_in_corresponding_train',
                identifier_resolution='sample_id_exact',pdb_id=record_by_id[r['sample_id']]['biological_system'][4:].lower(),
                original_PDB_chain_IDs='',official_split_family='original_record_biological_group_five_fold',same_PDB_split='',
                original_fold=r['original_fold'],family_fold=r['family_fold'],
                own_same_family_in_train=int(r['same_family_other_train_records'])>0,snapshot_utc=timestamp))
    assert len(membership) == 18378
    csv_write('training_exposure_membership_frozen.csv', membership)
    protocol = dict(generated_utc=timestamp, activity='read_only_local_exact_identifier_audit',
        records=3063, chains=4973, no_torch_import=True, checkpoint_tensors_not_loaded=True,
        no_network=True, no_inference=True, no_training=True, no_sequence_similarity_or_clustering=True,
        cohort_and_main_ranking_unchanged=True,
        match_rule='Normalize only PDB accession case and dot/underscore separator; author-chain case retained. Both PDB and original author chain required for exact match. MD chain labels are not silently treated as original PDB chains.',
        ATLAS_resolution='sample_id encodes original PDB+chain; frozen structure chain _ is not used as original author chain.',
        non_ATLAS_resolution='Original PDB from frozen biological_system; frozen MD chain labels retained separately. Existing cohort audit warns source author chains may be renamed; no certified inverse map located. PDB-only/current-label matches do not count as exact exposure.',
        exposure_scope='Documented training-list exposure is an ID overlap with published checkpoint-bound split; it does not prove historical SGD inclusion, coordinate/trajectory equality, or absence/presence of sequence-homology leakage.',
        sensitivity_membership_frozen_before_performance_read=True,
        performance_tables_read_in_this_audit=[],
        sensitivity_policy='Audit-derived sensitivity only; report all train/val/test/unlisted/unresolved/unknown categories. No performance-based selection, no changes to main cohort, predictions or primary rankings.',
        membership_schema=list(membership[0]),membership_rows=len(membership),
        limitations='Absence from one public ID list does not mean no training exposure via homologs, other data or pretraining. val is not train; test membership alone does not establish independence of this common-reference evaluation.',
        bindings=bindings, official_split_sizes={f:{s:len(ids) for s,ids in sets.items()} for f,sets in split_sets.items()},
        own_Full_Intra_OOF=own_audit,
        flexpert_NPY_vs_original_JSON=cross_check,
        atlas_base_list_not_in_three_splits=sorted(all_atlas-atlas_union),
        atlas_split_IDs_not_in_base=sorted(atlas_union-all_atlas),
        summaries=summary,
        source_commits=dict(dynaprot='e1bc8e9f5d80afd91f79ea6e8def8029121a8e6e',
                            backflip1='60b8c83511b965925bbe990bb2be578bf3529460',
                            backflip2='887c6a73b216e3fe2ff9f551ae6e2c983c21b292'),
        source_hashes=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(SOURCES)],
        script_sha256=sha(Path(__file__)),
        output_hashes={n:sha(HERE/n) for n in ['training_exposure_chain_details.csv','training_exposure_record_details.csv',
            'training_exposure_official_split_ids.csv','training_exposure_summary.csv',
            'training_exposure_own_oof_details.csv','training_exposure_own_family_details.csv',
            'training_exposure_membership_frozen.csv']})
    (HERE/'training_exposure_protocol_sources_sha256.json').write_text(json.dumps(protocol,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(dict(summaries=summary, source_count=len(SOURCES), chain_detail_rows=len(details),
                         record_detail_rows=len(record_details), flexpert_cross_check=cross_check,
                         own_Full_Intra_OOF=own_audit),indent=2),flush=True)


def audit_own_folds(records):
    core = Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
    export = core/'neural_exports/evaluation'
    by_sample = {r['sample_id']:r for r in records}
    family_groups = group_by(records, lambda r:r['family_component_id'])
    biological_groups = group_by(records, lambda r:r['biological_cluster_id'])
    family_rows = []
    for ident, items in family_groups.items():
        original = sorted({int(r['original_fold']) for r in items})
        family = sorted({int(r['family_fold']) for r in items})
        family_rows.append(dict(family_component_id=ident, records=len(items),
            datasets=';'.join(sorted({r['dataset'] for r in items})),
            original_folds=';'.join(map(str,original)),family_folds=';'.join(map(str,family)),
            crosses_original_folds=len(original)>1,crosses_family_folds=len(family)>1))
    csv_write('training_exposure_own_family_details.csv', family_rows)
    train_ids, test_ids = {}, {}
    for fold in range(5):
        for role, store in [('train',train_ids), ('test',test_ids)]:
            path = core/'training_payload/splits'/f'fold_{fold}_{role}.jsonl'
            SOURCES.add(path)
            with path.open(encoding='utf-8') as handle:
                ids = [json.loads(line)['sample_id'] for line in handle if line.strip()]
            assert len(ids)==len(set(ids))
            store[fold] = set(ids)
        assert not train_ids[fold]&test_ids[fold]
        assert train_ids[fold]|test_ids[fold] == set(by_sample)
        assert test_ids[fold] == {r['sample_id'] for r in records if int(r['original_fold'])==fold}
    path = export/'prediction_manifest.jsonl'
    SOURCES.add(path)
    details, seen, versions = [], set(), Counter()
    with path.open(encoding='utf-8') as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            variant, sample, fold = row['variant'], row['sample_id'], int(row['fold'])
            assert variant in ('full','intrachain') and (variant,sample) not in seen
            seen.add((variant,sample))
            r = by_sample[sample]
            metadata = row['record']['metadata']
            members = family_groups[r['family_component_id']]
            training_family = [x for x in members if x['sample_id'] in train_ids[fold]]
            training_biology = [x for x in biological_groups[r['biological_cluster_id']]
                                if x['sample_id'] in train_ids[fold]]
            versions[metadata.get('source_split_version','')] += 1
            assert int(metadata['fold']) == fold == int(r['original_fold'])
            assert sample in test_ids[fold] and sample not in train_ids[fold]
            details.append(dict(variant=variant,sample_id=sample,dataset=r['dataset'],
                group=r['main_group_candidate'],prediction_manifest_fold=fold,
                nested_metadata_fold=metadata['fold'],original_fold=r['original_fold'],
                family_fold=r['family_fold'],matches_original_fold=fold==int(r['original_fold']),
                matches_family_fold=fold==int(r['family_fold']),
                sample_in_corresponding_train_manifest=sample in train_ids[fold],
                sample_in_corresponding_test_manifest=sample in test_ids[fold],
                family_component_id=r['family_component_id'],
                same_family_other_train_records=len(training_family),
                same_biological_cluster_other_train_records=len(training_biology),
                checkpoint_sha256=row['checkpoint_sha256'],
                source_split_version=metadata.get('source_split_version',''),
                split_version=metadata.get('split_version',''),
                prediction_path=r['full_prediction_path'] if variant=='full' else r['intrachain_prediction_path']))
    assert len(details)==6126 and len(seen)==6126
    csv_write('training_exposure_own_oof_details.csv', details)
    evaluation = json.loads(text(export/'evaluation_audit.json'))
    strict = json.loads(text(core/'neural_exports/strict_external_evaluation/evaluation_audit.json'))
    per_variant = {}
    for variant, items in group_by(details, lambda x:x['variant']).items():
        per_variant[variant] = dict(records=len(items),matches_original_fold=sum(x['matches_original_fold'] for x in items),
            matches_family_fold=sum(x['matches_family_fold'] for x in items),
            sample_in_corresponding_train=sum(x['sample_in_corresponding_train_manifest'] for x in items),
            sample_in_corresponding_test=sum(x['sample_in_corresponding_test_manifest'] for x in items),
            records_with_same_family_in_train=sum(x['same_family_other_train_records']>0 for x in items),
            records_with_same_biological_cluster_in_train=sum(x['same_biological_cluster_other_train_records']>0 for x in items))
    return dict(manifest_rows=len(details),manifest_variants=dict(Counter(x['variant'] for x in details)),
        manifest_path=str(path),manifest_sha256=sha(path),per_variant=per_variant,
        families=len(family_rows),families_crossing_original_folds=sum(x['crosses_original_folds'] for x in family_rows),
        families_crossing_family_folds=sum(x['crosses_family_folds'] for x in family_rows),
        original_biological_clusters=len(biological_groups),
        biological_clusters_crossing_original_folds=sum(len({r['original_fold'] for r in items})>1 for items in biological_groups.values()),
        TCR_original_fold_records=dict(Counter(r['original_fold'] for r in records if r['dataset']=='itrace')),
        TCR_family_fold_records=dict(Counter(r['family_fold'] for r in records if r['dataset']=='itrace')),
        source_split_versions=dict(versions),evaluation_protocol=evaluation['protocol'],
        conclusion='Record/biological-cluster OOF; not sequence-family-held-out generalization. Family-cluster bootstrap describes dependence of scores and does not change training isolation.',
        strict_external_separate=dict(prespecified_records=strict['prespecified_records'],
            scored_records=strict['scored_records'],training_protocol=strict['external_prediction'],
            not_used_to_replace_current_3063_predictions=True))


def group_by(rows, key):
    groups = defaultdict(list)
    for row in rows:
        groups[key(row)].append(row)
    return groups


if __name__ == '__main__':
    main()
