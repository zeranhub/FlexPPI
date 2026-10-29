"""Freeze five folds whose complete chain-family components never cross folds.

Only metadata and sequence-search edges determine folds. Labels and predictions
are not read. Existing external sequences participate in the dependence graph
but are never assigned to the training/test folds.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


class Components:
    def __init__(self, items):
        self.parent = {item: item for item in items}

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def join(self, items):
        items = sorted(items)
        if items:
            first = self.find(items[0])
            for item in items[1:]:
                self.parent[self.find(item)] = first


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference-root', required=True)
    parser.add_argument('--training-root', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--supplemental-similarity', action='append', default=[])
    args = parser.parse_args()
    reference, output = Path(args.reference_root), Path(args.output)
    if (output / 'training_freeze.json').exists():
        raise RuntimeError('Frozen output already exists; use a new directory')
    chains_path = reference / 'sequence_audit/chain_inventory.json'
    similarity_path = reference / 'sequence_audit/all_chain_similarity.tsv'
    original_path = reference / 'training_payload/splits/all_ready.jsonl'
    fasta_path = reference / 'sequence_audit/long_unique.fasta'
    chains = json.loads(chains_path.read_text(encoding='utf-8'))
    records = [json.loads(line) for line in original_path.read_text(encoding='utf-8').splitlines()]
    assert len(records) == 3063 and len({r['sample_id'] for r in records}) == 3063
    training = {r['sample_id'] for r in records}
    assert {c['sample_id'] for c in chains if c['role'] == 'training'} == training
    assert len(chains) == 5068
    sequences, biological, by_record = defaultdict(set), defaultdict(set), defaultdict(list)
    sequence_text = {}
    for chain in chains:
        assert len(chain['sequence']) == chain['length'] and chain['length'] > 0
        sequence_id = chain['sequence_id']
        assert sequence_id == 'seq_' + hashlib.sha256(chain['sequence'].encode()).hexdigest()
        if sequence_id in sequence_text:
            assert sequence_text[sequence_id] == chain['sequence']
        sequence_text[sequence_id] = chain['sequence']
        sequences[sequence_id].add(chain['sample_id'])
        biological[(chain['dataset'], chain['biological_system'])].add(chain['sample_id'])
        by_record[chain['sample_id']].append(chain)
    graph = Components(by_record)
    for group in list(sequences.values()) + list(biological.values()):
        graph.join(group)
    fasta_sequences = {}
    fasta_id = None
    for line in fasta_path.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            fasta_id = line[1:].split()[0]
            assert fasta_id not in fasta_sequences
            fasta_sequences[fasta_id] = ''
        else:
            assert fasta_id is not None
            fasta_sequences[fasta_id] += line.strip()
    long_sequences = {key for key, text in sequence_text.items() if len(text) >= 20}
    assert set(fasta_sequences) == long_sequences
    assert all(fasta_sequences[key] == sequence_text[key] for key in long_sequences)
    accepted, self_queries = set(), set()
    search_paths = [similarity_path] + [Path(p) for p in args.supplemental_similarity]
    for search_path in search_paths:
      with search_path.open(encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip().split('\t')
            assert len(fields) == 14
            query, target = fields[:2]
            assert query in sequences and target in sequences
            assert len(sequence_text[query]) >= 20 and len(sequence_text[target]) >= 20
            assert int(fields[6]) == len(sequence_text[query])
            assert int(fields[9]) == len(sequence_text[target])
            if query == target:
                self_queries.add(query)
            qcoverage = (int(fields[5]) - int(fields[4]) + 1) / int(fields[6])
            tcoverage = (int(fields[8]) - int(fields[7]) + 1) / int(fields[9])
            if float(fields[2]) < .30 or min(qcoverage, tcoverage) < .70:
                continue
            if query != target:
                accepted.add(tuple(sorted((query, target))))
            graph.join(sequences[query] | sequences[target])
    assert self_queries == long_sequences, 'A long-chain query is absent from the search'
    groups = defaultdict(set)
    for record in by_record:
        groups[graph.find(record)].add(record)
    components = []
    lookup = {r['sample_id']: r for r in records}
    for members in groups.values():
        members_train = members & training
        if not members_train:
            continue
        component_id = 'family_' + hashlib.sha256('\n'.join(sorted(members)).encode()).hexdigest()[:16]
        components.append((component_id, members_train, members))
    # Deterministic largest-first allocation prioritizes record balance; source
    # balance breaks ties. Shared partners may conservatively join large groups.
    components.sort(key=lambda item: (-len(item[1]), item[0]))
    fold_counts, fold_source = [0] * 5, [Counter() for _ in range(5)]
    assignment, component_for_record, component_rows = {}, {}, []
    for component_id, members, all_members in components:
        sources = Counter(lookup[m]['dataset'] for m in members)
        fold = min(range(5), key=lambda f: (fold_counts[f], sum(fold_source[f][s] for s in sources), f))
        fold_counts[fold] += len(members)
        fold_source[fold].update(sources)
        for member in members:
            assignment[member] = fold
            component_for_record[member] = component_id
        component_rows.append(dict(component_id=component_id, fold=fold, training_records=len(members),
            external_records=len(all_members - training), sources=dict(sources), members=sorted(all_members)))
    assert set(assignment) == training
    def assert_no_cross(group):
        assert len({assignment[s] for s in group if s in assignment}) <= 1
    for group in list(sequences.values()) + list(biological.values()):
        assert_no_cross(group)
    for query, target in accepted:
        assert_no_cross(sequences[query] | sequences[target])
    output.mkdir(parents=True, exist_ok=True)
    (output / 'splits').mkdir(exist_ok=True)
    (output / 'configs').mkdir(exist_ok=True)
    rewritten = []
    for record in records:
        record = json.loads(json.dumps(record))
        sample = record['sample_id']
        for field, directory in [('structure_path', 'structures'), ('label_path', 'labels')]:
            filename = record[field].replace('\\', '/').split('/')[-1]
            record[field] = args.training_root.rstrip('/') + '/dataset/' + directory + '/' + filename
        record['metadata']['original_combination_fold'] = record['metadata']['fold']
        record['metadata']['fold'] = assignment[sample]
        record['metadata']['family_component_id'] = component_for_record[sample]
        record['metadata']['split_version'] = 'all_chain_family_disjoint_v14_20261010'
        rewritten.append(record)
    def write_jsonl(path, rows):
        path.write_text(''.join(json.dumps(row, ensure_ascii=False) + '\n' for row in rows), encoding='utf-8')
    write_jsonl(output / 'splits/all_ready.jsonl', rewritten)
    for fold in range(5):
        for role in ['train', 'test']:
            selected = [r for r in rewritten if (r['metadata']['fold'] == fold) == (role == 'test')]
            write_jsonl(output / f'splits/fold_{fold}_{role}.jsonl', selected)
        for variant in ['full', 'intrachain']:
            template = (reference / f'training_payload/configs/{variant}_fold0.yaml').read_text(encoding='utf-8')
            assert 'epochs: 40' in template and 'seed: 7' in template and 'init_checkpoint: null' in template
            template = template.replace(f'output_dir: runs/{variant}/fold_0', f'output_dir: runs/{variant}/fold_{fold}')
            template = template.replace('fold_0_train.jsonl', f'fold_{fold}_train.jsonl').replace('fold_0_test.jsonl', f'fold_{fold}_test.jsonl')
            (output / f'configs/{variant}_fold{fold}.yaml').write_text(template, encoding='utf-8')
    with (output / 'record_family_folds.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sample_id', 'dataset', 'biological_system', 'original_fold', 'fold', 'component_id', 'chains', 'short_chains'])
        writer.writeheader()
        for record in rewritten:
            sample = record['sample_id']
            writer.writerow(dict(sample_id=sample, dataset=record['dataset'], biological_system=record['metadata']['biological_system'],
                original_fold=record['metadata']['original_combination_fold'], fold=assignment[sample], component_id=component_for_record[sample],
                chains=len(by_record[sample]), short_chains=sum(c['length'] < 20 for c in by_record[sample])))
    (output / 'family_components.json').write_text(json.dumps(component_rows, indent=2), encoding='utf-8')
    report = dict(status='frozen', records=3063, training_chains=sum(c['role'] == 'training' for c in chains),
        all_graph_records=len(by_record), all_graph_chains=len(chains), training_family_components=len(components),
        accepted_nonself_sequence_edges=len(accepted), long_unique_sequences=len(long_sequences), short_unique_sequences=sum(len(s) < 20 for s in sequence_text.values()),
        short_training_chains=sum(c['length'] < 20 and c['role'] == 'training' for c in chains),
        largest_training_component=len(components[0][1]), largest_all_graph_component=max(len(m) for m in groups.values()),
        fold_records=fold_counts, fold_sources=[dict(c) for c in fold_source],
        fold_multichain_records=[sum(assignment[s] == f and len(by_record[s]) > 1 for s in training) for f in range(5)],
        all_members_exact_and_long_edge_no_cross=True, biological_system_no_cross=True,
        rule='Any exact sequence including chains shorter than20, any long-chain MMseqs identity>=0.30 and integer-span coverage>=0.70 in both directions, or same source-specific biological system; transitive complete-complex components.',
        graph_external_policy='32 external records participate only in family connectivity; none enter the3063 training/test records.',
        allocation='Deterministic largest-first record-count balancing; source counts break ties; no labels or predictions used.',
        seed=7, epochs=40, trained_from_scratch=True, fits=10, validation_selection=False,
        sequence_search_sha256=sha(similarity_path), supplemental_search_sha256={str(p): sha(p) for p in search_paths[1:]},
        long_fasta_sha256=sha(fasta_path), long_fasta_complete_mapping_verified=True,
        chain_inventory_sha256=sha(chains_path), original_manifest_sha256=sha(original_path),
        limitations='MMseqs search sensitivity and identity rounded to3 decimals follow the existing protocol. Shared short peptides/partners can conservatively connect unrelated domains; group size and source concentration are disclosed.',
        files={str(p.relative_to(output)).replace('\\', '/'): sha(p) for p in sorted(list(output.glob('splits/*.jsonl')) + list(output.glob('configs/*.yaml')))})
    (output / 'training_freeze.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'files'}, indent=2))


if __name__ == '__main__':
    main()
