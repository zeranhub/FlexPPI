"""ThermoPCD32 after complete-chain purge633; reuses family fold0 weights."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

EVALUATION_METRICS = ('rmsf', 'chi_entropy')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def statistics(target, prediction):
    result = dict(count=len(target), mae=float(np.mean(abs(target-prediction))), rmse=float(np.sqrt(np.mean((target-prediction)**2))), bias=float(np.mean(prediction-target)), spearman=np.nan, pearson=np.nan)
    if len(target) >= 4 and np.ptp(target) > 0 and np.ptp(prediction) > 0:
        result.update(spearman=float(spearmanr(target,prediction).statistic), pearson=float(pearsonr(target,prediction).statistic))
    return result


def main():
    started = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-root', required=True)
    parser.add_argument('--external-root', required=True)
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--static-output', help='Optional completed2430 static directory')
    parser.add_argument('--static-reference', help='Reference feature cache used for static fitting')
    args = parser.parse_args()
    sys.path.insert(0, args.source_root)
    import torch
    from flexppi.config import ModelConfig
    from flexppi.constants import METRIC_NAMES
    from flexppi.data.dataset import StateDataset, collate_states
    from flexppi.data.structures import parse_pdb
    from flexppi.model.information_ablation import build_information_ablation
    from flexppi.trainer import move_to_device
    torch.set_num_threads(2)
    root, external = Path(__file__).resolve().parent, Path(args.external_root)
    freeze = json.loads((root / 'training_freeze.json').read_text())
    for name, expected_hash in freeze['files'].items():
        assert sha(root / name) == expected_hash, 'Frozen training configuration/manifest changed: ' + name
    assert tuple(METRIC_NAMES) == ('rmsf', 'chi_entropy', 'sasa')
    label_indices = [tuple(METRIC_NAMES).index(name) for name in EVALUATION_METRICS]
    assert freeze['largest_training_component'] == 633 and freeze['fold_records'][0] == 633
    allowed = {json.loads(line)['sample_id'] for line in (root/'splits/fold_0_train.jsonl').read_text().splitlines()}
    assert len(allowed) == 2430
    assert sha(root/'splits/fold_0_train.jsonl') == freeze['files']['splits/fold_0_train.jsonl']
    expected = [json.loads(line) for line in (external/'ready_test.jsonl').read_text().splitlines()]
    assert len(expected) == 32
    component = json.loads((root/'family_components.json').read_text())[0]
    assert component['training_records'] == 633 and component['external_records'] == 32
    assert {r['sample_id'] for r in expected} <= set(component['members'])
    assert not allowed & set(component['members'])
    dataset = StateDataset(str(external/'ready_test.jsonl'), max_residues=1024)
    assert len(dataset) == 32
    assert [row.sample_id for row in dataset.records] == [row['sample_id'] for row in expected]
    expected_counts = np.zeros(len(EVALUATION_METRICS), dtype=np.int64)
    input_provenance = []
    for index, record in enumerate(dataset.records):
        example = dataset[index]
        assert tuple(example['metric_names']) == tuple(METRIC_NAMES)
        assert example['label_values'].shape == example['label_mask'].shape
        assert example['label_values'].shape[1] == len(METRIC_NAMES)
        valid = (example['label_mask'][:,label_indices] & example['target_mask'][:,None]).numpy()
        expected_counts += valid.sum(0)
        input_provenance.append(dict(sample_id=record.sample_id, structure_sha256=sha(record.structure_path),
            label_sha256=sha(record.label_path), residues=len(valid), expected_valid_counts=valid.sum(0).tolist()))
    assert expected_counts.tolist() == [18748,16265], 'Frozen32-record two-target mask totals changed'
    output = root/'external_2430_evaluation'
    output.mkdir(exist_ok=True)
    records, residue_tables, provenance, paired_truth = [], [], [], {}
    valid_counts = {}
    for variant in ['full','intrachain']:
        checkpoint_path = root/f'runs/{variant}/fold_0/epoch_040.pt'
        checkpoint = torch.load(checkpoint_path, map_location='cpu', weights_only=False)
        config = checkpoint['config']
        assert checkpoint['epoch'] == 40 and config['seed'] == 7
        assert config['model']['information_ablation'] == variant
        assert tuple(config['model']['metrics']) == EVALUATION_METRICS, 'Checkpoint output channels differ from evaluation targets'
        assert config['init_checkpoint'] is None and config['train']['epochs'] == 40
        training_manifest = Path(config['data']['state_manifest'])
        if not training_manifest.is_absolute():
            training_manifest = root / training_manifest
        assert training_manifest.resolve() == (root/'splits/fold_0_train.jsonl').resolve()
        model = build_information_ablation(ModelConfig(**config['model']), variant).to(args.device).eval()
        model.load_state_dict(checkpoint['model'], strict=True)
        folder = output/variant
        folder.mkdir(exist_ok=True)
        counts = np.zeros(len(EVALUATION_METRICS), dtype=np.int64)
        for index, record in enumerate(dataset.records):
            example = dataset[index]
            structure = parse_pdb(record.structure_path)
            batch = move_to_device(collate_states([example]), torch.device(args.device))
            with torch.inference_mode():
                prediction = model.forward_state(batch)
            mean = prediction['mean'][0].cpu().numpy()
            variance = prediction['log_variance'][0].exp().cpu().numpy()
            assert tuple(example['metric_names']) == tuple(METRIC_NAMES)
            source_indices = [tuple(example['metric_names']).index(name) for name in config['model']['metrics']]
            assert source_indices == label_indices
            target = example['label_values'][:,source_indices].numpy()
            valid = (example['label_mask'][:,source_indices] & example['target_mask'][:,None] & prediction['mask'][0].cpu()[:,None]).numpy()
            assert mean.shape == variance.shape == target.shape == valid.shape == (len(structure.sequence),len(EVALUATION_METRICS))
            assert np.isfinite(mean).all() and np.isfinite(variance).all()
            assert np.isfinite(target[valid]).all()
            if variant == 'full':
                paired_truth[record.sample_id] = (target.copy(), valid.copy())
            else:
                old_target, old_valid = paired_truth[record.sample_id]
                assert np.array_equal(valid, old_valid) and np.array_equal(target, old_target, equal_nan=True)
            counts += valid.sum(0)
            np.savez_compressed(folder/(record.sample_id+'.npz'), mean=mean, model_variance=variance, target=target, valid=valid,
                output_metric_names=np.asarray(EVALUATION_METRICS), label_metric_names=np.asarray(METRIC_NAMES), label_source_indices=np.asarray(source_indices),
                amino_acid=np.asarray(list(structure.sequence)), chain=np.asarray(structure.chain_ids), author_residue_id=np.asarray(structure.author_residue_ids))
            row = dict(sample_id=record.sample_id, variant=variant, residue_index=np.arange(len(mean)), chain=structure.chain_ids, author_residue_id=structure.author_residue_ids)
            for m,metric in enumerate(EVALUATION_METRICS):
                records.append(dict(sample_id=record.sample_id, variant=variant, metric=metric, **statistics(target[valid[:,m],m],mean[valid[:,m],m])))
                row.update({f'target_{metric}':target[:,m], f'pred_{metric}':mean[:,m], f'valid_{metric}':valid[:,m]})
            residue_tables.append(pd.DataFrame(row))
        valid_counts[variant] = counts.tolist()
        assert counts.tolist() == expected_counts.tolist()
        provenance.append(dict(variant=variant, checkpoint_sha256=sha(checkpoint_path), training_manifest_sha256=sha(root/'splits/fold_0_train.jsonl')))
        del model
        if args.device.startswith('cuda'):
            torch.cuda.empty_cache()
        print('Completed',variant,32,counts.tolist(),flush=True)
    assert valid_counts['full'] == valid_counts['intrachain'] == [18748,16265]
    if args.static_output:
        assert args.static_reference
        static_output, reference = Path(args.static_output), Path(args.static_reference)
        static_audit = json.loads((static_output/'run_audit.json').read_text())
        assert static_audit['status'] == 'complete' and static_audit['training_records'] == 2430
        assert static_audit['training_manifest_sha256'] == sha(root/'splits/fold_0_train.jsonl')
        schema = json.loads((reference/'external_schema.json').read_text())
        for name in ['external_features.npz','external_schema.json']:
            assert sha(reference/name) == static_audit['reference_hashes'][name], 'Frozen static external reference changed'
        with np.load(reference/'external_features.npz') as cache:
            target, valid = cache['target'], cache['valid']
        assert target.shape == valid.shape == (18748,len(EVALUATION_METRICS))
        assert valid.sum(0).tolist() == expected_counts.tolist() == static_audit['valid_external']
        assert [r['sample_id'] for r in schema['records']] == [r.sample_id for r in dataset.records]
        for record in schema['records']:
            ids = np.arange(record['start'],record['stop'])
            neural_target, neural_valid = paired_truth[record['sample_id']]
            assert np.array_equal(valid[ids], neural_valid), 'Neural/static valid positions differ'
            assert np.array_equal(target[ids], neural_target, equal_nan=True), 'Neural/static label values differ'
            input_row = next(row for row in input_provenance if row['sample_id'] == record['sample_id'])
            assert record['input_sha256'] == input_row['structure_sha256'] and record['label_sha256'] == input_row['label_sha256']
        for variant in ['static_sequence_complex','aa_only']:
            with np.load(static_output/(variant+'_external.npz')) as cache:
                prediction = cache['prediction']
            assert prediction.shape == target.shape and np.isfinite(prediction[valid]).all()
            for record in schema['records']:
                ids = np.arange(record['start'],record['stop'])
                for m,metric in enumerate(EVALUATION_METRICS):
                    chosen = ids[valid[ids,m]]
                    records.append(dict(sample_id=record['sample_id'],variant=variant,metric=metric,**statistics(target[chosen,m],prediction[chosen,m])))
    per = pd.DataFrame(records)
    expected_variants = ['full','intrachain'] + (['static_sequence_complex','aa_only'] if args.static_output else [])
    assert len(per) == len(expected_variants)*32*len(EVALUATION_METRICS)
    assert not per.duplicated(['sample_id','variant','metric']).any()
    for variant in expected_variants:
        assert set(per[per.variant == variant].sample_id) == {row.sample_id for row in dataset.records}
        assert [int(per[(per.variant == variant)&(per.metric == metric)]['count'].sum()) for metric in EVALUATION_METRICS] == expected_counts.tolist()
    per.to_csv(output/'per_record_metrics.csv',index=False)
    per.groupby(['variant','metric'])[['mae','rmse','bias','spearman','pearson']].mean().to_csv(output/'record_macro_metrics.csv')
    pd.concat(residue_tables,ignore_index=True).to_csv(output/'per_residue.csv.gz',index=False,compression='gzip')
    audit = dict(status='complete', training_records=2430, purged_records=633, external_records=32, external_family_components=1,
        checkpoint_protocol='Reuse fresh family-fold0 full/intrachain models, seed7 epoch40, one model per arm',
        no_external_tuning=True, family_freeze_sha256=sha(root/'training_freeze.json'), external_manifest_sha256=sha(external/'ready_test.jsonl'),
        code_sha256=sha(__file__), checkpoints=provenance, valid_counts=valid_counts,
        evaluation_metric_names=list(EVALUATION_METRICS), source_label_metric_names=list(METRIC_NAMES),
        label_source_indices=label_indices, expected_valid_counts=expected_counts.tolist(), record_order_exact=True,
        static_and_neural_labels_masks_identical=bool(args.static_output), inputs=input_provenance,
        completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), seconds=time.perf_counter()-started,
        uncertainty='One connected external family; no independent-family CI or p-value. Source MD protocol and released comparator pretraining overlap remain distinct limitations.')
    (output/'evaluation_audit.json').write_text(json.dumps(audit,indent=2))
    print('Complete external evaluation32 records, targets',list(EVALUATION_METRICS),'valid',expected_counts.tolist(),'seconds',audit['seconds'],flush=True)


if __name__ == '__main__':
    main()
