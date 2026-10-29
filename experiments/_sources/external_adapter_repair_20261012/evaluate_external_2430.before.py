"""ThermoPCD32 after complete-chain purge633; reuses family fold0 weights."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def statistics(target, prediction):
    result = dict(count=len(target), mae=float(np.mean(abs(target-prediction))), rmse=float(np.sqrt(np.mean((target-prediction)**2))), bias=float(np.mean(prediction-target)), spearman=np.nan, pearson=np.nan)
    if len(target) >= 4 and np.ptp(target) > 0 and np.ptp(prediction) > 0:
        result.update(spearman=float(spearmanr(target,prediction).statistic), pearson=float(pearsonr(target,prediction).statistic))
    return result


def main():
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
    from flexppi.data.dataset import StateDataset, collate_states
    from flexppi.data.structures import parse_pdb
    from flexppi.model.information_ablation import build_information_ablation
    from flexppi.trainer import move_to_device
    torch.set_num_threads(2)
    root, external = Path(__file__).resolve().parent, Path(args.external_root)
    freeze = json.loads((root / 'training_freeze.json').read_text())
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
    output = root/'external_2430_evaluation'
    output.mkdir(exist_ok=True)
    records, residue_tables, provenance = [], [], []
    valid_counts = {}
    for variant in ['full','intrachain']:
        checkpoint_path = root/f'runs/{variant}/fold_0/epoch_040.pt'
        checkpoint = torch.load(checkpoint_path, map_location='cpu', weights_only=False)
        config = checkpoint['config']
        assert checkpoint['epoch'] == 40 and config['seed'] == 7
        assert config['model']['information_ablation'] == variant
        model = build_information_ablation(ModelConfig(**config['model']), variant).to(args.device).eval()
        model.load_state_dict(checkpoint['model'], strict=True)
        folder = output/variant
        folder.mkdir(exist_ok=True)
        counts = np.zeros(2, dtype=np.int64)
        for index, record in enumerate(dataset.records):
            example = dataset[index]
            structure = parse_pdb(record.structure_path)
            batch = move_to_device(collate_states([example]), torch.device(args.device))
            with torch.inference_mode():
                prediction = model.forward_state(batch)
            mean = prediction['mean'][0].cpu().numpy()
            variance = prediction['log_variance'][0].exp().cpu().numpy()
            target = example['label_values'].numpy()
            valid = (example['label_mask'] & example['target_mask'][:,None] & prediction['mask'][0].cpu()[:,None]).numpy()
            assert np.isfinite(mean).all() and np.isfinite(variance).all()
            counts += valid.sum(0)
            np.savez_compressed(folder/(record.sample_id+'.npz'), mean=mean, model_variance=variance, target=target, valid=valid,
                amino_acid=np.asarray(list(structure.sequence)), chain=np.asarray(structure.chain_ids), author_residue_id=np.asarray(structure.author_residue_ids))
            row = dict(sample_id=record.sample_id, variant=variant, residue_index=np.arange(len(mean)), chain=structure.chain_ids, author_residue_id=structure.author_residue_ids)
            for m,metric in enumerate(['rmsf','chi_entropy']):
                records.append(dict(sample_id=record.sample_id, variant=variant, metric=metric, **statistics(target[valid[:,m],m],mean[valid[:,m],m])))
                row.update({f'target_{metric}':target[:,m], f'pred_{metric}':mean[:,m], f'valid_{metric}':valid[:,m]})
            residue_tables.append(pd.DataFrame(row))
        valid_counts[variant] = counts.tolist()
        provenance.append(dict(variant=variant, checkpoint_sha256=sha(checkpoint_path), training_manifest_sha256=sha(root/'splits/fold_0_train.jsonl')))
        del model
        if args.device.startswith('cuda'):
            torch.cuda.empty_cache()
        print('Completed',variant,32,counts.tolist(),flush=True)
    assert valid_counts['full'] == valid_counts['intrachain'] and valid_counts['full'][0] == 18748
    if args.static_output:
        assert args.static_reference
        static_output, reference = Path(args.static_output), Path(args.static_reference)
        static_audit = json.loads((static_output/'run_audit.json').read_text())
        assert static_audit['status'] == 'complete' and static_audit['training_records'] == 2430
        assert static_audit['training_manifest_sha256'] == sha(root/'splits/fold_0_train.jsonl')
        schema = json.loads((reference/'external_schema.json').read_text())
        with np.load(reference/'external_features.npz') as cache:
            target, valid = cache['target'], cache['valid']
        assert [r['sample_id'] for r in schema['records']] == [r.sample_id for r in dataset.records]
        for variant in ['static_sequence_complex','aa_only']:
            with np.load(static_output/(variant+'_external.npz')) as cache:
                prediction = cache['prediction']
            for record in schema['records']:
                ids = np.arange(record['start'],record['stop'])
                for m,metric in enumerate(['rmsf','chi_entropy']):
                    chosen = ids[valid[ids,m]]
                    records.append(dict(sample_id=record['sample_id'],variant=variant,metric=metric,**statistics(target[chosen,m],prediction[chosen,m])))
    per = pd.DataFrame(records)
    per.to_csv(output/'per_record_metrics.csv',index=False)
    per.groupby(['variant','metric'])[['mae','rmse','bias','spearman','pearson']].mean().to_csv(output/'record_macro_metrics.csv')
    pd.concat(residue_tables,ignore_index=True).to_csv(output/'per_residue.csv.gz',index=False,compression='gzip')
    audit = dict(status='complete', training_records=2430, purged_records=633, external_records=32, external_family_components=1,
        checkpoint_protocol='Reuse fresh family-fold0 full/intrachain models, seed7 epoch40, one model per arm',
        no_external_tuning=True, family_freeze_sha256=sha(root/'training_freeze.json'), external_manifest_sha256=sha(external/'ready_test.jsonl'),
        code_sha256=sha(__file__), checkpoints=provenance, valid_counts=valid_counts,
        uncertainty='One connected external family; no independent-family CI or p-value. Source MD protocol and released comparator pretraining overlap remain distinct limitations.')
    (output/'evaluation_audit.json').write_text(json.dumps(audit,indent=2))


if __name__ == '__main__':
    main()
