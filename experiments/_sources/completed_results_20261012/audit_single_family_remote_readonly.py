"""Read-only CPU comparison of existing target caches; no model imports or writes.

Run remotely through standard input and save stdout locally. Existing inference
certificates are never regenerated. Full caches already exist on the host.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np

ROOT=Path('/root/flexppi_completion_20261010')
OLD=Path('/root/flexppi_revision_20261002/corrected_training_v1/evaluation')
SINGLE=ROOT/'single_target'


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def require(condition,message):
    if not condition:
        raise ValueError(message)


def jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def cache(path,row):
    local=path/row['variant']/(row['sample_id']+'.npz')
    return local if local.is_file() else Path(row['cache_path'])


def main():
    started=time.monotonic()
    dual=jsonl(OLD/'prediction_manifest.jsonl')
    dual={row['sample_id']:row for row in dual if row['variant']=='full'}
    single=jsonl(SINGLE/'evaluation/prediction_manifest.jsonl')
    require(len(dual)==3063 and len(single)==6126,'Wrong manifest sizes')
    require(len({(row['variant'],row['sample_id']) for row in single})==6126,'Duplicate single-target predictions')
    report=json.loads((SINGLE/'analysis/analysis_audit.json').read_text())
    require(report['status']=='complete','Analysis not complete')
    require(sha(OLD/'prediction_manifest.jsonl')==report['inputs']['dual_manifest_sha256'],'Dual manifest changed since analysis')
    require(sha(SINGLE/'evaluation/prediction_manifest.jsonl')==report['inputs']['single_manifest_sha256'],'Single manifest changed since analysis')
    counts={'rmsf_only':0,'chi_entropy_only':0}
    pairs={'rmsf_only':0,'chi_entropy_only':0}
    labels=['rmsf','chi_entropy','sasa']
    for row in single:
        arm=row['variant']
        require(arm in counts,'Unexpected single-target arm')
        reference=dual[row['sample_id']]
        for field in ['fold','dataset','input_sha256','label_sha256']:
            require(row[field]==reference[field],'Pair differs in '+field+': '+row['sample_id'])
        metric='rmsf' if arm=='rmsf_only' else 'chi_entropy'
        with np.load(cache(SINGLE/'evaluation',row),allow_pickle=False) as a,np.load(cache(OLD,reference),allow_pickle=False) as b:
            ai=list(a['label_names']) if 'label_names' in a.files else labels
            bi=list(b['label_names']) if 'label_names' in b.files else labels
            acol,bcol=ai.index(metric),bi.index(metric)
            amask,bmask=a['valid'][:,acol],b['valid'][:,bcol]
            require(amask.dtype==bmask.dtype==np.dtype(bool),'Nonboolean frozen masks')
            require(np.array_equal(amask,bmask),'Actual saved masks differ: '+row['sample_id'])
            require(np.array_equal(a['target'][:,acol][amask],b['target'][:,bcol][bmask]),'Actual valid targets differ: '+row['sample_id'])
            anames=list(a['metric_names']) if 'metric_names' in a.files else [metric]
            bnames=list(b['metric_names']) if 'metric_names' in b.files else ['rmsf','chi_entropy']
            apred,bpred=a['mean'][:,anames.index(metric)],b['mean'][:,bnames.index(metric)]
            require(np.isfinite(apred[amask]).all() and np.isfinite(bpred[bmask]).all(),'Nonfinite existing prediction')
            require(len(apred)==len(bpred)==len(amask),'Residue shape mismatch')
            # Residue identity is preserved by the same frozen input hash. If
            # both caches save explicit identity arrays, check those too.
            for field in ['chain','author_residue_id','amino_acid']:
                if field in a.files and field in b.files:
                    require(np.array_equal(a[field],b[field]),'Saved residue identity differs: '+row['sample_id'])
            counts[arm]+=int(amask.sum())
            pairs[arm]+=1
    require(counts=={'rmsf_only':961321,'chi_entropy_only':824988},'Frozen effective counts changed')
    require(pairs=={'rmsf_only':3063,'chi_entropy_only':3063},'Incomplete pair coverage')
    result=dict(status='passed_readonly_comparison',no_training=True,no_inference=True,no_gpu=True,no_existing_audit_modified=True,
        records_per_arm=pairs,effective_residues_per_arm=counts,
        checks=['Exact corresponding heldout fold and source','Identical frozen input and label hashes',
            'Every saved effective mask equal','Every valid target exactly equal','Finite same-shape predictions','Saved residue identities equal where present'],
        manifest_hashes=report['inputs'],analysis_audit_sha256=sha(SINGLE/'analysis/analysis_audit.json'),
        seconds=time.monotonic()-started,checked_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
