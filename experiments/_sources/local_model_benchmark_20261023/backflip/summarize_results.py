"""Summarize every available new result, with explicit coverage and failures."""
from pathlib import Path
import json
import argparse
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output-subdir',default='results')
    args=ap.parse_args()
    rows=[]
    failures=[]
    for folder in (ROOT/args.output_subdir).glob('*'):
        path=folder/'status.jsonl'
        if not path.exists():
            continue
        statuses={}
        for line in path.read_text(encoding='utf-8').splitlines():
            status=json.loads(line)
            statuses[status['sample_id'],status['context']]=status
        for status in statuses.values():
            if status['status']!='complete':
                failures.append(status)
                continue
            metric=status['metrics']
            rows.append(dict(tag=status['tag'],context=status['context'],role=status['role'],
                dataset=status['dataset'],sample_id=status['sample_id'],
                predicted=status['predicted'],inference_seconds=status['seconds'],
                reused_fresh_context=status.get('reused_fresh_context',''),**metric))
    if not rows:
        print('No completed inference available')
        return
    table=pd.DataFrame(rows)
    table.to_csv(ROOT/(args.output_subdir+'_per_record_metrics.csv'),index=False)
    aggregates=[]
    for keys,group in table.groupby(['tag','context','role']):
        for dataset,pick in [('ALL',group)]+list(group.groupby('dataset')):
            aggregates.append(dict(zip(['tag','context','role'],keys),dataset=dataset,
                completed_records=len(pick),valid_residues=int(pick['count'].sum()),
                actually_executed_records=int((pick.reused_fresh_context=='').sum()),
                shared_current_single_chain_forwards=int((pick.reused_fresh_context!='').sum()),
                inference_seconds=float(pick.inference_seconds.sum()),
                **{name:float(pick[name].mean()) for name in ['mae','rmse','bias','spearman','pearson']}))
    summary=pd.DataFrame(aggregates)
    summary.to_csv(ROOT/(args.output_subdir+'_macro_metrics.csv'),index=False)
    (ROOT/(args.output_subdir+'_failures.json')).write_text(json.dumps(failures,indent=2),encoding='utf-8')
    print(summary.to_string(index=False))
    print('FAILURES',len(failures),flush=True)


if __name__=='__main__':
    main()
