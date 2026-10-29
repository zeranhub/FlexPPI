"""Concrete ranked candidate lists, using the first prespecified sampling seed."""
from pathlib import Path
import sys,json,csv
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'evaluation'))
import evaluate_multimode as ev

def main():
    panels=ev.read_json(ROOT/'features/primary__contracts.json')['panels'];records=[];summaries=[]
    for panel in panels:
        ids,y,hits,kernels,settings,zero,bindings,hitinfo=ev.load_panel(panel)
        n=len(ids);initial=np.sort(np.random.default_rng(20261022).choice(n,size=18,replace=False)).tolist()
        hidden=sorted(set(range(n))-set(initial));lp={r['candidate_id']:r for r in ev.read_csv(panel['labels_csv'])}
        actual=ev.read_csv(ROOT/'evaluation/results_primary'/panel['panel_id']/'selected_candidates.csv')
        for method in list(ev.MODEL_BLOCKS)+['published_HERMES_fixed','published_ProteinMPNN']:
            score=ev.ridge_predict(kernels[method],initial,y[initial],hidden) if method in ev.MODEL_BLOCKS else zero[hidden,0 if method=='published_HERMES_fixed' else 1]
            indexed={i:float(p) for i,p in zip(hidden,score)}
            order=sorted(hidden,key=lambda i:(-indexed[i],ids[i]))
            archived=next(r for r in actual if r['mode']=='mode2' and r['strategy']==method and r['seed']=='20261022')
            assert [ids[i] for i in order[:30]]==json.loads(archived['selected_candidate_ids'])
            for rank,i in enumerate(order,1):
                row=lp[ids[i]]
                records.append(dict(panel_id=panel['panel_id'],strategy=method,seed=20261022,rank=rank,
                    candidate_id=ids[i],mutation=row['mutation'],sequence=row['sequence'],selected_in_30_assay_budget=rank<=30,
                    auxiliary_predicted_functional_score=indexed[i],source_measured_highest_dose_response=float(y[i]),
                    operational_response_retained=bool(hits[i]),source_outcome_used_for_ranking=False,
                    initial_candidate_ids=json.dumps([ids[j] for j in initial]),initial_n=18))
            summaries.append(dict(panel_id=panel['panel_id'],strategy=method,new_hits_30=int(hits[order[:30]].sum())))
    ev.dump_csv(ROOT/'evaluation/candidate_rankings_seed20261022.csv',records)
    ev.dump_csv(ROOT/'evaluation/ranking_example_yields_seed20261022.csv',summaries)
    print(json.dumps([r for r in summaries if r['panel_id']=='MSKCC_TCR6']))

if __name__=='__main__':main()
