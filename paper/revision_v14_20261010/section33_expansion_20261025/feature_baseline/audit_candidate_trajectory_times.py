"""Read actual local XTC timestamps, never infer duration from metadata labels."""
import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import time
import numpy as np
import pandas as pd
from mdtraj.formats import XTCTrajectoryFile

BASE=Path(__file__).resolve().parent
OUT=BASE/'case_candidates'
RAW=Path('D:/data/dynarepo_non_nucleic_full')


def audit(sid):
    path=RAW/sid/'trajectory.xtc'
    row=dict(sample_id=sid,trajectory_path=str(path),topology_path=str(RAW/sid/'structure.pdb'))
    try:
        with XTCTrajectoryFile(str(path),'r') as xtc:
            n=len(xtc)
            first=xtc.read(n_frames=min(2,n),atom_indices=np.array([0]))
            xtc.seek(n-1)
            last=xtc.read(n_frames=1,atom_indices=np.array([0]))
        row.update(status='ok',n_frames=n,first_time_ps=float(first[1][0]),
                   second_time_ps=float(first[1][1]) if n>1 else np.nan,last_time_ps=float(last[1][0]),
                   first_step=int(first[2][0]),second_step=int(first[2][1]) if n>1 else -1,last_step=int(last[2][0]),
                   dt_ps=float(first[1][1]-first[1][0]) if n>1 else np.nan,
                   duration_ns=float((last[1][0]-first[1][0])/1000),
                   bytes=path.stat().st_size)
        row['at_least_50ns']=row['duration_ns']>=49.9
        row['at_least_100ns']=row['duration_ns']>=99.9
    except Exception as error:
        row.update(status='error',error=repr(error),at_least_50ns=False,at_least_100ns=False)
    return row


def main():
    ranking=pd.read_csv(OUT/'non_tcr_candidate_chain_ranking.csv')
    previous=pd.read_csv(BASE.parent.parent.parent.parent/'experiments/non_tcr_case_screen_20261025/all_616_record_screen.csv') if False else pd.read_csv('E:/Project/FlexPPI/experiments/non_tcr_case_screen_20261025/all_616_record_screen.csv')
    # Explicitly exclude genuine TCR descriptions from the legacy non-ITrace screen.
    previous=previous[~previous.description.str.contains(r'T[ -]?cell receptor|TCR',case=False,na=False)]
    ordered=list(dict.fromkeys(list(ranking.sample_id)+list(previous.sample_id)))
    file=OUT/'actual_local_trajectory_time_audit.csv'
    cache={}
    if file.exists():
        cache={r['sample_id']:r for r in pd.read_csv(file).to_dict('records')}
    results=[];tick=time.monotonic()
    with ProcessPoolExecutor(max_workers=3) as pool:
        for start,stop in [(0,50),(50,len(ordered))]:
            batch=ordered[start:stop]
            todo=[sid for sid in batch if sid not in cache]
            results.extend(cache[sid] for sid in batch if sid in cache)
            futuremap={pool.submit(audit,sid):sid for sid in todo}
            for future in as_completed(futuremap):
                result=future.result();results.append(result)
                if len(results)%25==0:
                    print('TIME_AUDIT',len(results),len(ordered),round(time.monotonic()-tick,1),flush=True)
            pd.DataFrame(results).to_csv(file,index=False)
            annotated=ranking.merge(pd.DataFrame(results),on='sample_id',how='left',validate='many_to_one')
            annotated.to_csv(OUT/'candidate_ranking_with_actual_duration.csv',index=False)
            eligible=annotated[annotated.at_least_50ns.eq(True)].copy()
            eligible.to_csv(OUT/'candidates_minimum_50ns.csv',index=False)
            print('DURATION_BATCH_READY',stop,len(results),'eligible_chains',len(eligible),flush=True)
    summary=dict(status='complete',records=len(results),passed50ns=sum(r.get('at_least_50ns',False) for r in results),
                 passed100ns=sum(r.get('at_least_100ns',False) for r in results),
                 timestamp_protocol='mdtraj XTCTrajectoryFile read-only first,second,last times and steps; duration is last minus first, ps/1000; no metadata rescaling',
                 caveat='A declared 100ns source with reindexed 0..999ps timestamps is flagged, not silently promoted to 100ns. Independent repeats are not inferred.',
                 code_modifies_original_files=False)
    (OUT/'trajectory_time_audit_contract.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
    if len(eligible):
        print(eligible.head(20)[['sample_id','pdb_id','chain','duration_ns','n_disjoint_superior_pairs','full_mae','aa_sasa_contacts_mae','static47_sasa_mae','full_typed_rho','static47_sasa_typed_rho','interface_mae_gain_vs_intra','family_clean_original_fold','selected_residue_pairs']].to_string(index=False))


if __name__=='__main__':main()
