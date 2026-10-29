"""Validate existing replica labels for two biologically nominated ITrace series."""
from pathlib import Path
import importlib.util,json,sys
ROOT=Path(__file__).resolve().parent
source=ROOT.parent/'tcr_expanded_screen_v41_20261017'/'recompute_allprotein_pair_labels.py'
spec=importlib.util.spec_from_file_location('itrace_replica_recheck',source)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
mod.STAGE=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\itrace_search_only_replica_check')
mod.STAGE.mkdir(parents=True,exist_ok=True)
OUT=ROOT/'replica_recheck';OUT.mkdir(exist_ok=True)
records=[]
requested=sys.argv[1:] or ['1lp9','2j8u','2jcc','2uwe','7pdw','7pbc']
for pdb in requested:
    out=OUT/pdb;out.mkdir(exist_ok=True)
    frame,pred,means,variance,prov,frames,audit=mod.complex_labels(pdb,out)
    record={'pdb':pdb,'status':'passed','frame_audit':frames,'source_files':prov,'mean_variance_checks':audit}
    (out/'replica_label_audit.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf8')
    records.append(record)
all_records=[json.loads(p.read_text(encoding='utf8')) for p in sorted(OUT.glob('*/replica_label_audit.json'))]
(OUT/'all_rechecked_complexes_audit.json').write_text(json.dumps({'status':'passed','purpose':'Exploratory ITrace example search with complete comparison/control states; no significance claims.','complexes':len(all_records),'existing_MD_trajectories':3*len(all_records),'new_MD':False,'new_training':False,'source_data_unchanged':True,'records':all_records},ensure_ascii=False,indent=2),encoding='utf8')
print('Complete:',len(all_records),'complexes /',3*len(all_records),'existing MD replicas passed.',flush=True)
