"""Inventory public/local paired MD without reusing legacy response-model scores."""
import json,hashlib
from pathlib import Path
import pandas as pd

ROOT=Path(r'E:\Project\FlexPPI');OUT=Path(__file__).resolve().parent
core=ROOT/'experiments/core_corrected_pipeline'
def read(path):return [json.loads(x) for x in path.read_text().splitlines() if x.strip()]
ready={r['sample_id']:r for r in read(core/'training_payload/splits/all_ready.jsonl')}
state=read(ROOT/'artifacts/zenodo_6638504/data/zenodo_6638504_states.jsonl')
zenodo=[]
for r in state:
    trajectories=r.get('trajectory_paths',[])
    zenodo.append(dict(sample_id=r['sample_id'],pdb_id=r['metadata']['pdb_id'],state=r['source_group'],
                       structure_path=r['structure_path'],structure_exists=Path(r['structure_path']).exists(),
                       trajectory_paths=';'.join(trajectories),local_trajectory_count=sum(Path(p).exists() for p in trajectories),
                       frames=r['metadata']['frames'],production_ns=r['metadata']['production_ns'],
                       in_corrected3063=r['sample_id'] in ready,
                       current_full_prediction_present=(core/'neural_exports/evaluation/full'/(r['sample_id']+'.npz')).exists(),
                       current_chemical_label_provenance_verified=False))
pd.DataFrame(zenodo).to_csv(OUT/'zenodo24_state_inventory.csv',index=False)
summaries=[dict(dataset='DynaRepo_test53',available_states=73,true_pairs=37,independent_systems=36,
                local_independent_repeats_per_state=1,window_unit='five100ns blocks withinone500ns trajectory',
                current_corrected_oof_pairs=36,current_corrected_oof_systems=35,
                eligibility='36same-heldout-fold pairs recomputed after canonical sequence remapping; oneholoA00BV absentcorrected corpus',
                public_source='https://dynarepo.org/',
                next_step='Independent random-seed repeats require new trajectories or original producer raw replicates; declaredmetadata3 does not establish local3.'),
           dict(dataset='Zenodo6638504',available_states=24,true_pairs=16,independent_systems=8,
                local_independent_repeats_per_state=1,window_unit='threeor five100ns blocks inone300or500ns trajectory',
                current_corrected_oof_pairs=0,current_corrected_oof_systems=0,
                eligibility='RawstrictapoA/apoB/holoavailable; no corrected3063 FullOOF cache. Existing evaluation belongs tolegacyresponse experimentandcannotbe substituted.',
                public_source='https://zenodo.org/records/6638504',
                next_step='Beforefuture inference: remake chemical-correct labels, audit full-sequence familyoverlap, then choose a prespecified frozencheckpoint. Notrunhere.'),
           dict(dataset='KIX_local_four_states',available_states=4,true_pairs=5,independent_systems=1,
                local_independent_repeats_per_state=1,window_unit='single100nsproductionperstate; checkpoints/windowsnotreplicates',
                current_corrected_oof_pairs=5,current_corrected_oof_systems=1,
                eligibility='Correctedfold1RMSFsourceexport; exploratory singlefamily case, not independent-family performance',
                public_source='https://doi.org/10.1371/journal.pcbi.1009977',
                next_step='Noavailable3independentlocalrepeats. Currentfold1Hchiexportnotpresentinpaperbundle; notfilledfromoldmodels.'),
           dict(dataset='ATLAS_corrected_absolute',available_states=1929,true_pairs=0,independent_systems=1929,
                local_independent_repeats_per_state=3,window_unit='independentreplicasusedforabsoluteRMSF/Hchi',
                current_corrected_oof_pairs=0,current_corrected_oof_systems=0,
                eligibility='Canassessabsolute-stateandreplicatereliability; not a matchedapo/holo set',
                public_source='https://www.dsimb.inserm.fr/ATLAS/',
                next_step='Dousethree-replicateabsolutevalidation; do notcall3foldsorwindowsthreebindingexperiments.'),
           dict(dataset='KIX_Peiffer2022_publication',available_states=11,true_pairs=0,independent_systems=1,
                local_independent_repeats_per_state=0,window_unit='PublicationSIstatesfiveindependent>=100nstrialspercomplex',
                current_corrected_oof_pairs=0,current_corrected_oof_systems=0,
                eligibility='Public methodsandaggregateEPSfiguresareavailable; no downloadablepertrialrawtrajectorylinkidentifiedonarticleorSI',
                public_source='https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009977',
                next_step='Rawfive-trialtrajectorieswouldenable reuse aftertopology/protocolvalidation; authorcontactwouldrequireexplicitinstruction, notsent.')]
pd.DataFrame(summaries).to_csv(OUT/'candidate_dataset_inventory.csv',index=False)
print(pd.DataFrame(zenodo)[['local_trajectory_count','structure_exists','in_corrected3063','current_full_prediction_present']].value_counts().to_string())
