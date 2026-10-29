from pathlib import Path
import json, csv, re, hashlib

ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
PROJECT=Path('E:/Project/FlexPPI')
OLD=PROJECT/'output/tcr_external_validation_search_20260917'
ANN=EXP/'functional_annotations_complete/record_annotations.csv'
with ANN.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
pat=re.compile(r't.?cell receptor|TCR[ -]|histocompatibility|HLA[ -]|MHC[ -]',re.I)
hits=[r for r in rows if pat.search(r['protein_descriptions']+' '+r['pdb_classification'])]
with (ROOT/'atlas_dynarepo_immune_annotation_hits.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows([r for r in hits if r['dataset']!='itrace'])

manifest_path=Path('D:/data/dynarepo_non_nucleic_manifest.json')
dyn=json.loads(manifest_path.read_text())
def walk(obj):
    if isinstance(obj,dict):
        if 'accession' in obj and any(k in obj for k in ['name','description']): yield obj
        else:
            for v in obj.values(): yield from walk(v)
    elif isinstance(obj,list):
        for v in obj: yield from walk(v)
dhits=[r for r in walk(dyn) if pat.search(str(r.get('name',''))+' '+str(r.get('description','')))]
(ROOT/'dynarepo_manifest_keyword_hits.json').write_text(json.dumps(dhits,indent=2),encoding='utf-8')

with Path('E:/data/protein_md_datasets/thermopcd_index.csv').open(encoding='utf-8-sig',newline='') as f: thermo=list(csv.DictReader(f))
thits=[r for r in thermo if pat.search(r['name'])]
with (ROOT/'thermopcd_index_keyword_hits.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(thermo[0]));w.writeheader();w.writerows(thits)

with (OLD/'catalogues/DynaTPH_pdb_overlap.csv').open(newline='',encoding='utf-8-sig') as f: dyntph=list(csv.DictReader(f))
new=[r for r in dyntph if r['in_itrace_source'].lower()!='true']
with (ROOT/'DynaTPH_59_outside_ITrace.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(dyntph[0]));w.writeheader();w.writerows(new)

train0=[json.loads(l) for l in (EXP/'family_folds/splits/fold_0_train.jsonl').read_text().splitlines()]
test0=[json.loads(l) for l in (EXP/'family_folds/splits/fold_0_test.jsonl').read_text().splitlines()]
summary={
    'scope':'Local inventory, not new online availability verification.',
    'annotation_counts':{d:sum(r['dataset']==d for r in rows) for d in sorted({r['dataset'] for r in rows})},
    'atlas_dynarepo_immune_annotation_hits':[r for r in hits if r['dataset']!='itrace'],
    'dynarepo_local_manifest_keyword_hits':len(dhits),
    'thermopcd_titles_scanned':len(thermo),'thermopcd_keyword_hits':thits,
    'DynaTPH':{'catalogue_date':'2026-09-17','total_pdb':len(dyntph),'outside_ITrace':len(new),
       'full_trajectories_local':False,'available_local':'descriptor.csv and peptide-only mean RMSF table (3 replicas), no individual-repeat arrays or chi entropy',
       'units_caution':'rmsf.csv header declares nm but values around0.3–0.8; verify fit protocol and units against authors before comparing absolute predictions.',
       'old_README_access':'Full trajectories restricted for peer review at Zenodo21971877. Re-check if needed; not presumed currently inaccessible.'},
    '3KXF':{'old_exclusion_reason':'Exact PDB overlap with DynaRepo A005I in previous general training pool, not a trajectory integrity exclusion.',
      'in_current_family_fold0_train':any(r['sample_id']=='A005I' for r in train0),
      'in_current_family_fold0_test':any(r['sample_id']=='A005I' for r in test0),
      'ITrace_md_paths':[str(p) for p in sorted(Path('E:/data/ITrace/trajectories').glob('3kxf*'))],
      'dynarepo_structure':str(Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline/dataset/structures/A005I.pdb')),
      'use':'Can evaluate with same family-fold0 checkpoint; must independently generate ITrace targets and preserve source-specific chain mapping. Repeated structure/second MD protocol does not add independent PDB.'},
    'knapp2019':{'root':str(PROJECT/'artifacts/tcr_pmhc_reassessment_20260915/knapp2019'),
      'local_A6_partial_size':(PROJECT/'artifacts/tcr_pmhc_reassessment_20260915/knapp2019/A6.tar.gz.partial').stat().st_size,
      'expected_A6_bytes':399223244,'full_archive_checksum_verified':False,
      'other_sources':'1G4, JM22 and LC13 public Figshare URLs and exact archive sizes in local literature/knapp2019_*_figshare.json; not present as complete local trajectories.',
      'A6_construct_caution':'Three internal beta-chain regions removed in prepared simulation construct; do not treat as intact natural A6.'},
    'verified_current_train0_size':len(train0),
    'source_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ANN,manifest_path,OLD/'catalogues/DynaTPH_pdb_overlap.csv']}
}
(ROOT/'additional_sources_inventory.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:summary[k] for k in ['annotation_counts','atlas_dynarepo_immune_annotation_hits','dynarepo_local_manifest_keyword_hits','thermopcd_titles_scanned','thermopcd_keyword_hits','3KXF']},ensure_ascii=False,indent=2))
