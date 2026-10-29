"""Frozen-OOF all-residue SASA features and train-only nonlinear controls.

No original model, split, label, structure or manuscript is modified.
Run with the local_model_benchmark_20261023/runtime_py312 Python runtime.
"""
from __future__ import annotations

import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import argparse
import csv
import hashlib
import importlib.util
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
from pathlib import Path
import warnings

import numpy as np

ROOT = Path('E:/Project/FlexPPI')
CORE = ROOT / 'experiments/core_corrected_pipeline'
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(CORE / 'training_payload/src'))
from flexppi.constants import THREE_TO_ONE, CHI_ATOMS
from flexppi.data.structures import parse_pdb

SPHERE_POINTS = 960
METRICS = ['rmsf', 'hchi']
SASA_COLUMNS = [
    'sidechain_sasa_complex_A2', 'residue_sasa_complex_A2',
    'sidechain_sasa_isolated_A2', 'residue_sasa_isolated_A2',
    'delta_sidechain_sasa_A2', 'delta_residue_sasa_A2',
]
EXPECTED = {
    'ALA':'CB', 'ARG':'CB CG CD NE CZ NH1 NH2', 'ASN':'CB CG OD1 ND2',
    'ASP':'CB CG OD1 OD2', 'CYS':'CB SG', 'GLN':'CB CG CD OE1 NE2',
    'GLU':'CB CG CD OE1 OE2', 'GLY':'', 'HIS':'CB CG ND1 CD2 CE1 NE2',
    'ILE':'CB CG1 CG2 CD1', 'LEU':'CB CG CD1 CD2', 'LYS':'CB CG CD CE NZ',
    'MET':'CB CG SD CE', 'PHE':'CB CG CD1 CD2 CE1 CE2 CZ',
    'PRO':'CB CG CD', 'SER':'CB OG', 'THR':'CB OG1 CG2',
    'TRP':'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR':'CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL':'CB CG1 CG2',
}
EXPECTED = {name:set(atoms.split()) for name, atoms in EXPECTED.items()}


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576), b''):
            result.update(block)
    return result.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf-8')
    temporary.replace(path)


def load_schema():
    return json.loads((CORE / 'static_corrected/training_schema.json').read_text(encoding='utf-8'))


def priority_ids():
    import pandas as pd
    path = ROOT / 'experiments/_sources/model_argument_v24_20261013/experiment_predictive_occupancy/replica1_per_residue.csv.gz'
    return set(pd.read_csv(path, usecols=['sample_id'])['sample_id'].unique())


def prepare():
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    OUT.mkdir(parents=True, exist_ok=True)
    schema = load_schema()
    records = schema['records']
    ready = [json.loads(line) for line in (CORE / 'training_payload/splits/all_ready.jsonl').read_text(encoding='utf-8').splitlines()]
    assert len(records) == len(ready) == 3063
    original = {record['sample_id']:record for record in ready}
    table = []
    for record in records:
        sid = record['sample_id']
        source = original[sid]
        assert record['fold'] == source['metadata']['fold']
        table.append(dict(record, structure_path=str(CORE / 'training_payload/dataset/structures' / (sid+'.pdb')),
                          label_path=str(CORE / 'training_payload/dataset/labels' / (sid+'.npz')),
                          full_oof_path=str(CORE / 'neural_exports/evaluation/full' / (sid+'.npz')),
                          intra_oof_path=str(CORE / 'neural_exports/evaluation/intrachain' / (sid+'.npz')),
                          source_group=source['source_group'],
                          pdb_id=source['metadata'].get('pdb_id',''),
                          priority83=sid in priority_ids_cached))
    pd.DataFrame(table).to_csv(OUT / 'cohort_records.csv', index=False)
    write_json(OUT / 'cohort_records.json', table)
    write_json(OUT / 'feature_contract.json', dict(
        records=3063, rows=961321, valid_counts=schema['valid_counts'],
        key=['sample_id','residue_index'], residue_index='0-based within frozen parsed/OOF record; never author numbering',
        author_key=['sample_id','chain','author_residue_id'],
        sasacolumns=SASA_COLUMNS, delta_definition='isolated-chain area minus complete-complex area, same coordinates and atom point grids',
        atoms='All heavy atoms in all protein residues recognized by the frozen parser; other hetero ligands/water/ions excluded',
        radius_probe_A=1.4,n_sphere_points=SPHERE_POINTS,
        sasa_input_richness='Full protein sidechains supply additional input beyond FlexPPI N/CA/C backbone; explanatory stronger-input baseline, not a matched-information ablation',
        static_feature_names=schema['names'], source_schema_sha256=sha(CORE / 'static_corrected/training_schema.json'),
        source_features_sha256=sha(CORE / 'static_corrected/training_features.npz'),
        schema_version='section33_sasa_full_cohort_v1'))
    if (OUT / 'base_residues.parquet').exists():
        print('PREPARE_REUSE base_residues.parquet', flush=True)
        return
    arrays = np.load(CORE / 'static_corrected/training_features.npz')
    x, target, valid = arrays['X'], arrays['target'], arrays['valid']
    static = np.load(CORE / 'static_corrected/static_sequence_complex_oof.npz')['prediction']
    aaonly = np.load(CORE / 'static_corrected/aa_only_oof.npz')['prediction']
    temp = OUT / 'base_residues.parquet.tmp'
    writer = None
    start_time = time.monotonic()
    for i, record in enumerate(table):
        lo, hi = record['start'], record['stop']
        full = np.load(record['full_oof_path'], allow_pickle=False)
        intra = np.load(record['intra_oof_path'], allow_pickle=False)
        assert hi-lo == len(full['mean'])
        assert np.array_equal(full['aatype'], intra['aatype'])
        assert np.array_equal(full['chain'], intra['chain'])
        assert np.array_equal(full['author_residue_id'], intra['author_residue_id'])
        assert np.array_equal(full['valid'][:,:2], valid[lo:hi])
        assert np.allclose(full['target'][:,:2][valid[lo:hi]], target[lo:hi][valid[lo:hi]], rtol=0, atol=1e-6)
        n = hi-lo
        block = dict(sample_id=np.repeat(record['sample_id'], n),record_index=np.full(n,i,np.int32),
                     residue_index=np.arange(n,dtype=np.int32),dataset=np.repeat(record['dataset'],n),
                     biological_system=np.repeat(record['biological_system'],n),fold=np.full(n,record['fold'],np.int8),
                     chain=full['chain'],author_residue_id=full['author_residue_id'],amino_acid=full['amino_acid'],
                     aatype=full['aatype'].astype(np.int8))
        for j, metric in enumerate(METRICS):
            block['valid_'+metric] = valid[lo:hi,j]
            block['md_'+metric] = target[lo:hi,j]
            block['pred_full_'+metric] = full['mean'][:,j]
            block['pred_intra_'+metric] = intra['mean'][:,j]
            block['pred_static_'+metric] = static[lo:hi,j]
            block['pred_aa_only_'+metric] = aaonly[lo:hi,j]
        for j,name in enumerate(schema['names']):
            block[name] = x[lo:hi,j]
        block['chain_length'] = np.rint(np.expm1(block['log_chain_length'])).astype(np.int32)
        block['distance_to_chain_end'] = np.rint(np.expm1(block['log_distance_to_chain_end'])).astype(np.int32)
        for radius in [6,8,10,12]:
            for which in ['intra','partner']:
                block[f'{which}_ca_contacts_{radius}A'] = np.rint(np.expm1(block[f'{which}_contacts_{radius}A'])).astype(np.int16)
        block['interface_ca8'] = block['partner_ca_contacts_8A'] > 0
        arrow = pa.Table.from_pydict(block)
        if writer is None:
            writer = pq.ParquetWriter(temp, arrow.schema, compression='zstd')
        writer.write_table(arrow)
        full.close();intra.close()
        if (i+1)%300 == 0:
            print('PREPARE', i+1, len(table), round(time.monotonic()-start_time,1), flush=True)
    writer.close()
    temp.replace(OUT / 'base_residues.parquet')
    print('PREPARE_DONE',len(table), round(time.monotonic()-start_time,1), flush=True)


def heavy_structure(path):
    """Mirror frozen parse_pdb residue/altloc rules, retain protein heavy atoms."""
    import mdtraj as md
    residues = {}
    alt_discard = 0
    with Path(path).open(encoding='utf-8', errors='replace') as stream:
        for line in stream:
            tag = line[:6].strip()
            if tag == 'ENDMDL':
                break
            if tag not in ['ATOM','HETATM'] or len(line)<54:
                continue
            name = line[17:20].strip().upper()
            if name not in THREE_TO_ONE:
                continue
            key = (line[21:22].strip() or '_', line[22:26].strip(), line[26:27].strip())
            if key not in residues:
                residues[key] = dict(name=name, atoms={}, altloc_discarded=0, duplicate_discarded=0)
            if line[16:17] not in [' ','A','1']:
                residues[key]['altloc_discarded'] += 1
                alt_discard += 1
                continue
            atom = line[12:16].strip().upper()
            symbol = line[76:78].strip().upper() if len(line)>=78 else ''
            if not symbol:
                symbol = 'SE' if atom=='SE' else next((x for x in atom if x.isalpha()), 'C')
            if symbol in ['H','D']:
                continue
            xyz = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
            if atom in residues[key]['atoms']:
                residues[key]['duplicate_discarded'] += 1
            else:
                residues[key]['atoms'][atom] = (xyz,symbol)
    # Frozen parser creates a residue only for accepted altlocs.
    residues = {key:r for key,r in residues.items() if r['atoms']}
    topology = md.Topology()
    chains = {}
    positions, atom_to_residue, atom_is_sidechain = [], [], []
    completeness, chifull, expected_n, observed_n, nonstandard, discarded = [],[],[],[],[],[]
    for ri,(key,row) in enumerate(residues.items()):
        ch,number,icode = key
        if ch not in chains:
            chains[ch] = topology.add_chain(chain_id=ch)
        residue = topology.add_residue(row['name'],chains[ch],resSeq=int(number))
        canonical = {'HID':'HIS','HIE':'HIS','HIP':'HIS','MSE':'MET'}.get(row['name'],row['name'])
        expected = EXPECTED[canonical].copy()
        if row['name']=='MSE':
            expected.discard('SD');expected.add('SE')
        observed = set(row['atoms'])
        completeness.append(expected <= observed)
        definitions = CHI_ATOMS.get(canonical, [])
        chi_atoms = set(atom for quad in definitions for atom in quad)
        if row['name']=='MSE':
            chi_atoms.discard('SD');chi_atoms.add('SE')
        chifull.append(chi_atoms <= observed)
        expected_n.append(len(expected))
        observed_n.append(len(observed - {'N','CA','C','O','OXT','OT1','OT2'}))
        nonstandard.append(row['name'] not in EXPECTED)
        discarded.append(row['altloc_discarded'])
        for atom,(xyz,symbol) in row['atoms'].items():
            element = md.element.get_by_symbol(symbol)
            topology.add_atom(atom,element,residue)
            positions.append(xyz);atom_to_residue.append(ri)
            atom_is_sidechain.append(atom not in {'N','CA','C','O','OXT','OT1','OT2'})
    trajectory = md.Trajectory(np.asarray(positions,np.float32)[None,:,:]/10, topology)
    return trajectory, list(residues), np.array(atom_to_residue), np.array(atom_is_sidechain), dict(
        sidechain_complete=np.asarray(completeness,bool),chi_input_complete=np.asarray(chifull,bool),
        expected_sidechain_heavy_atoms=np.asarray(expected_n,np.int16),observed_sidechain_heavy_atoms=np.asarray(observed_n,np.int16),
        nonstandard_residue=np.asarray(nonstandard,bool),discarded_altloc_atoms=np.asarray(discarded,np.int16))


def extract_one(record):
    import mdtraj as md
    from scipy.spatial import cKDTree
    sid = record['sample_id']
    cache = OUT / 'cache' / (sid+'.npz')
    if cache.exists():
        with np.load(cache) as values:
            assert str(values['input_sha256']) == record['input_sha256']
            assert int(values['sphere_points']) == SPHERE_POINTS
            return dict(sample_id=sid,status='cached',n=len(values['sidechain_complete']))
    tick = time.monotonic()
    assert sha(record['structure_path']) == record['input_sha256']
    structure = parse_pdb(record['structure_path'])
    trajectory,keys,residue_index,sidechain,quality = heavy_structure(record['structure_path'])
    assert len(keys) == structure.num_residues == record['stop']-record['start']
    assert [key[0] for key in keys] == structure.chain_ids
    assert [key[1]+key[2] for key in keys] == structure.author_residue_ids
    area = md.shrake_rupley(trajectory,probe_radius=.14,n_sphere_points=SPHERE_POINTS,mode='atom')[0]*100
    n = len(keys)
    complex_whole = np.bincount(residue_index,weights=area,minlength=n)
    complex_sc = np.bincount(residue_index[sidechain],weights=area[sidechain],minlength=n)
    isolated_whole,isolated_sc = complex_whole.copy(),complex_sc.copy()
    if len(trajectory.topology._chains)>1:
        for chain in trajectory.topology.chains:
            atoms = np.array([a.index for a in chain.atoms],int)
            chainarea = md.shrake_rupley(trajectory.atom_slice(atoms),probe_radius=.14,n_sphere_points=SPHERE_POINTS,mode='atom')[0]*100
            rids = residue_index[atoms];sc = sidechain[atoms]
            unique = np.unique(rids)
            whole = np.bincount(rids,weights=chainarea,minlength=n)
            scarea = np.bincount(rids[sc],weights=chainarea[sc],minlength=n)
            isolated_whole[unique] = whole[unique];isolated_sc[unique] = scarea[unique]
    deltas = [isolated_sc-complex_sc, isolated_whole-complex_whole]
    negative_count = [int(np.sum(d < -1e-5)) for d in deltas]
    negative_min = [float(d.min()) for d in deltas]
    assert min(negative_min)>=-1e-3, (sid, negative_min)
    mindist = np.full(n,np.inf,np.float32)
    for chain in np.unique(structure.chain_index):
        target = (structure.chain_index==chain)&structure.frame_mask
        partner = (structure.chain_index!=chain)&structure.frame_mask
        if partner.any():
            mindist[target] = cKDTree(structure.backbone[partner,1]).query(structure.backbone[target,1])[0]
    values = dict(zip(SASA_COLUMNS,[complex_sc,complex_whole,isolated_sc,isolated_whole,np.maximum(deltas[0],0),np.maximum(deltas[1],0)]))
    values.update(quality,min_partner_ca_distance_A=mindist,input_sha256=np.asarray(record['input_sha256']),
                  sphere_points=np.asarray(SPHERE_POINTS),heavy_atom_count=np.asarray(trajectory.n_atoms),
                  negative_delta_count=np.asarray(negative_count),negative_delta_min=np.asarray(negative_min))
    temp=cache.with_suffix('.npz.tmp')
    with temp.open('wb') as stream:
        np.savez_compressed(stream,**values)
    temp.replace(cache)
    return dict(sample_id=sid,status='new',n=n,heavy_atoms=trajectory.n_atoms,
                sidechain_incomplete=int((~quality['sidechain_complete']).sum()),seconds=time.monotonic()-tick)


def extract(workers=3):
    (OUT/'cache').mkdir(exist_ok=True)
    records = json.loads((OUT/'cohort_records.json').read_text(encoding='utf-8'))
    batches = [([r for r in records if r['priority83']],'priority83'),
               ([r for r in records if not r['priority83']],'remaining')]
    start = time.monotonic()
    progress=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for batch,label in batches:
            futures={pool.submit(extract_one,r):r for r in batch}
            for index,future in enumerate(as_completed(futures)):
                record=futures[future]
                try:
                    result=future.result()
                except Exception as error:
                    result=dict(sample_id=record['sample_id'],status='error',error=repr(error))
                progress.append(result)
                if result['status']=='error' or (index+1)%25==0:
                    print('SASA',label,index+1,len(batch),round(time.monotonic()-start,1),result,flush=True)
                if (index+1)%25==0 or index+1==len(batch):
                    write_json(OUT/'sasa_progress.json',dict(status='running',elapsed_seconds=time.monotonic()-start,results=progress))
            if label=='priority83':
                assemble(priority=True)
                print('PRIORITY83_READY',flush=True)
    write_json(OUT/'sasa_progress.json',dict(status='complete',elapsed_seconds=time.monotonic()-start,results=progress))
    print('SASA_ALL_DONE',round(time.monotonic()-start,1),flush=True)


def assemble(priority=False):
    import pyarrow as pa
    import pyarrow.parquet as pq
    import pandas as pd
    records=json.loads((OUT/'cohort_records.json').read_text(encoding='utf-8'))
    base=pq.ParquetFile(OUT/'base_residues.parquet')
    name='features_priority83.parquet' if priority else 'features_all.parquet'
    temp=OUT/(name+'.tmp')
    writer=None
    quality=[]
    for record in records:
        if priority and not record['priority83']:
            continue
        cache=OUT/'cache'/(record['sample_id']+'.npz')
        if not cache.exists():
            quality.append(dict(sample_id=record['sample_id'],status='missing'))
            continue
        block=base.read_row_group(record['record_index'])
        with np.load(cache) as z:
            for column in SASA_COLUMNS+['sidechain_complete','chi_input_complete','expected_sidechain_heavy_atoms',
                                      'observed_sidechain_heavy_atoms','nonstandard_residue','discarded_altloc_atoms','min_partner_ca_distance_A']:
                block=block.append_column(column,pa.array(z[column]))
            valid_hchi=block['valid_hchi'].to_numpy()
            quality.append(dict(sample_id=record['sample_id'],status='ok',n=len(block),dataset=record['dataset'],
                                incomplete_sidechains=int((~z['sidechain_complete']).sum()),
                                incomplete_valid_hchi=int((~z['sidechain_complete']&valid_hchi).sum()),
                                chi_incomplete_valid_hchi=int((~z['chi_input_complete']&valid_hchi).sum()),
                                nonstandard=int(z['nonstandard_residue'].sum()),
                                discarded_altloc_atoms=int(z['discarded_altloc_atoms'].sum()),
                                heavy_atoms=int(z['heavy_atom_count']),
                                negative_delta_count=int(z['negative_delta_count'].sum()),
                                negative_delta_min=float(z['negative_delta_min'].min())))
        if writer is None:
            writer=pq.ParquetWriter(temp,block.schema,compression='zstd')
        writer.write_table(block)
    if writer:
        writer.close();temp.replace(OUT/name)
    pd.DataFrame(quality).to_csv(OUT/('quality_priority83.csv' if priority else 'quality_all.csv'),index=False)
    write_json(OUT/('priority83_ready.json' if priority else 'features_all_ready.json'),dict(
        file=name,records=sum(r['status']=='ok' for r in quality),missing=[r['sample_id'] for r in quality if r['status']!='ok'],
        sha256=sha(OUT/name),n_sphere_points=SPHERE_POINTS))


def balance_weights(records,record_ids,train):
    counts=np.bincount(record_ids[train],minlength=len(records))
    active=[r for r in records if counts[r['record_index']]>0]
    systems=Counter((r['dataset'],r['biological_system']) for r in active)
    sources=Counter(source for source,system in systems)
    w=np.zeros(len(records),float)
    for r in active:
        index=r['record_index']
        w[index]=1/(len(sources)*sources[r['dataset']]*systems[(r['dataset'],r['biological_system'])]*counts[index])
    result=w[record_ids[train]]
    return result/result.mean()


def train(threads=3, variants=None):
    import pyarrow.parquet as pq
    import joblib
    from sklearn.ensemble import HistGradientBoostingRegressor
    from threadpoolctl import threadpool_limits
    from scipy.stats import spearmanr
    import pandas as pd
    import sklearn
    schema=load_schema()
    records=json.loads((OUT/'cohort_records.json').read_text(encoding='utf-8'))
    names=schema['names']+SASA_COLUMNS
    columns=names+['record_index','residue_index','fold','valid_rmsf','valid_hchi','md_rmsf','md_hchi','sidechain_complete']
    frame=pq.read_table(OUT/'features_all.parquet',columns=columns).to_pandas()
    aa=['aa_offset_+0']
    contacts=[name for name in schema['names'] if '_contacts_' in name]
    sets={'aa_sasa':aa+SASA_COLUMNS,'aa_sasa_contacts':aa+SASA_COLUMNS+contacts,
          'static47_sasa':schema['names']+SASA_COLUMNS}
    if variants:
        sets={name:sets[name] for name in variants}
    config=dict(loss='squared_error',learning_rate=.07,max_iter=120,max_leaf_nodes=31,max_depth=None,
                min_samples_leaf=50,l2_regularization=1.,max_bins=63,early_stopping=False,random_state=7)
    # All structure rows remain in their original fold; no feature-dependent record selection.
    # Complete-sidechain sensitivity is evaluated separately, never used to tune a model.
    contract=dict(configuration=config,feature_sets=sets,n_sphere_points=SPHERE_POINTS,threads=threads,
                  sklearn_version=sklearn.__version__,fold_protocol='original main fivefold frozen OOF',
                  all_training_records=3063,feature_input='full protein heavy atoms; more information than FlexPPI',
                  weights='equal source, equal biological system within source, equal record within system, equal valid residue within record',
                  no_test_tuning=True,missing_features='native histogram-boosting missing-value handling; all incomplete sidechains flagged')
    write_json(OUT/'baseline_contract.json',contract)
    jobs=OUT/'models';jobs.mkdir(exist_ok=True)
    fold=frame['fold'].to_numpy()
    record_ids=frame['record_index'].to_numpy()
    rowids=np.array([r['start'] for r in records],np.int64)[record_ids]+frame['residue_index'].to_numpy()
    valid=np.stack([frame['valid_'+m].to_numpy() for m in METRICS],axis=1)
    target=np.stack([frame['md_'+m].to_numpy() for m in METRICS],axis=1)
    progress=[]
    with threadpool_limits(limits=threads):
        for variant,cols in sets.items():
            x=frame[cols].to_numpy(dtype=np.float32)
            categorical=[name.startswith('aa_offset_') for name in cols]
            prediction=np.full(target.shape,np.nan,np.float32)
            for f in range(5):
                for mi,metric in enumerate(METRICS):
                    stem=f'{variant}_fold{f}_{metric}'
                    jobfile=jobs/(stem+'.npz')
                    trainmask=(fold!=f)&valid[:,mi]
                    testmask=(fold==f)&valid[:,mi]
                    if jobfile.exists():
                        with np.load(jobfile) as z:
                            assert np.array_equal(z['rowids'],rowids[testmask])
                            pred=z['prediction']
                        elapsed=0
                    else:
                        tick=time.monotonic()
                        model=HistGradientBoostingRegressor(**config,categorical_features=categorical)
                        model.fit(x[trainmask],target[trainmask,mi],sample_weight=balance_weights(records,record_ids,trainmask))
                        pred=model.predict(x[testmask]).astype(np.float32)
                        assert np.isfinite(pred).all()
                        joblib.dump(model,jobs/(stem+'.joblib'),compress=3)
                        np.savez_compressed(jobfile,rowids=rowids[testmask],prediction=pred)
                        elapsed=time.monotonic()-tick
                    prediction[testmask,mi]=pred
                    result=dict(variant=variant,fold=f,metric=metric,n_train=int(trainmask.sum()),n_test=int(testmask.sum()),
                                mae=float(np.abs(pred-target[testmask,mi]).mean()),seconds=elapsed)
                    progress.append(result)
                    write_json(OUT/'training_progress.json',dict(status='running',jobs=progress))
                    print('BASELINE',result,flush=True)
            np.savez_compressed(OUT/(variant+'_oof.npz'),prediction=prediction,rowids=rowids)
    write_json(OUT/'training_progress.json',dict(status='complete',jobs=progress))


def summarize():
    import pandas as pd
    import pyarrow.parquet as pq
    from scipy.stats import spearmanr
    columns=['sample_id','record_index','residue_index','dataset','biological_system','fold','sidechain_complete','interface_ca8']
    for metric in METRICS:
        columns += ['valid_'+metric,'md_'+metric]+['pred_'+v+'_'+metric for v in ['full','intra','static','aa_only']]
    frame=pq.read_table(OUT/'features_all.parquet',columns=columns).to_pandas()
    for variant in ['aa_sasa','aa_sasa_contacts','static47_sasa']:
        path=OUT/(variant+'_oof.npz')
        if path.exists():
            prediction=np.load(path)['prediction']
            assert len(prediction)==len(frame)
            for i,metric in enumerate(METRICS):
                frame['pred_'+variant+'_'+metric]=prediction[:,i]
    frame.to_parquet(OUT/'oof_predictions_all.parquet',index=False,compression='zstd')
    rows=[]
    for (sid,record_id),part in frame.groupby(['sample_id','record_index'],sort=False):
        for metric in METRICS:
            for subset,keep in [('all',np.ones(len(part),bool)),('complete_sidechains',part['sidechain_complete'].to_numpy()),
                                ('interface_ca8',part['interface_ca8'].to_numpy()),('noninterface_ca8',~part['interface_ca8'].to_numpy())]:
                mask=part['valid_'+metric].to_numpy()&keep
                if not mask.any():
                    continue
                y=part['md_'+metric].to_numpy()[mask]
                for name in part.columns:
                    if name.startswith('pred_') and name.endswith('_'+metric):
                        p=part[name].to_numpy()[mask]
                        if not np.isfinite(p).all():
                            continue
                        rho=float(spearmanr(y,p).statistic) if len(y)>=4 and np.ptp(y)>0 and np.ptp(p)>0 else np.nan
                        rows.append(dict(sample_id=sid,record_index=record_id,dataset=part['dataset'].iloc[0],
                                         biological_system=part['biological_system'].iloc[0],fold=int(part['fold'].iloc[0]),
                                         variant=name[len('pred_'):-len('_'+metric)],metric=metric,subset=subset,n=int(mask.sum()),
                                         mae=float(np.abs(p-y).mean()),rmse=float(np.sqrt(np.mean((p-y)**2))),rho=rho,bias=float((p-y).mean())))
    results=pd.DataFrame(rows)
    results.to_csv(OUT/'per_record_metrics.csv',index=False)
    results.groupby(['variant','metric','subset','dataset'])[['mae','rmse','rho','bias']].mean().to_csv(OUT/'record_macro_by_source.csv')
    results.groupby(['variant','metric','subset'])[['mae','rmse','rho','bias']].mean().to_csv(OUT/'record_macro_all.csv')
    print('SUMMARIZE_DONE',len(results),flush=True)


priority_ids_cached=set()
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('stage',choices=['prepare','extract','assemble','train','summarize','all'])
    ap.add_argument('--workers',type=int,default=3)
    ap.add_argument('--threads',type=int,default=3)
    ap.add_argument('--variants',nargs='+')
    args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    priority_ids_cached=priority_ids()
    if args.stage in ['prepare','all']:prepare()
    if args.stage in ['extract','all']:extract(args.workers)
    if args.stage in ['assemble','all']:assemble()
    if args.stage in ['train','all']:train(args.threads,args.variants)
    if args.stage in ['summarize','all']:summarize()
