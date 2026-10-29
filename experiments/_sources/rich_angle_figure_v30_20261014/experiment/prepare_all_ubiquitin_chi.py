"""Existing MD, all 66 valid ubiquitin Hchi sites and all defined chi angles.

No result-based selection, model inference or new MD. Only this directory is
written; missing angles remain NaN, never zero-filled. Source files are read-only.
"""
from __future__ import annotations
import sys
sys.dont_write_bytecode = True
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
from pathlib import Path
from datetime import datetime,timezone
import ast
import csv
import json
import hashlib
import shutil
import tempfile
from collections import Counter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
CORRECTED = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
PREDICTION = CORRECTED/'neural_exports/evaluation/full/A003R.npz'
CONSTANTS = CORRECTED/'training_payload/src/flexppi/constants.py'
LABEL_CODE = CORRECTED/'training_payload/src/flexppi/data/labels.py'
PROVENANCE = Path(r'E:\Project\FlexPPI\experiments\_case_validation_sources\ubiquitin\provenance.json')
PREVIOUS = BASE/'free_energy_curves_v27_20261014/case_curves'
ORIGINS = (0.,1.25,2.5,3.75)


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''): result.update(block)
    return result.hexdigest()


def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [clean(v) for v in x]
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x


def write_json(path,x):
    Path(path).write_text(json.dumps(clean(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')


def definitions():
    namespace={}
    exec(compile(CONSTANTS.read_text(encoding='utf8'),str(CONSTANTS),'exec'),namespace)
    tree=ast.parse(LABEL_CODE.read_text(encoding='utf8'))
    node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='_dihedral_degrees')
    angle={'np':np}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(LABEL_CODE),'exec'),angle)
    return namespace['CHI_ATOMS'],namespace['THREE_TO_ONE'],angle['_dihedral_degrees']


def freeze():
    HERE.mkdir(parents=True,exist_ok=True)
    prov=json.loads(PROVENANCE.read_text(encoding='utf8'))
    sources={str(PREDICTION):'b8a764739993e3bce18dedb3cff493c50a1a52613a08baac16f1803c4edabb5a',
        str(CONSTANTS):sha(CONSTANTS),str(LABEL_CODE):'b4bd55066faffeebe76f3fcf8c3ad1f6541ab592445a98b58b18277cfbcd7f10',
        str(PROVENANCE):sha(PROVENANCE),prov['source_structure']:prov['source_structure_sha256'],prov['source_trajectory']:prov['source_trajectory_sha256']}
    previous=json.loads((PREVIOUS/'extraction_audit.json').read_text(encoding='utf8'))
    assert previous['status']=='complete'
    for name in ('angles_501frames_all6GLN.csv','histograms_501frames.csv'):
        sources[str(PREVIOUS/name)]=previous['outputs'][name]['sha256']
    for path,digest in sources.items():assert sha(path)==digest,path
    prediction=np.load(PREDICTION,allow_pickle=False)
    indices=np.flatnonzero((prediction['chain']=='A')&prediction['valid'][:,1])
    assert len(indices)==66
    assert int((prediction['chain']=='A').sum())==72
    inventory=[dict(residue_index=int(i),chain=str(prediction['chain'][i]),author_residue_id=str(prediction['author_residue_id'][i]),
        amino_acid=str(prediction['amino_acid'][i]),residue=str(prediction['amino_acid'][i])+str(prediction['author_residue_id'][i]),
        predicted_Hchi_nats=float(prediction['mean'][i,1]),frozen_MD_joint_Hchi_nats=float(prediction['target'][i,1]),
        target_valid=True,is_proline=bool(prediction['amino_acid'][i]=='P')) for i in indices]
    pro=[r for r in inventory if r['is_proline']]
    assert len(pro)==3 and all(r['frozen_MD_joint_Hchi_nats']==0 for r in pro)
    plan=dict(protocol='All66_ubiquitin_valid_Hchi_per_angle_MD_occupancy_v1',created_utc=datetime.now(timezone.utc).isoformat(),
        analysis_code_sha256=sha(__file__),source_hash_files=sources,inventory=inventory,
        eligibility='Exactly every chain A residue with frozen Full A003R valid[:,1], including all three Pro with original zero joint Hchi; 66 sites, no prediction or MD outcome exclusion',
        trajectory=dict(start_ns=0.,end_ns=100.,frames=501,step_ps=200.,source_frame_indices=list(range(0,1001,2))),
        angles='Original unmodified float32 _dihedral_degrees and CHI_ATOMS from frozen corrected pipeline; original atom order retained; ILE CD canonicalized to CD1 in task copy only; HIS aliases canonicalized in task copy only',
        histogram='72 equal periodic 5degree bins; origin0; wrap(theta+180)%360; no smoothing, pseudocount or interpolation',
        Fstar='-ln(p_bin/pmax); empty observed bins have probability0 and FstarNaN; chemically absent chi slots have NaN for every measure including count/probability',
        per_angle_metrics='Per-angle occupied-bin Shannon entropy(natural log); total measure of all occupied bins Fstar<=1 including disconnected bins; report origin0 and arithmetic mean over origins0,1.25,2.5,3.75; do not interpret absent chi as zero',
        residue_metrics='Mean metrics over actual chemically defined angles only; model predicts original joint Hchi, neither per-angle histograms nor these width/entropy metrics',
        outputs='Recomputable float32 angles NPZ(501,66,4) with NaN padding and complete identity/model reference; all264 angle slots in CSV metric/bin tables and all132264 frame/site/angle combinations in long angle table',
        historical_checks='All six GLN 501frame angles must match v27 within1e-5 degrees circular tolerance; all432 two-selected-GLN origin0 histogram/probability/Fstar rows must match v27',
        raw_access='Verify raw source full SHA256, copy bytes to owned task staging, MDAnalysis reads task copy only; after close auto-remove temporary copy; verify sources unchanged after analysis',
        no_new_experiment='Existing MD analysis only; no new trajectory/training/model inference/server work')
    path=HERE/'analysis_plan.json'
    if path.exists():
        old=json.loads(path.read_text(encoding='utf8'))
        assert old['analysis_code_sha256']==sha(__file__) and old['source_hash_files']==sources and old['inventory']==inventory
        plan=old
    else:write_json(path,plan)
    return plan,prov


def canonical_copy(source,target):
    lines=Path(source).read_text(encoding='utf8').splitlines(True)
    changes=Counter()
    names={}
    for line in lines:
        if line.startswith(('ATOM  ','HETATM')) and len(line)>=54 and line[17:20]=='ILE':
            names.setdefault(line[21:27],set()).add(line[12:16].strip())
    assert not any({'CD','CD1'}.issubset(x) for x in names.values())
    output=[]
    for line in lines:
        if line.startswith(('ATOM  ','HETATM')) and len(line)>=54:
            old=line[17:20]
            if old in ('HSD','HSE','HSP','HID','HIE','HIP'):
                line=line[:17]+'HIS'+line[20:]
                changes[old+'->HIS']+=1
            if old=='ILE' and line[12:16].strip()=='CD':
                line=line[:12]+' CD1'+line[16:]
                changes['ILE_CD->CD1']+=1
        output.append(line)
    target.write_text(''.join(output),encoding='utf8')
    return dict(changes)


def bins_for(values,origin):
    indices=np.floor(((values.astype(float)+180-origin)%360)/5).astype(int)
    assert indices.min()>=0 and indices.max()<72
    counts=np.bincount(indices,minlength=72)
    assert counts.sum()==501
    p=counts.astype(float)/501
    f=np.full(72,np.nan,float)
    visited=counts>0
    f[visited]=-np.log(p[visited]/p.max())
    entropy=float(-np.dot(p[visited],np.log(p[visited])))
    keep=visited&(f<=1.)
    np.testing.assert_array_equal(keep,visited&(counts>=counts.max()*np.exp(-1.)))
    return counts,p,f,float(5*keep.sum()),entropy


def run(plan,prov):
    sys.path.insert(0,str(BASE/'label_sensitivity/runtime'))
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import calc_dihedrals
    chi_atoms,three_to_one,angle_function=definitions()
    inventory=plan['inventory']
    n=len(inventory)
    angles=np.full((501,n,4),np.nan,np.float32)
    available=np.zeros((n,4),bool)
    qcounts=np.zeros(n,np.int8)
    times=np.empty(501,float)
    source_stat={p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in plan['source_hash_files']}
    mappings=[]
    max_independent_error=0.
    max_bond_length=0.
    staging=HERE/'_staging'
    staging.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='all66_existing_',dir=staging) as folder:
        folder=Path(folder)
        assert folder.resolve().is_relative_to(staging.resolve())
        top,traj=folder/'canonical.pdb',folder/'trajectory.xtc'
        chemistry=canonical_copy(prov['source_structure'],top)
        shutil.copy2(prov['source_trajectory'],traj)
        assert sha(traj)==prov['source_trajectory_sha256']
        universe=mda.Universe(str(top),str(traj))
        try:
            assert len(universe.trajectory)==5001
            quads,locations=[],[]
            for si,row in enumerate(inventory):
                selected=universe.select_atoms(f"chainID A and resid {row['author_residue_id']}").residues
                assert len(selected)==1
                residue=selected[0]
                assert residue.ix==row['residue_index'] and three_to_one[residue.resname]==row['amino_acid']
                lookup={a.name:int(a.index) for a in residue.atoms}
                assert len(lookup)==len(residue.atoms)
                definitions_here=chi_atoms[residue.resname]
                assert 1<=len(definitions_here)<=4
                qcounts[si]=len(definitions_here)
                for qi,quad in enumerate(definitions_here):
                    assert all(name in lookup for name in quad),(row,quad)
                    quads.append([lookup[name] for name in quad])
                    locations.append((si,qi))
                    available[si,qi]=True
                mappings.append(row|dict(mda_resindex=int(residue.ix),canonical_resname=str(residue.resname),chi_count=int(qcounts[si]),
                    chi_atom_names=definitions_here,chi_atom_indices=[[lookup[name] for name in quad] for quad in definitions_here]))
            quads=np.asarray(quads,np.int64)
            for j,ts in enumerate(universe.trajectory[0:1001:2]):
                assert abs(ts.time-j*200)<=.05
                times[j]=ts.time
                points=ts.positions[quads]
                values=angle_function(points)
                assert np.isfinite(values).all()
                reference=np.degrees(calc_dihedrals(points[:,0].astype(float),points[:,1].astype(float),points[:,2].astype(float),points[:,3].astype(float)))
                error=abs((values.astype(float)-reference+180)%360-180)
                max_independent_error=max(max_independent_error,float(error.max()))
                max_bond_length=max(max_bond_length,float(np.linalg.norm(np.diff(points,axis=1),axis=-1).max()))
                for k,(si,qi) in enumerate(locations):angles[j,si,qi]=values[k]
        finally:universe.trajectory.close()
    assert max_independent_error<1e-3 and max_bond_length<5.
    assert np.isfinite(angles[:,available]).all() and np.isnan(angles[:,~available]).all()
    assert np.array_equal(available.sum(axis=1),qcounts)
    oldangles=pd.read_csv(PREVIOUS/'angles_501frames_all6GLN.csv')
    max_GLN_error=0.
    for site,frame in oldangles.groupby('residue',sort=False):
        si=next(i for i,row in enumerate(inventory) if row['residue']==site)
        assert len(frame)==501 and qcounts[si]==3
        np.testing.assert_array_equal(frame.frame_zero_based.to_numpy(),np.arange(0,1001,2))
        observed=angles[:,si,:3].astype(float)
        original=frame[[f'chi{q}_degrees' for q in (1,2,3)]].to_numpy(float)
        error=abs((observed-original+180)%360-180)
        max_GLN_error=max(max_GLN_error,float(error.max()))
        np.testing.assert_allclose(observed,original,atol=1e-5,rtol=0)
    counts=np.full((n,4,72),np.nan,float)
    probability=np.full_like(counts,np.nan)
    fstar=np.full_like(counts,np.nan)
    width_origin=np.full((n,4,4),np.nan,float)
    entropy_origin=np.full((n,4,4),np.nan,float)
    per_angle_rows,bin_rows,angle_rows=[],[],[]
    for si,row in enumerate(inventory):
        for qi in range(4):
            defined=bool(available[si,qi])
            meta=row|dict(chi=f'chi{qi+1}',chi_index=qi,chi_defined=defined,chi_count=int(qcounts[si]),n_frames=501)
            if defined:
                for oi,origin in enumerate(ORIGINS):
                    c,p,f,w,h=bins_for(angles[:,si,qi],origin)
                    width_origin[si,qi,oi]=w
                    entropy_origin[si,qi,oi]=h
                    if oi==0:counts[si,qi],probability[si,qi],fstar[si,qi]=c,p,f
            per_angle_rows.append(meta|dict(width_origin0_degrees=float(width_origin[si,qi,0]),
                width_four_origins_mean_degrees=float(width_origin[si,qi].mean()),
                fine_entropy_origin0_nats=float(entropy_origin[si,qi,0]),
                fine_entropy_four_origins_mean_nats=float(entropy_origin[si,qi].mean()),
                **{f'width_origin_{origin:g}_degrees':float(width_origin[si,qi,oi]) for oi,origin in enumerate(ORIGINS)},
                **{f'fine_entropy_origin_{origin:g}_nats':float(entropy_origin[si,qi,oi]) for oi,origin in enumerate(ORIGINS)}))
            for b in range(72):
                bin_rows.append(meta|dict(angle_lower_degrees=-180+5*b,angle_upper_degrees=-175+5*b,angle_center_degrees=-177.5+5*b,
                    count=float(counts[si,qi,b]),probability=float(probability[si,qi,b]),density_per_degree=float(probability[si,qi,b]/5),
                    F_star=float(fstar[si,qi,b]),low_Fstar_le1=bool(counts[si,qi,b]>0 and fstar[si,qi,b]<=1.) if defined else None))
            for j in range(501):
                angle_rows.append(meta|dict(frame_zero_based=j*2,time_ps=float(times[j]),time_ns=float(times[j]/1000),angle_degrees=float(angles[j,si,qi])))
    bins=pd.DataFrame(bin_rows)
    historic=pd.read_csv(PREVIOUS/'histograms_501frames.csv')
    matches=0
    for (site,chi),part in historic.groupby(['residue','chi'],sort=False):
        current=bins[(bins.residue==site)&(bins.chi==chi)]
        assert len(current)==len(part)==72
        for column in ('angle_lower_degrees','angle_upper_degrees','angle_center_degrees','count','probability','density_per_degree','F_star'):
            np.testing.assert_allclose(current[column].to_numpy(),part[column].to_numpy(),atol=1e-14,rtol=0,equal_nan=True)
        np.testing.assert_array_equal(current.low_Fstar_le1.to_numpy(),part.low_Fstar_le1.to_numpy())
        matches+=len(part)
    assert matches==432
    residue_rows=[]
    for si,row in enumerate(inventory):
        residue_rows.append(row|dict(chi_count=int(qcounts[si]),all_defined_chi_valid=True,
            mean_per_chi_width_origin0_degrees=float(np.nanmean(width_origin[si,:,0])),
            mean_per_chi_width_four_origins_degrees=float(np.nanmean(width_origin[si])),
            mean_per_chi_fine_entropy_origin0_nats=float(np.nanmean(entropy_origin[si,:,0])),
            mean_per_chi_fine_entropy_four_origins_nats=float(np.nanmean(entropy_origin[si]))))
    npz=HERE/'ubiquitin_all66_chi_501frames.npz'
    np.savez_compressed(npz,chi_degrees=angles,chi_defined=available,chi_count=qcounts,times_ps=times,source_frame_indices=np.arange(0,1001,2),
        residue_index=np.array([r['residue_index'] for r in inventory],np.int64),chain=np.array([r['chain'] for r in inventory]),
        author_residue_id=np.array([r['author_residue_id'] for r in inventory]),amino_acid=np.array([r['amino_acid'] for r in inventory]),
        residue=np.array([r['residue'] for r in inventory]),predicted_Hchi_nats=np.array([r['predicted_Hchi_nats'] for r in inventory]),
        frozen_MD_joint_Hchi_nats=np.array([r['frozen_MD_joint_Hchi_nats'] for r in inventory]),
        counts_origin0=counts,probability_origin0=probability,Fstar_origin0=fstar,width_per_origin_degrees=width_origin,
        fine_entropy_per_origin_nats=entropy_origin,bin_origins_degrees=np.array(ORIGINS),bin_centers_degrees=np.arange(-177.5,180,5),
        plan_sha256=np.array(sha(HERE/'analysis_plan.json')),prediction_source_sha256=np.array(plan['source_hash_files'][str(PREDICTION)]))
    outputs={npz.name:dict(sha256=sha(npz),bytes=npz.stat().st_size)}
    for name,data in [('all66_per_angle_metrics.csv',per_angle_rows),('all66_histograms_origin0.csv',bin_rows),
        ('all66_angles_long_501frames.csv',angle_rows),('all66_residue_summary.csv',residue_rows)]:
        path=HERE/name
        pd.DataFrame(data).to_csv(path,index=False)
        outputs[name]=dict(sha256=sha(path),bytes=path.stat().st_size,rows=len(data))
    for path,digest in plan['source_hash_files'].items():
        assert sha(path)==digest
        assert (Path(path).stat().st_size,Path(path).stat().st_mtime_ns)==source_stat[path]
    audit=dict(status='complete',plan_sha256=sha(HERE/'analysis_plan.json'),code_sha256=sha(__file__),
        source_files_sha256_size_and_mtime_unchanged=True,prediction_source=str(PREDICTION),
        MDAnalysis_version=mda.__version__,numpy_version=np.__version__,topology_canonicalization=chemistry,
        raw_source_access='SHA256/read copy only; MDAnalysis opens task staging copy only; staging auto-deleted',
        chainA_original_positions=72,chainA_Hchi_valid_positions=66,all66_have_all_chemically_defined_angles=True,
        chi_count_distribution={str(q):int((qcounts==q).sum()) for q in range(1,5)},defined_per_angle_positions=int(available.sum()),
        absent_angle_slots=int((~available).sum()),all_missing_angles_remain_NaN=True,
        angle_array_shape=list(angles.shape),frames=501,time_range_ns=[0.,100.],step_ps=200.,
        original_zero_reference_Pro_sites=[r['residue'] for r in inventory if r['is_proline']],
        independent_float64_dihedral_max_circular_error_degrees=max_independent_error,maximum_adjacent_atom_distance_angstrom=max_bond_length,
        six_GLN_501frame_angle_values_checked=501*6*3,six_GLN_max_circular_error_degrees=max_GLN_error,
        prior_two_site_histogram_rows_checked=matches,prior_histogram_probability_Fstar_match_atol=1e-14,
        original_target_and_prediction_unchanged=True,atom_mapping=mappings,outputs=outputs)
    write_json(HERE/'extraction_audit.json',audit)
    print(json.dumps(clean({k:v for k,v in audit.items() if k not in ('atom_mapping','outputs')}),ensure_ascii=False,indent=2))


if __name__=='__main__':
    plan,provenance=freeze()
    run(plan,provenance)
