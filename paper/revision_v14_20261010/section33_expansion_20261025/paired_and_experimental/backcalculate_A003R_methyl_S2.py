"""Read-only geometric methyl-axis S2 from the bound A003R trajectory.

This supplies a correctly defined observable, not condition-matched NMR or
model validation. Raw XTC reader is used to avoid writing index sidecars.
"""
import json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.formats.libmdaxdr import XTCFile

OUT=Path(__file__).resolve().parent
PDB=Path(r'D:\data\dynarepo_non_nucleic_full\A003R\structure.pdb')
XTC=PDB.with_name('trajectory.xtc')

def rotation(mobile,reference):
    a=mobile-mobile.mean(axis=0);b=reference-reference.mean(axis=0)
    u,s,vt=np.linalg.svd(a.T@b)
    sign=np.linalg.det(u@vt)
    u[:,-1]*=sign
    return u@vt

def s2(vectors):
    # vectors: frames x methyl axes x xyz, normalized and in molecule frame.
    moment=np.einsum('fmi,fmj->mij',vectors,vectors)/len(vectors)
    return 1.5*np.square(moment).sum(axis=(1,2))-.5

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    # Meaningful checks for the geometric observable and the fit convention.
    rng=np.random.default_rng(7)
    vector=np.array([1.,0.,0.])
    assert np.allclose(s2(np.tile(vector,(30,1,1))),1)
    sample=rng.normal(size=(100000,1,3));sample/=np.linalg.norm(sample,axis=-1,keepdims=True)
    assert abs(s2(sample)[0])<.001
    reference=rng.normal(size=(30,3))
    u,_,vt=np.linalg.svd(rng.normal(size=(3,3)));r=u@vt
    if np.linalg.det(r)<0:r[:,0]*=-1
    rr=rotation(reference@r+3.,reference)
    assert np.allclose(r@rr,np.eye(3),atol=1e-9)
    universe=mda.Universe(str(PDB))
    data=pd.read_csv(OUT/'Fares2009_to_A003R_all50_mapping.csv')
    data=data[data.A003R_present].copy()
    target=universe.select_atoms('chainID A')
    ca=target.select_atoms('name CA').indices
    assert len(ca)==72
    atom_pairs=[];labels=[];actual_names=[]
    for _,row in data.iterrows():
        res=target.select_atoms('resid '+str(int(row.residue_number)))
        first=res.select_atoms('name '+row.atom1)
        second=res.select_atoms('name '+row.atom2)
        # ILE aliases in raw force-field topologies are explicit and audited.
        if row.residue_name=='ILE' and row.atom2=='CD1' and not len(second):
            second=res.select_atoms('name CD')
        assert len(first)==len(second)==1,(row.residue_number,row.atom1,row.atom2)
        atom_pairs.append([first.indices[0],second.indices[0]])
        actual_names.append((first.names[0],second.names[0]))
        labels.append(row.methyl_number)
    atom_pairs=np.asarray(atom_pairs)
    vectors=[];times=[];reference=None
    axis_min=np.inf;axis_max=0.;ca_max=0.
    with XTCFile(str(XTC),'r') as trajectory:
        for frame in trajectory:
            xyz=np.asarray(frame.x,dtype=np.float64)*10. # native XTC nm -> A
            if reference is None:reference=xyz[ca].copy()
            rr=rotation(xyz[ca],reference)
            delta=xyz[atom_pairs[:,1]]-xyz[atom_pairs[:,0]]
            lengths=np.linalg.norm(delta,axis=1)
            axis_min=min(axis_min,float(lengths.min()));axis_max=max(axis_max,float(lengths.max()))
            ca_max=max(ca_max,float(np.linalg.norm(np.diff(xyz[ca],axis=0),axis=1).max()))
            if lengths.min()<1.1 or lengths.max()>2.2 or ca_max>4.5:
                raise ValueError('Broken covalent geometry/PBC; no silent repair')
            vectors.append((delta/lengths[:,None])@rr)
            times.append(frame.time)
    vectors=np.asarray(vectors);times=np.asarray(times)
    assert vectors.shape==(5001,48,3)
    assert np.allclose(np.diff(times),100.)
    # Match original label frame coverage: first5000 of5001, split five100ns.
    values=vectors[:5000]
    full=s2(values)
    blocks=np.asarray([s2(b) for b in np.array_split(values,5)])
    out=data.copy();out['md_bound_S2_500ns']=full
    out['actual_topology_atom1']=[x[0] for x in actual_names]
    out['actual_topology_atom2']=[x[1] for x in actual_names]
    out['md_bound_block_mean_S2']=blocks.mean(axis=0)
    out['md_bound_block_SD_S2']=blocks.std(axis=0,ddof=1)
    for i in range(5):out['md_bound_S2_block'+str(i+1)]=blocks[i]
    out['analysis_interpretation']='same_observable_backcalculation_only_not_state_matched_NMR_validation'
    out.to_csv(OUT/'A003R_bound_methyl_S2_with_unbound_NMR_reference.csv',index=False)
    np.savez_compressed(OUT/'A003R_bound_aligned_methyl_unit_vectors.npz',vectors=values,
                        time_ps=times[:5000],methyl_number=np.array(labels),atom_pairs=atom_pairs)
    protocol=dict(sample_id='A003R',state='2QHO_ubiquitin_EDD1_UBA_bound',temperature_K=310,
                  atom_selection='48 methyl axes present in72residue chainA; includes axes without identifiable NMR fit',
                  source_NMR_valid_axes_present=35,
                  alignment='Independent proper Kabsch fit of72 chainA CA to first trajectory frame, then rotate every methyl axial unit vector by same rotation.',
                  estimator='S2=1.5*sum_ij(mean_t(u_i*u_j)^2)-0.5; axis methyl rotation excluded',
                  time_sampling='5000frames;100ps spacing;0-499.9ns;five nonoverlapping100ns blocks; blocks not replicates',
                  interpretation='Ensemble geometrical order parameter.100ps saving cannot resolve fast correlation functions; finite500ns coverage does not match all RDC averaging times. Bound vs free mismatch prevents direct validation; no correlation or fitted calibration is performed.',
                  PBC='Raw coordinates, no wrapping/unwrapping; every axis bond1.1-2.2A and all consecutive targetCA <=4.5A checked for every frame.',
                  bond_min_angstrom=axis_min,bond_max_angstrom=axis_max,consecutive_CA_max_angstrom=ca_max,
                  checks='rigid-axisS2=1; isotropicS2~0; proper Kabsch recovery ofknown rigid rotation; exact48axis5001frame shape;100ps spacing; atom-name uniqueness',
                  source_files=[dict(path=str(p),sha256=sha(p)) for p in [PDB,XTC]])
    (OUT/'A003R_methyl_S2_protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    print(json.dumps(dict(axes=len(data),frames_used=len(values),valid_experimental_axes=35,
                         mean_bound_S2=float(full.mean()),bond_range_A=[axis_min,axis_max])))

if __name__=='__main__':main()
