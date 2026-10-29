"""Exploratory refinement after complete-residue contact analysis.

Separate sidechain from backbone contacts for every H3 residue; check fixed
three-state joint chi occupancies and four consecutive blocks of ONE trace.
"""
from pathlib import Path
import csv, importlib.util, json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

ROOT=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
OUT=ROOT/'hcd3_functional_analysis_20261029/md_audit'
PDB=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
XTC=Path(r'D:\data\dynarepo_non_nucleic_full\A00LZ\trajectory.xtc')
CP=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\src\flexppi\constants.py')
spec=importlib.util.spec_from_file_location('frozen_constants',CP)
constants=importlib.util.module_from_spec(spec);spec.loader.exec_module(constants)
u=mda.Universe(str(PDB),str(XTC))
h3=[r for r in u.residues if r.segid=='C' and 95<=r.resid<=109]
heavy=u.atoms[np.char.upper(u.atoms.elements.astype(str))!='H']
partners={chain:heavy[heavy.segids==chain] for chain in ['A','B']}
backbone={'N','CA','C','O','OXT'}
groups={}
quads={}
for r in h3:
    atoms=r.atoms[np.char.upper(r.atoms.elements.astype(str))!='H']
    groups[r.resid]={
        'sidechain':atoms[np.array([n not in backbone for n in atoms.names])],
        'backbone':atoms[np.array([n in backbone for n in atoms.names])],
    }
    names={a.name:a.index for a in r.atoms}
    chi=constants.CHI_ATOMS.get(r.resname,[])
    if chi:
        quads[r.resid]=np.array([[names[n] for n in q] for q in chi])

def dihedral(points):
    b0=-(points[:,1]-points[:,0]);b1=points[:,2]-points[:,1];b2=points[:,3]-points[:,2]
    b1=b1/np.maximum(np.linalg.norm(b1,axis=-1,keepdims=True),1e-8)
    v=b0-np.sum(b0*b1,axis=-1,keepdims=True)*b1
    w=b2-np.sum(b2*b1,axis=-1,keepdims=True)*b1
    return np.degrees(np.arctan2(np.sum(np.cross(b1,v)*w,axis=-1),np.sum(v*w,axis=-1)))

traces={r.resid:{f'{part}_{chain}':[] for part in groups[r.resid] for chain in partners} for r in h3}
codes={r.resid:[] for r in h3 if r.resid in quads}
for ts in u.trajectory:
    for r in h3:
        for part,atoms in groups[r.resid].items():
            for chain,partner in partners.items():
                traces[r.resid][f'{part}_{chain}'].append(float(distance_array(atoms.positions,partner.positions,box=ts.dimensions).min()))
        if r.resid in quads:
            angles=dihedral(u.atoms.positions[quads[r.resid]])
            bins=np.clip(np.floor((angles+180)/120).astype(int),0,2)
            codes[r.resid].append(int(sum(int(b)*(3**i) for i,b in enumerate(bins))))
all_rows=[];states=[]
labels={int(r['author']):r for r in csv.DictReader((ROOT/'6A77_all_positions_predictions_and_MD.csv').open(encoding='utf-8-sig')) if r['variant']=='full' and r['chain']=='C'}
max_entropy_error=0.
for r in h3:
    for block,start,stop in [('all',0,1000),('block1',0,250),('block2',250,500),('block3',500,750),('block4',750,1000)]:
        row={'author':int(r.resid),'resname':r.resname,'block':block,'start_frame':start,'stop_frame_exclusive':stop}
        for name,trace in traces[r.resid].items():
            a=np.array(trace[start:stop])
            for cutoff in [4.,4.5,5.]:
                row[f'{name}_fraction_{cutoff:g}A']=float(np.mean(a<=cutoff))
        for part in ['sidechain','backbone']:
            for cutoff in [4.,4.5,5.]:
                a=np.array(traces[r.resid][f'{part}_A'][start:stop])<=cutoff
                b=np.array(traces[r.resid][f'{part}_B'][start:stop])<=cutoff
                row[f'{part}_both_A_and_B_fraction_{cutoff:g}A']=float(np.mean(a&b))
        if r.resid in codes:
            c=np.array(codes[r.resid][start:stop]);values,counts=np.unique(c,return_counts=True)
            probs=counts/counts.sum()
            row['chi_joint_entropy']=float(-np.sum(probs*np.log(probs)))
            row['chi_joint_observed_states']=len(counts)
            row['chi_joint_dominant_fraction']=float(probs.max())
            if block=='all':
                max_entropy_error=max(max_entropy_error,abs(row['chi_joint_entropy']-float(labels[r.resid]['chi_entropy_MD'])))
            for code,count in zip(values,counts):
                states.append({'author':int(r.resid),'resname':r.resname,'block':block,'joint_code':int(code),'frames':int(count),'fraction':float(count/(stop-start))})
        all_rows.append(row)
fields=list(dict.fromkeys(k for row in all_rows for k in row))
for name,rows,fieldnames in [('all_h3_sidechain_backbone_and_blocks.csv',all_rows,fields),('all_h3_joint_state_occupancies.csv',states,list(states[0]))]:
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fieldnames);w.writeheader();w.writerows(rows)
audit={'status':'complete','max_joint_entropy_difference_from_frozen_reference':max_entropy_error,
       'scope':'post hoc exploratory refinement of all 15 H3 positions, separating sidechain/backbone; consecutive blocks of one trajectory are not independent replicates',
       'sidechain_definition':'all non-hydrogen atoms except N/CA/C/O/OXT; partner includes all heavy atoms',
       'proline_caution':'zero coarse chi entropy is not evidence of zero backbone motion or energetic importance'}
(OUT/'sidechain_refinement_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
assert max_entropy_error<1e-5
print(json.dumps(audit,indent=2))
for row in all_rows:
    if row['author'] in [96,97,101,104,105,109]:
        print({k:row.get(k) for k in ['author','block','sidechain_A_fraction_4A','sidechain_B_fraction_4A','sidechain_both_A_and_B_fraction_4A','backbone_B_fraction_4A','chi_joint_entropy','chi_joint_observed_states','chi_joint_dominant_fraction']})
