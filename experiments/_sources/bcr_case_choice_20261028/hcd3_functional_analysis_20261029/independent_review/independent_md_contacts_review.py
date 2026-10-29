"""Read-only independent verification of the existing H3 contact analysis."""
from pathlib import Path
import csv, json, hashlib
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array, calc_dihedrals

R = Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
D = R/'hcd3_functional_analysis_20261029/md_audit'
O = R/'hcd3_functional_analysis_20261029/independent_review'
O.mkdir(parents=True,exist_ok=True)
P=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
T=Path(r'D:\data\dynarepo_non_nucleic_full\A00LZ\trajectory.xtc')
u=mda.Universe(str(P),str(T))
stored=np.load(D/'all_h3_partner_minimum_distances.npz')
res=list(u.residues)
h3=[r for r in res if r.segid=='C' and 95<=r.resid<=109]
heavy=u.atoms[u.atoms.elements!='H']
assert np.array_equal(u.atoms.chainIDs,u.atoms.segids)
assert set(u.atoms.elements)=={'H','C','N','O','S'}
starts=np.flatnonzero(np.r_[True,np.diff(heavy.resindices)!=0])
sample_frames=[0,111,249,250,499,500,749,750,888,999]
max_minimum_distance_error=0.
for fr in sample_frames:
    u.trajectory[fr]
    assert np.allclose(u.dimensions[3:],90.)
    box=u.dimensions[:3].astype(float)
    target=heavy.positions.astype(float)
    for j,r in enumerate(h3):
        atoms=r.atoms[r.atoms.elements!='H']
        delta=atoms.positions.astype(float)[:,None,:]-target[None,:,:]
        delta-=box*np.rint(delta/box)
        minimum=np.minimum.reduceat(np.sqrt(np.sum(delta*delta,axis=2)).min(axis=0),starts)
        max_minimum_distance_error=max(max_minimum_distance_error,float(np.max(np.abs(minimum-stored['minimum_distances_A'][fr,j]))))
assert max_minimum_distance_error<2e-5

original=list(csv.DictReader((D/'all_h3_sidechain_backbone_and_blocks.csv').open(encoding='utf-8-sig')))
original_index={(int(r['author']),r['block']):r for r in original}
all_context=list(csv.DictReader((D/'all_h3_dynamic_context.csv').open(encoding='utf-8-sig')))
context={int(r['author']):r for r in all_context}
focus=[96,97,101,104,105,109]
backbone={'N','CA','C','O','OXT'}
partner={c:heavy[heavy.chainIDs==c] for c in ['A','B']}
groups={}
quadruplets={
    'ARG': [('N','CA','CB','CG'),('CA','CB','CG','CD'),('CB','CG','CD','NE'),('CG','CD','NE','CZ')],
    'GLU': [('N','CA','CB','CG'),('CA','CB','CG','CD'),('CB','CG','CD','OE1')],
    'TYR': [('N','CA','CB','CG'),('CA','CB','CG','CD1')],
}
selected={int(r.resid):r for r in h3 if r.resid in focus}
angles={n:[] for n in focus}
traces={}
for n,r in selected.items():
    at=r.atoms[r.atoms.elements!='H']
    groups[n]={
        'sidechain':at[np.array([name not in backbone for name in at.names])],
        'backbone':at[np.array([name in backbone for name in at.names])],
    }
    assert len(groups[n]['sidechain'])+len(groups[n]['backbone'])==len(at)
    assert len(set(groups[n]['sidechain'].indices)&set(groups[n]['backbone'].indices))==0
    for part in groups[n]:
        for chain in partner:
            traces[n,part,chain]=[]
for ts in u.trajectory:
    for n,r in selected.items():
        for part,at in groups[n].items():
            for chain,pa in partner.items():
                traces[n,part,chain].append(float(distance_array(at.positions,pa.positions,box=ts.dimensions).min()))
        names={at.name:at for at in r.atoms}
        quad=quadruplets[r.resname]
        arrays=[np.array([names[q[col]].position for q in quad]) for col in range(4)]
        angles[n].append(np.degrees(calc_dihedrals(*arrays,box=ts.dimensions)))

max_occupancy_error=0.
max_entropy_error=0.
focus_summaries=[]
blocks=[('all',0,1000),('block1',0,250),('block2',250,500),('block3',500,750),('block4',750,1000)]
for n in focus:
    bins=np.clip(np.floor((np.array(angles[n])+180.)/120.).astype(int),0,2)
    joint=(bins*(3**np.arange(bins.shape[1]))).sum(axis=1)
    for block,start,stop in blocks:
        origin=original_index[n,block]
        for part in ['sidechain','backbone']:
            for cutoff in [4.,4.5,5.]:
                by_chain={}
                for chain in ['A','B']:
                    events=np.array(traces[n,part,chain][start:stop])<=cutoff
                    by_chain[chain]=events
                    observed=float(events.mean())
                    expected=float(origin[f'{part}_{chain}_fraction_{cutoff:g}A'])
                    max_occupancy_error=max(max_occupancy_error,abs(observed-expected))
                both=float(np.mean(by_chain['A']&by_chain['B']))
                expected=float(origin[f'{part}_both_A_and_B_fraction_{cutoff:g}A'])
                max_occupancy_error=max(max_occupancy_error,abs(both-expected))
        _,counts=np.unique(joint[start:stop],return_counts=True)
        probs=counts/counts.sum()
        entropy=float(-np.sum(probs*np.log(probs)))
        max_entropy_error=max(max_entropy_error,abs(entropy-float(origin['chi_joint_entropy'])))
    origin=original_index[n,'all']
    focus_summaries.append({k:origin[k] for k in ['author','resname','sidechain_A_fraction_4A','sidechain_B_fraction_4A','sidechain_both_A_and_B_fraction_4A','backbone_B_fraction_4A','chi_joint_entropy']})
assert max_occupancy_error<1e-12
assert max_entropy_error<1e-7
review={
    'status':'passed_independent_sample_and_focused_full_trace_recomputation',
    'full_H3_sampled_frames':sample_frames,
    'sampled_PBC_method':'manual minimum image over orthorhombic box, rather than MDAnalysis distance_array',
    'max_stored_pair_minimum_distance_error_A':max_minimum_distance_error,
    'full_trace_focus_authors':focus,
    'max_sidechain_backbone_and_joint_contact_fraction_error':max_occupancy_error,
    'max_joint_chi_entropy_difference_nats':max_entropy_error,
    'chi_method':'independent MDAnalysis calc_dihedrals + same specified 120-degree three-state bins',
    'focus_values':focus_summaries,
    'code_findings':[
        'All protein heavy atoms retained; all elements explicitly H/C/N/O/S; chainIDs equal segids.',
        'Residue heavy atoms are contiguous; minimum.reduceat correctly reduces atom minima into 519 residue minima.',
        'Backbone N/CA/C/O/OXT and remaining sidechain heavy atoms partition residues correctly.',
        'Both-chain contact fractions require simultaneous events in the same frame.',
        'Whole-heavy-chain fit reproduces label; other fits define different observables and do not validate alternative labels.'
    ],
    'claim_limits':[
        'One correlated trajectory only. Four contiguous blocks do not constitute independent replicate evidence.',
        'A contact is a geometric proximity event. It does not itself define a hydrogen bond, interaction energy, hotspot, or causal binding requirement.',
        'Low coarse joint Hchi means concentrated occupancy under fixed bins; it is not zero motion or thermodynamic binding entropy.',
        'R96/R101 are a useful same-residue-type contrast; they do not isolate Full-model benefit because Intra also recovers the entropy contrast.',
        'Y104/Y105 may illustrate simultaneous cross-chain contact and relatively concentrated sidechain occupancy in the sampled bound state, but not an accurately recovered complete same-Tyr ranking or affinity mechanism.'
    ],
    'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P,T,Path(__file__)]},
}
(O/'independent_review.json').write_text(json.dumps(review,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(review,ensure_ascii=False,indent=2))
