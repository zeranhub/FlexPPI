"""Read-only independent model/coordinate checks for the KRAS case.

Only the output audit is written. Existing predictions and structural inputs are
not changed. Contacts are geometric distances in the single model-input frame.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

EXP = Path(r'E:\Project\FlexPPI\experiments\_sources')
OUT = Path(__file__).resolve().parent
OLD = EXP / 'kras_case_v41_20261017'
CSV = OLD / 'figure6/source_all10_peptide_model_and_MD_changes.csv'
PRED = EXP / 'tcr_case_v40_20261017/family_predictions/full'
STATES = {'WT': '7ow5', 'G12D': '7ow6'}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def pdb_atoms(path):
    atoms = []
    for l in path.read_text().splitlines():
        if not l.startswith('ATOM  '):
            continue
        name = l[12:16].strip()
        element = l[76:78].strip() or name[0]
        if element in ('H', 'D'):
            continue
        atoms.append({'name': name, 'aa': l[17:20].strip(),
                      'chain': l[21], 'residue': int(l[22:26]),
                      'xyz': np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])})
    return atoms

def fmt_atom(a):
    return f"{a['chain']}:{a['aa']}{a['residue']}:{a['name']}"

def residue_contacts(pep, partner, cutoff=4.5):
    x = np.array([a['xyz'] for a in pep])
    y = np.array([a['xyz'] for a in partner])
    distances = np.sqrt(((x[:, None, :] - y[None, :, :]) ** 2).sum(2))
    minimum = np.unravel_index(np.argmin(distances), distances.shape)
    result = {}
    for j, p in enumerate(partner):
        i = int(np.argmin(distances[:, j]))
        d = float(distances[i, j])
        if d <= cutoff:
            key = f"{p['chain']}:{p['aa']}{p['residue']}"
            if key not in result or d < result[key]['min_distance_A']:
                result[key] = {'partner_residue': key, 'min_distance_A': d,
                               'peptide_atom': fmt_atom(pep[i]), 'partner_atom': fmt_atom(p)}
    return {'nearest_distance_A': float(distances[minimum]),
            'nearest_pair': [fmt_atom(pep[minimum[0]]), fmt_atom(partner[minimum[1]])],
            'heavy_atom_pair_count_within_cutoff': int((distances <= cutoff).sum()),
            'partner_residues_within_cutoff': sorted(result.values(), key=lambda d: d['min_distance_A'])}

def polar_pairs(pep, partner, cutoff=3.5):
    result=[]
    for a in pep:
        if a['name'][0] not in 'NO':
            continue
        for b in partner:
            if b['name'][0] not in 'NO':
                continue
            d=float(np.linalg.norm(a['xyz']-b['xyz']))
            if d<=cutoff:
                result.append({'peptide_atom':fmt_atom(a),'partner_atom':fmt_atom(b),'distance_A':d})
    return sorted(result,key=lambda r:r['distance_A'])

rows=list(csv.DictReader(CSV.open(encoding='utf8',newline='')))
assert len(rows)==10 and [int(r['position_a']) for r in rows]==list(range(1,11))
out={'status':'passed','policy':{'prediction_model':'same original Full family-fold-0 joint model trained without ITrace',
      'conditions':['WT 7OW5 complex-bound input','G12D 7OW6 complex-bound input'],
      'units':'Angstrom for RMSF and distances','all_10_positions_retained':True,
      'coordinate_contact_cutoff_A':4.5,'polar_atom_proximity_cutoff_A':3.5,
      'contacts_are_single_frame_geometry_not_hydrogen_bond_or_energy_assignments':True},
     'source_files':[],'all_10_model_values':[],'contacts':{}}
predictions={}
for state,pdb in STATES.items():
    p=PRED/f'itrace_{pdb}.npz'
    f=np.load(p)
    indices=np.flatnonzero(f['chain']=='C')
    assert len(indices)==10
    predictions[state]=f['mean'][indices,0].astype(float)
    assert f['author_residue_id'][indices].tolist()==[str(i) for i in range(1,11)]
    assert f['amino_acid'][indices].tolist()==list('VVVGAGGVGK' if state=='WT' else 'VVVGADGVGK')
    out['source_files'].append({'path':str(p),'sha256':sha(p)})
    atoms_path=OLD/f'case_sources/{pdb}_model_input.pdb'
    out['source_files'].append({'path':str(atoms_path),'sha256':sha(atoms_path)})
    atoms=pdb_atoms(atoms_path)
    out['contacts'][state]=[]
    for i in range(1,11):
        pep=[a for a in atoms if a['chain']=='C' and a['residue']==i]
        mhc=[a for a in atoms if a['chain']=='A']
        tcr=[a for a in atoms if a['chain'] in ('D','E')]
        out['contacts'][state].append({'position':i,'model_input_peptide_residue':pep[0]['aa'],
            'MHC':residue_contacts(pep,mhc),'TCR':residue_contacts(pep,tcr),
            'polar_atom_proximities_to_MHC':polar_pairs(pep,mhc),
            'polar_atom_proximities_to_TCR':polar_pairs(pep,tcr)})
for i,row in enumerate(rows):
    a=float(predictions['WT'][i]); b=float(predictions['G12D'][i]); delta=b-a
    assert abs(a-float(row['rmsf_full_a']))<1e-14
    assert abs(b-float(row['rmsf_full_b']))<1e-14
    assert abs(delta-float(row['rmsf_full_delta']))<1e-14
    out['all_10_model_values'].append({'position':i+1,'WT_AA':row['aa_a'],'G12D_AA':row['aa_b'],
       'WT_RMSF_A':a,'G12D_RMSF_A':b,'G12D_minus_WT_A':delta,
       'fraction_decrease':(a-b)/a})
out['source_files'].append({'path':str(CSV),'sha256':sha(CSV)})
out['regional_summary']={}
for name,lo,hi in [('N_terminal_P1_3',1,3),('response_region_P4_10',4,10),('central_P5_8',5,8),('central_P6_7',6,7)]:
    a=predictions['WT'][lo-1:hi];b=predictions['G12D'][lo-1:hi]
    out['regional_summary'][name]={'positions':list(range(lo,hi+1)),'WT_mean_RMSF_A':float(a.mean()),
      'G12D_mean_RMSF_A':float(b.mean()),'mean_delta_A':float((b-a).mean()),
      'fractional_mean_decrease':float((a.mean()-b.mean())/a.mean())}
mutant_atoms=pdb_atoms(OLD/'case_sources/7ow6_model_input.pdb')
carboxyl=[a for a in mutant_atoms if a['chain']=='C' and a['residue']==6 and a['name'] in ('OD1','OD2')]
tcr_atoms=[a for a in mutant_atoms if a['chain'] in ('D','E')]
arg114=[a for a in mutant_atoms if a['chain']=='A' and a['residue']==114 and a['name'] in ('NE','NH1','NH2')]
asp116=[a for a in mutant_atoms if a['chain']=='A' and a['residue']==116 and a['name'] in ('OD1','OD2')]
out['local_network_geometry']={
  'G12D_D6_carboxyl_to_TCR':residue_contacts(carboxyl,tcr_atoms),
  'G12D_MHC_Arg114_to_Asp116_polar_atoms':residue_contacts(arg114,asp116),
  'input_TCR_numbering_note':'Input chain D=TCR alpha and E=TCR beta. Beta HIS101/ASN102/GLY99 correspond to the article HIS100/ASN101/GLY98; the input is renumbered and should not silently be quoted as article numbering.'}
out['scientific_boundaries']=[
 'The RMSF quantities describe backbone motion in the bound TCR-pMHC complex, not unbound peptide preorganization.',
 'Lower bound-state RMSF is not itself a binding affinity, free energy or binding entropy calculation.',
 'The largest predicted reductions occur at P6 and the unchanged P7, and extend through unchanged P4/P5/P8/P9/P10.',
 'The WT high-amplitude segment P6-P10 remains higher in the mutant than P1-P5; the mutant is not globally flat or motion-free.',
 'Model-input coordinate distances are not original-crystal distances; use original study for named interaction and mechanistic claims.',
 'No new training, model inference, molecular dynamics or frame sampling was performed.']
OUT.mkdir(parents=True,exist_ok=True)
(OUT/'model_phenomena.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':out['status'],'values':out['all_10_model_values'],'regions':out['regional_summary'],
 'WT_P6_MHC':out['contacts']['WT'][5]['MHC'], 'G12D_P6_MHC':out['contacts']['G12D'][5]['MHC'],
 'G12D_P6_TCR':out['contacts']['G12D'][5]['TCR']},ensure_ascii=False,indent=2))
