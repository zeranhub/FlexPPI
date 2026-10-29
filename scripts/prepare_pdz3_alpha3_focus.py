"""No-ligand PDZ3 alpha3 perturbation protocol, all author IDs retained.

This builds fixed-backbone computational inputs, not experimental mutants.
Never creates a deletion-plus-ligand structure. Use pvb_cpu Python.
"""
from pathlib import Path
import json
import hashlib
from prepare_mechanistic_case_studies import atoms, residues, AA3

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/pdz3_alpha3_focus_20260914'
STATES = []

def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding='utf-8')

def build(background, context, mutations=()):
    pdb = {'apo1bfe':'1BFE','strip5heb':'5HEB','strip1be9':'1BE9'}[background]
    low = 306 if background == 'apo1bfe' else 303
    if context in ('dN','dNd8'):
        low = 311
    high = {'full':402,'d7':395,'d8':394,'dN':402,'dNd8':394}[context]
    mutname = '_'.join(mutations) or 'WT'
    sid = f'{background}_{context}_{mutname}'
    if any(x['id'] == sid for x in STATES):
        return sid
    changes = {int(m[1:-1]):m[-1] for m in mutations}
    source = atoms(pdb)
    ids, seq = residues(source, 'A')
    source_seq = dict(zip(map(int, ids), seq))
    canonical = dict(source_seq)
    canonical[328] = 'V'
    for m in mutations:
        pos = int(m[1:-1])
        assert canonical[pos] == m[0], (sid,m,canonical[pos])
        assert low <= pos <= high, (sid,m,low,high)
    changes = {328:'V', **changes}
    lines = []
    for line in source:
        pos = int(line[22:26])
        if line[21] != 'A' or not low <= pos <= high:
            continue
        if pos in changes and source_seq[pos] != changes[pos]:
            if line[12:16].strip() not in ('N','CA','C','O','OXT'):
                continue
            line = line[:17] + AA3[changes[pos]] + line[20:]
        lines.append(line)
    observed, aaseq = residues(lines, 'A')
    assert list(map(int, observed)) == list(range(low,high+1)), (sid,observed)
    path = BASE/'states'/f'{sid}.pdb'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines)+'\nEND\n', encoding='ascii')
    STATES.append(dict(id=sid,path=str(path),chains=['A'],background=background,context=context,
        mutations=list(mutations),mutant=mutname,pdb_id=pdb,source_model=1,
        author_low=low,author_high=high,n_residues=len(observed),sequence=aaseq,
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        provenance='Experimental apo backbone, sequence I328V standardized' if pdb=='1BFE'
          else 'Ligand-stripped crystal backbone; not an experimentally observed apo ensemble',
        intervention='Sequence-only substitutions and/or coordinate-preserving terminal truncation; no relaxation'))
    return sid

def main():
    BASE.mkdir(exist_ok=True,parents=True)
    keymut = ['V328A','V328I','F337A','F337Y','E334A','E334Q','E334K',
              'K355A','K355Q','K355E','Q358A','L323A','H372A',
              'R312A','D357A','I336A','I341A','E373A']
    helixmut = ['E396A','Y397A','Y397F','S398A','R399A','R399Q','R399E','R399K',
                'F400A','F400Y','E401A','E401Q','E401K']
    doubles = [('F337A','F400A'),('F337A','Y397A'),('E334Q','R399Q'),('K355Q','E401Q')]
    # Protocol locked before new inference. All tested mutations are retained.
    dump(BASE/'protocol.json',dict(question='Does alpha3 constrain specific residues, and which contacts/identities matter?',
        prohibited='No deletion-ligand inputs or deletion-bound rescue conclusion.',
        primary_backgrounds=['apo1bfe','strip5heb'],sensitivity_background='strip1be9',
        primary_deletion='396-402 (d7)',boundary_sensitivity='395-402 (d8, 2026 paper construct)',
        no_inference_on_bound_states=True,metrics=['rmsf','chi_entropy'],
        predeclared_regions={'carboxylate_loop':[318,323],'beta2_beta3_loop':[330,336],'middle_region':[341,356]},
        key_mutations=keymut,helix_mutations=helixmut,double_mutations=doubles,
        systematic_screen='Every non-A residue in312-377 -> A, full/d8, both primary backbones; all19F337substitutions.',
        mutation_readout='Mutated residue excluded from per-residue comparisons and output-based interactions.',
        hchi_mask='A/G have no chi; P flagged near-constant, retained raw but not primary Hchi plots.',
        boundary_control='Report core <=390 separately from391-395; compare d7 and d8.',
        epistemic_limits=['Fixed backbones do not sample relaxation or dynamics.',
            'Hchi != S2 != folding/binding free energy.',
            'Five trained models are model spread, not independent biological replicates.',
            'Static heavy-atom contacts support structural plausibility, not a dynamic transmission pathway.',
            'Learned AA embeddings plus backbone geometry; no explicit side-chain geometry/energies.']))
    ids, seq = residues(atoms('5HEB'),'A')
    amap = dict(zip(map(int,ids),seq))
    alanines = [f'{amap[i]}{i}A' for i in range(312,378) if amap[i] != 'A']
    for bg in ('apo1bfe','strip5heb'):
        for ctx in ('full','d7','d8'):
            build(bg,ctx)
            for m in keymut:
                build(bg,ctx,(m,))
        for ctx in ('full','d8'):
            for m in alanines:
                build(bg,ctx,(m,))
            for aa in AA3:
                if aa!='F':
                    build(bg,ctx,(f'F337{aa}',))
        for m in helixmut:
            build(bg,'full',(m,))
        for ms in doubles:
            build(bg,'full',ms)
    for ctx in ('full','d7','d8'):
        build('strip1be9',ctx)
    for ctx in ('dN','dNd8'):
        for m in ('','F337A','V328A','V328I','R312A','D357A','H372A'):
            build('strip5heb',ctx,(m,) if m else ())
    dump(BASE/'states.json',STATES)
    print(json.dumps({'states':len(STATES),'checkpoint_evaluations':5*len(STATES),'no_ligand':True}))

if __name__=='__main__':
    main()
