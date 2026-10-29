"""Create explicit, coordinate-preserving mechanistic screens (not simulated apo).

Run with D:/Anaconda/envs/pvb_cpu/python.exe (Biopython is installed there).
All generated PDBs retain author numbering. Sequence-only reversion removes
incompatible side-chain atoms; no missing loops or atoms are fabricated.
"""
from pathlib import Path
import hashlib
import json
import urllib.request
import numpy as np
from Bio.Align import PairwiseAligner

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/mechanistic_case_studies_20260913'
SRC = ROOT / 'artifacts/apo_holo_mechanistic_literature_v2_20260913/sources'
AA3 = dict(zip('ARNDCQEGHILKMFPSTWYV', 'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()))
assert len(AA3) == 20 and AA3['V'] == 'VAL' and AA3['P'] == 'PRO'
AA1 = {v:k for k,v in AA3.items()}
STATES, CONTRASTS, AUDIT = [], [], []

def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding='utf-8')

def source(pdb):
    p = SRC / (pdb + '.pdb')
    if not p.exists():
        p = BASE / 'sources' / (pdb + '.pdb')
        p.parent.mkdir(parents=True, exist_ok=True)
        if not p.exists():
            data = urllib.request.urlopen('https://files.rcsb.org/download/' + pdb + '.pdb', timeout=45).read()
            if b'ATOM  ' not in data:
                raise ValueError('Not a PDB response')
            p.write_bytes(data)
    return p

def atoms(pdb, model=1):
    lines, active = [], model == 1
    for line in source(pdb).read_text().splitlines():
        if line.startswith('MODEL '):
            active = int(line[10:14]) == model
        elif line.startswith('ENDMDL') and active:
            break
        elif active and line.startswith('ATOM  ') and line[17:20] in AA1 and line[16] in ' A':
            lines.append(line)
    if not lines:
        raise ValueError((pdb, model, 'empty model'))
    return lines

def residues(lines, chain):
    d = {}
    for l in lines:
        if l[21] == chain:
            d.setdefault(l[22:27].strip(), AA1[l[17:20]])
    return list(d), ''.join(d.values())

def revert_map(lines, chain, reference, refchain):
    ids, seq = residues(lines, chain)
    rid, rseq = residues(reference, refchain)
    aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-3, extend_gap_score=-.5)
    aln = aligner.align(seq, rseq)[0]
    changes = {}
    for (s, e), (rs, re) in zip(*aln.aligned):
        for i, j in zip(range(s, e), range(rs, re)):
            if seq[i] != rseq[j]:
                changes[(chain, ids[i])] = rseq[j]
    return changes

def make(case, name, pdb, chains, ranges=None, changes=None, model=1, note=''):
    lines = atoms(pdb, model)
    ranges, changes = ranges or {}, changes or {}
    result, altered = [], set()
    for l in lines:
        ch, author, num = l[21], l[22:27].strip(), int(l[22:26])
        if ch not in chains:
            continue
        if ch in ranges and not any(a <= num <= b for a, b in ranges[ch]):
            continue
        key = (ch, author)
        if key in changes:
            if l[12:16].strip() not in ('N', 'CA', 'C', 'O', 'OXT'):
                continue
            altered.add((ch, author, l[17:20], AA3[changes[key]]))
            l = l[:17] + AA3[changes[key]] + l[20:]
        result.append(l)
    path = BASE / case / 'states' / (name + '.pdb')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(result) + '\nEND\n', encoding='ascii')
    stats = {ch:dict(n=len(residues(result,ch)[0]), sequence=residues(result,ch)[1], author_ids=residues(result,ch)[0]) for ch in chains}
    spec = dict(id=name, case=case, path=str(path), pdb_id=pdb, chains=list(chains), source_model=model,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(), note=note, residues=stats,
                sequence_only_changes=[list(x) for x in sorted(altered)])
    STATES.append(spec)
    return name

def comparison(case, name, before, after, chainmap, label, kind='fixed-coordinate context', regions=None):
    CONTRASTS.append(dict(id=name, case=case, before=before, after=after, chainmap=chainmap,
                         label=label, kind=kind, regions=regions or []))

def region(name, chain, low, high, expected=0):
    return dict(name=name, chain=chain, low=low, high=high, expected_rmsf_sign=expected)

def prepare():
    # Exact observed protein303-402; deletion removes396-402, not missing-coordinate 6QJD.
    for cut in (False, True):
        for bound in (False, True):
            name = 'pdz_' + ('d7' if cut else 'wt') + ('_bound' if bound else '_free')
            make('pdz3',name,'5HEB','AB' if bound else 'A', {'A':[(303,395 if cut else 402)]},
                 note='5HEB observed 7/9 CRIPT residues; free is partner-stripped, not experimental apo; deletion is coordinate-preserving.')
        make('pdz3','pdz_' + ('d7' if cut else 'wt') + '_apo1bfe','1BFE','A', {'A':[(306,395 if cut else 402)]},
             {('A','328'):'V'}, note='Independent apo sensitivity; I328V sequence-only reversion; missing303-305 not modeled.')
    preg = [region('beta2-beta3','A',318,323), region('alpha1 loop','A',330,336),region('core loop','A',341,356),region('shared core','A',303,395)]
    for name, before, after, label in [
        ('pdz_deletion_free','pdz_wt_free','pdz_d7_free','Deletion, free'),
        ('pdz_deletion_bound','pdz_wt_bound','pdz_d7_bound','Deletion, bound'),
        ('pdz_binding_wt','pdz_wt_free','pdz_wt_bound','Binding, WT'),
        ('pdz_binding_d7','pdz_d7_free','pdz_d7_bound','Binding, deletion'),
        ('pdz_deletion_apo','pdz_wt_apo1bfe','pdz_d7_apo1bfe','Deletion, independent apo')]:
        comparison('pdz3',name,before,after,{'A':'A'},label,regions=preg)
    # TCR: biological heterodimer only; revert4FTV sequence at aligned variable sites.
    rev = {}
    for ch, rc in [('D','A'),('E','B')]:
        rev.update(revert_map(atoms('4FTV'),ch,atoms('3QH3'),rc))
    for suffix, changes in [('wtseq',rev),('engineered',{})]:
        make('tcr','tcr_bound_'+suffix,'4FTV','ABCDE',changes=changes,
                 note='A6c134 crystallographic backbone; aligned TCR substitutions matched to3QH3 sequence, including engineered constant-domain cysteines. Not native wild-type.' if changes else 'A6c134 engineered affinity variant sensitivity.')
        make('tcr','tcr_strip_'+suffix,'4FTV','DE',changes=changes)
    make('tcr','tcr_apo_ab','3QH3','AB',note='Independent apo crystallographic heterodimer AB.')
    make('tcr','tcr_apo_cd','3QH3','CD',note='Independent apo crystallographic copy CD; missing residues remain missing.')
    treg = [region('V alpha','D',1,115),region('C alpha','D',116,206),region('V beta','E',1,115),region('C beta','E',116,246),region('C beta FG','E',218,231,-1)]
    comparison('tcr','tcr_context','tcr_strip_wtseq','tcr_bound_wtseq',{'D':'D','E':'E'},'Partner context',regions=treg)
    comparison('tcr','tcr_apo_holo','tcr_apo_ab','tcr_bound_wtseq',{'A':'D','B':'E'},'Apo to holo',kind='Independent crystals; holo backbone engineered',regions=treg)
    comparison('tcr','tcr_geometry','tcr_apo_ab','tcr_strip_wtseq',{'A':'D','B':'E'},'Backbone/construct',kind='Independent stripped crystals',regions=treg)
    comparison('tcr','tcr_copy_cd','tcr_apo_cd','tcr_bound_wtseq',{'C':'D','D':'E'},'Apo copy CD sensitivity',kind='Independent crystals; holo backbone engineered',regions=treg)
    comparison('tcr','tcr_engineered','tcr_strip_engineered','tcr_bound_engineered',{'D':'D','E':'E'},'Engineered sequence sensitivity',regions=treg)
    # 1DX5 contains four copies. A/M/I is one contacting biological complex.
    # MD sequential91 corresponds to PDB60I, NOT author91 (SEQADV validated).
    thchanges = {('I','364'):'N',('I','456'):'R',('I','457'):'H',('M','60I'):'T'}
    threg=[region('60s loop','M',60,64,1),region('Exosite I','M',70,80),region('170s loop','M',170,179),region('220s loop','M',220,225)]
    for name, ch, ran in [('th_free','AM',{}),('th_tm56','AMI',{'I':[(388,462)]}),('th_tm456','AMI',{})]:
        make('thrombin',name,'1DX5',ch,ran,thchanges,note='Protein-only input; same inhibitor-bound crystal backbone; inhibitor and Ca not encoded. TM4 author345-387.')
    for n,b,a,l in [('th_tm456_effect','th_free','th_tm456','TM456 binding'),('th_tm56_effect','th_free','th_tm56','TM56 binding'),('th_tm4_increment','th_tm56','th_tm456','Adding TM4')]:
        comparison('thrombin',n,b,a,{'A':'A','M':'M'},l,regions=threg)
    # All deposited NMR conformers, no cherry-picking of one model.
    models = [int(l[10:14]) for l in source('2AGH').read_text().splitlines() if l.startswith('MODEL ')] or [1]
    for m in models:
        for suffix,chs in [('free','B'),('myb','AB'),('mll','BC'),('both','ABC')]:
            make('kix',f'kix_m{m:02d}_{suffix}','2AGH',chs,ranges={'C':[(840,858)]},model=m,
                 note='Same deposited ternary NMR conformer; MLL840-858 matches Peiffer2022 peptide; stripped states are not measured apo ensembles.')
        for suffix,label in [('myb','cMyb binding'),('mll','MLL binding'),('both','Dual binding')]:
            comparison('kix',f'kix_m{m:02d}_{suffix}_effect',f'kix_m{m:02d}_free',f'kix_m{m:02d}_{suffix}',{'B':'B'},label,
                       regions=[region('L12 loop','B',614,621),region('all KIX','B',586,672)])
    make('kix','kix_m01_mll_full31','2AGH','BC',note='Full deposited MLL839-869 construct sensitivity.')
    comparison('kix','kix_m01_mll_full31_effect','kix_m01_free','kix_m01_mll_full31',{'B':'B'},'MLL31 sensitivity',regions=[region('L12 loop','B',614,621)])
    for b,ch in [(False,'A'),(True,'AC')]:
        make('giv','giv_'+('bound' if b else 'free'),'6MHF',ch,note='GDP not encoded; only18observed GIV residues; crystal context screen.')
    comparison('giv','giv_binding','giv_free','giv_bound',{'A':'A'},'GIV binding',regions=[region('Switch II','A',201,216,-1)])
    # BAP1 incomplete UCH crystal construct, covalent ubiquitin probe.
    for mut in (False,True):
        for bound in (False,True):
            name='bap_'+('l49v' if mut else 'wt')+('_bound' if bound else '_free')
            changes={('L','91'):'C'}
            if mut:
                changes[('L','49')]='V'
            make('bap1',name,'8H1T','LM' if bound else 'L',ranges={'L':[(3,238)]},changes=changes,
                 note='UCH3-238 from cryo-EM BAP1-ASXL1/nucleosome complex; C91S reverted by sequence only; L49V on same backbone. Covalent/chemical context is not encoded.')
    for suffix in ('free','bound'):
        comparison('bap1','bap_mutation_'+suffix,'bap_wt_'+suffix,'bap_l49v_'+suffix,{'L':'L'},'L49V, '+suffix,
                   kind='Sequence-only mutation',regions=[region('L49 neighborhood','L',45,53)])
    dump(BASE/'states.json',STATES)
    dump(BASE/'contrasts.json',CONTRASTS)
    for st in STATES:
        for c,v in st['residues'].items():
            assert v['n'] > 0, (st['id'],c)
    print(json.dumps(dict(states=len(STATES),contrasts=len(CONTRASTS),kix_models=len(models),tcr_reversions=[list(k)+[v] for k,v in rev.items()]),indent=2))

if __name__ == '__main__':
    prepare()
