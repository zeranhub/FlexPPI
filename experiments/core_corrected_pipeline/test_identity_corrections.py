"""Targeted regression tests for the actual observed data corruption mechanisms."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from prepare_corrected_labels import canonical_pdb, parse_pdb, replica_map

def line(atom,residue,number,x):
    return f"ATOM  {number:5d} {atom:>4s} {residue:3s} A{number:4d}    {x:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           C\n"

def test_his_alias_is_retained_and_ile_cd_coordinates_preserved(tmp_path):
    text=line("CA","HSD",1,1.2)+line("CA","ILE",2,3.4)+line("CD","ILE",2,5.6)
    p,q=tmp_path/"old.pdb",tmp_path/"corrected.pdb"
    p.write_text(text)
    canonical_pdb(p,q)
    s=parse_pdb(str(q))
    assert s.sequence=="HI" and s.author_residue_ids==["1","2"]
    corrected=q.read_text().splitlines()
    assert corrected[-1][12:16].strip()=="CD1"
    assert corrected[-1][30:54]==text.splitlines()[-1][30:54]
    assert p.read_text()==text

def test_conflicting_ile_terminal_atom_names_do_not_silently_choose(tmp_path):
    p,q=tmp_path/"ambiguous.pdb",tmp_path/"output.pdb"
    p.write_text(line("CA","ILE",2,1)+line("CD","ILE",2,2)+line("CD1","ILE",2,3))
    with pytest.raises(ValueError,match="Ambiguous"):
        canonical_pdb(p,q)

class ToyStructure:
    def __init__(self,chains):
        self.chains=chains
    def unique_chains(self):return list(self.chains)
    def chain_sequence(self,c):return self.chains[c]
    def indices_for_chain(self,c):
        keys=list(self.chains)
        start=sum(len(self.chains[k]) for k in keys[:keys.index(c)])
        return np.arange(start,start+len(self.chains[c]))

def test_swapped_chain_names_follow_sequence_identity_not_letter():
    a=ToyStructure({"C":"ITDQVPFSV","D":"CASSLGQETQYF"})
    b=ToyStructure({"D":"ITDQVPFSV","C":"CASSLGQETQYF"})
    mapping,stats=replica_map(a,b)
    assert {(r["canonical_chain"],r["replica_chain"]) for r in stats}=={("C","D"),("D","C")}
    assert len(mapping)==21 and len({p[0] for p in mapping})==21

def test_sequence_substitution_is_rejected_as_replica_identity():
    a=ToyStructure({"A":"CASSLGQETQYF"})
    b=ToyStructure({"A":"CASSAGQETQYF"})
    with pytest.raises(ValueError,match="exact-identity"):
        replica_map(a,b)
