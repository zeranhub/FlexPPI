"""Physical sanity checks for the exact SASA implementation; synthetic atoms only."""
import importlib.util,argparse,math,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--runtime-path',action='append',default=[]);a=p.parse_args()
spec=importlib.util.spec_from_file_location('features_static_test',HERE/'extract_features.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.dependencies(a.runtime_path)
np=m.np
one={('A','1'):{'aa':'A','resname':'ALA','atoms':{'CA':{'xyz':np.zeros(3),'element':'C'}}}}
area=m.sasa(one,{'A'});expected=4*math.pi*(1.7+1.4)**2
assert abs(area-expected)<1e-10,(area,expected)
# Widely separated components have exactly zero buried surface under this rule.
two={**one,('D','1'):{'aa':'A','resname':'ALA','atoms':{'CA':{'xyz':np.array([100.,0.,0.]),'element':'C'}}}}
bsa=(m.sasa(two,{'A'})+m.sasa(two,{'D'})-m.sasa(two,{'A','D'}))/2
assert abs(bsa)<1e-10,bsa
# Two atom-level contacts from the same two residues count as one pair.
roles={'peptide':'C','tcr_alpha':'D','tcr_beta':'E'}
contact={('C','1'):{'aa':'A','resname':'ALA','atoms':{'N':{'xyz':np.zeros(3),'element':'N'},'CA':{'xyz':np.array([1.,0.,0.]),'element':'C'}}},('D','5'):{'aa':'A','resname':'ALA','atoms':{'N':{'xyz':np.array([2.,0.,0.]),'element':'N'},'CA':{'xyz':np.array([3.,0.,0.]),'element':'C'}}}}
assert len(m.contact_pairs(contact,roles))==1
report={'status':'passed','isolated_carbon_SASA_A2':area,'analytic_SASA_A2':expected,'noninteracting_BSA_A2':bsa,'multiple_atom_contacts_unique_residue_pair_count':1,'probe_A':1.4,'sphere_points':960,'no_real_candidate_or_label_data_used':True}
m.dump(HERE/'static_calculation_selftest.json',report);print(json.dumps(report,indent=2))
