"""Independent standard-AA and exact real-input parsing checks; no predictions."""
from pathlib import Path
import argparse,importlib.util,json,runpy,hashlib,datetime

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
parser=argparse.ArgumentParser()
parser.add_argument('--runtime-path',action='append',default=[])
parser.add_argument('--static-manifest',default=str(ROOT/'structures/static_complete_inputs.json'))
parser.add_argument('--all-inputs',action='store_true')
args=parser.parse_args()
spec=importlib.util.spec_from_file_location('feature_impl',HERE/'extract_features.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
mod.dependencies(args.runtime_path)
source=Path('E:/Project/FlexPPI/src/flexppi/constants.py')
constant=runpy.run_path(str(source))
standard={k:v for k,v in constant['THREE_TO_ONE'].items() if k in {'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'}}
assert len(mod.ONE_TO_THREE)==len(mod.THREE_TO_ONE)==20
assert mod.THREE_TO_ONE==standard
assert {mod.ONE_TO_THREE[a] for a in constant['RESTYPES']}==set(standard)
assert set(mod.CANONICAL)==set(constant['RESTYPES'])
for aa,name in mod.ONE_TO_THREE.items():
    assert mod.CHI_COUNTS[aa]==len(constant['CHI_ATOMS'].get(name,[]))
# Import just structures.py independently: it has no model/checkpoint dependency.
import sys,types
# The package __init__ eagerly imports Torch. Load the unchanged parser source
# with only its real constants dependency so this check needs no model runtime.
pkg=types.ModuleType('flexppi');pkg.__path__=[];sys.modules['flexppi']=pkg
cm=types.ModuleType('flexppi.constants');cm.__dict__.update(constant);sys.modules['flexppi.constants']=cm
parser_path=Path('E:/Project/FlexPPI/src/flexppi/data/structures.py')
parser_spec=importlib.util.spec_from_file_location('official_structure_parser',parser_path)
parser_module=importlib.util.module_from_spec(parser_spec);sys.modules['official_structure_parser']=parser_module
parser_spec.loader.exec_module(parser_module)
parse_pdb=parser_module.parse_pdb
fixed=mod.records(mod.read_json(ROOT/'structures/fixed_backbone_manifest.json'))
static={mod.input_id(r):r for r in mod.records(mod.read_json(args.static_manifest))}
checks=[]
for r in fixed:
    if not args.all_inputs and not mod.true(r.get('is_panel_parent')):continue
    p=mod.input_path(r);s=mod.static_path(static[mod.input_id(r)])
    assert s.exists(),('Complete parent file not ready',s)
    f,fi=mod.read_structure(p);c,ci=mod.read_structure(s)
    authoritative=parse_pdb(str(p))
    expected=dict(zip(zip(authoritative.chain_ids,authoritative.author_residue_ids),authoritative.sequence))
    assert list(f)==list(expected)
    assert {k:v['aa'] for k,v in f.items()}==expected
    assert set(c)==set(expected)
    assert {k:v['aa'] for k,v in c.items()}==expected
    roles=mod.role_map(r)
    assert mod.sequence(f,roles['peptide'])==r['peptide_sequence']==mod.sequence(c,roles['peptide'])
    assert not ci,('Unexpected nonprotein heavy atoms in complete protein structure',r['template_pdb'],len(ci),ci[:3])
    compare=mod.backbone_comparison(f,c);assert compare['passed']
    missing=[k for k,v in c.items() if not mod.complete(v)]
    assert not missing,(r['template_pdb'],'missing canonical residues',len(missing),missing[:3])
    checks.append({'template_pdb':r['template_pdb'],'parent_input_id':mod.input_id(r),'residues_checked':len(expected),'all_AA_and_author_ids_match_model_parser':True,'all_complete_PDB_residues_match_fixed_PDB':True,'all_backbone_coordinates_unchanged':True,'all_canonical_heavy_atom_sets_complete':True,'fixed_nonprotein_atoms_excluded_like_model_parser_count':len(fi),'fixed_nonprotein_atoms_first3':fi[:3],'fixed_sha256':mod.digest(p),'static_complete_sha256':mod.digest(s)})
assert len(checks)==(76 if args.all_inputs else 5)
report={'status':'passed','checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'standard_AA_mapping_checked':20,'canonical_chi_count_checked':20,'authoritative_constants_path':str(source),'authoritative_constants_sha256':mod.digest(source),'extract_script_sha256':mod.digest(HERE/'extract_features.py'),'scope':'all76' if args.all_inputs else 'five_template_parents','checked_inputs':checks,'model_predictions_or_SPR_read':False}
mod.dump(HERE/'input_contract_test.json',report)
print(json.dumps({'status':report['status'],'standard_AA_mapping_checked':20,'inputs_checked':len(checks),'residues_checked':sum(r['residues_checked'] for r in checks)},indent=2))
