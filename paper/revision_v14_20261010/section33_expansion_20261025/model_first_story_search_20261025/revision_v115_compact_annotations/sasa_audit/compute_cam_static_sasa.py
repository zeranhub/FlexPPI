"""Static solvent exposure with the original ubiquitin PyMOL 3.0 protocol.

Writes only beside this script; does not edit source structures or run a model.
"""
from pathlib import Path
import csv
import json
import hashlib
import shutil
import sys
import pymol
from pymol import cmd

OUT = Path(__file__).resolve().parent
S = OUT.parents[1]
CAM_SOURCE = Path(r'E:\Project\FlexPPI\experiments\_sources\baselines\cohort\structures\record_aba5e1f3db118b0aba8895de.pdb')
CAM_TABLE = S / 'revision_v108_dual_axis_case_profiles/figure3/source_data/cam_all25.csv'
UB_SOURCE = S / 'revision_v80_redesign/pymol/inputs/A003R_input_structure.pdb'
UB_TABLE = S / 'revision_v93_dual_case_entropy/ubiquitin_audit/all119_input_environment.csv'
BB = {'N', 'CA', 'C', 'O', 'OXT', 'OT1', 'OT2'}
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
AA.update({'HSE': 'H', 'HSD': 'H', 'HSP': 'H', 'HIE': 'H', 'HID': 'H', 'HIP': 'H'})


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(name, rows):
    with (OUT / name).open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(path, name):
    cmd.load(str(path), name)
    cmd.remove(name + ' and hydro')
    cmd.set('dot_solvent', 1)
    cmd.set('solvent_radius', 1.4)
    cmd.set('dot_density', 4)
    total = cmd.get_area(name, load_b=1)
    atoms = cmd.get_model(name).atom
    groups = {}
    atom_rows = []
    for atom in atoms:
        groups.setdefault((atom.chain, int(atom.resi), atom.resn), []).append(atom)
        atom_rows.append(dict(chain=atom.chain, input_residue_id=int(atom.resi),
                              resname=atom.resn, atom_name=atom.name, element=atom.symbol,
                              x_A=atom.coord[0], y_A=atom.coord[1], z_A=atom.coord[2],
                              atom_SASA_A2=atom.b, is_sidechain_heavy=atom.name not in BB))
    residues = []
    for (chain, residue, resname), group in groups.items():
        side = [atom for atom in group if atom.name not in BB]
        residues.append(dict(chain=chain, input_residue_id=residue, resname=resname,
                             amino_acid=AA[resname], whole_residue_SASA_A2=sum(atom.b for atom in group),
                             sidechain_SASA_A2=sum(atom.b for atom in side),
                             n_sidechain_heavy_atoms=len(side)))
    return total, atom_rows, residues


pymol.finish_launching(['pymol', '-cq'])
assert cmd.get_version()[0] == '3.0.0', cmd.get_version()
source_hash = sha(CAM_SOURCE)
frozen = OUT / 'CaMKI_exact_input_complex.pdb'
shutil.copy2(CAM_SOURCE, frozen)
assert sha(frozen) == source_hash
with CAM_TABLE.open(encoding='utf-8-sig', newline='') as handle:
    cam_reference = list(csv.DictReader(handle))
assert len(cam_reference) == 25
total, atom_rows, residues = analyze(frozen, 'cam_complex')
selected = [dict(row) for row in residues if row['chain'] == 'B']
assert len(selected) == 25
selected.sort(key=lambda row: row['input_residue_id'])
for row, reference in zip(selected, cam_reference):
    native = row['input_residue_id'] + 149
    assert native == int(reference['CaMKI_native_position'])
    assert row['amino_acid'] == reference['amino_acid']
    row['CaMKI_native_position'] = native
    row['sample_id'] = 'A001K'
    row['pdb_id'] = '1MXE'
assert [row['CaMKI_native_position'] for row in selected] == list(range(294, 319))
write_csv('CaMKI_all25_PyMOL_sidechain_SASA.csv', selected)
write_csv('CaMKI_all169_complex_residue_SASA.csv', residues)
write_csv('CaMKI_all_complex_atom_SASA.csv', atom_rows)

# Check the same runtime/protocol against the frozen ubiquitin values without
# changing those data or their source structure.
ub_total, ub_atoms, ub_rows = analyze(UB_SOURCE, 'ub_complex')
with UB_TABLE.open(encoding='utf-8-sig', newline='') as handle:
    ub_reference = list(csv.DictReader(handle))
lookup = {(row['chain'], int(row['input_residue_id'])): row for row in ub_reference}
assert len(lookup) == len(ub_rows) == 119
differences = []
for row in ub_rows:
    ref = lookup[(row['chain'], row['input_residue_id'])]
    for key in ['whole_residue_SASA_A2', 'sidechain_SASA_A2']:
        differences.append(abs(row[key] - float(ref[key])))
max_difference = max(differences)
assert max_difference <= 1e-9, max_difference

manifest = dict(
    status='passed', geometry_only=True, no_model_or_MD_run=True,
    runtime_executable=sys.executable, pymol_version=list(cmd.get_version()),
    source_structure=str(CAM_SOURCE), source_sha256=source_hash,
    frozen_structure=str(frozen), frozen_sha256=sha(frozen),
    script_sha256=sha(__file__),
    source_identity_table=str(CAM_TABLE), source_identity_table_sha256=sha(CAM_TABLE),
    n_heavy_atoms=len(atom_rows), n_complex_residues=len(residues), n_peptide_positions=len(selected),
    total_complex_SASA_A2=total, probe_radius_A=1.4, dot_solvent=1, dot_density=4,
    sidechain_definition='All residue heavy atoms except N,CA,C,O,OXT,OT1,OT2; all complex heavy atoms provide shielding',
    residue_mapping='chain B author 145..169 -> native CaMKI 294..318 by author+149, all 25 amino acids matched to existing case table',
    sidechain_SASA_min_A2=min(row['sidechain_SASA_A2'] for row in selected),
    sidechain_SASA_max_A2=max(row['sidechain_SASA_A2'] for row in selected),
    all_positions_including_Ala_retained=True,
    original_sources_unchanged=(sha(CAM_SOURCE) == source_hash),
    ubiquitin_protocol_validation=dict(source_structure=str(UB_SOURCE), source_sha256=sha(UB_SOURCE),
                                      reference_table=str(UB_TABLE), reference_sha256=sha(UB_TABLE),
                                      n_residues=119, n_ubiquitin_positions=72,
                                      max_absolute_SASA_difference_A2=max_difference,
                                      total_complex_SASA_A2=ub_total),
    scope='Exact single starting complex, static solvent exposure in Angstrom squared; not MD averaged or a model-predicted quantity')
(OUT / 'CaMKI_static_SASA_provenance.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(manifest, ensure_ascii=False))
cmd.quit()
