"""Read-only RCSB discovery and coordinate audit for the PDZ3 alpha3 study.

This downloads original database files only. It creates no model input states.
"""
import concurrent.futures
import datetime
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/pdz3_alpha3_focus_20260914'
SRC = OUT / 'sources'


def get_json(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'FlexPPI-structure-audit/1.0'})
    with urllib.request.urlopen(req, timeout=45) as response:
        return json.load(response)


def one_entity(item):
    identifier = item['identifier']
    entry_id, entity_id = identifier.split('_')
    entity = get_json(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}')
    entry = get_json(f'https://data.rcsb.org/rest/v1/core/entry/{entry_id}')
    return {
        'identifier': identifier,
        'title': entry.get('struct', {}).get('title'),
        'method': entry.get('exptl'),
        'resolution': entry.get('rcsb_entry_info', {}).get('resolution_combined'),
        'length': entity.get('entity_poly', {}).get('rcsb_sample_sequence_length'),
        'sequence': entity.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can'),
        'chains': entity.get('rcsb_polymer_entity_container_identifiers', {}).get('auth_asym_ids'),
        'alignments': entity.get('rcsb_polymer_entity_align'),
        'mutations': entity.get('rcsb_polymer_entity', {}).get('pdbx_mutation'),
        'entity': entity,
        'entry': entry,
    }


def cif_audit(pdb_id):
    paths = []
    for extension in ('pdb', 'cif'):
        url = f'https://files.rcsb.org/download/{pdb_id}.{extension}'
        path = SRC / f'{pdb_id}.{extension}'
        if not path.exists():
            with urllib.request.urlopen(url, timeout=45) as response:
                path.write_bytes(response.read())
        paths.append({'path': str(path), 'source_url': url, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    cif = MMCIF2Dict(str(SRC / f'{pdb_id}.cif'))
    def rows(prefix):
        keys = [key for key in cif if key.startswith(prefix + '.')]
        if not keys:
            return []
        return [dict(zip(keys, vals)) for vals in zip(*(cif[key] for key in keys))]
    polymer_chains = {}
    for atom in rows('_atom_site'):
        if atom['_atom_site.label_atom_id'] != 'CA' or atom['_atom_site.label_seq_id'] == '.':
            continue
        chain = atom['_atom_site.auth_asym_id']
        key = (atom['_atom_site.auth_seq_id'], atom['_atom_site.pdbx_PDB_ins_code'])
        polymer_chains.setdefault(chain, {})[key] = {
            'author_residue_number': atom['_atom_site.auth_seq_id'],
            'insertion_code': atom['_atom_site.pdbx_PDB_ins_code'],
            'resname': atom['_atom_site.auth_comp_id'],
            'label_seq_id': atom['_atom_site.label_seq_id'],
        }
    return {
        'pdb_id': pdb_id,
        'files': paths,
        'coordinate_chains': {chain: list(residues.values()) for chain, residues in polymer_chains.items()},
        'entity_poly': rows('_entity_poly'),
        'missing_residues': rows('_pdbx_unobs_or_zero_occ_residues'),
        'modified_residues': rows('_pdbx_struct_mod_residue'),
        'sequence_differences': rows('_struct_ref_seq_dif'),
        'nonpolymer': rows('_pdbx_nonpoly_scheme'),
        'database_sequence_mapping': rows('_struct_ref_seq'),
    }


def main():
    SRC.mkdir(parents=True, exist_ok=True)
    query = {
        'query': {'type': 'terminal', 'service': 'text', 'parameters': {
            'attribute': 'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
            'operator': 'in', 'value': ['P31016', 'P78352']}},
        'return_type': 'polymer_entity',
        'request_options': {'paginate': {'start': 0, 'rows': 500}},
    }
    response = get_json('https://search.rcsb.org/rcsbsearch/v2/query?json=' + urllib.parse.quote(json.dumps(query)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(one_entity, response.get('result_set', [])))
    (OUT / 'psd95_rcsb_inventory.json').write_text(json.dumps({
        'accessed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope': 'UniProt P31016 rat and P78352 human annotated polymer entities; not an exhaustive sequence-only homology search',
        'query': query, 'response': response, 'entities': results,
    }, indent=2), encoding='utf-8')
    for entity in results:
        print(entity['identifier'], entity['length'], entity['method'], entity['title'])
    audits = [cif_audit(pdb_id) for pdb_id in ('1BFE', '3I4W', '3K82')]
    (OUT / 'pdz3_experimental_structure_audit.json').write_text(json.dumps(audits, indent=2), encoding='utf-8')
    for audit in audits:
        print('AUDIT', audit['pdb_id'], 'MOD', audit['modified_residues'], 'DIFF', audit['sequence_differences'])
        for chain, residues in audit['coordinate_chains'].items():
            print('CHAIN', chain, 'N', len(residues), 'FIRST', residues[0], 'LAST', residues[-1])


if __name__ == '__main__':
    main()
