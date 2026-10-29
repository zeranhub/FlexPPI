"""Download only small public catalogues and audit exact-PDB source overlap."""
from pathlib import Path
import json
import hashlib
import datetime
import math
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/tcr_external_validation_search_20260917'
DEST = OUT / 'catalogues'
COMMIT = 'fc2837a91e3978832f9e74c7483d7ed7c5ff8b91'


def fetch(url, path, **kwargs):
    response = requests.get(url, timeout=45, **kwargs)
    response.raise_for_status()
    path.write_bytes(response.content)
    return response


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    raw = f'https://raw.githubusercontent.com/leifugroup/DynaTPH/{COMMIT}'
    for filename in ['descriptor.csv', 'rmsf.csv', 'README.md', 'Dynamic_Data/README.md']:
        fetch(raw + '/' + filename, DEST / ('DynaTPH_' + filename.replace('/', '_')))
    # PDB IDs and header are ASCII; other numerical text in this CSV uses
    # a non-UTF8 encoding, so do not interpret those columns in this audit.
    descriptor = pd.read_csv(DEST / 'DynaTPH_descriptor.csv', encoding='latin1')
    ids = set(descriptor.PDBID.str.upper())
    assert len(ids) == len(descriptor) == 256
    registry = pd.read_csv(OUT / 'overlap_agent/pdb_exclusion_registry.csv').set_index('pdb_id')

    # The API defaults to pagination. Never call a partial page the full catalogue.
    response = fetch('https://dynarepo.inria.fr/api/rest/current/projects',
                     DEST / 'dynarepo_live_page_1.json', params={'limit': 100, 'page': 1})
    live = response.json()
    live_projects = live['projects']
    def get_page(page):
        return fetch('https://dynarepo.inria.fr/api/rest/current/projects',
                     DEST / f'dynarepo_live_page_{page}.json', params={'limit': 100, 'page': page}).json()
    with ThreadPoolExecutor(max_workers=3) as pool:
        for page in pool.map(get_page, range(2, math.ceil(live['filteredCount']/100)+1)):
            assert page['filteredCount'] == live['filteredCount'], 'Catalogue changed while paginating'
            live_projects.extend(page['projects'])
    live_complete = len({r['accession'] for r in live_projects}) == len(live_projects) == live['filteredCount']
    assert live_complete, 'Incomplete or duplicate pagination; do not claim a complete live audit'
    (DEST / 'dynarepo_live_projects.json').write_text(json.dumps(live), encoding='utf-8')
    live_ids = {p.upper() for r in live_projects for p in r['metadata'].get('PDBIDS', [])}
    local_ids = set(registry[registry.in_dynarepo_source].index)
    rows = []
    for pdb in sorted(ids):
        rec = dict(pdb_id=pdb)
        for column in ['in_itrace_source', 'in_dynarepo_source', 'in_atlas_source', 'in_actual_training_pool']:
            rec[column] = bool(registry.at[pdb, column]) if pdb in registry.index else False
        rec['in_dynarepo_live_returned_page'] = pdb in live_ids
        rec['passes_local_pdb_exclusion'] = not any(rec[c] for c in
            ['in_itrace_source', 'in_dynarepo_source', 'in_atlas_source', 'in_actual_training_pool'])
        rec['passes_pdb_exclusion_including_live'] = (rec['passes_local_pdb_exclusion'] and
            pdb not in live_ids) if live_complete else None
        rows.append(rec)
    table = pd.DataFrame(rows)
    table.to_csv(DEST / 'DynaTPH_pdb_overlap.csv', index=False)
    candidates = pd.read_csv(OUT / 'overlap_agent/candidate_pdb_screen.csv')
    candidates['in_dynarepo_live'] = candidates.pdb_id.isin(live_ids) if live_complete else None
    candidates.to_csv(DEST / 'individual_candidates_with_live_check.csv', index=False)
    summary = dict(retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        DynaTPH_commit=COMMIT, DynaTPH_systems=len(ids),
        DynaTPH_local_overlap_counts={c: int(table[c].sum()) for c in table.columns if c.startswith('in_')},
        DynaTPH_passes_local_pdb_exclusion=int(table.passes_local_pdb_exclusion.sum()),
        DynaTPH_passes_pdb_exclusion_including_live=int(table.passes_pdb_exclusion_including_live.sum()) if live_complete else None,
        dynarepo_live_complete=live_complete, dynarepo_live_reported_projects=live['filteredCount'],
        dynarepo_live_returned_projects=len(live_projects), dynarepo_live_returned_unique_pdb=len(live_ids),
        dynarepo_live_added_pdb=sorted(live_ids-local_ids) if live_complete else None,
        dynarepo_live_removed_pdb=sorted(local_ids-live_ids) if live_complete else None,
        sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in DEST.iterdir() if p.is_file()},
        limitations=['Exact PDB exclusion is not clonotype or CDR3 independence.',
          'DynaTPH full trajectories are restricted for peer review per the pinned author README.',
          'The public RMSF table describes peptide residues, not TCR CDR residues.',
          'ITrace screening uses the complete local source manifest, not a verified live inventory.'])
    (DEST / 'catalogue_audit.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='sha256'}, indent=2))


if __name__ == '__main__':
    main()
