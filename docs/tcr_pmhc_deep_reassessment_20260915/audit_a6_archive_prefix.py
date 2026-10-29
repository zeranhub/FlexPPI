"""Read only complete members from the existing, incomplete A6 download.

The partial archive is not checksum-verified. Nothing here is a full dataset audit.
No archive path is extracted, and archive scripts are never executed.
"""
from pathlib import Path, PurePosixPath
from collections import Counter
import hashlib
import json
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/tcr_pmhc_reassessment_20260915/knapp2019'
ARCHIVE = OUT / 'A6.tar.gz.partial'
assert ARCHIVE.resolve().parent == OUT.resolve()
members, pdbs = [], []
failure = None
try:
    with tarfile.open(ARCHIVE, 'r|gz') as archive:
        for member in archive:
            name = PurePosixPath(member.name)
            record = {'name': member.name, 'declared_size': member.size,
                      'regular_file': member.isfile(), 'payload_read': False}
            members.append(record)
            safe = not name.is_absolute() and '..' not in name.parts and '\\' not in member.name and ':' not in member.name
            if member.isfile() and safe and name.suffix == '.pdb' and member.size < 1_000_000:
                data = archive.extractfile(member).read(member.size + 1)
                if len(data) != member.size:
                    raise ValueError('Incomplete member; not written as PDB')
                record['payload_read'] = True
                target = OUT / ('prefix_firstframe_' + hashlib.sha256(member.name.encode()).hexdigest()[:10] + '.pdb')
                target.write_bytes(data)
                text = data.decode('ascii', errors='replace')
                ca = []
                for line in text.splitlines():
                    if line.startswith(('ATOM  ', 'HETATM')) and line[12:16].strip() == 'CA':
                        ca.append({'chain': line[21], 'resid': line[22:27].strip(), 'resname': line[17:20].strip()})
                pdbs.append({'member': member.name, 'local': str(target.relative_to(ROOT)),
                             'sha256': hashlib.sha256(data).hexdigest(), 'CA_count': len(ca),
                             'chain_counts': dict(Counter(r['chain'] for r in ca)), 'CA_residues': ca})
except (EOFError, tarfile.ReadError, OSError, ValueError) as error:
    failure = str(error)
table = json.loads((ROOT / 'artifacts/tcr_pmhc_reassessment_20260915/literature/knapp2019_s005_tables.json').read_text())
rmsf_rows = [r for t in table for r in t if any('RMSF' in v.upper() for v in r[:1])]
report = {'archive': str(ARCHIVE.relative_to(ROOT)), 'cached_bytes': ARCHIVE.stat().st_size,
          'expected_bytes': 399223244, 'checksum_verified': False,
          'scope': 'Partial archive prefix only. Member headers do not prove complete trajectory payloads.',
          'stop_reason': failure, 'seen_headers': members, 'complete_PDB_members': pdbs,
          'published_s2_RMSF_rows': rmsf_rows}
(OUT / 'prefix_inventory.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps({'bytes': report['cached_bytes'], 'stop_reason': failure,
                  'seen_headers': members, 'pdbs': [{k:v for k,v in p.items() if k != 'CA_residues'} for p in pdbs],
                  'RMSF_rows': rmsf_rows}, indent=2))
