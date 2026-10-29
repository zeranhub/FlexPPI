"""Archive bounded primary-source evidence and inspect real structural controls.

This creates a research record, not predictions or simulated experimental data.
Run with .venv-md/Scripts/python.exe from the FlexPPI project root.
"""
from pathlib import Path
import hashlib
import io
import json
import re
import urllib.request
import urllib.parse
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'artifacts/tcr_pmhc_reassessment_20260915'
DEST = BASE / 'literature'
DEST.mkdir(parents=True, exist_ok=True)
RECORDS = []


def fetch(url, name):
    path = DEST / name
    if path.exists():
        data = path.read_bytes()
    else:
        req = urllib.request.Request(url, headers={'User-Agent': 'FlexPPI-research-audit/1.0'})
        with urllib.request.urlopen(req, timeout=40) as response:
            data = response.read(12_000_001)
        if len(data) > 12_000_000:
            raise ValueError('Bounded retrieval exceeded 12 MB; no partial file saved')
        path.write_bytes(data)
    RECORDS.append({'url': url, 'file': str(path.relative_to(ROOT)),
                    'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
    return data


PAPERS = {
    'rangarajan2018': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6187629/',
    'he2020': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8356565/',
    'natarajan2017': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5440810/',
    'hawse2012': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3375328/',
    'chen2026': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC13004167/',
    'notti2025': 'https://www.nature.com/articles/s41467-025-66939-7',
}
errors = []
for name, url in PAPERS.items():
    try:
        data = fetch(url, name + '.html')
        soup = BeautifulSoup(data, 'html.parser')
        article = soup.find('article') or soup.find('main') or soup
        for junk in article.find_all(['script', 'style', 'nav']):
            junk.decompose()
        text = article.get_text('\n', strip=True)
        if len(text) < 3000:
            raise ValueError('Article response too short; may be an access-check page')
        (DEST / (name + '.txt')).write_text(text, encoding='utf-8')
        links = sorted(set(urllib.parse.urljoin(url, a['href'])
                           for a in soup.find_all('a', href=True)
                           if re.search(r'\.(docx|xlsx)(\?|$)', a['href'])))
        (DEST / (name + '_data_links.json')).write_text(json.dumps(links, indent=2), encoding='utf-8')
        # Small native data attachments are inspected; large MD archives are not fetched here.
        for index, link in enumerate(links):
            extension = '.docx' if '.docx' in link else '.xlsx'
            try:
                prior_access_response = DEST / (name + '_attachment_' + str(index) + '.access_response.html')
                if prior_access_response.exists():
                    raise ValueError('Known access-check response archived; no automatic repeated request')
                attachment = fetch(link, name + '_attachment_' + str(index) + extension)
                if not attachment.startswith(b'PK'):
                    invalid = DEST / (name + '_attachment_' + str(index) + extension)
                    archived = invalid.with_suffix('.access_response.html')
                    if archived.exists():
                        raise ValueError('Previously archived access-check response exists; inspect rather than overwrite')
                    invalid.rename(archived)
                    RECORDS[-1]['file'] = str(archived.relative_to(ROOT))
                    RECORDS[-1]['validation'] = 'invalid_office_container_access_response'
                    raise ValueError('Attachment is not a valid Office ZIP container')
                if extension == '.docx':
                    with zipfile.ZipFile(io.BytesIO(attachment)) as archive:
                        root = ET.fromstring(archive.read('word/document.xml'))
                    ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
                    paragraphs = [''.join(t.text or '' for t in p.findall('.//w:t', ns))
                                  for p in root.findall('.//w:p', ns)]
                    (DEST / (name + '_attachment_' + str(index) + '.txt')).write_text('\n'.join(paragraphs), encoding='utf-8')
                    tables = []
                    for table in root.findall('.//w:tbl', ns):
                        tables.append([[''.join(t.text or '' for t in cell.findall('.//w:t', ns))
                                        for cell in row.findall('w:tc', ns)]
                                       for row in table.findall('w:tr', ns)])
                    (DEST / (name + '_attachment_' + str(index) + '_tables.json')).write_text(json.dumps(tables, indent=2), encoding='utf-8')
                else:
                    with zipfile.ZipFile(io.BytesIO(attachment)) as archive:
                        ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
                        shared = []
                        if 'xl/sharedStrings.xml' in archive.namelist():
                            strings = ET.fromstring(archive.read('xl/sharedStrings.xml'))
                            shared = [''.join(t.text or '' for t in si.findall('.//s:t', ns))
                                      for si in strings.findall('s:si', ns)]
                        summary = {}
                        for member in archive.namelist():
                            if not re.match(r'xl/worksheets/sheet\d+\.xml$', member):
                                continue
                            sheet = ET.fromstring(archive.read(member))
                            rows = []
                            for row in sheet.findall('.//s:row', ns):
                                cells = {}
                                for cell in row.findall('s:c', ns):
                                    value = cell.find('s:v', ns)
                                    value = value.text if value is not None else ''.join(t.text or '' for t in cell.findall('.//s:t', ns))
                                    if cell.get('t') == 's' and value:
                                        value = shared[int(value)]
                                    cells[cell.get('r')] = value
                                rows.append(cells)
                            summary[member] = {'rows': len(rows), 'first_rows': rows[:15]}
                    (DEST / (name + '_attachment_' + str(index) + '_inventory.json')).write_text(json.dumps(summary, indent=2, default=str), encoding='utf-8')
            except Exception as error:
                errors.append({'url': link, 'error': str(error)})
        print(name, 'readable characters', len(text), 'data attachments', len(links), flush=True)
    except Exception as error:
        errors.append({'url': url, 'error': str(error)})
        print(name, 'retrieval failed:', str(error), flush=True)

for url, name in [
    ('https://zenodo.org/api/records/17065125', 'chen2026_zenodo_record.json'),
    ('https://zenodo.org/api/records/17065125/files/readme.txt/content', 'chen2026_zenodo_readme.txt'),
]:
    try:
        fetch(url, name)
    except Exception as error:
        errors.append({'url': url, 'error': str(error)})

(DEST / 'retrieval_manifest.json').write_text(json.dumps({'sources': RECORDS, 'errors': errors}, indent=2), encoding='utf-8')


def fit_rmsd(left, right):
    left = np.array(left)
    right = np.array(right)
    a = left - left.mean(axis=0)
    b = right - right.mean(axis=0)
    u, _, vh = np.linalg.svd(a.T @ b)
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vh)
    moved = a @ (u @ correction @ vh)
    return float(np.sqrt(np.mean(np.sum((moved-b)**2, axis=1))))


structures = json.loads((BASE / 'structure_sequence_audit.json').read_text())
codes = ['8TW4', '8TW6', '9BBC', '9C3E']
pairs = []
for chain in ['A', 'B']:
    common = set.intersection(*[set(structures[code]['chains'][chain]) for code in codes])
    for key in sorted(common, key=int):
        residues = [structures[code]['chains'][chain][key] for code in codes]
        complete = all(all(atom in structures[code]['backbone'][chain][key]
                           for atom in ['N', 'CA', 'C', 'O']) for code in codes)
        if len(set(residues)) == 1 and complete:
            pairs.append((chain, key))
results = []
for before in codes[:-1]:
    left = [structures[before]['backbone'][chain][key]['CA'] for chain, key in pairs]
    right = [structures['9C3E']['backbone'][chain][key]['CA'] for chain, key in pairs]
    result = {'before': before, 'after': '9C3E', 'n_common_all_four_structures': len(pairs),
              'fitted_CA_RMSD_A': fit_rmsd(left, right)}
    results.append(result)
    print('ALL-FOUR common coordinates', json.dumps(result), flush=True)
assert len(pairs) > 300
(BASE / 'notti_all_four_common_geometry.json').write_text(json.dumps({
    'description': 'Same exact residues, complete N/CA/C/O in all four structures; joint TCR alpha-beta CA rigid alignment. Static RMSD, not RMSF.',
    'positions': pairs, 'results': results}, indent=2), encoding='utf-8')
