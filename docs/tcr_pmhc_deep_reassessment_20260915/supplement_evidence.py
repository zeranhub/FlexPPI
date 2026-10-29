"""Collect small original-source metadata and audit cached evidence.

This script does not download large MD archives or calculate model predictions.
Use the existing .venv-md Python interpreter from the project root.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import io
import json
import urllib.request
import urllib.parse
import zipfile
import xml.etree.ElementTree as ET

from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'artifacts/tcr_pmhc_reassessment_20260915/literature'


def get_small(url):
    request = urllib.request.Request(url, headers={'User-Agent': 'FlexPPI-evidence-audit/1.0'})
    with urllib.request.urlopen(request, timeout=40) as response:
        data = response.read(12_000_001)
    if len(data) > 12_000_000:
        raise ValueError('Response exceeds 12 MB; not saved')
    return data


def figshare(item):
    label, article_id = item
    url = 'https://api.figshare.com/v2/articles/' + str(article_id)
    data = get_small(url)
    record = json.loads(data)
    path = DEST / ('knapp2019_' + label + '_figshare.json')
    path.write_bytes(data)
    return {'system': label, 'url': url, 'file': str(path.relative_to(ROOT)),
            'sha256': hashlib.sha256(data).hexdigest(), 'doi': record['doi'],
            'files': record['files']}


with ThreadPoolExecutor(max_workers=4) as pool:
    records = list(pool.map(figshare, [('A6', 8067827), ('1G4', 8067743),
                                     ('LC13', 8067746), ('JM22', 8067803)]))

url = 'https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007338'
data = get_small(url)
(DEST / 'knapp2019.html').write_bytes(data)
soup = BeautifulSoup(data, 'html.parser')
article = soup.find('article') or soup.find('main') or soup
for junk in article.find_all(['script', 'style', 'nav']):
    junk.decompose()
(DEST / 'knapp2019.txt').write_text(article.get_text('\n', strip=True), encoding='utf-8')
links = sorted(set(urllib.parse.urljoin(url, tag['href']) for tag in soup.find_all('a', href=True)
                   if 'article/file?' in tag['href'] or 'suppinfo' in tag['href']))
(DEST / 'knapp2019_supplement_links.json').write_text(json.dumps(links, indent=2), encoding='utf-8')
records.append({'system': 'Knapp2019_primary', 'url': url,
                'file': str((DEST / 'knapp2019.html').relative_to(ROOT)),
                'sha256': hashlib.sha256(data).hexdigest()})
for link in links:
    if not link.endswith(('.s004', '.s005')):
        continue
    name = 'knapp2019_' + link.rsplit('.', 1)[1]
    raw = get_small(link)
    if not raw.startswith(b'PK'):
        (DEST / (name + '.access_response.html')).write_bytes(raw)
        records.append({'url': link, 'validation': 'not_office_zip'})
        continue
    (DEST / (name + '.docx')).write_bytes(raw)
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        doc = ET.fromstring(archive.read('word/document.xml'))
    word_ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    tables = []
    for table in doc.findall('.//w:tbl', word_ns):
        tables.append([[''.join(t.text or '' for t in cell.findall('.//w:t', word_ns))
                        for cell in row.findall('w:tc', word_ns)]
                       for row in table.findall('w:tr', word_ns)])
    (DEST / (name + '_tables.json')).write_text(json.dumps(tables, indent=2), encoding='utf-8')
    records.append({'url': link, 'file': str((DEST / (name + '.docx')).relative_to(ROOT)),
                    'sha256': hashlib.sha256(raw).hexdigest(), 'tables': len(tables)})

# Preserve failed responses, with honest extensions; change only named, validated targets.
manifest_path = DEST / 'retrieval_manifest.json'
manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
for name in ['rangarajan2018_attachment_0.docx', 'he2020_attachment_0.docx',
             'he2020_attachment_1.docx']:
    invalid = (DEST / name).resolve()
    assert invalid.parent == DEST.resolve()
    archived = invalid.with_suffix('.access_response.html')
    if invalid.exists() and not invalid.read_bytes().startswith(b'PK'):
        if archived.exists():
            raise ValueError('Refusing to overwrite an archived response')
        invalid.rename(archived)
    if archived.exists():
        for entry in manifest['sources']:
            if Path(entry['file']).name == name:
                entry['file'] = str(archived.relative_to(ROOT))
                entry['validation'] = 'invalid_office_container_access_response'
manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')

# Inspect original XLSX XML directly; preserve source values and sheet names.
with zipfile.ZipFile(DEST / 'notti2025_attachment_0.xlsx') as archive:
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    strings = ET.fromstring(archive.read('xl/sharedStrings.xml'))
    shared = [''.join(t.text or '' for t in si.findall('.//s:t', ns))
              for si in strings.findall('s:si', ns)]
    sheets = {}
    for number in [1, 2, 3, 4, 5]:
        root = ET.fromstring(archive.read('xl/worksheets/sheet%d.xml' % number))
        cells = {}
        for cell in root.findall('.//s:c', ns):
            value = cell.find('s:v', ns)
            if value is not None:
                value = value.text
                value = shared[int(value)] if cell.get('t') == 's' else value
                cells[cell.get('r')] = value
        sheets[str(number)] = cells


def mean_cells(sheet, addresses):
    values = [float(sheets[str(sheet)][address]) for address in addresses]
    return values, sum(values) / len(values)


functional = []
for label, sheet, wt_cells, mut_cells in [
    ('tetramer_peak_activation', 1, ['B3', 'B4', 'B5'], ['C3', 'C4', 'C5']),
    ('COS7_A2_10uM_peptide_percentage_ZsGreen', 3, ['C13', 'C14', 'C15'], ['D13', 'D14', 'D15']),
    ('PMA_ionomycin_bypass_percentage_ZsGreen', 3, ['C3', 'C4', 'C5'], ['D3', 'D4', 'D5']),
]:
    wt, wt_mean = mean_cells(sheet, wt_cells)
    mutant, mutant_mean = mean_cells(sheet, mut_cells)
    functional.append({'assay': label, 'source_sheet_xml': sheet, 'WT_cells': wt_cells,
                       'mutant_cells': mut_cells, 'WT': wt, 'mutant_S104C_V182C': mutant,
                       'WT_mean': wt_mean, 'mutant_mean': mutant_mean,
                       'relative_reduction_percent': 100 * (1 - mutant_mean / wt_mean)})
report = {'source': 'notti2025_attachment_0.xlsx',
          'note': 'Recomputed from source triplicates; these are experimental values, not FlexPPI outputs. '
                  'Sheet 1 concentration header says mg/mL, whereas primary Methods says micrograms/mL; '
                  'do not use the spreadsheet header as a corrected unit.',
          'functional_checks': functional, 'original_cells': sheets}
(DEST / 'notti2025_functional_source_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
(DEST / 'supplement_retrieval_manifest.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
print(json.dumps({'retrievals': records, 'functional_checks': functional,
                  'supplement_links': links}, indent=2))
