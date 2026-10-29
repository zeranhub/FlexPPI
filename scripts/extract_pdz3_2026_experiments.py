"""Read-only XLSX extraction with bundled Python; no workbook is edited."""
from pathlib import Path
import json
import hashlib
import csv
from openpyxl import load_workbook

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/pdz3_alpha3_focus_20260914'
records=[]
for n in (2,3):
    p=BASE/'sources'/f'2026_supplementary_data{n}.xlsx'
    wb=load_workbook(p,read_only=True,data_only=True)
    sheet=wb['Sheet1'];it=iter(sheet.values);header=list(next(it))
    keep=[i for i,x in enumerate(header) if x is not None]
    rows=[[r[i] if i<len(r) else None for i in keep] for r in it]
    out=BASE/f'experiment_2026_data{n}.csv'
    with out.open('w',encoding='utf-8',newline='') as f:
        writer=csv.writer(f);writer.writerow([header[i] for i in keep]);writer.writerows(rows)
    notes=[[str(v) for v in row if v is not None] for row in wb['README'].values]
    notes=[r for r in notes if r]
    records.append(dict(file=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),rows=len(rows),
        columns=[header[i] for i in keep],sheet='Sheet1',header_row=1,readme=notes,
        source=f'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-69673-w/MediaObjects/41467_2026_69673_MOESM{n+1}_ESM.xlsx'))
    wb.close()
(BASE/'experiment_extraction_provenance.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
print(json.dumps(records,indent=2))
