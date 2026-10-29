"""Archive official PDB classifications and entity descriptions; preserve unknowns."""
from pathlib import Path
import argparse,json,time,urllib.request,urllib.parse,csv,hashlib,re
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument('--manifest',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--limit',type=int,default=0)
args=parser.parse_args()
records=[json.loads(line) for line in args.manifest.read_text(encoding='utf-8').splitlines() if line.strip()]
entries=set()
for record in records:
    metadata=record.get('metadata',{})
    ids=metadata.get('pdb_ids',[]) or [metadata.get('pdb_id','')]
    entries.update(str(value).upper() for value in ids if re.fullmatch(r'[0-9][A-Za-z0-9]{3}',str(value)))
entries=sorted(entries)
if args.limit:entries=entries[:args.limit]
args.output.mkdir(parents=True,exist_ok=True)
query='''query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title} struct_keywords{pdbx_keywords text} polymer_entities{rcsb_id rcsb_polymer_entity{pdbx_description pdbx_ec} rcsb_polymer_entity_container_identifiers{auth_asym_ids uniprot_ids}}}}'''
for offset in range(0,len(entries),50):
    ids=entries[offset:offset+50]
    archive=args.output/f'batch_{offset:04d}.json'
    if archive.exists():continue
    encoded=urllib.parse.urlencode({'query':query,'variables':json.dumps({'ids':ids})})
    request=urllib.request.Request('https://data.rcsb.org/graphql?'+encoded,headers={'User-Agent':'FlexPPI-manuscript-revision/1.0'})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request,timeout=45) as response:result=json.load(response)
            if result.get('errors'):raise RuntimeError(result['errors'])
            archive.write_text(json.dumps(dict(retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),requested_ids=ids,response=result),ensure_ascii=False,indent=2),encoding='utf-8')
            print(offset+len(ids),len(entries),flush=True)
            break
        except Exception as error:
            if attempt==2:raise
            print(type(error).__name__,str(error)[:160],flush=True)
            time.sleep(2*(attempt+1))
annotations={}
for path in sorted(args.output.glob('batch_*.json')):
    data=json.loads(path.read_text(encoding='utf-8'))
    for entry in data['response']['data']['entries']:
        if entry:annotations[entry['rcsb_id'].upper()]=entry
rows=[]
for record in records:
    metadata=record.get('metadata',{})
    ids=metadata.get('pdb_ids',[]) or [metadata.get('pdb_id','')]
    ids=[str(value).upper() for value in ids if re.fullmatch(r'[0-9][A-Za-z0-9]{3}',str(value))]
    found=[annotations[value] for value in ids if value in annotations]
    keywords=sorted({entry['struct_keywords']['pdbx_keywords'] for entry in found if entry.get('struct_keywords') and entry['struct_keywords'].get('pdbx_keywords')})
    descriptions=sorted({entity['rcsb_polymer_entity']['pdbx_description'] for entry in found for entity in entry.get('polymer_entities') or [] if entity.get('rcsb_polymer_entity') and entity['rcsb_polymer_entity'].get('pdbx_description')})
    enzyme_codes=sorted({entity['rcsb_polymer_entity']['pdbx_ec'] for entry in found for entity in entry.get('polymer_entities') or [] if entity.get('rcsb_polymer_entity') and entity['rcsb_polymer_entity'].get('pdbx_ec')})
    uniprot=sorted({value for entry in found for entity in entry.get('polymer_entities') or [] for value in (entity.get('rcsb_polymer_entity_container_identifiers') or {}).get('uniprot_ids') or []})
    rows.append(dict(sample_id=record['sample_id'],dataset=record['dataset'],pdb_ids=';'.join(ids),
        annotation_status='official_metadata' if found else 'unresolved',pdb_classification=';'.join(keywords),
        protein_descriptions=';'.join(descriptions),enzyme_codes=';'.join(enzyme_codes),uniprot_ids=';'.join(uniprot)))
with (args.output/'record_annotations.csv').open('w',newline='',encoding='utf-8-sig') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
counts=Counter((row['dataset'],row['pdb_classification'] or 'unresolved') for row in rows)
with (args.output/'classification_counts.csv').open('w',newline='',encoding='utf-8-sig') as handle:
    writer=csv.writer(handle);writer.writerow(['dataset','official_pdb_classification','records'])
    writer.writerows((*key,value) for key,value in sorted(counts.items()))
audit=dict(records=len(records),requested_entries=len(entries),returned_entries=len(annotations),resolved_records=sum(row['annotation_status']=='official_metadata' for row in rows),
    method='Official entry-level PDB classification and entity descriptions, not manually validated biological functional categories; original unknowns retained.',
    documentation='https://data.rcsb.org/index.html',manifest_sha256=hashlib.sha256(args.manifest.read_bytes()).hexdigest())
(args.output/'annotation_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2))
