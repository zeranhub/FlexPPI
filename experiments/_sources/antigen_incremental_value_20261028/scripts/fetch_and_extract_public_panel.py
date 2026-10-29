from pathlib import Path
import concurrent.futures, hashlib, json, urllib.request
import pandas as pd

R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
S=R/'sources';(S/'pdb').mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
urls={p+'.pdb':f'https://files.rcsb.org/download/{p}.pdb' for p in ['2BNQ','2BNR','3MV7','4PRP','5BRZ','5BS0']}
def fetch(item):
 name,url=item;p=S/'pdb'/name
 if not p.exists():
  with urllib.request.urlopen(url,timeout=45) as f:data=f.read()
  assert len(data)>20000
  p.write_bytes(data)
 return {'path':str(p),'url':url,'sha256':sha(p),'bytes':p.stat().st_size}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:downloads=list(pool.map(fetch,urls.items()))
rows=[];controls=[]
for sheet,system,tcr,hla,primary,secondary in [('NYES0-designs','nyeso','1G4','HLA-A*02:01','2BNQ','2BNR'),('EBV-designs','ebv','TK3','HLA-B*35:01','3MV7','4PRP'),('MAGE-designs','mage','engineered_MAGE_TCR','HLA-A*01:01','5BRZ','5BS0')]:
 d=pd.read_excel(S/'All-designs.xlsx',sheet_name=sheet)
 assert len(d)==96
 for x in d.to_dict('records'):
  seq=x['sequence'];is_design=pd.notna(x['protocol'])
  row={'candidate_id':x['name'],'system':system,'tcr':tcr,'hla':hla,'sequence':seq if pd.notna(seq) else '',
       'is_design':is_design,'functional_positive':bool(x['resp_el4' if system=='mage' else 'resp']),
       'functional_response':float(x['r_el4_mean' if system=='mage' else 'r_mean']),
       'primary_template':primary,'secondary_template':secondary,'source_sheet':sheet,
       'generation_protocol':x.get('protocol',''),'wt_for_design':x.get('wt_for_design',''),
       'tcrdock_pae':x.get('tcrdock_pae'),'af3_interface_pae':x.get('af3_tcr_pmhc_pae_score'),
       'mhc_el_rank_binder_0p5':x.get('is_binder_by_netmhcpan_at_el_rank_0p5'),
       'mhc_el_rank_binder_2p0':x.get('is_binder_by_netmhcpan_at_el_rank_2p0'),
       'min_hamming_distance_from_wildtypes':x.get('min_hamming_distance_from_wildtypes')}
  (rows if is_design else controls).append(row)
 assert sum(r['system']==system for r in rows)==93
pd.DataFrame(rows).to_csv(R/'functional_candidates.csv',index=False)
pd.DataFrame(controls).to_csv(R/'functional_controls_not_ranked.csv',index=False)
audit={'dataset':'HERMES 2025, peptide_designs/All-designs.xlsx','url':'https://raw.githubusercontent.com/StatPhysBio/tcr_antigen_design/main/peptide_designs/All-designs.xlsx',
       'xlsx_sha256':sha(S/'All-designs.xlsx'),'design_candidates':len(rows),'controls_not_ranked':len(controls),
       'endpoint':'Author-supplied functional resp/resp_el4, not direct physical binding or affinity.',
       'preserve_failed_candidates':True,'no_filter_on_prediction_or_presentation_score':True,
       'counts':pd.DataFrame(rows).groupby('system').functional_positive.agg(['size','sum']).to_dict('index'),
       'pdb_downloads':downloads,'note':'Public response labels were inspected to establish dataset suitability; no new FlexPPI or baseline outputs viewed at protocol selection.'}
(R/'audits/public_dataset_extraction.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2))
