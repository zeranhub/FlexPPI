"""Acquire Fares2009 SI and extract all50 methyl axes without treating blanks as0."""
import re,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import requests
import fitz
from bs4 import BeautifulSoup

OUT=Path(__file__).resolve().parent
SRC=OUT/'sources';SRC.mkdir(exist_ok=True)
ARTICLE='https://link.springer.com/article/10.1007/s10858-009-9354-7'
SUPPLEMENT='https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs10858-009-9354-7/MediaObjects/10858_2009_9354_MOESM1_ESM.pdf'
def download(url,name):
    f=SRC/name
    if not f.exists():
        r=requests.get(url,timeout=30);r.raise_for_status();f.write_bytes(r.content)
    return f

pdf=download(SUPPLEMENT,'Fares2009_SI.pdf');doc=fitz.open(pdf)
text='\n'.join('PAGE '+str(i+1)+'\n'+p.get_text(sort=True) for i,p in enumerate(doc))
(SRC/'Fares2009_SI.txt').write_text(text,encoding='utf-8')
page=doc[2];page.get_pixmap(matrix=fitz.Matrix(2,2)).save(SRC/'Fares2009_TableS2.png')
pattern=r'^\s*(\d+)\s+(\d+)\s+(MET|ILE|VAL|THR|LEU|ALA)\s+(\w+)\s+(\w+)\s+(n|[-\d.]+)\s+(n|[-\d.]+)\s+(n|[-\d.]+)\s*$'
rows=[]
for line in page.get_text(sort=True).splitlines():
    m=re.match(pattern,line)
    if not m:continue
    number,res,aa,atom1,atom2,s2,exp,model=m.groups()
    missing=s2=='n'
    # Two numerical zero triples are NOT measured zero mobility: S4/S5 mark
    # their MF parameters missing, and the main text explicitly excludes them.
    sentinel=int(res)==30 and atom2=='CD1' or int(res)==69 and atom2=='CD1'
    valid=not missing and not sentinel
    rows.append(dict(methyl_number=int(number),residue_number=int(res),residue_name=aa,atom1=atom1,atom2=atom2,
                     s2_raw=float(s2) if not missing else np.nan,
                     experimental_error=float(exp) if exp!='n' else np.nan,
                     combined_experimental_model_error=float(model) if model!='n' else np.nan,
                     valid_order_parameter=valid,
                     s2_physical_capped=min(1,float(s2)) if valid else np.nan,
                     exclusion_reason='unidentifiable_zero_sentinel_crosschecked_S4_S5_and_main_text' if sentinel else ('missing_in_TableS2' if missing else ''),
                     experiment_state='wild_type_free_human_ubiquitin',
                     observable='RDC-based methyl axial order parameter; fast methyl rotation removed',
                     time_window='motions faster than millisecond; includes ps-ns and supra-tumbling ns-us',
                     source_table='Fares2009_SI_TableS2_page3',source_url=SUPPLEMENT))
data=pd.DataFrame(rows)
assert len(data)==50 and data.valid_order_parameter.sum()==37
assert data.methyl_number.tolist()==list(range(1,51))
data.to_csv(OUT/'Fares2009_all50_methyl_order_parameters.csv',index=False)

# Map using CURRENT corrected chain author numbering, not an assumed offset.
cache=np.load(Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\neural_exports\evaluation\full\A003R.npz'))
chain=cache['chain'];author=cache['author_residue_id'];aa=cache['amino_acid']
three=dict(MET='M',ILE='I',VAL='V',THR='T',LEU='L',ALA='A')
maps=[]
for _,r in data.iterrows():
    ix=np.flatnonzero((chain=='A')&(author==str(r.residue_number)))
    present=len(ix)==1
    if present:assert aa[ix[0]]==three[r.residue_name]
    maps.append(dict(methyl_number=r.methyl_number,A003R_chain='A',A003R_author_residue=r.residue_number,
                     A003R_global_index=int(ix[0]) if present else np.nan,
                     A003R_present=present,chemical_identity_match=present,
                     A003R_md_state='ubiquitin_bound_to_EDD1_UBA_2QHO',
                     direct_condition_matched_experimental_validation=False,
                     mismatch='NMR unbound full-length ubiquitin; A003R bound truncated72residue ubiquitin,310K MD; exact NMR temperature unresolved from Fares2009 accessible main text.'))
data=data.merge(pd.DataFrame(maps),on='methyl_number',validate='one_to_one')
data.to_csv(OUT/'Fares2009_to_A003R_all50_mapping.csv',index=False)

html=download(ARTICLE,'Fares2009_article.html')
soup=BeautifulSoup(html.read_bytes(),'html.parser')
full=soup.get_text(' ',strip=True)
(SRC/'Fares2009_article_text.txt').write_text(full,encoding='utf-8')
table_html=download(ARTICLE+'/tables/1','Fares2009_table1_conditions.html')
tables=BeautifulSoup(table_html.read_bytes(),'html.parser').find_all('table')
if tables:
    conditions=pd.read_html(str(tables[0]))[0]
    conditions.to_csv(OUT/'Fares2009_alignment_conditions.csv',index=False)
provenance=dict(article_url=ARTICLE,supplement_url=SUPPLEMENT,raw_methyl_rows=50,valid_axes=37,
                usable_axes_in_A003R=int((data.valid_order_parameter&data.A003R_present).sum()),
                zero_sentinels='I30CG1-CD1 and L69CG-CD1 are missing fits, not measuredS2=0. Verified against main-text rejection list and SI S4/S5 missing orientation fields.',
                upper_bound='S2>1 raw values retained; cap at1 only in separate column following original paper.',
                data_specificity='Methyl-bearing residues only; glutamines Q41/Q62 have no measured methyl axis in this experiment.',
                cannot_validate_model_directly='The model predicts joint3-bin side-chain entropy, not methyl axialS2; available A003R trajectory is bound while NMR is free. No cross-condition correlation is reported as experimental validation.',
                files=[dict(path=str(x),sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in [pdf,html,table_html]])
(OUT/'Fares2009_provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:provenance[k] for k in ['raw_methyl_rows','valid_axes','usable_axes_in_A003R']}))

if __name__=='__main__':pass
