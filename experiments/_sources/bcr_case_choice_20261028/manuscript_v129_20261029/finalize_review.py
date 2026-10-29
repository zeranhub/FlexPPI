from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from zipfile import ZipFile
import csv,json
from PIL import Image
from lxml import etree

R=Path(__file__).resolve().parent
D=R.parent
P=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v129_20261029_中文修订_抗原识别与抗体动态.docx')
changes=json.loads((R/'revision_changes.json').read_text(encoding='utf-8'))
assert sha256(P.read_bytes()).hexdigest()==changes['target_sha256']
pages=json.loads((R/'render_final/page_text.json').read_text(encoding='utf-8'))
assert len(pages)==37
modified=[x['page'] for x in pages if not x['body_pixel_identical_to_previous_pages']]
assert modified==[6,10,12,14,16,18,37]
for n in [6,10,12]:
    assert (R/'render_final'/f'page-{n}.png').read_bytes()==(R/'render'/f'page-{n}.png').read_bytes()
for row in pages:
    # Word's PDF extraction order places its footer before the body text.
    assert row['text'].splitlines()[0].strip()==str(row['page'])

def rows(path):
    with path.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
pred=rows(D/'6A77_all_positions_predictions_and_MD.csv')
full={int(x['author']):x for x in pred if x['variant']=='full' and x['chain']=='C'}
intra={int(x['author']):x for x in pred if x['variant']=='intrachain' and x['chain']=='C'}
assert round(float(full[101]['rmsf_pred']),3)==1.149
assert round(float(full[101]['rmsf_MD']),3)==1.123
assert round(float(intra[101]['rmsf_pred']),3)==1.640
for author,prediction,reference in [(101,3.431,3.419),(96,1.969,1.801)]:
    assert round(float(full[author]['chi_entropy_pred']),3)==prediction
    assert round(float(full[author]['chi_entropy_MD']),3)==reference
contacts=rows(D/'hcd3_functional_analysis_20261029/md_audit/all_h3_sidechain_backbone_and_blocks.csv')
r101=next(x for x in contacts if x['author']=='101' and x['block']=='all')
assert round(float(r101['backbone_B_fraction_4A'])*100,1)==90.3
assert round(float(r101['sidechain_B_fraction_4A'])*100,1)==2.3
static=rows(D/'hcd3_functional_analysis_20261029/static_audit/complete_HCDR3_sidechain_contacts_SASA_sensitivity.csv')
for author,expected in [(96,16.20),(101,189.13)]:
    value=next(x for x in static if int(x['author'])==author)
    assert round(float(value['complex_sidechain_SASA_960']),2)==expected

with ZipFile(P) as z:
    tree=etree.fromstring(z.read('word/document.xml'))
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
text=''.join(tree.xpath('//w:t/text()',namespaces=ns))
assert '图4e' not in text and '图4f' not in text and '图4g' not in text
reference={'PMID':'30595454','DOI':'10.1016/j.str.2018.11.002',
           'title':'Affinity Improvement of a Cancer-Targeted Antibody through Alanine-Induced Adjustment of Antigen-Antibody Interface',
           'journal':'Structure','year':2019,'volume':27,'issue':3,'pages':'519–527.e5','authors':17,
           'primary_sources':['https://pubmed.ncbi.nlm.nih.gov/30595454/',str(D/'structures/6A77_raw.pdb')],
           'verification':'PubMed metadata and raw PDB JRNL authors/title/year/DOI compared with new reference [37]'}
(R/'reference37_verification.json').write_text(json.dumps(reference,ensure_ascii=False,indent=2),encoding='utf-8')
qa={'status':'passed','reviewed_utc':datetime.now(timezone.utc).isoformat(),'document':str(P),
    'document_sha256':changes['target_sha256'],'pages':37,
    'visually_reviewed_changed_pages':modified,'unchanged_body_pages':30,
    'previous_QA':r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\manuscript\v128\render_QA.json',
    'unchanged_page_verification':'30 page bodies pixel-identical to prior individually reviewed pages; automatic footer numbers verified 1..37',
    'renderer':'Unmodified packaged render_docx.py PDF-to-PNG and DPI pipeline, using existing read-only Word PDF-export adapter because no native LibreOffice is installed',
    'visual_checks':{'new_section_heading_and_paragraph_flow':'pass','new_6A77_first_paragraph_kept_together':'pass',
                     'Chinese_Greek_chi_and_Angstrom_glyphs':'pass','reference37_and_supplementH':'pass',
                     'clipping_overlap_missing_glyphs':'none observed',
                     'existing_full_page_figure4':'retained on page15; new figure panel preview intentionally remains outside manuscript'},
    'numeric_verification':'R101/R96 RMSF/Hchi, contact percentages, and sidechain SASA independently checked against source CSVs',
    'new_images_inserted':0,'original_media_binary_preserved':50,
    'existing_TCR_results_and_figure4_caption_unchanged':True,
    'page_reviews':[{'page':x['page'],'sha256':x['sha256'],'review':'visually inspected' if x['page'] in modified else 'body pixel-identical to previously reviewed page','previous_body_pages':x['body_pixel_identical_to_previous_pages']} for x in pages]}
(R/'render_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
changes['render_validation']='passed:37 pages; all7 changed pages visually reviewed;30 remaining page bodies pixel-identical to previously reviewed v128 pages; all footer numbers verified'
changes['numeric_validation']='passed against original full/intra predictions, atom-part contact source, SASA source'
changes['reference37_verification']=str(R/'reference37_verification.json')
(R/'revision_changes.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:qa[k] for k in ['status','pages','new_images_inserted','original_media_binary_preserved','document_sha256']},ensure_ascii=False))
