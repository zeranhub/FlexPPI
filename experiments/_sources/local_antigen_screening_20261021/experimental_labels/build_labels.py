"""Rebuild experimental labels from immutable, archived public source files only.

No FlexPPI predictions are read. A label is an experimental study mean, not an MD
replicate or structural model. Cross-study combination follows HERMES' log mean.
"""
from pathlib import Path
import csv,hashlib,json,math,re,statistics,datetime,xml.etree.ElementTree as ET
from collections import defaultdict

ROOT=Path(__file__).parent
RAW=ROOT/'raw'
META=json.loads((ROOT/'source_candidates.json').read_text())
COMMIT=META['commit']
REPO='StatPhysBio/tcr_antigen_design'
BASE=f'https://raw.githubusercontent.com/{REPO}/{COMMIT}/'
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
PANEL={
 '1G4_NYESO':dict(tcr='1G4',antigen='NY-ESO-1',hla='HLA-A*02:01',panel_reference='SLLMWITQV',native='SLLMWITQC',stem='nyeso',templates=['2bnq','2bnr']),
 'A6_TAX':dict(tcr='A6',antigen='HTLV-1 Tax',hla='HLA-A*02:01',panel_reference='LLFGYPVYV',native='LLFGYPVYV',stem='tax',templates=['1ao7','1qsf']),
}
def readcsv(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writecsv(name,rows):
    with (ROOT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def mutations(a,b):return '|'.join(f'{x}{i+1}{y}' for i,(x,y) in enumerate(zip(a,b)) if x!=y)
def chains(path):
    found={};seen=set()
    for line in Path(path).read_text().splitlines():
        if line[:6].strip()=='ENDMDL':break
        if line[:6].strip()!='ATOM' or line[16] not in ' A1':continue
        key=(line[21],line[22:27]);name=line[17:20]
        if key in seen or name not in AA:continue
        seen.add(key);found.setdefault(key[0],[]).append((line[22:27].strip(),AA[name]))
    return {k:{'sequence':''.join(a for _,a in v),'author_ids':[i for i,_ in v]} for k,v in found.items()}

# Primary Pettmann files contain two header rows. Use the Bmax-constrained values,
# the measurement method used by HERMES; retain fitted Bmax as sensitivity data.
primary={};method_rows=[]
for pid,info in PANEL.items():
    p=RAW/('elife-67092-fig1-data2-v3.csv' if pid=='1G4_NYESO' else 'elife-67092-fig1-data3-v3.csv')
    data=list(csv.reader(p.open(encoding='utf-8-sig',newline='')))
    assert data[1][6]=='Mean (μM)'
    for line,row in enumerate(data[2:],start=3):
        if not row:continue
        name,seq=row[:2]
        record={'name':name,'peptide_sequence':seq,'KD_uM':float(row[6]),'error_value':float(row[7]) if row[7]!='N/A' else None,'error_type':'SD','n':int(row[9]),'source_file':str(p),'source_locator':f'Figure 1—source data {2 if pid=="1G4_NYESO" else 3}; CSV line {line}; columns 7-10 (Bmax constrained)','source_url':f'https://cdn.elifesciences.org/articles/67092/{p.name}'}
        primary[(pid,'pettmann_et_al',seq)]=record
        method_rows.append(dict(panel_id=pid,**record,KD_Bmax_fitted_uM=float(row[2]),KD_Bmax_fitted_SD_uM=float(row[3]) if row[3]!='N/A' else '',KD_Bmax_fitted_n=int(row[5])))

# Aleksic Table S1 is an original page-16 PDF table with 17 pMHC variants. Eleven
# alter only peptide on WT HLA; six HLA substitutions are retained separately.
alpdf=RAW/'primary_supplements/aleksic2010_supplementaryFiles/mmc1.pdf'
if not alpdf.with_suffix('.txt').exists():
    from pypdf import PdfReader
    alpdf.with_suffix('.txt').write_text('\n\n'.join('PAGE '+str(i+1)+'\n'+p.extract_text() for i,p in enumerate(PdfReader(alpdf).pages)),encoding='utf8')
altxt=alpdf.with_suffix('.txt').read_text(encoding='utf8')
block=altxt.split('PAGE 16\n',1)[1].split('peptide name peptide sequence',2)[1]
alrows=[]
for line in block.splitlines():
    m=re.match(r'^(ESO-\S+|A2-\S+) ([A-Z]{9}) ([0-9.]+) ± ([0-9.]+) ',line)
    if m:alrows.append(dict(zip(['name','peptide_sequence','KD_uM','KD_SEM_uM'],m.groups())))
assert len(alrows)==17,len(alrows)
for i,entry in enumerate(alrows,start=1):
    entry['KD_uM']=float(entry['KD_uM']);entry['KD_SEM_uM']=float(entry['KD_SEM_uM'])
    entry['HLA_background']='HLA-A*02:01 WT' if entry['name'].startswith('ESO-') else 'HLA-A*02:01 '+entry['name'][3:]+'→A (author naming)'
    entry['retained_in_fixed_HLA_antigen_panel']=entry['name'].startswith('ESO-')
    entry['source_locator']=f'Table S1, PDF page 16, upper table data row {i}, KD (μM) column'
    if entry['retained_in_fixed_HLA_antigen_panel']:
        primary[('1G4_NYESO','aleksic_et_al',entry['peptide_sequence'])]=dict(name=entry['name'],peptide_sequence=entry['peptide_sequence'],KD_uM=entry['KD_uM'],error_value=entry['KD_SEM_uM'],error_type='SEM',n=None,source_file=str(alpdf),source_locator=entry['source_locator']+'; mean±SEM; n≥3 independent experiments',source_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2862301/supplementaryFiles')

template_rows=[];template_by_pdb={}
for pid,info in PANEL.items():
    for pdb in info['templates']:
        path=RAW/'github/mutation_effects'/info['stem']/'pdbs'/f'{pdb}.pdb'
        cs=chains(path);pep=cs['C']['sequence']
        rec=dict(panel_id=pid,pdb_id=pdb.upper(),peptide_chain='C',peptide_sequence=pep,pdb_path=str(path),pdb_sha256=sha(path),source_url=BASE+f"mutation_effects/{info['stem']}/pdbs/{pdb}.pdb",chain_sequences=cs)
        template_rows.append(rec);template_by_pdb[pdb]=rec
assert template_by_pdb['2bnr']['peptide_sequence']=='SLLMWITQC'
assert template_by_pdb['2bnq']['peptide_sequence']=='SLLMWITQV'
assert template_by_pdb['1ao7']['peptide_sequence']=='LLFGYPVYV'
assert template_by_pdb['1qsf']['peptide_sequence']=='LLFGYPVAV'

study_rows=[];aggregated=[];checks=[]
for pid,info in PANEL.items():
    stem=info['stem'];relative=f'mutation_effects/{stem}/{stem}_peptide_kd.csv'
    raw=readcsv(RAW/'github'/relative)
    groups=defaultdict(list)
    for i,r in enumerate(raw,start=2):
        seq=r['sequence'];study=r['study'];label=primary[(pid,study,seq)]
        kd=float(r['Kd [microM]']);assert kd==label['KD_uM'],(pid,study,seq,kd,label['KD_uM'])
        assert abs(-math.log10(kd)-float(r['-log10(Kd)']))<1e-12
        rec=dict(observation_id=f'{pid}_{study}_{seq}',panel_id=pid,tcr=info['tcr'],hla=info['hla'],antigen=info['antigen'],peptide_sequence=seq,peptide_name=label['name'],assay='SPR',temperature_C=37,KD_value=kd,KD_unit='uM',KD_uM=kd,pKD_M=6-math.log10(kd),negative_log10_KD_uM=-math.log10(kd),KD_relation='=',censoring='none_reported',detection_limit_uM='',error_value=label['error_value'] if label['error_value'] is not None else '',error_type=label['error_type'],n=label['n'] if label['n'] is not None else '',n_min=3 if study=='aleksic_et_al' else label['n'],SPR_method='steady-state; Bmax constrained by W6/32 calibration' if study=='pettmann_et_al' else 'steady-state; conventional fitting',study_id=study,primary_doi='10.7554/eLife.67092' if study=='pettmann_et_al' else '10.1016/j.immuni.2009.11.013',primary_source=label['source_file'],primary_source_url=label['source_url'],primary_locator=label['source_locator'],primary_source_sha256=sha(label['source_file']),HERMES_source_url=BASE+relative,HERMES_source_locator=f'{relative}: CSV line {i}',HERMES_source_sha256=sha(RAW/'github'/relative))
        study_rows.append(rec);groups[seq].append(rec)
    closest=readcsv(RAW/'github/mutation_effects'/stem/f'{stem}_peptide_kd_closest.csv')
    expected_by_seq={x['sequence']:x for x in closest}
    for idx,(seq,obs) in enumerate(groups.items(),start=1):
        score=statistics.mean(x['negative_log10_KD_uM'] for x in obs)
        kd=10**(-score)
        temp=min(info['templates'],key=lambda pdb:sum(a!=b for a,b in zip(seq,template_by_pdb[pdb]['peptide_sequence'])))
        authors=expected_by_seq[seq]
        assert temp==authors['pdb']
        assert abs(score-float(authors['-log10(Kd)']))<1e-12
        checks.append(dict(panel_id=pid,peptide_sequence=seq,source_raw_KD_first=obs[0]['KD_uM'],HERMES_closest_raw_KD=float(authors.get('Kd [microM]',authors.get('Kd'))),recomputed_geometric_KD_uM=kd,log_mean_matches_HERMES=True,closest_template_matches=True))
        aggregated.append(dict(candidate_id=f'{pid}_{seq}',panel_id=pid,tcr=info['tcr'],hla=info['hla'],antigen=info['antigen'],peptide_sequence=seq,peptide_length=len(seq),panel_reference_peptide=info['panel_reference'],native_peptide=info['native'],mutations_vs_panel_reference=mutations(info['panel_reference'],seq),mutations_vs_native=mutations(info['native'],seq),assay='SPR',temperature_C=37,KD_uM=kd,KD_M=kd*1e-6,pKD_M=6+score,negative_log10_KD_uM=score,KD_relation='=',censoring='none_reported',detection_limit_uM='',analysis_role='HERMES-compatible cross-study aggregate; secondary sensitivity only',aggregation='equal-weight study means on -log10(KD/uM); KD is geometric mean',n_studies=len(obs),study_ids='|'.join(x['study_id'] for x in obs),study_observation_ids='|'.join(x['observation_id'] for x in obs),template_pdb=temp.upper(),template_peptide_chain='C',template_peptide_sequence=template_by_pdb[temp]['peptide_sequence'],template_hamming=sum(a!=b for a,b in zip(seq,template_by_pdb[temp]['peptide_sequence'])),template_selection='minimum sequence Hamming distance among PNAS templates; fixed listed order resolves ties; no labels used',template_pdb_path=template_by_pdb[temp]['pdb_path'],template_pdb_sha256=template_by_pdb[temp]['pdb_sha256'],primary_references='|'.join(sorted({x['primary_doi'] for x in obs})),HERMES_commit=COMMIT))

assert len(study_rows)==36 and len(aggregated)==30
assert sum(x['panel_id']=='1G4_NYESO' for x in aggregated)==14
assert sum(x['panel_id']=='A6_TAX' for x in aggregated)==16
assert len(primary)==36
by_key={(r['panel_id'],r['study_id'],r['peptide_sequence']):r for r in study_rows}
by_candidate={r['candidate_id']:r for r in aggregated}
for r in study_rows:
    info=PANEL[r['panel_id']];parent=by_key[(r['panel_id'],r['study_id'],info['panel_reference'])]
    cid=f"{r['panel_id']}_{r['peptide_sequence']}";candidate=by_candidate[cid]
    r.update(candidate_id=cid,study=r['study_id'],parent_sample=parent['observation_id'],parent_peptide=parent['peptide_sequence'],parent_KD_uM=parent['KD_uM'],parent_pKD_M=parent['pKD_M'],delta_pKD_vs_parent=r['pKD_M']-parent['pKD_M'],parent_source_locator=parent['primary_locator'],parent_same_study=True,parent_same_protocol_temperature=True,parent_same_measurement_batch='not_reported',pairing_level='within-study protocol-matched published means; not individually batch-paired measurements',batch_id='not_reported',is_parent=r['observation_id']==parent['observation_id'],template_pdb=candidate['template_pdb'],template_peptide_chain='C',template_peptide_sequence=candidate['template_peptide_sequence'],template_pdb_path=candidate['template_pdb_path'],analysis_role='primary: within-study delta-pKD relative to its own measured parent; absolute-pKD secondary',replication_unit='published peptide-study mean; raw independent measurements not provided here')
writecsv('SPR_labels_study_level.csv',study_rows)
writecsv('SPR_labels_observations.csv',study_rows)
writecsv('SPR_labels_candidates.csv',aggregated)
writecsv('SPR_labels_HERMES_aggregate_secondary.csv',aggregated)
unlabelled_columns=['candidate_id','panel_id','tcr','hla','antigen','peptide_sequence','peptide_length','panel_reference_peptide','native_peptide','mutations_vs_panel_reference','mutations_vs_native','template_pdb','template_peptide_chain','template_peptide_sequence','template_hamming','template_selection','template_pdb_path','template_pdb_sha256']
writecsv('SPR_candidates.csv',[{k:r[k] for k in unlabelled_columns} for r in aggregated])
writecsv('Pettmann2021_Bmax_methods.csv',method_rows)
writecsv('Aleksic2010_TableS1_complete_pMHC_panel.csv',alrows)
(ROOT/'templates.json').write_text(json.dumps(template_rows,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'SPR_labels.json').write_text(json.dumps(dict(schema_version='1.1',primary_endpoint='within-study delta-pKD relative to measured 9V or Tax-WT parent; same batch not reported',study_level=study_rows,observations=study_rows,HERMES_aggregate_secondary=aggregated,candidates_unlabelled=[{k:r[k] for k in unlabelled_columns} for r in aggregated],templates=template_rows),ensure_ascii=False,indent=2),encoding='utf8')

summary=dict(status='passed',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),repository=REPO,commit=COMMIT,source_selection='Complete 1G4 and A6 SPR panels from HERMES Fig2; selected before any FlexPPI inference',panels=[dict(panel_id=pid,n_candidates=sum(x['panel_id']==pid for x in aggregated),n_study_observations=sum(x['panel_id']==pid for x in study_rows),KD_min_uM=min(x['KD_uM'] for x in aggregated if x['panel_id']==pid),KD_max_uM=max(x['KD_uM'] for x in aggregated if x['panel_id']==pid),hla=info['hla'],tcr=info['tcr']) for pid,info in PANEL.items()],n_unique_candidates=30,n_study_observations=36,n_censored=0,primary_value_checks=36,log_mean_and_template_checks=checks,important_notes=[
 'HERMES closest CSV retains first-study Kd while averaging -log10(Kd); candidate KD_uM here is recomputed geometric mean, not that stale source column.',
 'pKD_M = -log10(KD/M) = 6 - log10(KD/uM); HERMES uses -log10(KD/uM), which differs by 6.',
 'Highest 1309 and 2321 uM values are published point estimates from constrained Bmax SPR, not thresholds/censoring. Large SD and some A6 n=1 remain visible at study level.',
 '2BNR is native 9C and 2BNQ is 9V in actual coordinates and PNAS supplement Table S1; PNAS main Fig2A prose reverses these identifiers.',
 'The primary Aleksic Table S1 contains 17 pMHCs; six HLA mutations are retained in a separate table and excluded from the fixed-HLA peptide panel.',
 'All candidates have measured point affinities; they are not an unbiased population of binders and nonbinders. No externally justified binary screening threshold is assigned.',
 'Two TCR series are only two independent biological groups; 30 peptide candidates do not establish broad generalization.',
 'Study means are retained separately; aggregate uncertainty is not manufactured by pooling SD and SEM.',
 'Primary analysis uses observations with within-study parent delta-pKD. HERMES cross-study average is secondary only. Same study/37C/method is verified, but exact same experimental batch IDs were not reported.',
],visual_primary_table_QA=dict(Aleksic_TableS1_page16='Individually rendered and opened; all 11 fixed-HLA peptide KD values and SEM checked against source CSV.',PNAS_TableS1_page14='Individually rendered and opened; 2BNR/2BNQ and 1AO7/1QSF peptide sequences and HLA checked.'),files={name:dict(path=str(ROOT/name),sha256=sha(ROOT/name)) for name in ['SPR_labels_observations.csv','SPR_candidates.csv','SPR_labels_study_level.csv','SPR_labels_candidates.csv','SPR_labels_HERMES_aggregate_secondary.csv','SPR_labels.json','templates.json','Pettmann2021_Bmax_methods.csv','Aleksic2010_TableS1_complete_pMHC_panel.csv']})
(ROOT/'label_validation_report.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
manifest=[]
known_urls={str(Path(x['path']).resolve()):x['url'] for name in ['downloaded_panel_files.json','primary_downloads.json'] for x in json.loads((ROOT/name).read_text()) if x.get('path') and 'default.jpg' not in x['path'] and Path(x['path']).exists()}
known_urls.update({str((RAW/name).resolve()):url for name,url in {'pnas2025_fullText.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12557503/fullTextXML','pettmann2021_fullText.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8219380/fullTextXML','aleksic2010_fullText.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2862301/fullTextXML','pettmann_elife_api.json':'https://api.elifesciences.org/articles/67092','github_README.md':BASE+'README.md'}.items()})
for p in sorted(RAW.rglob('*')):
    if p.is_file():
        valid='not_used_auxiliary' if p.name=='default.jpg' else 'blocked_HTML_not_used' if p.name=='pnas_PMC12557503.html' else 'archived_source_or_extracted_supplement'
        manifest.append(dict(path=str(p),relative_path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p),url=known_urls.get(str(p.resolve()),''),status=valid))
(ROOT/'source_manifest.json').write_text(json.dumps(dict(status='passed',commit=COMMIT,files=manifest),ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:summary[k] for k in ['status','n_unique_candidates','n_study_observations','primary_value_checks','n_censored','panels']},indent=2))
