"""Independent prediction-blind identity audit of the deposited 6A78/6A79 pair."""
from pathlib import Path
from collections import OrderedDict
import hashlib,json,csv,math

R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
parsed={};metadata={}
for pdb in ['6A78','6A79']:
    p=R/'structures'/(pdb+'_raw.pdb');raw=p.read_text().splitlines();seqres={};residues=OrderedDict();alternates={}
    for line in raw:
        if line.startswith('ENDMDL'):break
        if line.startswith('SEQRES'):seqres.setdefault(line[11],[]).extend(line[19:70].split())
        if not line.startswith('ATOM  ') or line[17:20] not in AA:continue
        key=(line[21],line[22:27].strip());name=line[12:16].strip();alt=line[16]
        score=(alt==' ',float(line[54:60]),alt=='A')
        if (key,name) in alternates and alternates[key,name]>=score:continue
        alternates[key,name]=score
        rec=residues.setdefault(key,{'aa':AA[line[17:20]],'resname':line[17:20],'atoms':{}})
        assert rec['aa']==AA[line[17:20]]
        rec['atoms'][name]=tuple(float(line[i:i+8]) for i in [30,38,46])
    chains={}
    for (ch,author),rec in residues.items():chains.setdefault(ch,OrderedDict())[author]=rec
    parsed[pdb]=chains
    metadata[pdb]={'path':str(p),'source_url':f'https://files.rcsb.org/download/{pdb}.pdb','sha256':sha(p),
        'chains':{ch:{'resolved_residues':len(rs),'author_ids':list(rs),'resolved_sequence':''.join(r['aa'] for r in rs.values()),
                      'SEQRES_length':len(seqres[ch]),'SEQRES_sequence':''.join(AA[a] for a in seqres[ch]),
                      'incomplete_backbone_frames':[au for au,r in rs.items() if not {'N','CA','C'}<=set(r['atoms'])]} for ch,rs in chains.items()},
        'roles_from_COMPND':{'antigen':['A','B'],'VL':['L','M'],'VH_and_SEQRES_linker':['H','I']},
        'biological_assemblies_from_REMARK350':[{'id':1,'chains':['A','L','H'],'operator':'identity'},{'id':2,'chains':['B','M','I'],'operator':'identity'}],
        'linker_scope':'SEQRES of H/I includes author121..136 linker/tail, none has resolved ATOM coordinates. No linker or missing residues fabricated.',
        'missing_residue_records':[l for l in raw if l.startswith('REMARK 465') and l[15:].strip() and l[15:18].strip() in AA]}

pairs=[]
for copy,triplet in [(1,['A','L','H']),(2,['B','M','I'])]:
    chain_pairs=[]
    for ch in triplet:
        wt=parsed['6A78'][ch];mut=parsed['6A79'][ch];common=[a for a in wt if a in mut]
        differences=[{'author':a,'WT':wt[a]['aa'],'mutant':mut[a]['aa']} for a in common if wt[a]['aa']!=mut[a]['aa']]
        chain_pairs.append({'chain':ch,'common_positions':len(common),'WT_only_author_ids':[a for a in wt if a not in mut],
                            'mutant_only_author_ids':[a for a in mut if a not in wt],'shared_identity_differences':differences,
                            'SEQRES_identity_differences':[(i+1,a,b) for i,(a,b) in enumerate(zip(metadata['6A78']['chains'][ch]['SEQRES_sequence'],metadata['6A79']['chains'][ch]['SEQRES_sequence'])) if a!=b]})
    changes=[(c['chain'],d['author'],d['WT'],d['mutant']) for c in chain_pairs for d in c['shared_identity_differences']]
    assert changes==[(triplet[2],'103','P','A')],changes
    pairs.append({'copy':copy,'primary':copy==1,'selection_basis':'Depositor biological assembly1 fixed before predictions; assembly2 is complete sensitivity, not best-result selection.',
                  'triplet':triplet,'WT_resolved_total':sum(len(parsed['6A78'][ch]) for ch in triplet),
                  'mutant_resolved_total':sum(len(parsed['6A79'][ch]) for ch in triplet),'chain_pairs':chain_pairs})

nr=json.loads((R/'audits/independent_pair_anarcii_raw.json').read_text())
domain_summary={};number_rows=[]
for key,x in nr.items():
    pdb,ch=key.split('_');chain=parsed[pdb][ch];seq=''.join(r['aa'] for r in chain.values());authors=list(chain)
    assert x['error'] is None and x['chain_type']==('K' if ch in ['L','M'] else 'H')
    cursor=x['query_start'];mapped=[]
    for (pos,ins),aa in x['numbering']:
        if aa=='-':continue
        assert seq[cursor]==aa
        row={'pdb':pdb,'chain':ch,'author':authors[cursor],'aa':aa,'imgt_number':pos,'imgt_insertion':ins.strip()}
        mapped.append(row);number_rows.append(row);cursor+=1
    domain_summary[key]={'chain_type':x['chain_type'],'variable_first':mapped[0],'variable_last':mapped[-1],
                         'named_sites':[r for r in mapped if r['author'] in ['31','97','103']]}
with (R/'audits/independent_pair_residue_numbering.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(number_rows[0]));w.writeheader();w.writerows(number_rows)

geometry=[]
for pdb in ['6A78','6A79']:
    for ag,h in [('A','H'),('B','I')]:
        arg=parsed[pdb][ag]['67'];assert arg['aa']=='R'
        distances={}
        for target in ['31','97']:
            acid=parsed[pdb][h][target];assert acid['aa']==('D' if target=='31' else 'E')
            n=[v for k,v in arg['atoms'].items() if k in ['NE','NH1','NH2']]
            o=[v for k,v in acid['atoms'].items() if k in ['OD1','OD2','OE1','OE2']]
            distances[h+':'+target]=min(math.dist(a,b) for a in n for b in o)
        geometry.append({'pdb':pdb,'antigen_chain':ag,'antigen_author':'67','antigen_aa':'R','minimum_ARG_N_to_carboxylate_O_Angstrom':distances,
                         'scope':'Static structure distance only, not MD occupancy, energetic strength or entropy.'})

sifts=json.loads((R/'reference/6a78_uniprot_mapping.json').read_text())
uniprot=''.join((R/'reference/Q9Y6N7.fasta').read_text().splitlines()[1:])
query=metadata['6A78']['chains']['A']['SEQRES_sequence'][2:]
assert uniprot.find(query)==454
assert uniprot[512]=='R' and uniprot[476]=='V'
report={'status':'passed_with_numbering_caveat','prediction_blind':True,'new_model_outputs_read':False,
    'inputs':metadata,'paired_comparisons':pairs,'domain_numbering':domain_summary,'static_named_network_geometry':geometry,
    'mutant_identity':'Shared positions differ only at heavy-chain author103, P->A; same substitution in both copies, IMGT112A. Do not use row103 as position103.',
    'recommended_primary':'Assembly1 A/L/H for WT and mutant; preserve as-deposited distinct chains; assembly2 B/M/I as sensitivity. All six asymmetric-unit chains must not be treated as one biological complex.',
    'recommended_alignment':'Join within each biological copy by role plus author_id, validate AA identity, retain all shared positions. Exclude changed103 from same-AA Hchi comparisons. Missing terminal H120/M108 are audited, never padded or matched by row index.',
    'scFv_construct_note':'Biological scFv has connected variable domains, but deposited modeled coordinates use separate VL and VH chains because linker is unresolved. Do not invent a continuous modeled linker or assert full physical covalent connectivity from ATOM chain IDs.',
    'network_numbering':{'structural_antigen_site':'A/B author67 ARG','current_canonical_UniProt':'Q9Y6N7 R513',
       'mapping_evidence':'Official PDBe/SIFTS: author9..97 equals UniProt455..543 at 100% identity; full sequence is identical to current UniProt455..543.',
       'author31':'VAL, not ARG','literature_R477':'Numbering basis unresolved; current Q9Y6N7 and isoform2 position477 is VAL. Never silently translate literature477 to PDB31 or assert isoform correspondence without evidence.'},
    'source_sha256':{str(p):sha(p) for p in [R/'structures/6A78_raw.pdb',R/'structures/6A79_raw.pdb',R/'audits/independent_pair_anarcii_raw.json',R/'reference/6a78_uniprot_mapping.json',R/'reference/Q9Y6N7.fasta',R/'reference/Q9Y6N7-2.fasta',Path(__file__)]}}
(R/'audits/paired_input_identity_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':report['status'],'pairs':pairs,'network_numbering':report['network_numbering']},ensure_ascii=False,indent=2))
