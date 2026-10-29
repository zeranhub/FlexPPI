"""Protocol-fixed antigen-panel features. No SPR table or model fitting is read.

Two stages are supported: ``freeze-regions`` archives heavy-atom 8A environments
from complete parent structures; ``extract`` reads all-residue physical-scale
NPZ means and complete static structures, subtracting each template's own parent.
All writes remain under this script's directory.
"""
from __future__ import annotations
import argparse,copy,csv,datetime,hashlib,json,math,sys
from collections import defaultdict
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
EXPECTED_CHECKPOINT='b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f'
BACKBONE={'N','CA','C','O'}
ONE_TO_THREE={
    'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE',
    'G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU',
    'M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG',
    'S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR',
}
THREE_TO_ONE={v:k for k,v in ONE_TO_THREE.items()}
SIDECHAINS={
 'A':'CB','C':'CB SG','D':'CB CG OD1 OD2','E':'CB CG CD OE1 OE2',
 'F':'CB CG CD1 CD2 CE1 CE2 CZ','G':'','H':'CB CG ND1 CD2 CE1 NE2',
 'I':'CB CG1 CG2 CD1','K':'CB CG CD CE NZ','L':'CB CG CD1 CD2',
 'M':'CB CG SD CE','N':'CB CG OD1 ND2','P':'CB CG CD',
 'Q':'CB CG CD OE1 NE2','R':'CB CG CD NE CZ NH1 NH2','S':'CB OG',
 'T':'CB OG1 CG2','V':'CB CG1 CG2','W':'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
 'Y':'CB CG CD1 CD2 CE1 CE2 CZ OH',
}
CANONICAL={aa:BACKBONE|set(names.split()) for aa,names in SIDECHAINS.items()}
CHI_COUNTS=dict(zip('ACDEFGHIKLMNPQRSTVWY',[0,1,2,3,2,0,2,2,4,2,3,2,2,3,4,1,1,1,2,2]))
# Explicit element radii are fixed across candidates and all three SASA calls.
RADII={'C':1.70,'N':1.55,'O':1.52,'S':1.80}
ROLE_NAMES=('hla_heavy','beta2m','peptide','tcr_alpha','tcr_beta')

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def dump(path,value):
    path=Path(path).resolve();assert path.is_relative_to(HERE)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
def dump_csv(path,rows,fields=None):
    path=Path(path).resolve();assert path.is_relative_to(HERE)
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(rows[0])
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def records(payload):
    if isinstance(payload,list):return payload
    for field in ['candidates','records','structures','results']:
        if isinstance(payload.get(field),list):return payload[field]
    raise ValueError('Manifest must contain candidates/records/structures/results list')
def value(record,*names):
    for name in names:
        if record.get(name) is not None:return record[name]
    return None
def biological_id(record):return value(record,'biological_candidate_id','candidate_id')
def input_id(record):return value(record,'inference_candidate_id','input_id','candidate_id')
def role_map(record):
    roles=record['roles'];assert set(ROLE_NAMES)<=set(roles)
    assert len(set(roles[k] for k in ROLE_NAMES))==5
    return roles
def input_path(record):
    p=value(record,'prepared_input_pdb','prepared_pdb','model_input_pdb','input_pdb')
    if p is None:raise ValueError('Missing fixed prepared PDB path')
    return Path(p)
def static_path(record):
    p=value(record,'static_complete_pdb','static_pdb','complete_pdb','output_pdb','input_pdb','prepared_pdb','prepared_input_pdb','pdb_path')
    if p is None:raise ValueError('Missing complete static PDB path')
    return Path(p)
def true(v):return v is True or str(v).lower() in {'true','1','yes','passed','success','complete','completed'}

def dependencies(paths):
    # Caller chooses an explicit existing dependency runtime; never install here.
    for p in reversed(paths):sys.path.insert(0,str(Path(p).resolve()))
    global np,cKDTree,Bio,ShrakeRupley,Structure,Model,Chain,Residue,Atom,substitution_matrices
    import numpy as np
    from scipy.spatial import cKDTree
    import Bio
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB import Structure,Model,Chain,Residue,Atom
    from Bio.Align import substitution_matrices

def read_structure(path):
    """One model, standard heavy protein atoms, deterministic altloc selection.

    Residue keys use author identifiers as in FlexPPI parse_pdb, never integer
    positions shared between templates. Only water/nonprotein atoms are ignored.
    """
    result={};ignored=[]
    for line in Path(path).read_text(encoding='utf8').splitlines():
        kind=line[:6].strip()
        if kind=='ENDMDL':break
        if kind not in {'ATOM','HETATM'} or len(line)<54:continue
        if line[16] not in {' ','A','1'}:continue
        name=line[17:20].strip().upper();atom=line[12:16].strip().upper()
        element=line[76:78].strip().upper() if len(line)>=78 else ''
        if not element:element=re_element(atom)
        if element in {'H','D'}:continue
        if name not in THREE_TO_ONE:
            if name not in {'HOH','WAT','DOD'}:ignored.append({'resname':name,'chain':line[21:22].strip(),'atom':atom})
            continue
        key=(line[21:22].strip() or '_',line[22:26].strip()+line[26:27].strip())
        aa=THREE_TO_ONE[name]
        if key not in result:result[key]={'aa':aa,'resname':name,'atoms':{}}
        if result[key]['aa']!=aa:raise ValueError(('Mixed amino-acid residue',path,key))
        if atom not in result[key]['atoms']:
            result[key]['atoms'][atom]={'xyz':np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])],dtype=np.float64),'element':element}
    if not result:raise ValueError('No protein residues: '+str(path))
    return result,ignored
def re_element(atom):
    # Standard protein atom names CA/CD etc denote carbon, never calcium.
    return next((x for x in atom if x.isalpha()),'').upper()
def sequence(st,chain):return ''.join(v['aa'] for (c,_),v in st.items() if c==chain)
def has_frame(residue):return {'N','CA','C'}<=set(residue['atoms'])
def complete(residue):return CANONICAL[residue['aa']]<=set(residue['atoms'])
def missing_atoms(residue):return sorted(CANONICAL[residue['aa']]-set(residue['atoms']))
def coords(st,chains):
    out=[];keys=[]
    for key,r in st.items():
        if key[0] in chains:
            for atom in r['atoms'].values():out.append(atom['xyz']);keys.append(key)
    return np.asarray(out,dtype=np.float64),keys
def structure_qc(st,roles,peptide):
    missing=[{'chain':k[0],'id':k[1],'aa':r['aa'],'missing':missing_atoms(r)} for k,r in st.items() if missing_atoms(r)]
    return {'peptide_sequence_matches':sequence(st,roles['peptide'])==peptide,'chain_roles_present':set(roles.values())=={k[0] for k in st},'canonical_missing':missing,'all_canonical_heavy_atoms_complete':not missing,'unsupported_atom_elements':sorted({a['element'] for r in st.values() for a in r['atoms'].values()}-set(RADII))}
def backbone_comparison(fixed,other):
    errors=[];max_delta=0.0
    if set(fixed)!=set(other):errors.append('residue_identity_keys_differ')
    for key,res in fixed.items():
        if key not in other:continue
        if res['aa']!=other[key]['aa']:errors.append('AA_changed:'+':'.join(key))
        for atom in BACKBONE:
            if atom in res['atoms']:
                if atom not in other[key]['atoms']:errors.append('backbone_atom_missing:'+':'.join(key)+':'+atom)
                else:max_delta=max(max_delta,float(np.max(np.abs(res['atoms'][atom]['xyz']-other[key]['atoms'][atom]['xyz']))))
    # PDB text is stored to 0.001A. Exact preservation is expected at this level.
    if max_delta>0.0001:errors.append('NCAO_backbone_moved')
    return {'passed':not errors,'max_NCAO_coordinate_difference_A':max_delta,'errors':errors}

def load_manifests(args):
    fixed_path=Path(args.fixed_manifest).resolve();static_manifest=Path(args.static_manifest).resolve()
    f=read_json(fixed_path);s=read_json(static_manifest)
    fr=records(f);sr=records(s);sindex={input_id(r):r for r in sr}
    assert len(sindex)==len(sr),'Duplicate static input ids'
    ids=[input_id(r) for r in fr];assert len(set(ids))==len(ids)
    assert len(fr)==76,('Expected complete 76-input panel',len(fr))
    for r in fr:
        assert input_id(r) in sindex,('Static manifest missing input',input_id(r))
        role_map(r)
    return f,s,fr,sindex,{'fixed_manifest':str(fixed_path),'fixed_manifest_sha256':digest(fixed_path),'static_manifest':str(static_manifest),'static_manifest_sha256':digest(static_manifest)}

def freeze_regions(args):
    protocol=read_json(args.protocol);f,s,fr,sindex,bindings=load_manifests(args)
    templates=defaultdict(list)
    for r in fr:templates[r['template_pdb'].upper()].append(r)
    payload={'status':'frozen_without_model_or_SPR_values','created_at_utc':utc(),'protocol_id':protocol['protocol_id'],'protocol_sha256':digest(args.protocol),'radius_A':8.0,'Hchi_rule':'unchanged non-Gly/non-Ala U with chemically defined chi, valid original NCA frame, complete canonical parent heavy atoms','fixed_input_bindings':bindings,'templates':{}}
    for template,rows in templates.items():
        parents=[r for r in rows if true(r.get('is_panel_parent'))];assert len(parents)==1,(template,'parent count')
        p=parents[0];roles=role_map(p);fixed,_=read_structure(input_path(p));sr=sindex[input_id(p)];sp=static_path(sr);st,ignored=read_structure(sp)
        qc=structure_qc(st,roles,p['peptide_sequence']);bb=backbone_comparison(fixed,st)
        assert qc['peptide_sequence_matches'] and qc['chain_roles_present'],(template,{'peptide_sequence_matches':qc['peptide_sequence_matches'],'chain_roles_present':qc['chain_roles_present'],'n_canonical_missing':len(qc['canonical_missing']),'missing_first3':qc['canonical_missing'][:3]})
        assert bb['passed'],(template,bb)
        assert not ignored,(template,'Unexpected nonprotein heavy atoms',ignored)
        # All parent side chains must be complete; deriving U from deleted mutant
        # side chains would make the static comparator and neighborhood invalid.
        assert qc['all_canonical_heavy_atoms_complete'],(template,{'n_canonical_missing':len(qc['canonical_missing']),'missing_first3':qc['canonical_missing'][:3]})
        pep,_=coords(st,{roles['peptide']});tree=cKDTree(pep)
        environment={roles['hla_heavy'],roles['tcr_alpha'],roles['tcr_beta']}
        role_rev={v:k for k,v in roles.items()};u=[]
        for key,r in st.items():
            if key[0] not in environment:continue
            xyz=np.asarray([a['xyz'] for a in r['atoms'].values()])
            minimum=float(np.min(tree.query(xyz,k=1)[0]))
            if minimum<=8.0:
                original=fixed.get(key)
                frame=original is not None and has_frame(original)
                assert original is not None and original['aa']==r['aa']
                hvalid=frame and CHI_COUNTS[r['aa']]>0 and complete(r)
                u.append({'chain':key[0],'author_residue_id':key[1],'amino_acid':r['aa'],'role':role_rev[key[0]],'reference_peptide_distance_min_A':minimum,'RMSF_valid':frame,'Hchi_valid':hvalid,'n_chi':CHI_COUNTS[r['aa']],'canonical_heavy_complete':complete(r)})
        assert u and any(x['RMSF_valid'] for x in u) and any(x['Hchi_valid'] for x in u)
        uid=hashlib.sha256(json.dumps(u,sort_keys=True).encode()).hexdigest()
        out={'feature_region_id':template+'_U8A_'+uid[:12],'template_pdb':template,'panel_id':p['panel_id'],'parent_input_id':input_id(p),'parent_biological_candidate_id':biological_id(p),'parent_peptide':p['peptide_sequence'],'roles':roles,'static_complete_parent_pdb':str(sp),'static_complete_parent_sha256':digest(sp),'fixed_parent_pdb':str(input_path(p)),'fixed_parent_sha256':digest(input_path(p)),'parent_backbone_comparison':bb,'n_U_all':len(u),'n_U_RMSF':sum(x['RMSF_valid'] for x in u),'n_U_Hchi':sum(x['Hchi_valid'] for x in u),'U':u}
        payload['templates'][template]=out
        dump_csv(HERE/'regions'/f'{template}_U8A.csv',u)
        dump(HERE/'regions'/f'{template}_U8A.json',out)
    assert set(templates)=={'2BNQ','2BNR','1AO7','1QSE','1QSF'}
    dump(HERE/'fixed_regions.json',payload)
    print(json.dumps({'status':payload['status'],'templates':{k:{x:v[x] for x in ['n_U_all','n_U_RMSF','n_U_Hchi']} for k,v in payload['templates'].items()}},indent=2))

def make_bio_structure(st,selected):
    structure=Structure.Structure('S');model=Model.Model(0);structure.add(model);chains={}
    for i,(key,res) in enumerate(st.items()):
        if key[0] not in selected:continue
        if key[0] not in chains:chains[key[0]]=Chain.Chain(key[0]);model.add(chains[key[0]])
        # Sequential internal residue ids avoid PDB author-id parsing assumptions.
        rr=Residue.Residue((' ',i+1,' '),res['resname'],' ');chains[key[0]].add(rr)
        for j,(name,a) in enumerate(res['atoms'].items()):
            assert a['element'] in RADII,a['element']
            rr.add(Atom.Atom(name,a['xyz'].copy(),0.,1.,' ',name.rjust(4),i*100+j,element=a['element']))
    return structure
def sasa(st,chains):
    entity=make_bio_structure(st,chains)
    engine=ShrakeRupley(probe_radius=1.4,n_points=960,radii_dict=RADII.copy())
    engine.compute(entity,level='A')
    return float(sum(atom.sasa for atom in entity.get_atoms()))
def contact_pairs(st,roles):
    tc={roles['tcr_alpha'],roles['tcr_beta']};xyz,keys=coords(st,tc);tree=cKDTree(xyz);pairs=set()
    for k,r in st.items():
        if k[0]!=roles['peptide']:continue
        for a in r['atoms'].values():
            for j in tree.query_ball_point(a['xyz'],r=4.5):pairs.add((k,keys[j]))
    return sorted(pairs)
def calculate_static(st,roles):
    tc={roles['tcr_alpha'],roles['tcr_beta']};pm={roles['hla_heavy'],roles['beta2m'],roles['peptide']}
    full=sasa(st,tc|pm);tcr=sasa(st,tc);pmhc=sasa(st,pm);pairs=contact_pairs(st,roles)
    return {'SASA_TCR_A2':tcr,'SASA_pMHC_A2':pmhc,'SASA_full_A2':full,'TCR_pMHC_BSA_A2':(tcr+pmhc-full)/2.,'peptide_TCR_contact_residue_pairs_4p5A':len(pairs),'contact_pairs':[[list(a),list(b)] for a,b in pairs]}

def load_predictions(npz_path,expected_structure):
    with np.load(npz_path,allow_pickle=False) as z:
        metrics=[str(x) for x in z['metrics']];assert metrics==['rmsf','chi_entropy'],metrics
        assert str(z['checkpoint_sha256'].item())==EXPECTED_CHECKPOINT
        means=np.array(z['mean'],dtype=np.float64);assert means.ndim==2 and means.shape[1]==2
        assert np.isfinite(means).all()
        chains=[str(x) for x in z['chain']];ids=[str(x) for x in z['author_residue_id']];aas=[str(x) for x in z['amino_acid']]
        frames=np.array(z['geometry_valid'],dtype=bool);target=np.array(z['target_mask'],dtype=bool);chem=np.array(z['hchi_chemically_defined'],dtype=bool)
    keys=list(zip(chains,ids));assert len(set(keys))==len(keys)==len(expected_structure)
    assert keys==list(expected_structure),('Prediction input residue order mismatch',npz_path)
    out={}
    for i,k in enumerate(keys):
        assert aas[i]==expected_structure[k]['aa']
        assert bool(frames[i])==has_frame(expected_structure[k]),('Frame mismatch',npz_path,k)
        assert bool(chem[i])==(CHI_COUNTS[aas[i]]>0)
        out[k]={'aa':aas[i],'geometry_valid':bool(frames[i]),'target':bool(target[i]),'hchi_defined':bool(chem[i]),'rmsf':float(means[i,0]),'hchi':float(means[i,1])}
    return out

def extract(args):
    protocol=read_json(args.protocol);f,s,fr,sindex,bindings=load_manifests(args)
    reg=read_json(args.regions)
    assert reg['protocol_sha256']==digest(args.protocol),'Protocol changed since U freeze; review/regenerate transparently before extraction'
    assert reg['fixed_input_bindings']==bindings,'Input manifests changed since U freeze'
    run=read_json(args.prediction_audit)
    assert run['status'] in {'complete','completed','passed'},'Wait until inference has completed'
    assert run['checkpoint_sha256']==EXPECTED_CHECKPOINT==protocol['base_model']['sha256']
    assert run['manifest_sha256']==bindings['fixed_manifest_sha256'],'Inference used a different fixed manifest'
    output_dir=Path(args.prediction_audit).resolve().parent
    predicted={input_id(r):r for r in records(run)}
    assert set(predicted)=={input_id(r) for r in fr},'Incomplete prediction coverage'
    matrix=substitution_matrices.load('BLOSUM62')
    cache_dir=HERE/'static_cache';cache_dir.mkdir(exist_ok=True)
    pdata={};sdata={};qdata={};inputs={}
    for r in fr:
        iid=input_id(r);roles=role_map(r);fixed,ignored=read_structure(input_path(r));inputs[iid]=fixed
        ppath=output_dir/iid/'prediction.npz';assert ppath.exists(),ppath
        assert predicted[iid]['input_sha256']==digest(input_path(r))
        pdata[iid]=load_predictions(ppath,fixed)
        sr=sindex[iid];sp=static_path(sr);st,extra=read_structure(sp)
        qc=structure_qc(st,roles,r['peptide_sequence']);bb=backbone_comparison(fixed,st)
        static_ok=qc['peptide_sequence_matches'] and qc['chain_roles_present'] and qc['all_canonical_heavy_atoms_complete'] and not qc['unsupported_atom_elements'] and not extra and bb['passed']
        qp={'static_complete_pdb':str(sp),'static_sha256':digest(sp),'fixed_input_sha256':digest(input_path(r)),'NPZ':str(ppath),'NPZ_sha256':digest(ppath),'static_structure_QC':qc,'backbone_comparison':bb,'fixed_input_nonprotein_atoms_excluded_like_model_parser_count':len(ignored),'fixed_input_nonprotein_atoms_first3':ignored[:3],'unexpected_nonprotein_heavy_atoms':extra,'static_eligible':static_ok}
        qdata[iid]=qp
        if static_ok:
            cp=cache_dir/(iid+'.json')
            expected_cache={'static_pdb_sha256':digest(sp),'radii_A':RADII,'SASA_probe_A':1.4,'SASA_n_points':960,'roles':roles,'contact_cutoff_A':4.5,'biopython_version':Bio.__version__}
            if cp.exists() and read_json(cp).get('settings')==expected_cache:calc=read_json(cp)['values']
            else:
                calc=calculate_static(st,roles);dump(cp,{'input_id':iid,'settings':expected_cache,'values':calc})
            sdata[iid]=calc
        print(json.dumps({'input':iid,'static_eligible':static_ok,'static_complete':iid in sdata}),flush=True)
    features=[];uroows=[];peptide_rows=[];role_summaries=[]
    for r in fr:
        iid=input_id(r);template=r['template_pdb'].upper();u=reg['templates'][template];parent=u['parent_input_id']
        assert parent in pdata
        roles=role_map(r);p=pdata[iid];ref=pdata[parent];fixed=inputs[iid];parent_fixed=inputs[parent]
        prep=[]
        if set(fixed)!=set(parent_fixed):prep.append('candidate_parent_residue_keys_differ')
        for key,res in parent_fixed.items():
            if key not in fixed:continue
            if key[0]!=roles['peptide'] and res['aa']!=fixed[key]['aa']:prep.append('environment_sequence_changed')
            for atom in BACKBONE:
                if atom in res['atoms'] and (atom not in fixed[key]['atoms'] or not np.array_equal(res['atoms'][atom]['xyz'],fixed[key]['atoms'][atom]['xyz'])):prep.append('candidate_parent_backbone_differ')
        if sequence(fixed,roles['peptide'])!=r['peptide_sequence']:prep.append('peptide_sequence_mismatch')
        pp=[k for k in fixed if k[0]==roles['peptide']]
        dynamic_errors=[]
        if len(pp)!=9 or not all(p[k]['geometry_valid'] and ref[k]['geometry_valid'] for k in pp):dynamic_errors.append('peptide_frames_incomplete')
        rvals=[];hvals=[];byrole=defaultdict(lambda:{'r':[],'h':[]})
        for item in u['U']:
            key=(item['chain'],item['author_residue_id'])
            if key not in p or key not in ref:dynamic_errors.append('U_residue_missing');continue
            if p[key]['aa']!=ref[key]['aa'] or p[key]['aa']!=item['amino_acid']:dynamic_errors.append('U_AA_changed');continue
            dr=p[key]['rmsf']-ref[key]['rmsf'];dh=p[key]['hchi']-ref[key]['hchi']
            rv=true(item['RMSF_valid']);hv=true(item['Hchi_valid'])
            if rv:
                if not (p[key]['geometry_valid'] and ref[key]['geometry_valid']):dynamic_errors.append('U_frame_missing')
                rvals.append(dr);byrole[item['role']]['r'].append(dr)
            if hv:
                if not (p[key]['hchi_defined'] and ref[key]['hchi_defined']):dynamic_errors.append('U_hchi_undefined')
                hvals.append(dh);byrole[item['role']]['h'].append(dh)
            uroows.append({'candidate_id':biological_id(r),'input_id':iid,'template_pdb':template,'feature_region_id':u['feature_region_id'],**item,'parent_rmsf_A':ref[key]['rmsf'],'candidate_rmsf_A':p[key]['rmsf'],'delta_rmsf_A':dr,'parent_hchi_nats':ref[key]['hchi'],'candidate_hchi_nats':p[key]['hchi'],'delta_hchi_nats':dh,'used_in_R':rv,'used_in_H':hv})
        assert len(rvals)==u['n_U_RMSF'] and len(hvals)==u['n_U_Hchi']
        for i,k in enumerate(pp,start=1):
            peptide_rows.append({'candidate_id':biological_id(r),'input_id':iid,'template_pdb':template,'peptide_position':i,'chain':k[0],'author_residue_id':k[1],'parent_aa':ref[k]['aa'],'candidate_aa':p[k]['aa'],'parent_rmsf_A':ref[k]['rmsf'],'candidate_rmsf_A':p[k]['rmsf'],'delta_rmsf_A':p[k]['rmsf']-ref[k]['rmsf'],'parent_hchi_nats':ref[k]['hchi'] if ref[k]['hchi_defined'] else '', 'candidate_hchi_nats':p[k]['hchi'] if p[k]['hchi_defined'] else '', 'same_aa_hchi_delta_nats':p[k]['hchi']-ref[k]['hchi'] if p[k]['aa']==ref[k]['aa'] and p[k]['hchi_defined'] else '', 'Hchi_peptide_is_archive_only_not_primary_feature':True})
        for role,vals in byrole.items():role_summaries.append({'candidate_id':biological_id(r),'input_id':iid,'template_pdb':template,'role':role,'n_R':len(vals['r']),'mean_delta_RMSF_A':float(np.mean(vals['r'])) if vals['r'] else '', 'n_H':len(vals['h']),'mean_delta_Hchi_nats':float(np.mean(vals['h'])) if vals['h'] else '', 'archive_only':True})
        parent_seq=r['panel_reference_peptide'];seq=r['peptide_sequence'];assert len(seq)==len(parent_seq)==9
        static_ok=qdata[iid]['static_eligible'] and qdata[parent]['static_eligible']
        dynamic_ok=not dynamic_errors and not prep
        R=float(np.mean(rvals));H=float(np.mean(hvals))
        nmut=sum(a!=b for a,b in zip(parent_seq,seq));blosum=float(sum(matrix[a,a]-matrix[a,b] for a,b in zip(parent_seq,seq)))
        dbsa=sdata[iid]['TCR_pMHC_BSA_A2']-sdata[parent]['TCR_pMHC_BSA_A2'] if static_ok else ''
        dcontact=sdata[iid]['peptide_TCR_contact_residue_pairs_4p5A']-sdata[parent]['peptide_TCR_contact_residue_pairs_4p5A'] if static_ok else ''
        if true(r.get('is_panel_parent')):assert R==H==0 and nmut==0 and blosum==0 and (not static_ok or dbsa==dcontact==0)
        features.append({'candidate_id':biological_id(r),'panel_id':r['panel_id'],'template_pdb':template,'input_id':iid,'peptide_sequence':seq,'parent_candidate_id':u['parent_biological_candidate_id'],'parent_input_id':parent,'template_is_primary':true(value(r,'is_primary_template','primary')),'is_panel_parent':true(r.get('is_panel_parent')),'sequence_n_mutations':nmut,'BLOSUM62_deficit':blosum,'delta_TCR_pMHC_BSA_A2':dbsa,'delta_peptide_TCR_contact_residue_pairs_4p5A':dcontact,'R':R,'H':H,'static_eligible':static_ok,'preparation_qc_pass':not prep,'dynamic_qc_pass':dynamic_ok,'feature_region_id':u['feature_region_id'],'n_U_RMSF':len(rvals),'n_U_Hchi':len(hvals),'n_U_all':u['n_U_all'],'checkpoint_sha256':EXPECTED_CHECKPOINT,'fixed_input_sha256':qdata[iid]['fixed_input_sha256'],'static_complete_sha256':qdata[iid]['static_sha256'],'prediction_npz_sha256':qdata[iid]['NPZ_sha256'],'QC_reasons':'|'.join(sorted(set(prep+dynamic_errors)))})
        dump(HERE/'input_QA'/(iid+'.json'),{'input_id':iid,'candidate_id':biological_id(r),'template_pdb':template,'static':qdata[iid],'preparation_errors':sorted(set(prep)),'dynamic_errors':sorted(set(dynamic_errors)),'region_id':u['feature_region_id']})
    assert len(features)==76 and sum(x['template_is_primary'] for x in features)==30
    assert len({(x['candidate_id'],x['template_pdb']) for x in features})==76
    dump_csv(HERE/'candidate_features.csv',features)
    dump_csv(HERE/'archive_U_residue_profiles.csv',uroows)
    dump_csv(HERE/'archive_peptide_profiles.csv',peptide_rows)
    dump_csv(HERE/'archive_U_role_summaries.csv',role_summaries)
    report={'status':'passed' if all(x['preparation_qc_pass'] and x['dynamic_qc_pass'] and x['static_eligible'] for x in features) else 'complete_with_recorded_QC_exclusions','created_at_utc':utc(),'checkpoint_sha256':EXPECTED_CHECKPOINT,'protocol_id':protocol['protocol_id'],'protocol_sha256':digest(args.protocol),'fixed_regions_sha256':digest(args.regions),'prediction_audit':str(Path(args.prediction_audit).resolve()),'prediction_audit_sha256':digest(args.prediction_audit),'input_bindings':bindings,'n_inputs':len(features),'n_primary_candidates':30,'per_template_counts':{k:sum(x['template_pdb']==k for x in features) for k in reg['templates']},'n_static_eligible':sum(x['static_eligible'] for x in features),'n_dynamic_eligible':sum(x['dynamic_qc_pass'] for x in features),'all_parent_features_zero':True,'SASA_settings':{'method':'Bio.PDB.SASA.ShrakeRupley','biopython_version':Bio.__version__,'probe_radius_A':1.4,'n_points':960,'element_radii_A':RADII,'definition':'[SASA(TCR)+SASA(pMHC incl beta2m)-SASA(full)]/2'},'contact_definition':'unique peptide/TCR heavy-atom residue pairs with distance<=4.5A','dynamic_definition':'mean of candidate-minus-same-template-parent per-residue physical-scale means in fixed U; no expm1 or variance weighting','static_source_policy':'Complete canonical side-chain PDBs only; deleted-sidechain inference inputs never used as SASA/contact structures','labels_or_outcomes_read':False,'scripts':{'extract_features.py':digest(__file__)},'files':{n:{'path':str(HERE/n),'sha256':digest(HERE/n)} for n in ['candidate_features.csv','fixed_regions.json','archive_U_residue_profiles.csv','archive_peptide_profiles.csv','archive_U_role_summaries.csv']}}
    dump(HERE/'feature_manifest.json',report)
    print(json.dumps({k:report[k] for k in ['status','n_inputs','n_primary_candidates','n_static_eligible','n_dynamic_eligible','per_template_counts']},indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['freeze-regions','extract'])
    p.add_argument('--protocol',default=str(ROOT/'evaluation_review/fixed_analysis_protocol.json'))
    p.add_argument('--fixed-manifest',default=str(ROOT/'structures/fixed_backbone_manifest.json'))
    p.add_argument('--static-manifest',default=str(ROOT/'structures/static_complete_manifest.json'))
    p.add_argument('--regions',default=str(HERE/'fixed_regions.json'))
    p.add_argument('--prediction-audit',default=str(ROOT/'predictions/fixed_backbone/run_audit.json'))
    p.add_argument('--runtime-path',action='append',default=[],help='Existing dependency package directory; repeatable, no installs.')
    args=p.parse_args();dependencies(args.runtime_path)
    if args.stage=='freeze-regions':freeze_regions(args)
    else:extract(args)

if __name__=='__main__':main()
