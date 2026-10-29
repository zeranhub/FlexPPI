"""Independent heavy-atom/QC review; never reads prediction arrays or SPR values."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import argparse,hashlib,json,sys
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,r'E:\Project\FlexPPI\src')
from flexppi.data.structures import parse_pdb
from flexppi.constants import THREE_TO_ONE as MODEL_THREE_TO_ONE

SIDECHAINS={
 'ALA':'CB','ARG':'CB CG CD NE CZ NH1 NH2','ASN':'CB CG OD1 ND2','ASP':'CB CG OD1 OD2',
 'CYS':'CB SG','GLN':'CB CG CD OE1 NE2','GLU':'CB CG CD OE1 OE2','GLY':'',
 'HIS':'CB CG ND1 CD2 CE1 NE2','ILE':'CB CG1 CG2 CD1','LEU':'CB CG CD1 CD2',
 'LYS':'CB CG CD CE NZ','MET':'CB CG SD CE','PHE':'CB CG CD1 CD2 CE1 CE2 CZ','PRO':'CB CG CD',
 'SER':'CB OG','THR':'CB OG1 CG2','TRP':'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
 'TYR':'CB CG CD1 CD2 CE1 CE2 CZ OH','VAL':'CB CG1 CG2'}
AA3=dict(zip(SIDECHAINS, ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']))
CHI=dict(zip('ARNDCQEGHILKMFPSTWYV',[0,4,2,2,1,3,3,0,2,2,2,4,3,2,2,1,1,2,2,1]))
assert len(SIDECHAINS)==len(AA3)==len(CHI)==20
assert AA3=={name:MODEL_THREE_TO_ONE[name] for name in SIDECHAINS}
assert AA3['LEU']=='L' and CHI['L']==2
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))

def heavy_atoms(path, allow_ignored_nonprotein=False):
    residues={}
    for line in Path(path).read_text().splitlines():
        if line[:6].strip()=='ENDMDL':break
        if not line.startswith(('ATOM  ','HETATM')) or len(line)<54 or line[16:17] not in {' ','A','1'}:continue
        name=line[17:20].strip();atom=line[12:16].strip()
        element=line[76:78].strip().upper() if len(line)>=78 else ''
        if not element:element=next((c for c in atom if c.isalpha()),'').upper()
        if element in {'H','D'}:continue
        if name not in SIDECHAINS:
            if allow_ignored_nonprotein or name in {'HOH','WAT','DOD'}:continue
            raise ValueError(('Unknown heavy residue',str(path),name))
        key=(line[21:22].strip() or '_',line[22:26].strip()+line[26:27].strip())
        residues.setdefault(key,{'resname':name,'aa':AA3[name],'atoms':{}})
        assert residues[key]['resname']==name
        residues[key]['atoms'].setdefault(atom,np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])],dtype=np.float64))
    return residues

def expected_atoms(residue):
    return {'N','CA','C','O'}|set(SIDECHAINS[residue['resname']].split())

def independent_U(parent,roles):
    pep=np.array([v for k,r in parent.items() if k[0]==roles['peptide'] for v in r['atoms'].values()])
    environment={roles['hla_heavy'],roles['tcr_alpha'],roles['tcr_beta']}
    rev={v:k for k,v in roles.items()}
    u=[]
    for key,r in parent.items():
        if key[0] not in environment:continue
        xyz=np.array(list(r['atoms'].values()))
        minimum=float(np.sqrt(np.min(np.sum((xyz[:,None,:]-pep[None,:,:])**2,axis=-1))))
        if minimum<=8.0:
            frame={'N','CA','C'}<=set(r['atoms'])
            complete=expected_atoms(r)<=set(r['atoms'])
            u.append(dict(chain=key[0],author_residue_id=key[1],amino_acid=r['aa'],role=rev[key[0]],
                reference_peptide_distance_min_A=minimum,RMSF_valid=frame,Hchi_valid=frame and complete and CHI[r['aa']]>0,
                n_chi=CHI[r['aa']],canonical_heavy_complete=complete))
    return u

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',default=str(ROOT/'structures/static_complete_manifest.json'))
    parser.add_argument('--allow-partial',action='store_true')
    parser.add_argument('--regions',help='Final frozen feature_extraction/fixed_regions.json to independently recompute and verify')
    parser.add_argument('--output')
    args=parser.parse_args()
    fixedpath=ROOT/'structures/fixed_backbone_manifest.json'
    fixed=read(fixedpath);frows={r['candidate_id']:r for r in fixed['candidates']}
    live=read(args.manifest)
    complete=live['status']=='passed' and len(live['candidates'])==76
    if not args.allow_partial:assert complete,'Static preparation not yet finalized'
    assert live['input_manifest_sha256']==sha(fixedpath)
    assert live['preparation_script_sha256']==sha(ROOT/'prepare_complete_structures.py')
    rule=read(ROOT/'structure_preparation_protocol.json')
    assert rule['method_fixed_before_candidate_outputs'] and rule['local_only'] and not rule['new_MD']
    assert set(rule['primary_templates'])=={'2BNQ','1AO7'}
    assert set(rule['sensitivity_templates'])=={'2BNR','1QSE','1QSF'}
    rows=[];parents={};rebuilt=0
    for row in live['candidates']:
        iid=row['candidate_id'];old=frows[iid]
        beforepath=Path(old['prepared_input_pdb']);afterpath=Path(row['static_pdb'])
        assert Path(row['input_pdb'])==afterpath
        audit=read(row['preparation_audit'])
        assert audit['status']=='passed' and audit['seed']==17 and audit['OpenMM']=='8.6.1' and not audit['MD']
        assert audit['backbone_rebuilt']==False and audit['all_observed_backbone_values_bit_identical']
        assert sha(beforepath)==old['prepared_input_sha256']==audit['input_sha256']
        assert sha(afterpath)==audit['output_sha256'] and Path(audit['output'])==afterpath
        a,b=parse_pdb(str(beforepath)),parse_pdb(str(afterpath))
        assert a.sequence==b.sequence and a.chain_ids==b.chain_ids and a.author_residue_ids==b.author_residue_ids
        assert np.array_equal(a.backbone,b.backbone) and np.array_equal(a.backbone_mask,b.backbone_mask)
        assert b.frame_mask.all()
        oldatoms,newatoms=heavy_atoms(beforepath,allow_ignored_nonprotein=True),heavy_atoms(afterpath)
        assert list(oldatoms)==list(newatoms)
        old_heavy_count=0;new_heavy_count=0;added=0
        for key,residue in oldatoms.items():
            other=newatoms[key]
            assert residue['aa']==other['aa']
            assert expected_atoms(other)<=set(other['atoms']),(iid,key,expected_atoms(other)-set(other['atoms']))
            assert set(residue['atoms'])<=set(other['atoms'])
            for atom,xyz in residue['atoms'].items():
                assert np.array_equal(xyz,other['atoms'][atom]),(iid,key,atom,'Existing heavy atom moved')
            old_heavy_count+=len(residue['atoms']);new_heavy_count+=len(other['atoms'])
            added+=len(set(other['atoms'])-set(residue['atoms']))
        rebuilt+=added
        if row['is_panel_parent']:
            assert row['template_pdb'] not in parents
            parents[row['template_pdb']]=dict(row=row,residues=newatoms)
        rows.append(dict(input_id=iid,biological_candidate_id=row['biological_candidate_id'],template_pdb=row['template_pdb'],
            status='passed',input_path=str(beforepath),input_sha256=sha(beforepath),static_path=str(afterpath),static_sha256=sha(afterpath),
            preparation_audit_sha256=sha(row['preparation_audit']),residues=b.num_residues,
            source_heavy_atoms=old_heavy_count,complete_heavy_atoms=new_heavy_count,added_heavy_atoms=added,
            all_existing_protein_heavy_atom_coordinates_exactly_preserved=True,all_NCAO_masks_positions_and_author_ids_identical=True,
            canonical_heavy_atoms_complete=True,missing_residue_loops_added=False,seed=17,openmm_version='8.6.1'))
    assert len({r['input_id'] for r in rows})==len(rows)
    if complete:assert set(frows)=={r['input_id'] for r in rows} and len(parents)==5
    u_review=[]
    if args.regions:
        assert complete,'Do not final-review U before all76 structural inputs finalized'
        frozen=read(args.regions)
        assert frozen['fixed_input_bindings']['fixed_manifest_sha256']==sha(fixedpath)
        bound_static=Path(frozen['fixed_input_bindings']['static_manifest'])
        assert frozen['fixed_input_bindings']['static_manifest_sha256']==sha(bound_static)
        bound=read(bound_static)
        bound_rows={r['candidate_id']:r for r in bound['candidates']}
        actual_rows={r['candidate_id']:r for r in live['candidates']}
        assert set(bound_rows)==set(actual_rows)
        for iid,entry in bound_rows.items():
            explicit=entry.get('static_complete_pdb') or entry.get('static_pdb') or entry.get('complete_pdb') or entry.get('output_pdb')
            assert explicit and Path(explicit)==Path(actual_rows[iid]['static_pdb'])
        # U may bind the immutable explicit-path index rather than the live
        # completion journal; both must resolve to the same independently
        # verified76 final files, never to old deleted-sidechain inputs.
        assert frozen['protocol_sha256']==sha(HERE/'fixed_analysis_protocol.json')
        for template,stored in frozen['templates'].items():
            parent=parents[template];row=parent['row']
            assert stored['parent_input_id']==row['candidate_id']
            assert stored['static_complete_parent_sha256']==sha(row['static_pdb'])
            independently=independent_U(parent['residues'],row['roles'])
            assert len(independently)==len(stored['U'])
            for own,given in zip(independently,stored['U']):
                for key in ['chain','author_residue_id','amino_acid','role','RMSF_valid','Hchi_valid','n_chi','canonical_heavy_complete']:
                    assert own[key]==given[key],(template,key,own,given)
                assert abs(own['reference_peptide_distance_min_A']-given['reference_peptide_distance_min_A'])<1e-12
            assert stored['n_U_all']==len(independently)
            assert stored['n_U_RMSF']==sum(r['RMSF_valid'] for r in independently)
            assert stored['n_U_Hchi']==sum(r['Hchi_valid'] for r in independently)
            u_review.append(dict(template_pdb=template,status='passed',parent_input_id=row['candidate_id'],
                complete_parent_sha256=sha(row['static_pdb']),n_U_all=len(independently),n_U_RMSF=stored['n_U_RMSF'],n_U_Hchi=stored['n_U_Hchi'],
                role_counts=dict(Counter(r['role'] for r in independently)),all_U_positions_masks_and_distances_independently_recomputed=True))
        assert len(u_review)==5
    status='passed_all76_static_inputs_and_frozen_U' if complete and args.regions else 'passed_all76_static_inputs_U_pending' if complete else 'partial_inputs_reviewed_pending_completion'
    result=dict(status=status,reviewed_at_utc=datetime.now(timezone.utc).isoformat(),manifest_path=args.manifest,manifest_sha256=sha(args.manifest),
        preparation_protocol_sha256=sha(ROOT/'structure_preparation_protocol.json'),preparation_script_sha256=sha(ROOT/'prepare_complete_structures.py'),
        pdbfixer_source_sha256=sha(ROOT/'tool_runtime/pdbfixer/pdbfixer.py'),
        original_fixed_manifest_sha256=sha(fixedpath),n_reviewed=len(rows),n_templates_in_review=len({r['template_pdb'] for r in rows}),
        n_added_heavy_atoms=rebuilt,n_completed_parent_templates=len(parents),all_existing_protein_heavy_coordinates_exactly_unchanged=True,
        all_canonical_residue_heavy_atoms_complete=True,model_input_AA_NCAO_and_masks_unchanged=True,
        all20_amino_acid_names_independently_matched_to_actual_model_parser=True,
        no_SPR_values_or_model_predictions_read=True,no_inference_training_or_remote_access=True,
        technical_method_note='PDBFixer freezes existing atoms and optimizes only new atoms. Its clash repair may internally integrate new atoms for up to10x200 steps; no protein-production MD, no backbone relaxation, no saved dynamical trajectory.',
        U_review=u_review,regions_path=args.regions,regions_sha256=sha(args.regions) if args.regions else None,rows=rows)
    if args.regions:
        result['U_bound_static_manifest_path']=str(bound_static)
        result['U_bound_static_manifest_sha256']=sha(bound_static)
        result['immutable_static_index_and_completed_journal_resolve_to_same76_files']=True
    output=Path(args.output) if args.output else HERE/('static_complete_inputs_independent_audit.json' if complete else 'static_complete_inputs_partial_review.json')
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:result[k] for k in ['status','n_reviewed','n_templates_in_review','n_added_heavy_atoms','n_completed_parent_templates']},ensure_ascii=False))

if __name__=='__main__':main()
