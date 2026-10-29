"""Independent read-only fixed-backbone/sequence/identity audit; no inference."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import argparse,csv,hashlib,json,sys
import numpy as np

HERE=Path(__file__).resolve().parent
SRC=Path(r'E:\Project\FlexPPI\src')
sys.path.insert(0,str(SRC))
from flexppi.data.structures import parse_pdb

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',default=str(HERE.parent/'structures/fixed_backbone_manifest.json'))
    parser.add_argument('--output',default=str(HERE/'fixed_backbone_inputs_independent_audit.json'))
    args=parser.parse_args()
    manifest=read(args.manifest)
    assert manifest['candidate_count']==76 and manifest['biological_candidate_count']==30
    protocol=read(HERE/'fixed_analysis_protocol.json')
    gate=read(HERE/'protocol_gate_audit.json')
    cross={r['label_candidate_id']:r for r in gate['candidate_identity_crosswalk']}
    templateinfo={r['template_pdb']:r for r in manifest['templates']}
    raw_audit=read(HERE/'raw_template_overlap_audit.json')
    raw_sha={r['pdb'].upper():r['source_sha256'] for r in raw_audit['rows']}
    rows=[];seen=set();parents=Counter()
    role_mapping={'hla_heavy':'A','beta2m':'B','peptide':'C','tcr_alpha':'D','tcr_beta':'E'}
    for entry in manifest['candidates']:
        cid=entry['biological_candidate_id'];t=entry['template_pdb']
        info=cross[cid];key=cid,t
        assert key not in seen;seen.add(key)
        assert t in {'nyeso':['2BNQ','2BNR'],'tax':['1AO7','1QSE','1QSF']}[info['system']]
        assert entry['roles']==role_mapping
        assert entry['peptide_sequence']==info['sequence']
        assert entry['labels_or_closest_template_used_for_preparation']==False
        raw=Path(entry['template_input_pdb']);prepared=Path(entry['prepared_input_pdb'])
        assert sha(raw)==entry['template_input_sha256']==raw_sha[t]
        assert sha(prepared)==entry['prepared_input_sha256']
        a,b=parse_pdb(str(raw)),parse_pdb(str(prepared))
        assert a.chain_ids==b.chain_ids and a.author_residue_ids==b.author_residue_ids
        assert np.array_equal(a.backbone,b.backbone) and np.array_equal(a.backbone_mask,b.backbone_mask)
        assert np.array_equal(a.frame_mask,b.frame_mask)
        assert b.frame_mask.all()
        assert b.num_residues==entry['prepared_input_residues']==entry['geometry_valid_count']
        assert b.chain_sequence(entry['peptide_chain'])==entry['peptide_sequence']
        assert len(b.chain_sequence(entry['peptide_chain']))==9
        for chain in b.unique_chains():
            if chain!=entry['peptide_chain']:
                assert a.chain_sequence(chain)==b.chain_sequence(chain)
        changed=[i for i,(aa,bb) in enumerate(zip(a.sequence,b.sequence)) if aa!=bb]
        assert all(b.chain_ids[i]==entry['peptide_chain'] for i in changed)
        assert all(aa!='X' for aa in b.sequence)
        parentid=entry['parent_biological_candidate_id']
        assert parentid in cross and cross[parentid]['analysis_set']=='reference_only'
        assert entry['parent_inference_candidate_id']==t+'__'+parentid
        if entry['is_panel_parent']:
            assert cid==parentid;parents[t]+=1
        rows.append(dict(input_id=entry['candidate_id'],biological_candidate_id=cid,template_pdb=t,
            prepared_input_pdb=str(prepared),prepared_input_sha256=sha(prepared),raw_input_sha256=sha(raw),
            status='passed',residues=b.num_residues,peptide_residues=9,valid_N_CA_C_frames=int(b.frame_mask.sum()),
            chain_roles_verified=True,all_original_backbone_coordinates_masks_and_residue_order_identical=True,
            no_added_loops_or_residues=True,no_nonpeptide_sequence_substitutions=True,
            peptide_sequence_matches_frozen_candidate=True,labels_or_closest_template_used=False,
            parent_input_id=entry['parent_inference_candidate_id'],is_parent=entry['is_panel_parent'],
            sequence_train_overlap_scope='Unchanged nonpeptide chains inherit actual raw-template direct-overlap audit; changed9peptides inherit exact candidate-vs-actual-training check, including addedA00IL.' ))
    assert len(seen)==76 and len({c for c,t in seen})==30
    assert parents==Counter({'2BNQ':1,'2BNR':1,'1AO7':1,'1QSE':1,'1QSF':1})
    actual={(c,t) for c,t in seen}
    expected={(cid,t) for cid,info in cross.items() for t in {'nyeso':['2BNQ','2BNR'],'tax':['1AO7','1QSE','1QSF']}[info['system']]}
    assert actual==expected
    checkpoint=Path(protocol['base_model']['local_run_path'])
    assert sha(checkpoint)==protocol['base_model']['sha256']
    result=dict(status='passed_fixed_backbone_inputs_static_completion_pending',reviewed_at_utc=datetime.now(timezone.utc).isoformat(),
        manifest_path=args.manifest,manifest_sha256=sha(args.manifest),protocol_sha256=sha(HERE/'fixed_analysis_protocol.json'),
        execution_checkpoint_path=str(checkpoint),execution_checkpoint_sha256=sha(checkpoint),
        n_inputs=76,n_biological_candidates=30,n_templates=5,parents_prepared=5,
        all76_actual_parsed_sequences_backbones_and_masks_verified=True,
        exact_input_representations_audited=True,complete_template_expansion_verified=True,
        inherited_raw_template_audit_sha256=sha(HERE/'raw_template_overlap_audit.json'),
        inherited_training_audit_sha256=sha(HERE/'training_overlap_audit.json'),
        actual_training_transitive_family_overlap_records=1,
        not_completely_family_disjoint=True,
        sidechain_static_completion_final_audit_pending=True,
        fixed_environment_U_final_definition_pending=True,
        no_candidate_predictions_inspected=True,no_training_or_inference_performed=True,
        rows=rows)
    Path(args.output).write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:result[k] for k in ['status','n_inputs','n_biological_candidates','n_templates','execution_checkpoint_sha256','sidechain_static_completion_final_audit_pending']},ensure_ascii=False))

if __name__=='__main__':main()
