from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
OUT=Path(__file__).resolve().parent;BASE=OUT.parent.parent
OLD=BASE/'local_antigen_screening_20261021'
MAN=OLD/'structures/fixed_backbone_manifest.json'
manifest=json.loads(MAN.read_text(encoding='utf8'))
labels=pd.read_csv(OLD/'experimental_labels/SPR_labels_candidates_reported_hybrid.csv').set_index('candidate_id')
sources=[];cache={};rows=[]
def load(cid,ch):
    if cid not in cache:
        p=OLD/'predictions/fixed_backbone'/cid/'prediction.npz';z=np.load(p,allow_pickle=False)
        assert z['metrics'].tolist()==['rmsf','chi_entropy']
        assert str(z['checkpoint_sha256'])=='b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f'
        mask=z['chain']==ch;ids=np.flatnonzero(mask);assert len(ids)==9
        cache[cid]={k:z[k][ids] for k in ['mean','chain','author_residue_id','amino_acid','geometry_valid','hchi_chemically_defined']}
        sources.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    return cache[cid]
for c in manifest['candidates']:
    a=load(c['parent_inference_candidate_id'],c['peptide_chain']);b=load(c['candidate_id'],c['peptide_chain'])
    assert ''.join(b['amino_acid'])==c['peptide_sequence']
    assert np.array_equal(a['author_residue_id'],b['author_residue_id'])
    for i in range(9):
        same=str(a['amino_acid'][i])==str(b['amino_acid'][i])
        hv=bool(a['hchi_chemically_defined'][i] and b['hchi_chemically_defined'][i])
        rows.append({'inference_candidate_id':c['candidate_id'],'candidate_id':c['biological_candidate_id'],'panel_id':c['panel_id'],
                     'template':c['template_pdb'],'primary_template':c['primary'],'mutation':c['mutations_vs_panel_reference'],
                     'peptide_position':i+1,'peptide_author_id':str(b['author_residue_id'][i]),'parent_aa':str(a['amino_acid'][i]),'candidate_aa':str(b['amino_acid'][i]),
                     'same_AA':same,'hchi_defined_both':hv,'hchi_environment_comparable':hv and same,
                     'parent_RMSF_A':float(a['mean'][i,0]),'candidate_RMSF_A':float(b['mean'][i,0]),'RMSF_drop_A':float(a['mean'][i,0])-float(b['mean'][i,0]),
                     'parent_Hchi_nats':float(a['mean'][i,1]) if a['hchi_chemically_defined'][i] else np.nan,
                     'candidate_Hchi_nats':float(b['mean'][i,1]) if b['hchi_chemically_defined'][i] else np.nan,
                     'same_AA_Hchi_drop_nats':float(a['mean'][i,1])-float(b['mean'][i,1]) if hv and same else np.nan,
                     'delta_pKD_reported_hybrid':float(labels.loc[c['biological_candidate_id'],'delta_pKD_vs_parent']),
                     'checkpoint_sha256':'b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f',
                     'input_type':'Fixed parent N CA C O backbone; changed residue identity; no physical mutant relaxation'})
pd.DataFrame(rows).to_csv(OUT/'fixed_backbone_76_predictions_all_peptide_RMSF_Hchi.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
(OUT/'fixed_backbone_peptide_export_audit.json').write_text(json.dumps({'status':'passed','source_manifest':str(MAN),'source_manifest_sha256':hashlib.sha256(MAN.read_bytes()).hexdigest(),'prediction_n':len(cache),'rows':len(rows),'new_inference':0,'sources':sources},ensure_ascii=False,indent=2),encoding='utf8')
print('PASS',len(cache),len(rows))
