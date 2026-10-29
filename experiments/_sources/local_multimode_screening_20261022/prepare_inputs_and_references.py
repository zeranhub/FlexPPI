"""Build label-free fixed-backbone inputs, static features and training ECDFs."""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import argparse,csv,json,hashlib,sys,shutil
import numpy as np
sys.path.insert(0,r'E:\Project\FlexPPI\src')
from flexppi.data.structures import parse_pdb
from flexppi.constants import THREE_TO_ONE,CHI_ATOMS

ROOT=Path(__file__).resolve().parent
REF=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
def readcsv(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def csvout(p,rr):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)

def atoms(path):
    result=defaultdict(list)
    for line in path.read_text().splitlines():
        if line[:6].strip() not in {'ATOM','HETATM'} or len(line)<78:continue
        if line[16] not in {' ','A','1'} or line[76:78].strip() in {'H','D'}:continue
        result[line[21],line[22:26].strip()+line[26].strip()].append((line[12:16].strip(),[float(line[s:s+8]) for s in (30,38,46)]))
    return result

def prepare():
    assert not list((ROOT/'predictions').glob('**/prediction.npz'))
    candidates=readcsv(ROOT/'antigen_labels/curated/curated_candidates.csv')
    downloads=json.loads((ROOT/'antigen_labels/downloads.json').read_text())['files']
    temps=[];entries=[];audit=[];static=[]
    for i in range(1,8):
        panel=f'MSKCC_TCR{i}';parentid=panel+'__parent';tid=panel+'__AF3'
        rr=[r for r in candidates if r['panel_id']==panel];assert len(rr)==171
        f=next(x for x in downloads if x['repository_path'].endswith(f'mskcc_tcr{i}_ec50_sat_mut_af3_no_template.csv'))
        src=readcsv(f['path']);p=ROOT/'antigen_labels/raw/pdbs'/f"{src[0]['wt_pdb']}.pdb"
        out=ROOT/'structures'/p.name;out.parent.mkdir(exist_ok=True);shutil.copy2(p,out)
        st=parse_pdb(str(out));assert st.unique_chains()==list('ABCD') and st.chain_sequence('B')==src[0]['wt_seq']
        assert st.backbone_mask.all() and set(st.sequence)<=set('ACDEFGHIKLMNPQRSTVWY')
        # These source AF3 inputs contain MHC peptide-binding domain, peptide,
        # TCR alpha variable domain and beta variable domain; no beta2m.
        roles=dict(hla_heavy='A',peptide='B',tcr_alpha='C',tcr_beta='D')
        aa=atoms(out);pep=np.asarray([xyz for key,values in aa.items() if key[0]=='B' for name,xyz in values])
        env=[];region=[]
        for j,(c,rid,a) in enumerate(zip(st.chain_ids,st.author_residue_ids,st.sequence)):
            coords=np.asarray([xyz for name,xyz in aa[c,rid]])
            distance=float(np.sqrt(np.min(np.sum((coords[:,None,:]-pep[None,:,:])**2,axis=2))))
            if c!='B' and distance<=8.:
                env.append(dict(chain=c,author_residue_id=rid))
                region.append(dict(chain=c,author_residue_id=rid,amino_acid=a,minimum_parent_heavy_atom_distance_A=distance,
                    R_included=True,H_included=a not in {'G','A'}))
        t=dict(template_id=tid,input_pdb=str(out),input_sha256=sha(out),roles=roles,parent_candidate_id=parentid,
            environment_residues=env,environment_rule='all unchanged nonpeptide spectator residues with any parent heavy atom within8A of peptide; frozen without labels or predictions')
        temps.append(t)
        entries.append(dict(candidate_id=parentid,template_id=tid,chain_sequences={'B':st.chain_sequence('B')},
            parent_candidate_id=parentid,is_parent=True,biological_candidate_id='parent',panel_id=panel,mode='parent'))
        for r in rr:
            entries.append(dict(candidate_id=r['candidate_id'],template_id=tid,chain_sequences={'B':r['sequence']},
                parent_candidate_id=parentid,is_parent=False,biological_candidate_id=r['mutation'],panel_id=panel,mode='fixed_backbone_antigen_mutant'))
        removed=ROOT/'structures'/f'{panel}__TCR_removed.pdb'
        removed.write_text(''.join(line+'\n' for line in out.read_text().splitlines() if line[:6].strip() in {'ATOM','HETATM','TER'} and line[21:22] in {'A','B'})+'END\n')
        rt=tid+'__TCR_removed';rp=panel+'__TCR_removed_parent'
        temps.append(dict(template_id=rt,input_pdb=str(removed),input_sha256=sha(removed),roles=dict(hla_heavy='A',peptide='B'),
            parent_candidate_id=rp,environment_residues=[e for e in env if e['chain']=='A'],
            computational_ablation='TCR alpha/beta chains removed; HLA and peptide geometry unchanged. Not physical apo.'))
        entries.append(dict(candidate_id=rp,template_id=rt,chain_sequences={'B':st.chain_sequence('B')},parent_candidate_id=rp,
            is_parent=True,biological_candidate_id='TCR_removed_parent',panel_id=panel,mode='computational_TCR_context_ablation'))
        dump(ROOT/'structures'/f'{panel}__regions.json',dict(panel_id=panel,template_id=tid,peptide_chain='B',spectator_regions=region,
            all_variable_component_positions_excluded=True,predictions_seen=False,endpoint_labels_used_for_regions=False))
        audit.append(dict(panel_id=panel,template=out.name,template_sha256=sha(out),residues=st.num_residues,
            roles=roles,chain_sequences={c:st.chain_sequence(c) for c in st.unique_chains()},complete_NCAO=True,
            environment_R=len(env),environment_H=sum(x['H_included'] for x in region),candidates=171,
            template_scope='Source AF3 truncated MHC+TCR variable domains; not full-length experimental complexes',
            role_evidence='MHC175aa peptide-binding domain; source9mer exact peptide; alpha variable IG/TRAV motifs inC and TRBV/CASS motifs inD.'))
        chosen={}
        for method in ['hermes_py_000','proteinmpnn_v_48_020_full_pep_masked','proteinmpnn_v_48_002_full_pep_masked']:
            pattern=f'mskcc_tcr{i}_ec50_sat_mut_af3_no_template-'
            ff=[x for x in downloads if '/results/'+method+'/' in x['repository_path'] and pattern in x['repository_path']]
            assert len(ff)==1
            ss=readcsv(ff[0]['path']);assert len(ss)==180 and all(x['wt_pdb']==src[0]['wt_pdb'] for x in ss)
            wt=[float(x['pnlogp']) for x in ss if x['is_wt']=='True'];assert len(wt)==9
            chosen[method]=(ss,float(np.mean(wt)),ff[0])
        for r in rr:
            result=dict(panel_id=panel,candidate_id=r['candidate_id'],sequence=r['sequence'])
            for key,method in [('HERMES_fixed','hermes_py_000'),('ProteinMPNN','proteinmpnn_v_48_020_full_pep_masked'),('ProteinMPNN002','proteinmpnn_v_48_002_full_pep_masked')]:
                ss,wt,ff=chosen[method];row=next(x for x in ss if x['mutant']==r['mutation']);assert row['sequence']==r['sequence']
                result[key+'_raw']=float(row['pnlogp']);result[key+'_parent_mean']=wt;result[key+'_delta']=float(row['pnlogp'])-wt
                result[key+'_source_sha256']=ff['sha256']
            static.append(result)
    dump(ROOT/'structures/input_manifest.json',dict(templates=temps,candidates=entries,created_at_utc=datetime.now(timezone.utc).isoformat(),
        labels_read_for_sequences_only=True,outcomes_not_used_in_inputs=True,predictions_seen=False))
    dump(ROOT/'structures/input_preparation_audit.json',dict(status='passed',panels=audit,total_variants=1197,total_main_parents=7,
        computational_context_ablation_parents=7,total_forward_entries=len(entries),no_backbone_rebuilding=True,
        static_choice='fixed before inference; HERMES000+ProteinMPNN020 primary,002 sensitivity; all AF3_no_template same input name',
        static_parent_processing='Mean9 stochastic WT scoring rows, not9 biological measurements. This constant subtraction does not affect raw-score ranking.'))
    csvout(ROOT/'structures/static_features.csv',static)
    print(json.dumps(dict(status='inputs_prepared',templates=len(temps),forward_entries=len(entries),variants=len(static))))

def references():
    out=ROOT/'runtime/training_reference';out.mkdir(exist_ok=True)
    manifest=REF/'strict_external_payload/splits/strict_external_train.jsonl'
    assert sha(manifest)=='e0c84c566dea09b436dd7d07e9e6fd8ce1a8cf0c72be20931179e6689c8b44bd'
    rr=[json.loads(line) for line in manifest.read_text().splitlines() if line.strip()];assert len(rr)==2431
    collected=defaultdict(list);bindings=[]
    for row in rr:
        p=REF/'dataset/labels'/Path(row['label_path']).name
        assert p.exists(),str(p)
        digest=sha(p)
        historical_digest=row['metadata'].get('label_sha256')
        if historical_digest is not None:assert digest==historical_digest,str(p)
        z=np.load(p,allow_pickle=False);metric=[str(x) for x in z['metric_names']]
        values=z['values'];mask=z['mask']
        st=parse_pdb(str(REF/'dataset/structures'/Path(row['structure_path']).name))
        names=list(st.sequence)
        if 'residue_names' in z.files:assert ''.join(THREE_TO_ONE[str(x)] for x in z['residue_names'])==st.sequence
        assert values.shape[0]==len(names)==row['metadata']['residue_count']
        target=np.array([not row['target_chains'] or c in row['target_chains'] for c in st.chain_ids])
        for aa in set(names):
            sel=np.array([x==aa for x in names])
            for m in ['rmsf','chi_entropy']:
                j=metric.index(m);valid=sel&mask[:,j]&np.isfinite(values[:,j])&target&st.frame_mask
                if m=='chi_entropy' and aa in {'G','A'}:valid[:]=False
                collected[aa,m].extend(values[valid,j].astype(float).tolist())
        bindings.append(dict(sample_id=row['sample_id'],label_path=str(p),sha256=digest,
            historical_digest_available=historical_digest is not None,structure_sha256=sha(REF/'dataset/structures'/Path(row['structure_path']).name)))
    arrays={};info={}
    for (aa,m),values in sorted(collected.items()):
        if not values:continue
        x=np.sort(np.asarray(values,dtype=np.float64));key=aa+'__'+m;arrays[key]=x
        scale=float(np.sqrt(np.mean(x*x)));assert scale>0
        info[key]=dict(n=len(x),RMS_scale=scale,min=float(x[0]),max=float(x[-1]),median=float(np.median(x)))
    np.savez_compressed(out/'same_AA_training_ECDF.npz',**arrays)
    dump(out/'summary.json',dict(status='passed2431_local_training_record_binding',training_manifest_sha256=sha(manifest),
        training_records=2431,physical_label_scale=True,values_transformed=False,one_record_residue_one_value=True,
        historical_label_digest_matched_records=sum(b['historical_digest_available'] for b in bindings),
        newly_hashed_local_DynaRepo_labels_without_historical_digest=sum(not b['historical_digest_available'] for b in bindings),
        exact_training_target_and_geometry_masks_applied=True,
        same_AA_stats=info,ECDF_sha256=sha(out/'same_AA_training_ECDF.npz'),label_bindings=bindings,
        new_functional_outcomes_used=False,model_predictions_used=False))
    print(json.dumps(dict(status='training_reference_complete',records=2431,groups=len(arrays))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','references']);a=p.parse_args()
    prepare() if a.mode=='prepare' else references()
