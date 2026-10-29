"""Fresh, label-free CaM inference with frozen corrected-pipeline code/weights.

All five folds are retained; no fold selection on NMR results. Only PDB files
are read for inference. The experimental NMR tables are used by a separate
validation program after this program finishes.
"""
import os
os.environ['OMP_NUM_THREADS'] = '2'
os.environ['MKL_NUM_THREADS'] = '2'
os.environ['OPENBLAS_NUM_THREADS'] = '2'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
import sys, json, csv, hashlib, time
from pathlib import Path
ROOT = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
OUT = Path(__file__).resolve().parent
ST = OUT.parent / 'figure_cam_partner_dynamics' / 'structure_data'
sys.path.insert(0, str(ROOT / 'training_payload' / 'src'))
import numpy as np
import torch
from flexppi.config import ModelConfig
from flexppi.constants import METRIC_NAMES
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb, target_mask
from flexppi.model.information_ablation import build_information_ablation
torch.set_num_threads(2)
torch.set_num_interop_threads(1)
torch.manual_seed(7)

CASES = [('CaMKK','1CKK'),('eNOS','1NIW'),('nNOS','2O60'),('smMLCK','1CDL')]
REF = 'ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK'
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value): path.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')
def write_csv(path,rows):
    if not rows: return
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def prepare(path):
    s=parse_pdb(str(path)); n=s.num_residues
    for chain,pos,aa in zip(s.chain_ids,s.author_residue_ids,s.sequence):
        if chain=='A': assert REF[int(pos)-1]==aa, (path,pos,aa)
    example=dict(sample_id=path.stem,dataset='cam_external_structures',source_group='cam_external_structures',metadata={},
        aatype=torch.from_numpy(s.aatype),residue_index=torch.from_numpy(s.residue_index),chain_index=torch.from_numpy(s.chain_index),
        backbone=torch.from_numpy(s.backbone),backbone_mask=torch.from_numpy(s.backbone_mask),geometry_mask=torch.from_numpy(s.frame_mask),
        target_mask=torch.from_numpy(target_mask(s,['A'])),condition=torch.zeros(3),
        label_values=torch.full((n,len(METRIC_NAMES)),float('nan')),label_variances=torch.zeros(n,len(METRIC_NAMES)),
        label_mask=torch.zeros(n,len(METRIC_NAMES),dtype=torch.bool))
    return s,collate_states([example])

def main():
    started=time.time(); (OUT/'inputs').mkdir(exist_ok=True); (OUT/'raw_predictions').mkdir(exist_ok=True)
    plan=dict(cases=CASES,variants=['full','intrachain'],folds=list(range(5)),epoch=40,
        primary_summary='Arithmetic mean of all five folds, every fold retained; no best-fold selection.',
        input='Matched deposited first models, CaM A and peptide P; mature-CaM numbering; no imputation.',
        controls='Retrained intrachain; full same-backbone with peptide deleted (input intervention only).',
        endpoints='All common methyl residues: rank association Hchi vs -S2 and partner-difference direction; NH vs RMSF only exploratory.',
        key_sites=[18,39,105,109,144],core_sites=[39,105,144],
        limits='S2 is not Hchi or RMSF. No NMR fit-error-to-prediction-error conversion. Five folds are not independent experimental repeats.',
        script_sha256=sha(__file__),source_root=str(ROOT/'training_payload'/'src'))
    dump(OUT/'inference_plan.json',plan)
    prepared=[]; inputs=[]
    for partner,pid in CASES:
        path=ST/(pid+'_CaM_A_peptide_P_clean.pdb')
        lines=path.read_text(encoding='utf-8').splitlines()
        deleted=OUT/'inputs'/(pid+'_CaM_A_partner_deleted.pdb')
        # Protein parser ignores Ca ions; retain exactly CaM coordinates.
        deleted.write_text('\n'.join(x for x in lines if x[:6].strip() in ('ATOM','HETATM') and x[21:22]=='A')+'\nEND\n',encoding='utf-8')
        for mode,p in [('complex',path),('partner_deleted',deleted)]:
            s,b=prepare(p); prepared.append((partner,pid,mode,s,b))
            inputs.append(dict(partner=partner,pdb=pid,input_mode=mode,path=str(p),sha256=sha(p),
                total_residues=s.num_residues,cam_residues=sum(c=='A' for c in s.chain_ids),
                cam_valid_frames=int(sum(c=='A' and v for c,v in zip(s.chain_ids,s.frame_mask))),
                cam_MSE_count=sum(c=='A' and r=='MSE' for c,r in zip(s.chain_ids,s.residue_names))))
    dump(OUT/'input_audit.json',inputs)
    print('PREPARED',len(prepared),'inputs; torch',torch.__version__,flush=True)
    rows=[]; checkpoints=[]
    for variant in ['full','intrachain']:
        for fold in range(5):
            path=ROOT/'trained_checkpoints'/variant/('fold_'+str(fold))/'epoch_040.pt'
            ckpt=torch.load(path,map_location='cpu',weights_only=False)
            cfg=ckpt['config']; mc=cfg['model']
            assert ckpt['epoch']==40
            assert mc['information_ablation']==variant
            assert tuple(mc['metrics'])==('rmsf','chi_entropy')
            assert not mc['use_condition_features'] and not mc['enable_response_head']
            model=build_information_ablation(ModelConfig(**mc),variant).cpu().eval()
            model.load_state_dict(ckpt['model'],strict=True); model.requires_grad_(False)
            checkpoints.append(dict(variant=variant,fold=fold,path=str(path),sha256=sha(path),epoch=ckpt['epoch'],config=cfg))
            for partner,pid,mode,s,b in prepared:
                # Deleted-partner intrachain is redundant but provides invariance QA.
                with torch.inference_mode(): result=model.forward_state(b)
                means=result['mean'][0].numpy(); lv=result['log_variance'][0].numpy(); mask=result['mask'][0].numpy()
                assert np.isfinite(means).all() and np.array_equal(mask,s.frame_mask)
                np.savez_compressed(OUT/'raw_predictions'/f'{variant}_fold{fold}_{pid}_{mode}.npz',
                    mean=means,log_variance=lv,mask=mask,chain=np.array(s.chain_ids),author_residue_id=np.array(s.author_residue_ids),
                    amino_acid=np.array(list(s.sequence)),metric_names=np.array(mc['metrics']))
                for i,(chain,pos,aa) in enumerate(zip(s.chain_ids,s.author_residue_ids,s.sequence)):
                    if chain=='A' and mask[i]:
                        rows.append(dict(variant=variant,fold=fold,partner=partner,pdb_id=pid,input_mode=mode,
                            seq_id=int(pos),amino_acid=aa,pred_rmsf=float(means[i,0]),pred_hchi=float(means[i,1]),
                            rmsf_head_variance=float(np.exp(lv[i,0])),hchi_head_variance=float(np.exp(lv[i,1]))))
                print('DONE',variant,fold,partner,mode,flush=True)
            del model,ckpt
            write_csv(OUT/'predictions_by_fold.csv',rows)
    dump(OUT/'checkpoint_manifest.json',checkpoints)
    grouped={}
    for row in rows:
        key=tuple(row[k] for k in ['variant','partner','pdb_id','input_mode','seq_id','amino_acid'])
        grouped.setdefault(key,[]).append(row)
    means=[]
    for key,group in grouped.items():
        assert sorted(r['fold'] for r in group)==list(range(5))
        row=dict(zip(['variant','partner','pdb_id','input_mode','seq_id','amino_acid'],key))
        for metric in ['rmsf','hchi']:
            v=np.array([r['pred_'+metric] for r in group]); row['pred_'+metric]=float(v.mean()); row[metric+'_between_fold_sd']=float(v.std(ddof=0))
        means.append(row)
    write_csv(OUT/'predictions_fivefold_mean.csv',means)
    dump(OUT/'completion.json',dict(status='complete',elapsed_seconds=time.time()-started,prediction_rows=len(rows),
        mean_rows=len(means),checkpoint_count=len(checkpoints),source_files={str(p.relative_to(ROOT/'training_payload'/'src')):sha(p) for p in (ROOT/'training_payload'/'src'/'flexppi').rglob('*.py')}))
    print('COMPLETE',len(rows),'rows',time.time()-started,'seconds',flush=True)
if __name__=='__main__': main()
