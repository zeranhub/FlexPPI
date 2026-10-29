"""Primary CaM inference: family-held-out fold1, full and retrained intrachain."""
import run_current_cam as base
import time
np,torch=base.np,base.torch
OUT=base.OUT

def main():
    start=time.time(); rows=[]; manifests=[]
    model_root=OUT/'family_heldout_model'
    expected={'full':'c1d7b1856190895ff5c5ce8442810c203ad0deb90b80e8fec029f82b616fdcf9',
              'intrachain':'8bd8657ef3a8f85a60e9e5802ccda90ac62ccf386e8d0b0a4d02ece2143aa68c'}
    (OUT/'family_raw_predictions').mkdir(exist_ok=True)
    base.dump(OUT/'family_inference_plan.json',dict(fold=1,variants=['full','intrachain'],
        family='family_1ad15d0d4afa1419',family_records=['A001K','A004C','A007U'],
        choice='Prespecified from family assignment before examining predictions; no NMR-based model selection.',
        script_sha256=base.sha(__file__),source_root=str(base.ROOT/'training_payload'/'src'),
        primary_pair='nNOS -> smMLCK',key_sites=[18,39,105,109,144],core_sites=[39,105,144]))
    for variant in ['full','intrachain']:
        path=model_root/'runs'/variant/'fold_1'/'epoch_040.pt'
        assert base.sha(path)==expected[variant]
        ckpt=torch.load(path,map_location='cpu',weights_only=False)
        cfg=ckpt['config']; mc=cfg['model']
        assert mc['information_ablation']==variant and tuple(mc['metrics'])==('rmsf','chi_entropy')
        assert not mc['use_condition_features'] and not mc['enable_response_head']
        assert 'fold_1' in cfg['output_dir']
        model=base.build_information_ablation(base.ModelConfig(**mc),variant).cpu().eval()
        model.load_state_dict(ckpt['model'],strict=True); model.requires_grad_(False)
        manifests.append(dict(variant=variant,fold=1,sha256=base.sha(path),path=str(path),config=cfg))
        for partner,pid in base.CASES:
            for mode,path_in in [('complex',base.ST/(pid+'_CaM_A_peptide_P_clean.pdb')),
                                 ('partner_deleted',OUT/'inputs'/(pid+'_CaM_A_partner_deleted.pdb'))]:
                s,b=base.prepare(path_in)
                with torch.inference_mode(): result=model.forward_state(b)
                means=result['mean'][0].numpy(); lv=result['log_variance'][0].numpy(); mask=result['mask'][0].numpy()
                assert np.isfinite(means).all() and np.array_equal(mask,s.frame_mask)
                np.savez_compressed(OUT/'family_raw_predictions'/f'{variant}_fold1_{pid}_{mode}.npz',
                    mean=means,log_variance=lv,mask=mask,chain=np.array(s.chain_ids),
                    author_residue_id=np.array(s.author_residue_ids),amino_acid=np.array(list(s.sequence)),metric_names=np.array(mc['metrics']))
                for i,(chain,pos,aa) in enumerate(zip(s.chain_ids,s.author_residue_ids,s.sequence)):
                    if chain=='A' and mask[i]: rows.append(dict(variant=variant,fold=1,partner=partner,pdb_id=pid,input_mode=mode,
                        seq_id=int(pos),amino_acid=aa,pred_rmsf=float(means[i,0]),pred_hchi=float(means[i,1]),
                        rmsf_head_variance=float(np.exp(lv[i,0])),hchi_head_variance=float(np.exp(lv[i,1]))))
                print('DONE family',variant,partner,mode,flush=True)
        del model,ckpt
    base.write_csv(OUT/'family_predictions.csv',rows)
    base.dump(OUT/'family_checkpoint_manifest.json',manifests)
    base.dump(OUT/'family_completion.json',dict(status='complete',rows=len(rows),checkpoints=2,elapsed_seconds=time.time()-start))
    print('COMPLETE FAMILY',len(rows),time.time()-start,flush=True)
if __name__=='__main__': main()
