"""Reconstruct matched CaM state contrast from NPZ independently of peer CSV.

The state contrast includes calcium occupancy and binding partner differences.
Frozen fold0 predictions are training-member predictions, as authorized by user.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OUT=Path(__file__).resolve().parent
PEER=OUT.parent/'partner_screen'
CORE=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
REF='ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK'

def met(y,p):
    return {'n':len(y),'mae':float(abs(y-p).mean()),'zero_mae':float(abs(y).mean()),'rho':float(spearmanr(y,p).statistic),'direction_all':float((np.sign(y)==np.sign(p)).mean()),'mean_md':float(y.mean()),'mean_pred':float(p.mean())}

def main():
    data={};mappings=[]
    for sid,chain in [('A004C','B'),('A007U','A')]:
        truth=np.load(CORE/f'neural_exports/evaluation/full/{sid}.npz')
        full=np.load(PEER/f'training_member_predictions/full_{sid}.npz')
        intra=np.load(PEER/f'training_member_predictions/intrachain_{sid}.npz')
        for field in ['chain','author_residue_id','amino_acid']:
            assert np.array_equal(truth[field],full[field]),(sid,field)
            assert np.array_equal(truth[field],intra[field]),(sid,field)
        ix=np.flatnonzero(truth['chain']==chain)
        seq=''.join(truth['amino_acid'][ix])
        assert seq.startswith(REF[3:]) and len(ix) in [145,146],(sid,seq)
        if len(ix)==146:assert seq==REF[3:]+'A'
        data[sid]=dict(truth=truth,full=full,intra=intra,ix=ix[:145])
        mappings.extend(dict(sample_id=sid,chain=chain,residue_index=int(i),input_author=str(truth['author_residue_id'][i]),mature_number=j+4,amino_acid=REF[j+3]) for j,i in enumerate(ix[:145]))
    rows=[]
    for j,pos in enumerate(range(4,149)):
        a,b=[data[s] for s in ['A004C','A007U']]
        ia,ib=a['ix'][j],b['ix'][j]
        for k,metric in [(0,'rmsf'),(1,'hchi')]:
            if not(a['truth']['valid'][ia,k] and b['truth']['valid'][ib,k]):continue
            if metric=='hchi' and REF[pos-1]=='P':continue
            row=dict(mature_cam_number=pos,amino_acid=REF[pos-1],metric=metric,index_a=int(ia),index_b=int(ib),region='N_lobe' if pos<=75 else 'linker' if pos<=81 else 'C_lobe')
            for state,sid,idx in [('a','A004C',ia),('b','A007U',ib)]:
                row['md_'+state]=float(data[sid]['truth']['target'][idx,k])
                for model in ['full','intra']:
                    row[model+'0_'+state]=float(data[sid][model]['mean'][idx,k])
            for var in ['md','full0','intra0']:
                row[var+'_delta']=row[var+'_b']-row[var+'_a']
            rows.append(row)
    tab=pd.DataFrame(rows)
    peer=pd.read_csv(PEER/'cam_md_state_contrast_all_common_residues.csv')
    audit=tab.merge(peer,on=['metric','mature_cam_number'],validate='one_to_one',suffixes=('_audit','_peer'))
    errors={col:float(abs(audit[col+'_audit']-audit[col+'_peer']).max()) for col in ['md_a','md_b','md_delta','full0_a','full0_b','full0_delta','intra0_a','intra0_b','intra0_delta']}
    assert len(audit)==len(tab)==len(peer)==267
    assert max(errors.values())<1e-6,errors
    pd.DataFrame(mappings).to_csv(OUT/'cam_state_independent_sequence_mapping.csv',index=False)
    tab.to_csv(OUT/'cam_state_independent_all_common_residues.csv',index=False)
    metrics=[];means=[]
    for metric,g in tab.groupby('metric'):
        for region,q in [('all',g)]+list(g.groupby('region'))+[('all_nine_methionines',g[g.amino_acid=='M'])]:
            for model in ['full0','intra0']:
                metrics.append(dict(metric=metric,region=region,model=model,**met(q.md_delta.to_numpy(),q[model+'_delta'].to_numpy())))
            means.append(dict(metric=metric,region=region,n=len(q),md_a=q.md_a.mean(),md_b=q.md_b.mean(),md_delta=q.md_delta.mean(),up=int((q.md_delta>0).sum()),down=int((q.md_delta<0).sum()),strong_up=int((q.md_delta>=.1).sum()),strong_down=int((q.md_delta<=-.1).sum())))
    pd.DataFrame(metrics).to_csv(OUT/'cam_state_independent_region_metrics.csv',index=False)
    pd.DataFrame(means).to_csv(OUT/'cam_state_independent_region_means.csv',index=False)
    meth=tab[tab.amino_acid=='M'].merge(peer[['metric','mature_cam_number','contact_a','contact_b','sasa_a','sasa_b']],on=['metric','mature_cam_number'],validate='one_to_one')
    meth.to_csv(OUT/'cam_state_all_nine_methionines.csv',index=False)
    # Quantify the actual joint pattern, without calling the whole backbone fixed.
    joint=tab[tab.metric=='rmsf'].merge(tab[tab.metric=='hchi'],on=['mature_cam_number','amino_acid','region'],suffixes=('_rmsf','_hchi'))
    jt=[]
    for region,q in [('all',joint)]+list(joint.groupby('region'))+[('all_nine_methionines',joint[joint.amino_acid=='M'])]:
        jt.append(dict(region=region,n=len(q),backbone_decrease_0p1=int((q.md_delta_rmsf<=-.1).sum()),backbone_decrease_with_hchi_increase_0p1=int(((q.md_delta_rmsf<=-.1)&(q.md_delta_hchi>=.1)).sum()),backbone_decrease_with_hchi_decrease_0p1=int(((q.md_delta_rmsf<=-.1)&(q.md_delta_hchi<=-.1)).sum()),near_unchanged_backbone_0p1=int((q.md_delta_rmsf.abs()<=.1).sum())))
    pd.DataFrame(jt).to_csv(OUT/'cam_state_joint_backbone_sidechain_patterns.csv',index=False)
    (OUT/'cam_state_independent_audit.json').write_text(json.dumps(dict(status='passed',n=267,absolute_field_error_max=errors,scope='Fixed corrected fold0 training-member predictions for both states. No inference of partner-only effect. Backbone main-chain RMSF decreases substantially; no fixed-backbone claim. Linker RMSF response prediction is wrong sign.'),indent=2))
    print(pd.DataFrame(metrics).to_string(index=False))
    print(pd.DataFrame(means).to_string(index=False))
    print(pd.DataFrame(jt).to_string(index=False))

if __name__=='__main__':main()
