"""Evaluate fixed fold-0 inference against original MD labels, not NMR."""
from pathlib import Path
import sys,json,itertools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from screen_exact_sequence_contexts import score
OUT=Path(__file__).resolve().parent;R=OUT.parent.parent
CORE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
sys.path.insert(0,str(CORE/'training_payload/src'))
from flexppi.data.pairing import global_sequence_alignment

def main():
    data=pd.read_csv(R/'paired_and_experimental/test53_corrected_full_oof_residue_deltas_with_family.csv.gz')
    inv=pd.read_csv(R/'paired_and_experimental/paired_inventory.csv')
    selected=inv[inv.pdb_id.isin(['2Y5B','4ZRJ','2XWT','2I3T','3QML'])]
    ids=set(selected.apo_sample)|set(selected.holo_sample)|{'A001K','A004C','A007U'}
    p=pd.read_parquet(R/'feature_baseline/oof_predictions_all.parquet',filters=[('sample_id','in',sorted(ids))])
    feat=pd.read_parquet(R/'feature_baseline/features_all.parquet',columns=['sample_id','residue_index','chain','author_residue_id','amino_acid','interface_ca8','partner_ca_contacts_8A','sidechain_sasa_complex_A2','sidechain_sasa_isolated_A2'],filters=[('sample_id','in',sorted(ids))])
    feat=feat.merge(p.drop(columns='interface_ca8'),on=['sample_id','residue_index'],validate='one_to_one').set_index(['sample_id','residue_index'])
    pred={};cache={}
    for sid in ids:
        cache[sid]=np.load(CORE/'neural_exports/evaluation/full'/(sid+'.npz'))
        for model in ['full','intrachain']:
            z=np.load(OUT/'training_member_predictions'/(model+'_'+sid+'.npz'))
            for name in ['amino_acid','author_residue_id','chain']:assert np.array_equal(z[name],cache[sid][name]),(sid,name)
            pred[sid,model]=z['mean']
    rows=[]
    for _,pair in selected.iterrows():
        g=data[data.pair_id.eq(pair.pair_id)]
        for _,row in g.iterrows():
            m='rmsf' if row.metric=='rmsf' else 'hchi';j=int(m=='hchi');x=int(row.apo_index);z=int(row.holo_index)
            a=feat.loc[(row.apo_sample,x)];b=feat.loc[(row.holo_sample,z)]
            new=dict(pair_id=pair.pair_id,pdb_a=pair.pdb_id,pdb_b=pair.pdb_id,sample_a=row.apo_sample,sample_b=row.holo_sample,
                     index_a=x,index_b=z,chain_a=a.chain,chain_b=b.chain,author_a=a.author_residue_id,author_b=b.author_residue_id,
                     amino_acid=row.amino_acid,metric=m,md_a=row.md_apo,md_b=row.md_holo,md_delta=row.md_delta,
                     full_oof_delta=row.model_delta,full0_delta=float(pred[row.holo_sample,'full'][z,j]-pred[row.apo_sample,'full'][x,j]),
                     intra0_delta=float(pred[row.holo_sample,'intrachain'][z,j]-pred[row.apo_sample,'intrachain'][x,j]),
                     interface=bool(a.interface_ca8 or b.interface_ca8),contact_a=a.partner_ca_contacts_8A,contact_b=b.partner_ca_contacts_8A,
                     sasa_a=a.sidechain_sasa_complex_A2,sasa_b=b.sidechain_sasa_complex_A2)
            for model in ['intra','aa_sasa_contacts','static47_sasa']:
                new[model+'_oof_delta']=b['pred_'+model+'_'+m]-a['pred_'+model+'_'+m]
            rows.append(new)
    # Three native CaM MD states; pairwise exact common AA after global sequence alignment.
    cams=['A001K','A004C','A007U'];mapping=[]
    for s1,s2 in itertools.combinations(cams,2):
        a,b=cache[s1],cache[s2]
        ca=max(set(a['chain']),key=lambda c:sum(a['chain']==c));cb=max(set(b['chain']),key=lambda c:sum(b['chain']==c))
        ia=np.flatnonzero(a['chain']==ca);ib=np.flatnonzero(b['chain']==cb)
        sa=''.join(a['amino_acid'][ia]);sb=''.join(b['amino_acid'][ib])
        assert 'LTEEQIA' in sa and 'LTEEQIA' in sb
        al=global_sequence_alignment(sa,sb)
        for x0,z0 in zip(al.indices_a,al.indices_b):
            x,z=int(ia[x0]),int(ib[z0]);aa=str(a['amino_acid'][x])
            if aa!=b['amino_acid'][z]:continue
            av,bv=feat.loc[(s1,x)],feat.loc[(s2,z)]
            mapping.append(dict(sample_a=s1,sample_b=s2,index_a=x,index_b=z,author_a=av.author_residue_id,author_b=bv.author_residue_id,aa=aa,alignment_identity=al.identity,alignment_coverage=al.coverage))
            for j,m in enumerate(['rmsf','hchi']):
                if not(a['valid'][x,j] and b['valid'][z,j]):continue
                if m=='hchi' and aa in 'AGP':continue
                row=dict(pair_id=s1+'__'+s2,pdb_a='CaM_MD',pdb_b='CaM_MD',sample_a=s1,sample_b=s2,index_a=x,index_b=z,
                         chain_a=ca,chain_b=cb,author_a=av.author_residue_id,author_b=bv.author_residue_id,amino_acid=aa,metric=m,
                         md_a=float(a['target'][x,j]),md_b=float(b['target'][z,j]),md_delta=float(b['target'][z,j]-a['target'][x,j]),
                         full_oof_delta=float(b['mean'][z,j]-a['mean'][x,j]),
                         full0_delta=float(pred[s2,'full'][z,j]-pred[s1,'full'][x,j]),
                         intra0_delta=float(pred[s2,'intrachain'][z,j]-pred[s1,'intrachain'][x,j]),
                         interface=bool(av.interface_ca8 or bv.interface_ca8),contact_a=av.partner_ca_contacts_8A,contact_b=bv.partner_ca_contacts_8A,
                         sasa_a=av.sidechain_sasa_complex_A2,sasa_b=bv.sidechain_sasa_complex_A2)
                for model in ['intra','aa_sasa_contacts','static47_sasa']:
                    row[model+'_oof_delta']=bv['pred_'+model+'_'+m]-av['pred_'+model+'_'+m]
                rows.append(row)
    pd.DataFrame(mapping).to_csv(OUT/'cam_md_sequence_mapping.csv',index=False)
    df=pd.DataFrame(rows);df.to_csv(OUT/'fixed_fold0_candidate_all_residue_deltas.csv.gz',index=False,compression='gzip')
    sums=[]
    for (pid,m),g in df.groupby(['pair_id','metric']):
        for region,q in [('all',g),('interface_union',g[g.interface]),('noninterface',g[~g.interface])]:
            if len(q)<5:continue
            for model in ['full0','intra0','full_oof','intra_oof','aa_sasa_contacts_oof','static47_sasa_oof']:
                if pid.startswith('A00') and model.endswith('oof'):continue # CaM states have different heldout folds.
                sums.append(dict(pair_id=pid,metric=m,region=region,model=model,**score(q.md_delta.to_numpy(),q[model+'_delta'].to_numpy())))
    sums=pd.DataFrame(sums);sums['gain_vs_zero']=sums.zero_mae-sums.mae
    sums.to_csv(OUT/'fixed_fold0_candidate_metrics.csv',index=False)
    print(sums[sums.region.eq('all')].to_string(index=False))
if __name__=='__main__':main()
