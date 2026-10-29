"""All-residue candidate diagnostics; selection is exploratory, never test performance.

Adds frozen OOF baselines, same-coordinate geometry and public-PDB numbering.
No subregion is selected using prediction agreement. Interface ranges are fixed.
"""
import json,sys,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT=Path(r'E:\Project\FlexPPI')
OUT=Path(__file__).resolve().parent/'case_candidates';OUT.mkdir(exist_ok=True)
PARENT=OUT.parent
CORE=ROOT/'experiments/core_corrected_pipeline'
sys.path.insert(0,str(CORE/'training_payload/src'))
from flexppi.constants import THREE_TO_ONE
from flexppi.data.structures import parse_pdb
from flexppi.data.pairing import global_sequence_alignment

CANDIDATES=['2I3T','3QML','3VV2','4ZRJ']

def raw_residues(path):
    rows={}
    for line in path.read_text(encoding='utf-8',errors='replace').splitlines():
        if line[:6].strip()=='ENDMDL':break
        if line[:6].strip() not in ['ATOM','HETATM'] or line[17:20].strip() not in THREE_TO_ONE:continue
        if line[16:17] not in [' ','A','1']:continue
        name=line[12:16].strip();element=line[76:78].strip().upper()
        if element in ['H','D'] or not element and name.lstrip('0123456789').startswith('H'):continue
        chain=line[21:22].strip() or '_';author=line[22:26].strip()+line[26:27].strip()
        key=(chain,author)
        if key not in rows:rows[key]=dict(chain=chain,author=author,aa=THREE_TO_ONE[line[17:20].strip()],atoms={})
        if name not in rows[key]['atoms']:rows[key]['atoms'][name]=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
    return list(rows.values())

def metrics(g,pred):
    y,p=g.md_delta.to_numpy(float),g[pred].to_numpy(float)
    mask=np.isfinite(y)&np.isfinite(p);y,p=y[mask],p[mask]
    strong=np.abs(y)>=.1;up=y>=.1;down=y<=-.1
    return dict(n=len(y),mae=float(np.abs(p-y).mean()) if len(y) else np.nan,
                zero_mae=float(np.abs(y).mean()) if len(y) else np.nan,
                gain_over_zero=float((np.abs(y)-np.abs(p-y)).mean()) if len(y) else np.nan,
                spearman=float(spearmanr(y,p).statistic) if len(y)>3 and np.std(p)>1e-9 and np.std(y)>1e-9 else np.nan,
                mean_truth=float(y.mean()) if len(y) else np.nan,mean_prediction=float(p.mean()) if len(y) else np.nan,
                strong_n=int(strong.sum()),strong_increase_n=int(up.sum()),strong_decrease_n=int(down.sum()),
                direction=float((np.sign(p[strong])==np.sign(y[strong])).mean()) if strong.any() else np.nan,
                always_decrease_direction=float(down[strong].mean()) if strong.any() else np.nan,
                increase_recall=float((p[up]>0).mean()) if up.any() else np.nan,
                decrease_recall=float((p[down]<0).mean()) if down.any() else np.nan)

def main():
    data=pd.read_csv(PARENT/'test53_corrected_full_oof_residue_deltas_with_family.csv.gz')
    inv=pd.read_csv(PARENT/'paired_inventory.csv');inv=inv[inv.pdb_id.isin(CANDIDATES)]
    ids=set(inv.apo_sample)|set(inv.holo_sample)
    cols=['sample_id','residue_index','chain','author_residue_id','amino_acid','md_rmsf','md_hchi',
          'pred_full_rmsf','pred_full_hchi','pred_intra_rmsf','pred_intra_hchi',
          'pred_static_rmsf','pred_static_hchi','pred_aa_only_rmsf','pred_aa_only_hchi',
          'intra_ca_contacts_8A','partner_ca_contacts_8A','interface_ca8',
          'sidechain_sasa_complex_A2','sidechain_sasa_isolated_A2','delta_sidechain_sasa_A2',
          'sidechain_complete','chi_input_complete','nonstandard_residue']
    features=pd.read_parquet(PARENT.parent/'feature_baseline/features_all.parquet',columns=cols,
                             filters=[('sample_id','in',sorted(ids))])
    lookup=features.set_index(['sample_id','residue_index'])
    frames=[];summaries=[];mapping_audit=[]
    for _,r in inv.iterrows():
        pid=r.pdb_id;sample=r.holo_sample
        g=data[data.pdb_id.eq(pid)].copy()
        pdb=OUT/(pid+'_public_structure.pdb')
        if not pdb.exists():
            response=requests.get('https://files.rcsb.org/download/'+pid+'.pdb',timeout=40)
            response.raise_for_status();pdb.write_bytes(response.content)
        public=raw_residues(pdb)
        local_path=CORE/'training_payload/dataset/structures'/(sample+'.pdb')
        local=raw_residues(local_path)
        cache=np.load(CORE/'neural_exports/evaluation/full'/(sample+'.npz'))
        assert [(v['chain'],v['author']) for v in local]==list(zip(cache['chain'],cache['author_residue_id']))
        chain=str(g.chain.iloc[0]);local_target=[v for v in local if v['chain']==chain]
        seq=''.join(v['aa'] for v in local_target)
        candidates=[]
        for pc in sorted(set(v['chain'] for v in public)):
            public_chain=[v for v in public if v['chain']==pc]
            alignment=global_sequence_alignment(seq,''.join(v['aa'] for v in public_chain))
            candidates.append((alignment.coverage,alignment.identity,alignment.score,pc,alignment,public_chain))
        best=max(candidates,key=lambda v:v[:3]);coverage,identity,score,pc,alignment,pc_res=best
        assert identity>.95 and coverage>.95,(pid,coverage,identity)
        map_original={local_target[i]['author']:(pc_res[j]['chain'],pc_res[j]['author']) for i,j in zip(alignment.indices_a,alignment.indices_b)}
        target_keys=set((v['chain'],v['author']) for v in local_target)
        partner=np.concatenate([list(v['atoms'].values()) for v in local if v['chain']!=chain])
        tree=cKDTree(partner)
        heavy_dist={v['author']:float(tree.query(np.array(list(v['atoms'].values())))[0].min()) for v in local_target}
        ca_partner=np.array([v['atoms']['CA'] for v in local if v['chain']!=chain and 'CA'in v['atoms']])
        ca_tree=cKDTree(ca_partner)
        ca_dist={v['author']:float(ca_tree.query(v['atoms']['CA'])[0]) for v in local_target}
        maprow=[]
        for author in sorted(g.author_residue_id.unique(),key=lambda x:int(x)):
            pc2,pn=map_original.get(str(author),('', ''))
            maprow.append(dict(pdb_id=pid,holo_sample=sample,local_chain=chain,local_author=author,
                               public_chain=pc2,public_author=pn,heavy_distance_A=heavy_dist[str(author)],
                               min_CA_distance_A=ca_dist[str(author)]))
        g=g.merge(pd.DataFrame(maprow),left_on=['pdb_id','holo_sample','chain','author_residue_id'],
                  right_on=['pdb_id','holo_sample','local_chain','local_author'],validate='many_to_one')
        for i,row in g.iterrows():
            a=lookup.loc[(row.apo_sample,row.apo_index)];h=lookup.loc[(row.holo_sample,row.holo_index)]
            assert a.amino_acid==h.amino_acid==row.amino_acid
            m='hchi' if row.metric=='chi_entropy' else 'rmsf'
            assert np.isclose(h['md_'+m]-a['md_'+m],row.md_delta,atol=1e-6)
            assert np.isclose(h['pred_full_'+m]-a['pred_full_'+m],row.model_delta,atol=1e-6)
            for model in ['intra','static','aa_only']:
                g.loc[i,model+'_delta']=h['pred_'+model+'_'+m]-a['pred_'+model+'_'+m]
            for name in ['intra_ca_contacts_8A','partner_ca_contacts_8A','sidechain_sasa_complex_A2',
                         'sidechain_sasa_isolated_A2','delta_sidechain_sasa_A2','sidechain_complete','chi_input_complete']:
                g.loc[i,'holo_'+name]=h[name];g.loc[i,'apo_'+name]=a[name]
        g['heavy_distance_region']=np.select([g.heavy_distance_A<=5,g.heavy_distance_A<=12],['interface_le5A','near_5to12A'],default='distal_gt12A')
        g['CA8_interface']=g.min_CA_distance_A<8
        g.to_csv(OUT/(pid+'_all_mapped_residues.csv'),index=False)
        mapping_audit.append(dict(pdb_id=pid,target_local_chain=chain,target_public_chain=pc,
                                  source_sequence_identity=identity,source_sequence_coverage=coverage,
                                  target_residues=len(local_target),mapped_pair_rmsf=int((g.metric=='rmsf').sum()),
                                  mapped_pair_hchi=int((g.metric=='chi_entropy').sum()),
                                  public_file_sha256=hashlib.sha256(pdb.read_bytes()).hexdigest(),
                                  caveat='All source numbering from exact sequence correspondences, notoffset assumption. Homomer public alternatives equivalent; first bestscore chain retained.'))
        for metric,gg in g.groupby('metric'):
            for region,sub in [('all_mapped',gg)]+list(gg.groupby('heavy_distance_region')):
                for model,pred in [('Full','model_delta'),('Intra','intra_delta'),('Static47','static_delta'),('AAonly','aa_only_delta')]:
                    summaries.append(dict(pdb_id=pid,metric=metric,region=region,model=model,
                                          family_clean=bool(r.family_clean),**metrics(sub,pred)))
        frames.append(g)
    combined=pd.concat(frames,ignore_index=True);combined.to_csv(OUT/'candidate_all_residues.csv.gz',index=False,compression='gzip')
    pd.DataFrame(summaries).to_csv(OUT/'candidate_region_baseline_metrics.csv',index=False)
    pd.DataFrame(mapping_audit).to_csv(OUT/'public_PDB_mapping_audit.csv',index=False)
    full=pd.DataFrame(summaries).query("region=='all_mapped'")
    print(full[['pdb_id','metric','model','mae','zero_mae','gain_over_zero','spearman','direction','increase_recall','decrease_recall']].to_string(index=False))
if __name__=='__main__':main()
