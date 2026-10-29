"""Exploratory screen of existing, correctly held-out predictions. No fitting."""
from pathlib import Path
from collections import defaultdict
import csv, json, hashlib
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
BASE = Path(r'E:/Project/FlexPPI/experiments/_sources')
FAMILY = BASE / 'family_folds'
COHORT = BASE / 'local_model_benchmark_20261023/cohort'
IDENTITY = BASE / 'local_model_benchmark_20261023/grouping_protein_complex_20261024'
AA = dict(zip(['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'], 'ARNDCQEGHILKMFPSTWYV'))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()

def rho(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float)
    good=np.isfinite(a)&np.isfinite(b)
    if good.sum()<4 or np.ptp(a[good])==0 or np.ptp(b[good])==0: return np.nan
    return float(pd.Series(a[good]).rank().corr(pd.Series(b[good]).rank()))

def centered(a,aa):
    a=np.asarray(a,float).copy()
    valid=np.zeros(len(a),bool)
    for typ in np.unique(aa):
        ix=(aa==typ)&np.isfinite(a)
        if ix.sum()>=4:
            a[ix]-=a[ix].mean(); valid|=ix
    a[~valid]=np.nan
    return a

def parse_structure(path):
    residues={}
    with Path(path).open() as f:
        for line in f:
            if line[:6] not in ['ATOM  ','HETATM'] or line[16:17] not in [' ','A']: continue
            name=line[17:20].strip()
            if name not in AA: continue
            key=(line[21:22].strip() or '_',line[22:26].strip()+line[26:27].strip())
            r=residues.setdefault(key,dict(chain=key[0],number=key[1],aa=AA[name],atoms=[]))
            atom=line[12:16].strip()
            xyz=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
            if atom=='CA': r['ca']=xyz
            if not atom.startswith('H'): r['atoms'].append(xyz)
    return [r for r in residues.values() if 'ca' in r]

def metrics(d,mask):
    result={}
    for metric in ['rmsf','chi_entropy']:
        good=mask&d['valid_'+metric].to_numpy(bool)
        y=d['target_'+metric].to_numpy(float)[good]
        f=d['pred_'+metric].to_numpy(float)[good]
        i=d['intra_'+metric].to_numpy(float)[good]
        a=d['amino_acid'].to_numpy(str)[good]
        q=dict(n=int(good.sum()))
        if len(y):
            q.update(full_mae=float(np.mean(abs(f-y))),intra_mae=float(np.mean(abs(i-y))),
                     full_rho=rho(y,f),intra_rho=rho(y,i),
                     full_type_centered_rho=rho(centered(y,a),centered(f,a)),
                     intra_type_centered_rho=rho(centered(y,a),centered(i,a)))
            q['mae_gain_fraction']=1-q['full_mae']/q['intra_mae'] if q['intra_mae']>0 else np.nan
            q['fraction_sites_full_lower_error']=float((abs(f-y)<abs(i-y)).mean())
        result.update({metric+'_'+k:v for k,v in q.items()})
    return result

def main():
    members=pd.read_csv(IDENTITY/'membership_frozen.csv').fillna('')
    folds=pd.read_csv(FAMILY/'record_family_folds.csv').set_index('sample_id')
    records=pd.read_csv(COHORT/'records.csv').set_index('sample_id')
    chains=pd.read_csv(IDENTITY/'chain_biological_identity_evidence.csv').fillna('')
    wanted=set(members.loc[(members.benchmark_group=='protein_complex')&~members.sample_id.str.contains('^\\d[a-z0-9]{3}$',regex=True),'sample_id'])
    met=pd.read_csv(FAMILY/'evaluation/per_record_metrics.csv')
    wanted&=set(met.loc[met.dataset!='itrace','sample_id'])
    assert len(wanted)==616
    paths=list((OUT/'archive_extract').rglob('oof_per_residue.csv.gz'))
    if len(paths)!=1: raise RuntimeError('Expected exactly one extracted family OOF table: '+str(paths))
    path=paths[0]
    pieces=[]
    for chunk in pd.read_csv(path,chunksize=200000):
        pieces.append(chunk.loc[chunk.sample_id.isin(wanted)].copy())
    data=pd.concat(pieces,ignore_index=True)
    intra=data.loc[data.variant=='intrachain',['sample_id','residue_index','pred_rmsf','pred_chi_entropy','target_rmsf','target_chi_entropy']].rename(columns={'pred_rmsf':'intra_rmsf','pred_chi_entropy':'intra_chi_entropy','target_rmsf':'intra_target_rmsf','target_chi_entropy':'intra_target_chi_entropy'})
    data=data.loc[data.variant=='full'].merge(intra,on=['sample_id','residue_index'],validate='one_to_one')
    for metric in ['rmsf','chi_entropy']:
        assert np.allclose(data['target_'+metric],data['intra_target_'+metric],equal_nan=True)
    summaries=[]; chain_rows=[]; spatial_rows=[]; audits=[]; all_profiles=[]
    member_map=members.set_index('sample_id')
    for n,(sid,d) in enumerate(data.groupby('sample_id',sort=True)):
        d=d.sort_values('residue_index').copy().reset_index(drop=True)
        assert set(d.fold)=={int(folds.loc[sid,'fold'])}
        actual=Path(records.loc[sid,'structure_path'])
        if not actual.exists(): actual=Path(r'E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/dataset/structures')/(sid+'.pdb')
        residues=parse_structure(actual)
        if len(residues)!=len(d) or ''.join(r['aa'] for r in residues)!=''.join(d.amino_acid):
            audits.append(dict(sample_id=sid,status='structure_mapping_failed',rows=len(d),ca_residues=len(residues))); continue
        identity=chains.loc[chains.sample_id==sid].drop_duplicates('input_chain').set_index('input_chain')
        d['partner_identity']=[identity.loc[r['chain'],'biological_identity_evidence'] if r['chain'] in identity.index else r['chain'] for r in residues]
        d['chain_description']=[identity.loc[r['chain'],'description'] if r['chain'] in identity.index else '' for r in residues]
        xyz=np.array([r['ca'] for r in residues])
        dist=np.sqrt(np.sum((xyz[:,None,:]-xyz[None,:,:])**2,axis=-1))
        np.fill_diagonal(dist,np.inf)
        identities=d.partner_identity.to_numpy(str)
        distinct=identities[:,None]!=identities[None,:]
        d['partner_ca_distance']=np.where(distinct,dist,np.inf).min(axis=1)
        d['interface_ca8']=d.partner_ca_distance<=8
        d['all_ca_neighbors8']=(dist<=8).sum(axis=1)
        d['partner_ca_neighbors8']=((dist<=8)&distinct).sum(axis=1)
        d['pdb_id']=member_map.loc[sid,'mapping_pdb_id'] or member_map.loc[sid,'pdb_id']
        d['structure_path']=str(actual)
        row=dict(sample_id=sid,pdb_id=d.pdb_id.iloc[0],fold=int(d.fold.iloc[0]),n_residues=len(d),
                 n_chains=d.chain.nunique(),n_protein_identities=d.partner_identity.nunique(),
                 min_chain_length=int(d.groupby('chain').size().min()),
                 description=' / '.join(dict.fromkeys(d.chain_description)),structure_path=str(actual),
                 **metrics(d,np.ones(len(d),bool)))
        for m in ['rmsf','chi_entropy']:
            good=d['valid_'+m].to_numpy(bool)
            y=d['target_'+m].to_numpy(float)[good]
            density=-np.log1p(d.all_ca_neighbors8.to_numpy(float)[good])
            row[m+'_density_rho']=rho(y,density)
            row[m+'_density_type_centered_rho']=rho(centered(y,d.amino_acid.to_numpy(str)[good]),centered(density,d.amino_acid.to_numpy(str)[good]))
        summaries.append(row)
        for chain,cd in d.groupby('chain',sort=False):
            chain_rows.append(dict(sample_id=sid,pdb_id=row['pdb_id'],fold=row['fold'],chain=chain,
                                   description=cd.chain_description.iloc[0],identity=cd.partner_identity.iloc[0],
                                   n_residues=len(cd),**metrics(cd,np.ones(len(cd),bool))))
        for region,mask in [('interface_ca8',d.interface_ca8.to_numpy(bool)),('outside_interface',~d.interface_ca8.to_numpy(bool)),('remote_ca12',(d.partner_ca_distance>12).to_numpy(bool))]:
            spatial_rows.append(dict(sample_id=sid,pdb_id=row['pdb_id'],region=region,**metrics(d,mask)))
        all_profiles.append(d)
        audits.append(dict(sample_id=sid,status='mapped',structure_sha256=sha(actual),n_residues=len(d)))
        if (n+1)%100==0: print('mapped',n+1,flush=True)
    s=pd.DataFrame(summaries)
    s['balanced_mae_gain']=np.minimum(s.rmsf_mae_gain_fraction,s.chi_entropy_mae_gain_fraction)
    s['exploratory_both_gain_and_rank']=(s.balanced_mae_gain>0)&(s.rmsf_full_rho>=.65)&(s.chi_entropy_full_rho>=.85)
    s=s.sort_values('balanced_mae_gain',ascending=False)
    s.to_csv(OUT/'all_616_record_screen.csv',index=False)
    pd.DataFrame(chain_rows).to_csv(OUT/'all_chain_metrics.csv',index=False)
    pd.DataFrame(spatial_rows).to_csv(OUT/'all_region_metrics.csv',index=False)
    pd.concat(all_profiles,ignore_index=True).to_csv(OUT/'all_non_tcr_profiles.csv.gz',index=False,compression='gzip')
    (OUT/'screen_audit.json').write_text(json.dumps(dict(status='complete',records_requested=616,records_mapped=len(s),
        exploratory_selection=True,models_retrained=False,new_inference=False,
        caveat='An exploratory case screen, not an unbiased aggregate performance evaluation. Family OOF Full/Intra compare separately trained architectures; this is not direct context deletion from one model.',
        input_files={str(p):sha(p) for p in [path,IDENTITY/'membership_frozen.csv',FAMILY/'evaluation/per_record_metrics.csv',FAMILY/'record_family_folds.csv']},structure_audits=audits),indent=2),encoding='utf-8')
    cols=['sample_id','pdb_id','fold','min_chain_length','rmsf_full_rho','chi_entropy_full_type_centered_rho','rmsf_mae_gain_fraction','chi_entropy_mae_gain_fraction','description']
    q=s.loc[s.exploratory_both_gain_and_rank&(s.min_chain_length>=35)].head(25)
    print(q[cols].to_string(index=False))

if __name__=='__main__':main()
