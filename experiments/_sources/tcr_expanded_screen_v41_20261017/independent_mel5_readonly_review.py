from pathlib import Path
import json,sys,hashlib
ROOT=Path(__file__).resolve().parent
MEL=ROOT/'mel5_validation'
sys.path.insert(0,str(ROOT.parent/'label_sensitivity/runtime'))
import numpy as np,pandas as pd
from Bio.Align import PairwiseAligner
from scipy.stats import spearmanr

SOURCES=ROOT.parent/'tcr_case_v40_20261017/family_predictions/full'
roles={'A':'mhc','B':'b2m','C':'peptide','D':'tcr_alpha','E':'tcr_beta'}
states={};checks=[]
for pdb in ['3hg1','4jff']:
    with np.load(SOURCES/f'itrace_{pdb}.npz') as z:p={k:z[k].copy() for k in z.files}
    with np.load(MEL/f'{pdb}_all_residue_replica_labels.npz') as z:v={k:z[k].copy() for k in z.files}
    for key in ['chain','author_residue_id','amino_acid']:np.testing.assert_array_equal(v[key],p[key])
    finite=np.isfinite(v['values']).any(axis=0)
    np.testing.assert_allclose(v['mean'][finite],np.nanmean(v['values'][:,finite],axis=0),atol=3e-6,rtol=0,equal_nan=True)
    for m in [0,1]:
        valid=p['valid'][:,m]
        np.testing.assert_allclose(v['mean'][valid,m],p['target'][valid,m],atol=3e-6,rtol=0)
        np.testing.assert_allclose(v['sample_variance'][valid,m],p['label_variance'][valid,m],atol=3e-6,rtol=0)
    cas={}
    for line in (MEL/f'{pdb}_model_input.pdb').read_text().splitlines():
        if line.startswith('ATOM  ') and line[12:16].strip()=='CA':
            key=(line[21],line[22:27].strip());assert key not in cas
            cas[key]=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
    xyz=np.array([cas[(str(c),str(r))] for c,r in zip(p['chain'],p['author_residue_id'])])
    pmhc=np.flatnonzero(np.isin(p['chain'],['A','B','C']));tcr=np.flatnonzero(np.isin(p['chain'],['D','E']))
    near=np.sum((xyz[pmhc,None]-xyz[None,tcr])**2,axis=-1)<=64
    interface=np.zeros(len(xyz),bool);interface[pmhc]=near.any(1);interface[tcr]=near.any(0)
    states[pdb]={'p':p,'v':v,'interface':interface}

al=PairwiseAligner(mode='global',match_score=2,mismatch_score=-2,open_gap_score=-6,extend_gap_score=-.5)
rows=[];alignments=[]
CDR={'D':[(26,31,'CDR1alpha'),(49,54,'CDR2alpha'),(89,97,'CDR3alpha')], 'E':[(27,32,'CDR1beta'),(50,54,'CDR2beta'),(92,104,'CDR3beta')]}
sa=states['3hg1'];sb=states['4jff']
for chain in ['A','B','C','D','E']:
    ia=np.flatnonzero(sa['p']['chain']==chain);ib=np.flatnonzero(sb['p']['chain']==chain)
    seq_a=''.join(sa['p']['amino_acid'][ia]);seq_b=''.join(sb['p']['amino_acid'][ib])
    result=al.align(seq_a,seq_b)[0]
    pairs=[]
    for ba,bb in zip(*result.aligned):pairs.extend(zip(range(*ba),range(*bb)))
    alignments.append({'chain':chain,'aligned':len(pairs),'mismatches':sum(seq_a[i]!=seq_b[j] for i,j in pairs),'blocks':result.aligned.tolist()})
    for i,j in pairs:
        a=ia[i];b=ib[j];ida=str(sa['p']['author_residue_id'][a]);idb=str(sb['p']['author_residue_id'][b])
        cdr=next((name for lo,hi,name in CDR.get(chain,[]) if lo<=int(ida)<=hi),'outside_CDR')
        row={'chain':chain,'id_a':ida,'id_b':idb,'role':roles[chain],'cdr':cdr,'aa_a':str(sa['p']['amino_acid'][a]),'aa_b':str(sb['p']['amino_acid'][b]),
             'same_aa':sa['p']['amino_acid'][a]==sb['p']['amino_acid'][b],'interface':bool(sa['interface'][a] or sb['interface'][b])}
        for m,label in enumerate(['rmsf','hchi']):
            valid=sa['p']['valid'][a,m] and sb['p']['valid'][b,m]
            row[label+'_valid']=bool(valid)
            row[label+'_md_delta']=float(sb['v']['mean'][b,m]-sa['v']['mean'][a,m])
            # Convert each float32 model value to float64 before differencing, as the source CSV does.
            row[label+'_model_delta']=float(sb['p']['mean'][b,m])-float(sa['p']['mean'][a,m])
            row[label+'_sem']=float(np.sqrt((sa['v']['sample_variance'][a,m]+sb['v']['sample_variance'][b,m])/3))
            va=sa['v']['values'][:,a,m];vb=sb['v']['values'][:,b,m]
            row[label+'_nonoverlap']=bool(va.max()<vb.min() or vb.max()<va.min())
        rows.append(row)
df=pd.DataFrame(rows)
source=pd.read_csv(MEL/'corrected_pair_all_aligned_positions.csv',dtype={'id_a':str,'id_b':str})
assert len(source)==len(df)==824
audit_source=df.merge(source,on=['role','id_a','id_b'],validate='one_to_one',suffixes=('_review','_source'))
np.testing.assert_array_equal(audit_source.interface_review,audit_source.interface_source)
for m in ['rmsf','hchi']:
    np.testing.assert_allclose(audit_source[m+'_md_delta_review'],audit_source[m+'_md_delta_source'],rtol=0,atol=1e-12,equal_nan=True)
    np.testing.assert_allclose(audit_source[m+'_model_delta'],audit_source[m+'_full_delta'],rtol=0,atol=1e-12,equal_nan=True)

def stats(data,m):
    z=data[data[m+'_valid']]
    x=z[m+'_md_delta'].to_numpy();y=z[m+'_model_delta'].to_numpy()
    robust=abs(x)>2*z[m+'_sem'].to_numpy()
    return {'n':len(z),'substituted_n':int((~z.same_aa).sum()),'rho':float(spearmanr(x,y).statistic) if len(x)>=4 else None,
      'MAE':float(np.mean(abs(x-y))),'zero_MAE':float(np.mean(abs(x))),'gain':float(1-np.mean(abs(x-y))/np.mean(abs(x))),
      'MD_mean_delta':float(np.mean(x)),'model_mean_delta':float(np.mean(y)),
      'direction_fraction':float(np.mean(x*y>0)),'robust_n':int(robust.sum()),'robust_direction_fraction':float(np.mean(x[robust]*y[robust]>0)) if robust.any() else None,
      'all3_nonoverlap_n':int(z[m+'_nonoverlap'].sum())}

cdr3=df.cdr.str.startswith('CDR3');pmhc=df.role.isin(['mhc','peptide'])&df.interface&df.same_aa
regions={'CDR3_all22':df[cdr3], 'CDR3_same14':df[cdr3&df.same_aa],'CDR3alpha_all9':df[df.cdr=='CDR3alpha'],'CDR3beta_all13':df[df.cdr=='CDR3beta'],
         'all6_CDR':df[df.cdr!='outside_CDR'],'pMHC_interface':df[pmhc],
         'pMHC_interface_without_MHC_R65':df[pmhc&~((df.role=='mhc')&(df.id_a=='65'))],
         'MHC_interface':df[df.interface&(df.role=='mhc')],'TCR_interface_sameAA':df[df.interface&df.role.str.startswith('tcr')&df.same_aa],
         'CDR3_sameAA_sidechain':df[cdr3&df.same_aa]}
stats_out={name:{m:stats(part,m) for m in ['rmsf','hchi']} for name,part in regions.items()}
plotcdr=pd.read_csv(MEL/'plot_sources/mel5_all6_CDR_per_position_source.csv',dtype={'id_a':str,'id_b':str})
plotpmhc=pd.read_csv(MEL/'plot_sources/mel5_natural_pMHC_interface_source.csv',dtype={'id_a':str,'id_b':str})
assert set(zip(plotcdr.role,plotcdr.id_a,plotcdr.id_b))==set(zip(df.loc[df.cdr!='outside_CDR','role'],df.loc[df.cdr!='outside_CDR','id_a'],df.loc[df.cdr!='outside_CDR','id_b']))
assert set(zip(plotpmhc.role,plotpmhc.id_a,plotpmhc.id_b))==set(zip(df.loc[pmhc,'role'],df.loc[pmhc,'id_a'],df.loc[pmhc,'id_b']))
assert len(df[cdr3])==22 and (~df[cdr3].same_aa).sum()==8
assert stats_out['pMHC_interface']['hchi']['n']==23
out={'status':'passed','independent_reconstruction':'Original family NPZ, new MD per-replica NPZ, direct CA8A interface union from both model-input PDBs, global full-chain alignment; original outputs never changed.',
 'aligned_residues':len(df),'all6_CDR_positions':int((df.cdr!='outside_CDR').sum()),'all_CDR3_positions':22,'CDR3_substitutions':8,
 'plot_source_completeness_passed':True,'interface_geometry_matches_annotation':True,'statistics':stats_out,'alignment':alignments,
 'method_warning':'abs(delta)>2SEM is a descriptive repeat-to-repeat threshold, not a calibrated95%CI or multiple-test-adjusted significance claim.',
 'source_SHA256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [MEL/'corrected_pair_all_aligned_positions.csv',MEL/'corrected_predefined_region_change_metrics.csv',MEL/'corrected_alignment_and_CDR_audit.json',MEL/'plot_sources/mel5_all6_CDR_per_position_source.csv',MEL/'plot_sources/mel5_natural_pMHC_interface_source.csv']}}
(ROOT/'independent_MEL5_review_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':out['status'],'statistics':stats_out},indent=2),flush=True)
