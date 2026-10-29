from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'mel5_validation';PLOT=OUT/'plot_sources';PLOT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT.parent/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
SOURCE=ROOT.parent/'tcr_case_v40_20261017'
annotation_path=Path(r'E:\Project\FlexPPI\output\itrace_residue_biology\raw\matched_residue_context_responses.csv')
cols=['pdb_id','condition','chain','author_residue_id','amino_acid_1','region','is_tcr_pmhc_interface_8A']
annotations=[]
for chunk in pd.read_csv(annotation_path,usecols=cols,dtype={'author_residue_id':str,'chain':str},chunksize=100000):
    z=chunk[chunk.pdb_id.isin(['3hg1','4jff'])&(chunk.condition=='full')]
    if len(z):annotations.append(z.drop(columns='condition'))
annotations=pd.concat(annotations,ignore_index=True);annotations.to_csv(OUT/'full_input_region_annotations.csv',index=False)
roles=json.loads(Path(r'E:\Project\FlexPPI\output\itrace_benchmark\raw\itrace_chain_roles.json').read_text())
states={}
for pdb in ['3hg1','4jff']:
    with np.load(SOURCE/f'family_predictions/full/itrace_{pdb}.npz',allow_pickle=False) as z:p={k:z[k].copy() for k in z.files}
    with np.load(OUT/f'{pdb}_all_residue_replica_labels.npz',allow_pickle=False) as z:v=z['values'].copy();means=z['mean'].copy();variances=z['sample_variance'].copy()
    f=pd.DataFrame({'pdb_id':pdb,'index':np.arange(len(p['chain'])),'chain':p['chain'],'author_residue_id':p['author_residue_id'],'aa':p['amino_acid']});f['role']=f.chain.map(roles[pdb]);f['position']=f.groupby('chain').cumcount()+1
    f=f.merge(annotations,on=['pdb_id','chain','author_residue_id'],how='left',validate='one_to_one');assert (f.aa==f.amino_acid_1).all()
    for j,m in enumerate(['rmsf','hchi']):
        f[m+'_MD']=means[:,j];f[m+'_SD']=np.sqrt(variances[:,j]);f[m+'_model']=p['mean'][:,j];f[m+'_valid']=p['valid'][:,j]
        for rep in (1,2,3):f[f'{m}_rep{rep}']=v[rep-1,:,j]
    states[pdb]=f

aligner=PairwiseAligner(mode='global',match_score=2,mismatch_score=-2,open_gap_score=-6,extend_gap_score=-.5)
frames=[];alignment_audit=[]
for role in ['mhc','b2m','peptide','tcr_alpha','tcr_beta']:
    a=states['3hg1'][states['3hg1'].role==role].reset_index(drop=True);b=states['4jff'][states['4jff'].role==role].reset_index(drop=True)
    sa=''.join(a.aa);sb=''.join(b.aa);aln=aligner.align(sa,sb)[0];pairs=[]
    for aa,bb in zip(*aln.aligned):pairs.extend(zip(range(*aa),range(*bb)))
    alignment_audit.append({'role':role,'WT_sequence':sa,'engineered_sequence':sb,'WT_length':len(sa),'engineered_length':len(sb),'aligned_positions':len(pairs),'mismatch_n':sum(sa[i]!=sb[j] for i,j in pairs),'aligned_blocks':aln.aligned.tolist(),'alignment_string':str(aln)})
    for i,j in pairs:
        ra,rb=a.iloc[i],b.iloc[j];region=str(ra.region)
        row={'a':'3hg1','b':'4jff','role':role,'chain_a':ra.chain,'chain_b':rb.chain,'id_a':ra.author_residue_id,'id_b':rb.author_residue_id,'position_a':ra.position,'position_b':rb.position,'aa_a':ra.aa,'aa_b':rb.aa,'same_aa':ra.aa==rb.aa,'region':region,'region_a_original':ra.region,'region_b_original':rb.region,'region_a':region,'region_b':region,'interface_a':bool(ra.is_tcr_pmhc_interface_8A),'interface_b':bool(rb.is_tcr_pmhc_interface_8A),'interface':bool(ra.is_tcr_pmhc_interface_8A or rb.is_tcr_pmhc_interface_8A)}
        for m in ['rmsf','hchi']:
            row[m+'_valid']=bool(ra[m+'_valid'] and rb[m+'_valid'])
            for side,r in [('a',ra),('b',rb)]:
                row[m+'_md_'+side]=float(r[m+'_MD']);row[m+'_sd_'+side]=float(r[m+'_SD']);row[m+'_full_'+side]=float(r[m+'_model'])
                for rep in (1,2,3):row[f'{m}_rep{rep}_{side}']=float(r[f'{m}_rep{rep}'])
            row[m+'_md_delta']=row[m+'_md_b']-row[m+'_md_a'];row[m+'_full_delta']=row[m+'_full_b']-row[m+'_full_a'];row[m+'_delta_sem']=np.sqrt((row[m+'_sd_a']**2+row[m+'_sd_b']**2)/3)
            av=np.array([row[f'{m}_rep{r}_a'] for r in (1,2,3)]);bv=np.array([row[f'{m}_rep{r}_b'] for r in (1,2,3)])
            row[m+'_all3_b_above_all3_a']=bool(bv.min()>av.max());row[m+'_all3_b_below_all3_a']=bool(bv.max()<av.min())
        frames.append(row)
frame=pd.DataFrame(frames);frame.to_csv(OUT/'corrected_pair_all_aligned_positions.csv',index=False)
same=frame.same_aa;iface=frame.interface;mhc=frame.role=='mhc';pep=frame.role=='peptide';tcr=frame.role.isin(['tcr_alpha','tcr_beta']);cdr=frame.region.str.startswith('CDR');cdr3=frame.region.str.contains('CDR3')
scopes={'all_common_AA':same,'all_interface_common_AA':same&iface,'MHC_interface_natural':mhc&iface,'MHC_interface_common_AA':same&mhc&iface,'peptide_all_natural':pep,'peptide_common_AA':same&pep,'MHC_plus_peptide_interface_natural':iface&(mhc|pep),'MHC_plus_peptide_interface_common_AA':same&iface&(mhc|pep),'MHC_groove_common_AA':same&mhc&(frame.position_a<=180),'TCR_interface_all_aligned':tcr&iface,'TCR_interface_common_AA':same&tcr&iface,'TCR_variable_common_AA':same&tcr&(frame.position_a<=125),'TCR_CDR_all_aligned':cdr,'TCR_CDR_common_AA':same&cdr,'CDR3_all_aligned':cdr3,'CDR3_common_AA':same&cdr3}
loops=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
for loop in loops:
    scopes[loop+'_all_aligned']=frame.region==loop;scopes[loop+'_common_AA']=same&(frame.region==loop)
for name in ['α1 helix','α2 helix']:scopes[name+'_common_AA']=same&(frame.region==name)

def rho(x,y):return float(pd.Series(x).rank().corr(pd.Series(y).rank())) if len(x)>=4 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
stats=[]
for scope,mask in scopes.items():
    for m in ['rmsf','hchi']:
        f=frame[mask&frame[m+'_valid']];x=f[m+'_md_delta'].to_numpy();y=f[m+'_full_delta'].to_numpy();sem=f[m+'_delta_sem'].to_numpy();robust=abs(x)>2*sem
        if not len(f):continue
        stats.append({'scope':scope,'metric':m,'n':len(f),'mutated_amino_acid_n':int((~f.same_aa).sum()),'rho_delta':rho(x,y),'MAE_delta':float(np.mean(abs(x-y))),'zero_change_MAE':float(np.mean(abs(x))),'relative_MAE_gain':float(1-np.mean(abs(x-y))/np.mean(abs(x))) if np.mean(abs(x))>0 else np.nan,'MD_mean_a':float(f[m+'_md_a'].mean()),'MD_mean_b':float(f[m+'_md_b'].mean()),'model_mean_a':float(f[m+'_full_a'].mean()),'model_mean_b':float(f[m+'_full_b'].mean()),'MD_mean_delta':float(x.mean()),'model_mean_delta':float(y.mean()),'robust_n':int(robust.sum()),'robust_direction_fraction':float(np.mean(x[robust]*y[robust]>0)) if robust.any() else np.nan,'all3_nonoverlap_n':int((f[m+'_all3_b_above_all3_a']|f[m+'_all3_b_below_all3_a']).sum())})
stat=pd.DataFrame(stats);stat.to_csv(OUT/'corrected_predefined_region_change_metrics.csv',index=False);stat[stat.scope.str.startswith('CDR')|stat.scope.str.startswith('TCR_CDR')].to_csv(PLOT/'mel5_all6_CDR_summary.csv',index=False)
frame[cdr].to_csv(PLOT/'mel5_all6_CDR_per_position_source.csv',index=False)
frame[iface&(mhc|pep)].to_csv(PLOT/'mel5_natural_pMHC_interface_source.csv',index=False)
definitions=[]
for loop in loops:
    f=frame[frame.region==loop]
    for side,pdb in [('a','3hg1'),('b','4jff')]:
        definitions.append({'loop':loop,'pdb':pdb,'chain':f['chain_'+side].iloc[0],'author_residue_ids':','.join(f['id_'+side].astype(str)),'author_range':f['id_'+side].iloc[0]+'–'+f['id_'+side].iloc[-1],'sequence':''.join(f['aa_'+side]),'length':len(f),'valid_rmsf_positions':int(f.rmsf_valid.sum()),'valid_hchi_in_both_states':int(f.hchi_valid.sum()),'common_AA_positions':int(f.same_aa.sum()),'substituted_positions':int((~f.same_aa).sum()),'definition':'WT ITrace boundary projected through globally aligned full chain; no internal gaps in either TCR chain'})
pd.DataFrame(definitions).to_csv(PLOT/'mel5_CDR_definitions.csv',index=False)
long=[]
for r in frame[cdr].itertuples():
    d=r._asdict()
    for m in ['rmsf','hchi']:
        for side,state in [('a','WT'),('b','Engineered')]:
            for rep in (1,2,3):long.append({'loop':d['region'],'WT_residue_id':d['id_a'],'engineered_residue_id':d['id_b'],'AA_WT':d['aa_a'],'AA_engineered':d['aa_b'],'same_aa':d['same_aa'],'state':state,'metric':m,'source':'MD','replica':rep,'value':d[f'{m}_rep{rep}_{side}']})
            long.append({'loop':d['region'],'WT_residue_id':d['id_a'],'engineered_residue_id':d['id_b'],'AA_WT':d['aa_a'],'AA_engineered':d['aa_b'],'same_aa':d['same_aa'],'state':state,'metric':m,'source':'Model','replica':0,'value':d[f'{m}_full_{side}']})
pd.DataFrame(long).to_csv(PLOT/'mel5_all6_CDR_state_values_long.csv',index=False)
k=frame[mhc&(frame.id_a=='68')].iloc[0];ks=[]
for side,pdb,label in [('a','3hg1','WT'),('b','4jff','Engineered')]:
    for rep in (1,2,3):ks.append({'pdb':pdb,'state':label,'source':'MD','replica':rep,'Hchi_nats':k[f'hchi_rep{rep}_{side}']})
    ks.append({'pdb':pdb,'state':label,'source':'MD_mean','replica':0,'Hchi_nats':k[f'hchi_md_{side}']});ks.append({'pdb':pdb,'state':label,'source':'Model','replica':0,'Hchi_nats':k[f'hchi_full_{side}']})
pd.DataFrame(ks).to_csv(PLOT/'mel5_MHC_K68_absolute_replica_source.csv',index=False)
contributions=[]
for m in ['rmsf','hchi']:
    f=frame[same&iface&frame[m+'_valid']].copy();x=f[m+'_md_delta'].to_numpy();y=f[m+'_full_delta'].to_numpy();rx=pd.Series(x).rank().to_numpy();ry=pd.Series(y).rank().to_numpy();cross=(rx-rx.mean())*(ry-ry.mean())
    for group,mask in [('MHC',f.role=='mhc'),('peptide',f.role=='peptide'),('TCR',f.role.isin(['tcr_alpha','tcr_beta']))]:
        z=f[mask];xx=z[m+'_md_delta'].to_numpy();yy=z[m+'_full_delta'].to_numpy();contributions.append({'metric':m,'component':group,'n':len(z),'sum_absolute_MD_change':float(np.sum(abs(xx))),'sum_absolute_prediction_error':float(np.sum(abs(xx-yy))),'sum_MAE_improvement_over_zero':float(np.sum(abs(xx)-abs(xx-yy))),'fraction_of_global_rank_covariance_numerator':float(cross[mask.to_numpy()].sum()/cross.sum())})
pd.DataFrame(contributions).to_csv(OUT/'corrected_interface_component_decomposition.csv',index=False)
old=pd.read_csv(OUT/'predefined_region_change_metrics.csv');joined=old.merge(stat,on=['scope','metric'],suffixes=('_old','_corrected'));joined.to_csv(OUT/'old_vs_corrected_region_metrics_audit.csv',index=False)
audit={'status':'complete','alignment':'Biopython PairwiseAligner global match2 mismatch−2 gap_open−6 gap_extend−0.5. Both TCR chains align WT residue i to engineered i+1 with no internal gaps.','boundary_rule':'All six WT ITrace CDR definitions are projected to engineered TCR through the same global alignment. No best-loop or same-AA selection is used for all_aligned RMSF. Common-AA Hchi is reported separately.','WT_CDR_author_ranges':{'alpha':{'CDR1':[26,31],'CDR2':[49,54],'CDR3':[89,97]},'beta':{'CDR1':[27,32],'CDR2':[50,54],'CDR3':[92,104]}},'engineered_CDR_author_ranges':{'alpha':{'CDR1':[27,32],'CDR2':[50,55],'CDR3':[90,98]},'beta':{'CDR1':[28,33],'CDR2':[51,55],'CDR3':[93,105]}},'old_issue':'SequenceMatcher skipped and misaligned engineered replacement blocks; 4JFF ITrace beta CDR annotations used WT author ranges despite one extra N-terminal residue. Earlier14-point CDR3 statistics are superseded, although common-AA count remains14.','all_aligned_CDR3_positions':22,'CDR3_substituted_positions':8,'CDR3_common_AA_positions':14,'all_aligned_residue_pairs':len(frame),'alignment_details':alignment_audit,'annotation_sha256':hashlib.sha256(annotation_path.read_bytes()).hexdigest(),'protein_trajectory_validation':'Existing raw_trajectory_audit.json passed; this is an alignment/region correction, no MD was rerun.'}
(OUT/'corrected_alignment_and_CDR_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(stat[stat.scope.isin(['all_interface_common_AA','MHC_plus_peptide_interface_common_AA','TCR_interface_common_AA','CDR3_all_aligned','CDR3_common_AA'])].to_string(index=False),flush=True)
