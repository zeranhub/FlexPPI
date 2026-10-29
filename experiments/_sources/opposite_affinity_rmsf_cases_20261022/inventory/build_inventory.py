"""Read existing case predictions; preserve complete regions and adverse contrasts."""
from pathlib import Path
import sys, json, hashlib, math
import zipfile,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

OUT=Path(__file__).resolve().parent
BASE=OUT.parent.parent
ATLAS=BASE/'candidate_case_atlas_20261021'
sys.path.insert(0,str(ATLAS/'data_and_general_figures'))
import prepare_all_data as prior
CAT=json.loads((ATLAS/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
SOURCE=[]
def source(p,kind):
    p=Path(p)
    x={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'kind':kind}
    if x not in SOURCE:SOURCE.append(x)
    return x
def dump(name,rows):
    df=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    df.to_csv(OUT/name,index=False,encoding='utf-8-sig',float_format='%.17g')
    return df
source(ATLAS/'catalog/candidate_inventory.json','prior exhaustive 30-series inventory')
paper=BASE.parent/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'
source(paper,'current manuscript read only')
with zipfile.ZipFile(paper) as archive:
    xml=ET.fromstring(archive.read('word/document.xml'))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    paragraphs=[''.join(p.itertext()) for p in xml.findall('.//w:p',ns)]
    contexts=[{'paragraph_1based':i+1,'text':text} for i,text in enumerate(paragraphs) if any(k in text for k in ['MEL5','α24β17','Q72','Gly54'])]
    (OUT/'current_manuscript_case_context.json').write_text(json.dumps({'path':str(paper),'sha256':hashlib.sha256(paper.read_bytes()).hexdigest(),'matching_paragraphs':contexts},ensure_ascii=False,indent=2),encoding='utf8')
source(prior.ANNOT,'prior per-residue biological annotation')
source(prior.MAN,'ITrace trajectory inventory') if isinstance(prior.MAN,Path) else source('E:/data/ITrace/tables/manifest.csv','ITrace trajectory inventory')

inv=[]
for c in CAT['cases']:
    for a,b in c['contrasts']:
        inv.append({'case_id':c['case_id'],'case_name':c['name'],'tier':c['tier'],'pdb_before':a,'pdb_after':b,
                    'background':c['background'],'existing_support':c['model_can_support'],
                    'existing_limitations':c['model_cannot_support_or_limitations'],
                    'structure_input_type':'Each observed experimental complex structure; not a fixed parent backbone',
                    'model_checkpoint_sha256':CAT['model']['checkpoint_sha256'],
                    'affinity_endpoint_verified':'See background and primary source; not all cases have direct KD',
                    'primary_sources_json':json.dumps(c['local_primary_sources'],ensure_ascii=False)})
dump('all_30_series_contrast_inventory.csv',inv)

PDBS=['3hg1','4jff','4jfd','4jfe']
KD={'3hg1':18e-6,'4jff':.6e-9,'4jfd':36e-6,'4jfe':31e-6}
STATE={'3hg1':'MEL5 + ELA','4jff':'alpha24beta17 + ELA','4jfd':'alpha24beta17 + ELA4A','4jfe':'alpha24beta17 + ELA7A'}
PAIRS=[('3hg1','4jff','TCR affinity maturation'),('4jff','4jfd','peptide G4A'),('4jff','4jfe','peptide I7A')]
frames={}; audits=[]; seqs=[]
for p in PDBS:
    f,a=prior.load(p)
    f=f.copy();f['pdb']=p;f['state']=STATE[p]
    z=np.load(a['prediction'],allow_pickle=False)
    assert z['mean'].shape[1]==2 and z['target'].shape[1]==3
    for k,m in enumerate(['rmsf','hchi']):
        assert np.array_equal(f[m+'_valid'].to_numpy(),z['valid'][:,k]&np.isfinite(z['target'][:,k]))
        f[m+'_model_variance']=z['model_variance'][:,k].astype(float)
    frames[p]=f
    source(a['prediction'],'existing 2430 joint Full fold0 prediction')
    source(a['structure'],'original ITrace run1 topology for static geometry')
    audits.append(a)
    for role,g in f.groupby('role',sort=False):
        seqs.append({'pdb':p,'state':STATE[p],'role':role,'chain':g.chain.iloc[0],'residue_count':len(g),'sequence':''.join(g.aa),'KD_M':KD[p],'KD_source':'Madura et al JBC 2013 main article/Table 1; parent MEL5 18uM discussed in text'})
dump('MEL5_four_states_all_residues.csv',pd.concat(frames.values(),ignore_index=True))
dump('MEL5_four_states_sequence_and_KD.csv',seqs)
dump('MEL5_four_states_complete_peptide.csv',pd.concat([f[f.role.eq('peptide')] for f in frames.values()],ignore_index=True))
manifest_path=BASE/'tcr_case_v40_20261017/family_predictions/prediction_manifest.jsonl'
source(manifest_path,'checkpoint and input provenance for all four predictions')
manifest_records=[json.loads(line) for line in manifest_path.read_text(encoding='utf8').splitlines() if line.strip()]
provenance=[]
for p in PDBS:
    records=[r for r in manifest_records if r['sample_id']=='itrace_'+p and r['variant']=='full']
    assert len(records)==1 and records[0]['checkpoint_sha256']==CAT['model']['checkpoint_sha256']
    provenance.extend(records)

# One shared interface scope from all four original annotations, projected to 3HG1.
# The original annotation uses nearest inter-partner C-alpha distance <= 8 Angstrom.
def alignment(a,b,role):
    al=PairwiseAligner();al.mode='global';al.match_score=2;al.mismatch_score=-1
    al.open_gap_score=-6;al.extend_gap_score=-.5
    if role!='peptide':al.target_end_gap_score=0;al.query_end_gap_score=0
    return al.align(''.join(a.aa),''.join(b.aa))[0]
ref=frames['3hg1']; union=set(); maps={}
for p,f in frames.items():
    mp={}
    for role in ['mhc','b2m','peptide','tcr_alpha','tcr_beta']:
        ra=ref[ref.role.eq(role)].reset_index(drop=True);rb=f[f.role.eq(role)].reset_index(drop=True)
        hit=alignment(ra,rb,role)
        for ia,ib in zip(*hit.indices):
            if ia>=0 and ib>=0:
                a=ra.iloc[ia];b=rb.iloc[ib];mp[int(b['index'])]=int(a['index'])
                if role in ['mhc','peptide'] and bool(b.interface):union.add(int(a['index']))
    maps[p]=mp

aligned=[]; summaries=[]; neighborhood=[]
for a,b,trigger in PAIRS:
    af,bf=frames[a],frames[b]
    for role in ['mhc','b2m','peptide','tcr_alpha','tcr_beta']:
        aa=af[af.role.eq(role)].reset_index(drop=True);bb=bf[bf.role.eq(role)].reset_index(drop=True)
        hit=alignment(aa,bb,role)
        for col,(i,j) in enumerate(zip(*hit.indices)):
            ar=aa.iloc[i] if i>=0 else None;br=bb.iloc[j] if j>=0 else None
            row={'pdb_before':a,'pdb_after':b,'trigger':trigger,'role':role,'alignment_column':col+1,
                 'paired':ar is not None and br is not None,'same_aa':bool(ar is not None and br is not None and ar.aa==br.aa),
                 'affinity_fold_increase':KD[a]/KD[b],'delta_pKD':math.log10(KD[a]/KD[b])}
            for suffix,r in [('before',ar),('after',br)]:
                for key in ['index','chain','id','aa','region','interface']:
                    row[key+'_'+suffix]=r[key] if r is not None else ''
                for m in ['rmsf','hchi']:
                    for key in ['full','md','var','valid','model_variance']:
                        row[m+'_'+key+'_'+suffix]=r[m+'_'+key] if r is not None else (False if key=='valid' else np.nan)
            ri=maps[a].get(int(ar['index'])) if ar is not None else None
            row['four_state_pMHC_interface_8A']=bool(role in ['mhc','peptide'] and ri in union)
            for m in ['rmsf','hchi']:
                row[m+'_valid']=bool(row[m+'_valid_before'] and row[m+'_valid_after'])
                row[m+'_change_after_minus_before']=row[m+'_full_after']-row[m+'_full_before']
                row[m+'_drop_before_minus_after']=-row[m+'_change_after_minus_before']
                row[m+'_MD_change_after_minus_before']=row[m+'_md_after']-row[m+'_md_before']
            aligned.append(row)
df=dump('MEL5_three_contrasts_all_aligned_residues.csv',aligned)
pep=dump('MEL5_three_contrasts_complete_peptide.csv',df[df.role.eq('peptide')])
cdr=dump('MEL5_three_contrasts_all_six_CDR.csv',df[df.role.isin(['tcr_alpha','tcr_beta'])&df.region_before.astype(str).str.startswith('CDR')])
iface=dump('MEL5_three_contrasts_four_state_union_pMHC_interface.csv',df[df.four_state_pMHC_interface_8A])
scopes=['complete_peptide','all_six_CDR','four_state_union_pMHC_interface','CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
for a,b,trigger in PAIRS:
    pair=df[(df.pdb_before==a)&(df.pdb_after==b)]
    for scope in scopes:
        if scope=='complete_peptide':q=pair[pair.role.eq('peptide')]
        elif scope=='all_six_CDR':q=pair[pair.role.isin(['tcr_alpha','tcr_beta'])&pair.region_before.astype(str).str.startswith('CDR')]
        elif scope=='four_state_union_pMHC_interface':q=pair[pair.four_state_pMHC_interface_8A]
        else:q=pair[pair.region_before.eq(scope)]
        for m in ['rmsf','hchi']:
            v=q[q[m+'_valid']&q.paired]
            if m=='hchi':v=v[v.same_aa]
            x=v[m+'_full_before'].to_numpy(float);y=v[m+'_full_after'].to_numpy(float)
            md=v[m+'_MD_change_after_minus_before'].to_numpy(float);pr=y-x
            imax=int(np.argmax(x-y)) if len(x) else None
            r=v.iloc[imax] if imax is not None else None
            summaries.append({'pdb_before':a,'pdb_after':b,'trigger':trigger,'affinity_fold_increase':KD[a]/KD[b],'scope':scope,'metric':m,
                              'n':len(v),'all_valid_Hchi_require_same_AA':m=='hchi','mean_before':np.mean(x) if len(x) else np.nan,'mean_after':np.mean(y) if len(x) else np.nan,
                              'mean_change_after_minus_before':np.mean(pr) if len(x) else np.nan,'n_decrease':int((y<x).sum()),'n_increase':int((y>x).sum()),
                              'maximum_drop':float((x-y)[imax]) if r is not None else np.nan,'maximum_drop_chain':r.chain_before if r is not None else '',
                              'maximum_drop_id':r.id_before if r is not None else '', 'maximum_drop_aa':r.aa_before if r is not None else '',
                              'change_Spearman_to_MD':prior.rho(pr,md),'change_MAE_to_MD':np.mean(abs(pr-md)) if len(x) else np.nan,
                              'zero_change_MAE_to_MD':np.mean(abs(md)) if len(x) else np.nan})
summ=dump('MEL5_three_contrasts_complete_scope_summaries.csv',summaries)

# All peptide-neighbor residues, without ranking on response: distances from original PDB C-alphas.
from Bio.PDB import PDBParser
contact_summaries=[]; heavy_contacts=[]
parser=PDBParser(QUIET=True)
for p in PDBS:
    f=frames[p];model=parser.get_structure(p,f'E:/data/ITrace/trajectories/{p}_run1_topology.pdb')[0]
    def residue(row):
        ch=model[str(row.chain)];rid=str(row.id)
        hits=[r for r in ch if str(r.id[1])+str(r.id[2]).strip()==rid]
        assert len(hits)==1,(p,row.chain,rid,len(hits))
        return hits[0]
    residues={int(r['index']):residue(r) for _,r in f.iterrows()}
    coords={i:r['CA'].coord.astype(float) for i,r in residues.items()}
    heavy={i:np.array([a.coord for a in r if a.element not in ['H','D']],float) for i,r in residues.items()}
    pp=f[f.role.eq('peptide')]
    for _,pr in pp.iterrows():
        tally={'pdb':p,'state':STATE[p],'peptide_id':pr.id,'peptide_aa':pr.aa}
        for label in ['HLA','TCR']:
            tally[label+'_CA8A_partner_residue_n']=0;tally[label+'_heavy4_5A_partner_residue_n']=0
            tally[label+'_heavy4_5A_atom_pair_n']=0;tally[label+'_minimum_heavy_distance_A']=float('inf')
        for _,r in f[f.role.isin(['mhc','tcr_alpha','tcr_beta'])].iterrows():
            d=float(np.linalg.norm(coords[int(pr['index'])]-coords[int(r['index'])]))
            label='HLA' if r.role=='mhc' else 'TCR'
            atomds=np.linalg.norm(heavy[int(pr['index'])][:,None,:]-heavy[int(r['index'])][None,:,:],axis=-1)
            hd=float(atomds.min());natom=int((atomds<=4.5).sum())
            tally[label+'_minimum_heavy_distance_A']=min(tally[label+'_minimum_heavy_distance_A'],hd)
            if natom:
                tally[label+'_heavy4_5A_partner_residue_n']+=1;tally[label+'_heavy4_5A_atom_pair_n']+=natom
                heavy_contacts.append({'pdb':p,'peptide_id':pr.id,'peptide_aa':pr.aa,'partner_chain':r.chain,'partner_role':r.role,'partner_id':r.id,'partner_aa':r.aa,
                                       'partner_region':r.region,'minimum_heavy_distance_A':hd,'atom_pair_n_within4_5A':natom,'CA_distance_A':d})
            if d<=8:
                tally[label+'_CA8A_partner_residue_n']+=1
                neighborhood.append({'pdb':p,'peptide_id':pr.id,'peptide_aa':pr.aa,'partner_chain':r.chain,'partner_role':r.role,'partner_id':r.id,'partner_aa':r.aa,
                                     'partner_region':r.region,'CA_distance_A':d,'partner_RMSF_A':r.rmsf_full,'partner_Hchi_nats':r.hchi_full,'partner_Hchi_valid':r.hchi_valid})
        contact_summaries.append(tally)
dump('MEL5_four_states_all_peptide_CA8A_neighbors.csv',neighborhood)
dump('MEL5_four_states_all10_peptide_contact_summary.csv',contact_summaries)
dump('MEL5_four_states_all_peptide_heavy4_5A_contacts.csv',heavy_contacts)

checks=[]
for p in ['3hg1','4jff']:
    path=BASE/f'tcr_expanded_screen_v41_20261017/mel5_validation/{p}_all_residue_replica_labels.npz'
    source(path,'previously independently recomputed three-replica MD labels')
    rep=np.load(path,allow_pickle=False);src=np.load(BASE/f'tcr_case_v40_20261017/family_predictions/full/itrace_{p}.npz')
    assert rep['metric_names'].tolist()==['rmsf','hchi']
    assert np.array_equal(rep['chain'],src['chain']) and np.array_equal(rep['author_residue_id'],src['author_residue_id'])
    errs=[]
    for j,m in enumerate(['rmsf','hchi']):
        mask=rep['valid'][:,j]&src['valid'][:,j]
        err=float(np.max(abs(rep['mean'][mask,j]-src['target'][mask,j])))
        assert err<5e-6,(p,m,err)
        errs.append({'metric':m,'target_column':j,'checked':int(mask.sum()),'max_abs_mean_error':err})
    checks.append({'pdb':p,'source_target_first2_columns_verified':True,'checks':errs})

source(BASE/'mel5_joint_case_v49_20261019/case_sources/primary_literature/mel5_thermodynamics_2013.md','locally archived original JBC article')
audit={'status':'passed','date':'2026-10-22','new_inference':0,'new_MD':0,'new_training':0,'main_manuscript_edited':False,
       'checkpoint_sha256':CAT['model']['checkpoint_sha256'],'training_records':2430,'training_ITrace':0,
       'fixed_pairs':PAIRS,'fixed_regions':scopes,'four_state_union_interface_ref_position_n':len(union),'MD_target_mapping_checks':checks,
       'target_column_meaning':{'0':'RMSF Angstrom','1':'joint sidechain Hchi nats','2':'unused third metric; not Hchi'},
       'metadata':audits,'original_prediction_manifest_records':provenance,'sources':SOURCE,'scope':'Exploratory case comparison from existing observed complex structures, not an affinity predictor validation.'}
(OUT/'inventory_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(summ[summ.scope.eq('complete_peptide')].to_string(index=False))
print('FOUR_STATE_INTERFACE_N',len(union),'sourceN',len(SOURCE))
