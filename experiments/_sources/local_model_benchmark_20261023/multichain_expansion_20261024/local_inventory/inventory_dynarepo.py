"""Bounded, read-only local inventory. Read small PDB/JSON only; stat trajectories."""
from pathlib import Path
import collections,datetime,hashlib,json,sys
import pandas as pd

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
RAW=Path('D:/data/dynarepo_non_nucleic_full')
MANIFEST=Path('D:/data/dynarepo_non_nucleic_manifest.json')
sys.path.insert(0,str(ROOT/'grouping_protein_complex_20261024/runtime'))
from Bio.Align import PairwiseAligner
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
AA.update({'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','CYX':'C','CYM':'C','LYN':'K','MSE':'M'})
AL=PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,open_gap_score=-5,extend_gap_score=-.1)
AL.query_left_gap_score=0;AL.query_right_gap_score=0

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stat(p):return {'path':str(p),'exists':p.is_file(),'bytes':p.stat().st_size if p.is_file() else 0}
def chains(p):
    result={};seen=set()
    for ln in p.read_text(encoding='utf-8',errors='replace').splitlines():
        if ln.startswith('ENDMDL'):break
        if ln[:6].strip() not in ['ATOM','HETATM'] or ln[12:16].strip()!='CA' or ln[17:20].strip() not in AA:continue
        key=(ln[21:22],ln[22:27])
        if key in seen:continue
        seen.add(key);result[key[0]]=result.get(key[0],'')+AA[ln[17:20].strip()]
    return result
def compare(first,second):
    a,b=(first,second) if len(first)>=len(second) else (second,first)
    alignment=AL.align(a,b)[0];co=alignment.coordinates;matches=0
    for j in range(co.shape[1]-1):
        t0,t1=map(int,co[0,j:j+2]);q0,q1=map(int,co[1,j:j+2])
        if t1>t0 and q1>q0:matches+=sum(u==v for u,v in zip(a[t0:t1],b[q0:q1]))
    return matches/min(len(a),len(b)),min(len(a),len(b))/max(len(a),len(b))
def as_csv(records,path):
    pd.DataFrame([{k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in records]).to_csv(path,index=False)

def main():
    cohort_path=ROOT/'cohort/records.csv';membership_path=ROOT/'grouping_protein_complex_20261024/membership_frozen.csv'
    cohort=pd.read_csv(cohort_path,keep_default_na=False);members=pd.read_csv(membership_path,keep_default_na=False)
    source=json.loads(MANIFEST.read_text(encoding='utf-8'));existing=set(cohort.sample_id)
    extra=[p for p in source['projects'] if p['accession'] not in existing]
    assert len(source['projects'])==922 and len(extra)==32
    inventory=[];chain_rows=[];hashes={str(MANIFEST):sha(MANIFEST),str(cohort_path):sha(cohort_path),str(membership_path):sha(membership_path)}
    for project in extra:
        sid=project['accession'];folder=RAW/sid;pdb=folder/'structure.pdb';seqs=chains(pdb)
        hashes[str(pdb)]=sha(pdb)
        first=next(iter(seqs.values()));pairs=[]
        for name,seq in seqs.items():
            identity,coverage=compare(first,seq)
            pairs.append({'chain':name,'length':len(seq),'matches_per_shorter_chain':identity,'length_coverage':coverage})
            chain_rows.append({'sample_id':sid,'pdb_id':';'.join(project['pdb_ids']),'input_chain':name,'sequence':seq,'residues':len(seq),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'structure_path':str(pdb)})
        candidate=len(seqs)>=2 and all(x['matches_per_shorter_chain']>=.97 and x['length_coverage']>=.8 for x in pairs)
        # Similarity is a candidate signal, not a substitute for official identity.
        classification='same_protein_homomer_candidate_needs_official_identity' if candidate else 'single_chain_exclude' if len(seqs)==1 else 'mixed_sequence_partners_not_a_homomer_candidate'
        inventory.append({'sample_id':sid,'pdb_ids':project['pdb_ids'],'name':project['name'],'source_description':project['description'],
            'directory':str(folder),'actual_protein_chains':len(seqs),'actual_residues':sum(map(len,seqs.values())),
            'chain_lengths':{k:len(v) for k,v in seqs.items()},'comparison_to_first_chain':pairs,
            'candidate_status':classification,'present_in_original_3063':False,
            'original_1024_residue_limit_exceeded':sum(map(len,seqs.values()))>1024,
            'source_declared_frames':project.get('frames'),'source_declared_temperature':project.get('temperature'),
            'source_declared_runs':project.get('mds'),
            'structure':stat(pdb),'trajectory':stat(folder/'trajectory.xtc'),'topology':stat(folder/'topology.tpr'),
            'precomputed_RMSF_label_path':'','label_status':'No standalone RMSF label identified in the project folder; real trajectory exists and would require a separately authorized matched-label analysis.',
            'independent_test_status':'Not one of the original 3063 exact records; family/sequence overlap and checkpoint exposure still need checking before claiming independent test.'})
    as_csv(inventory,OUT/'dynarepo_32_not_in_original_cohort.csv');pd.DataFrame(chain_rows).to_csv(OUT/'dynarepo_32_actual_chain_sequences.csv',index=False)
    (OUT/'dynarepo_32_not_in_original_cohort.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    # Existing ten are explicitly not new independent test proteins.
    joined=cohort.merge(members[['sample_id','benchmark_group']],on='sample_id',validate='one_to_one')
    old=joined[(joined.benchmark_group=='protein')&(joined.n_chains.astype(int)>1)]
    old_rows=[]
    for row in old.itertuples():
        folder=RAW/row.sample_id
        export=Path(row.full_prediction_path)
        corrected_structure=Path(row.structure_path)
        kind='source_homomer_or_near_identical_variant' if row.sample_id in ['A000J','A0046','A0088','A00C7'] else 'same_mature_protein_split_domains_or_engineered_fusion'
        old_rows.append({'sample_id':row.sample_id,'pdb_id':row.biological_system,'actual_input_n_chains':row.n_chains,'actual_input_n_residues':row.n_residues,
            'protein_interpretation':kind,'raw_structure':stat(folder/'structure.pdb'),'trajectory':stat(folder/'trajectory.xtc'),
            'topology':stat(folder/'topology.tpr'),'frozen_corrected_input_structure':stat(corrected_structure),
            'local_label_container':stat(export),'label_arrays':'target[:,0] with valid[:,0] in this previously frozen export; contains predictions too, not read in this inventory.',
            'original_remote_label_path':row.original_label_path,'original_label_sha256':row.original_label_sha256,
            'local_reference_table_path':row.reference_table_path,'reference_table_filter':row.sample_id,
            'original_fold':row.original_fold,'family_fold':row.family_fold,
            'present_in_original_3063':True,'independent_test_eligibility':'Not a new independent external test: this exact record is in the five-fold corpus and in four other folds training. Its corresponding original held-out model supplies record OOF only.'})
    assert len(old_rows)==10
    as_csv(old_rows,OUT/'existing_10_multichain_protein_local_assets.csv')
    status=collections.Counter(x['candidate_status'] for x in inventory)
    audit={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'read_scope':'32 previously omitted DynaRepo small topology PDBs plus metadata/stat of their files; existing ten cohort metadata/file stats; no trajectory contents read',
        'source_projects':922,'original_dynarepo_records':890,'not_in_original_records':32,'extra_records_actual_multichain':sum(x['actual_protein_chains']>=2 for x in inventory),
        'candidate_status_counts':dict(status),'promising_extra_homomer_candidates':[x['sample_id'] for x in inventory if x['candidate_status'].startswith('same_protein')],
        'already_in_original_multichain_proteins':10,'new_independent_test_homomers_confirmed':0,
        'network_used':False,'trajectory_opened_or_hashed':False,'new_labels_computed':False,'model_inference_run':False,
        'read_sources_sha256':hashes,'code_sha256':sha(Path(__file__))}
    (OUT/'local_inventory_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:audit[k] for k in ['source_projects','original_dynarepo_records','not_in_original_records','extra_records_actual_multichain','candidate_status_counts','promising_extra_homomer_candidates','new_independent_test_homomers_confirmed']},ensure_ascii=False))

if __name__=='__main__':main()
