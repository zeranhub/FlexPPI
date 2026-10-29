from pathlib import Path
import sys, json, csv, struct, hashlib
from collections import OrderedDict

ROOT = Path(__file__).resolve().parent
EXP = ROOT.parents[1]
sys.path.insert(0, str(EXP / 'label_sensitivity/runtime'))
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

SOURCE = Path('E:/data/CatchBond')
CORRECTED = Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
inventory = json.loads((CORRECTED / 'sequence_audit/chain_inventory.json').read_text())
train_records = [json.loads(l) for l in (EXP / 'family_folds/splits/fold_0_train.jsonl').read_text().splitlines()]
train_ids = {r['sample_id'] for r in train_records}
train_chains = [r for r in inventory if r['sample_id'] in train_ids]
aligner = PairwiseAligner(mode='local', match_score=2, mismatch_score=-1, open_gap_score=-3, extend_gap_score=-0.5)

def dcd_header(path):
    with path.open('rb') as f:
        start = f.read(8)
        endian = '<' if struct.unpack('<i',start[:4])[0] == 84 else '>'
        data = f.read(80)
        header = struct.unpack(endian+'20i', data)
        f.read(4)
        title_len=struct.unpack(endian+'i',f.read(4))[0]
        f.read(title_len+4)
        atom_record_len=struct.unpack(endian+'i',f.read(4))[0]
        n_atoms=struct.unpack(endian+'i',f.read(atom_record_len))[0]
    return {'path':str(path),'bytes':path.stat().st_size,'frames':header[0], 'n_atoms':n_atoms,
            'istart':header[1], 'nsavc':header[2], 'nominal_frame_interval_ns_from_readme':10,
            'analysis_start_index_zero_based':7,'maximum_analysis_stop_index_inclusive_from_readme':110,
            'analysis_frames_using_indices_7_to_110':max(0,min(header[0],111)-7)}

def chains_from_pdb(path, output):
    chains=OrderedDict()
    atom_lines=[]
    custom={'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H','CYX':'C','CYM':'C','ASH':'D','GLH':'E','LYN':'K'}
    for line in path.read_text().splitlines():
        if not line.startswith(('ATOM  ','HETATM')): continue
        rn=line[17:20].strip()
        one=seq1(rn,custom_map=custom,undef_code='X')
        if one=='X' or len(one)!=1: continue
        ch=line[21:22]
        if ch not in chains: chains[ch]=OrderedDict()
        resid=line[22:27].strip()
        chains[ch][resid]=one
        canon={'HID':'HIS','HIE':'HIS','HIP':'HIS','HSD':'HIS','HSE':'HIS','HSP':'HIS','CYX':'CYS','CYM':'CYS','ASH':'ASP','GLH':'GLU','LYN':'LYS'}.get(rn,rn)
        atom_lines.append(line[:17]+canon+line[20:])
    output.write_text('\n'.join(atom_lines)+'\nEND\n')
    return [{'chain':c,'length':len(rr),'sequence':''.join(rr.values()),'first_residue':next(iter(rr)), 'last_residue':next(reversed(rr))} for c,rr in chains.items()]

summaries=[]; chain_rows=[]
for condition in ['WT','S32M']:
    folder=SOURCE / condition / condition
    inp=folder / 'system_renum.pdb'
    clean=ROOT / f'catchbond_{condition}_prepared_protein.pdb'
    chains=chains_from_pdb(inp,clean)
    trajectories=[dcd_header(p) for p in sorted(folder.glob('*.dcd'))]
    for chain in chains:
        query=chain['sequence']
        hits=[]
        exact=[r for r in train_chains if r['sequence']==query]
        for record in train_chains:
            target=record['sequence']
            if min(len(query),len(target))/max(len(query),len(target)) < .5: continue
            al=aligner.align(query,target)[0]
            blocks_q,blocks_t=al.aligned
            aligned=sum(int(e-s) for s,e in blocks_q)
            if aligned==0: continue
            matches=sum(sum(a==b for a,b in zip(query[int(sq):int(eq)],target[int(st):int(et)]))
                        for (sq,eq),(st,et) in zip(blocks_q,blocks_t))
            cov=aligned/len(query); identity=matches/aligned
            if cov >= .7: hits.append((identity,cov,record['sample_id'],record['chain'],record['biological_system']))
        hits.sort(reverse=True)
        chain.update({'condition':condition,'exact_fold0_train_chains':exact,'top10_local_alignment_train_hits':hits[:10]})
        chain_rows.append(chain)
    summaries.append({'condition':condition,'original_system_pdb':str(inp),'prmtop':str(folder/'system.prmtop'),
        'protein_only_canonical_prepared_pdb':str(clean),'prepared_pdb_sha256':hashlib.sha256(clean.read_bytes()).hexdigest(),
        'trajectory_count':len(trajectories),'trajectories':trajectories,'chains':chains,
        'caveat':'Prepared structure accompanies MD; use deposited experimental 9NMU/9NMX for primary structure-based inference. DCD stored every 10 ns and source readme reports water-oxygen Lennard-Jones parameter assignment error.'})

out={'source':'https://zenodo.org/records/17065125','paper_doi':'10.1126/science.adx3162',
     'paper_fulltext_local':'E:/Project/FlexPPI/artifacts/tcr_pmhc_reassessment_20260915/literature/chen2026.txt',
     'tcr':'TCR156','peptide':'TLMSAMTNL','MHC':'HLA-A*02:01','experimental_PDBs':{'WT':'9NMU','S32M':'9NMX'},
     'fold0_train_records':len(train_records),'fold0_train_chains':len(train_chains),'conditions':summaries}
(ROOT/'catchbond_inventory.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
with (ROOT/'catchbond_chain_train_overlap.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=['condition','chain','length','sequence','first_residue','last_residue','exact_fold0_train_count','nearest_train_identity','nearest_train_coverage','nearest_train_sample','nearest_train_chain','nearest_train_system'])
    w.writeheader()
    for row in chain_rows:
        hit=row['top10_local_alignment_train_hits'][0] if row['top10_local_alignment_train_hits'] else [None]*5
        w.writerow({**{k:row[k] for k in ['condition','chain','length','sequence','first_residue','last_residue']},'exact_fold0_train_count':len(row['exact_fold0_train_chains']),
                    **dict(zip(['nearest_train_identity','nearest_train_coverage','nearest_train_sample','nearest_train_chain','nearest_train_system'],hit))})
print(json.dumps({'conditions':[{'condition':r['condition'],'trajectory_count':r['trajectory_count'], 'frame_counts':[x['frames'] for x in r['trajectories']], 'prepared_pdb':r['protein_only_canonical_prepared_pdb'],'chain_lengths':[(c['chain'],c['length']) for c in r['chains']]} for r in summaries], 'train_overlap_csv':str(ROOT/'catchbond_chain_train_overlap.csv')},indent=2))
