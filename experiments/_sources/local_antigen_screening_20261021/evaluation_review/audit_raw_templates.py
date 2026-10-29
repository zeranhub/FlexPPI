"""Audit actual downloaded X-ray PDB chains, independently of ITrace topologies."""
from pathlib import Path
from collections import OrderedDict, defaultdict
import csv,hashlib,json,sys
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
REF=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-1,open_gap_score=-3,extend_gap_score=-.5)
records=[json.loads(x) for x in (EXP/'family_folds/splits/fold_0_train.jsonl').read_text().splitlines() if x.strip()]
ids={x['sample_id'] for x in records}
inventory=json.loads((REF/'sequence_audit/chain_inventory.json').read_text())
trainseq=defaultdict(list)
for c in inventory:
    if c['sample_id'] in ids: trainseq[c['sequence']].append({k:c[k] for k in ['sample_id','chain','biological_system','length']})
roles={'A':'mhc','B':'b2m','C':'peptide','D':'tcr_alpha','E':'tcr_beta'}
rows=[]
paths=sorted((ROOT.parent/'experimental_labels/raw/github/mutation_effects').glob('*/pdbs/*.pdb'))
for p in paths:
    chains=OrderedDict()
    for line in p.read_text().splitlines():
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM  ','HETATM')) or line[16:17] not in {' ','A','1'}: continue
        aa=seq1(line[17:20].strip(),custom_map={'MSE':'M'},undef_code='X')
        if aa=='X':continue
        key=(line[22:26].strip(),line[26:27])
        chains.setdefault(line[21:22],OrderedDict()).setdefault(key,aa)
    for ch,seqdict in chains.items():
        seq=''.join(seqdict.values())
        role=roles.get(ch,'unassigned_duplicate')
        for scope,query in [('full_chain',seq)]+([('first_120_residues_variable_domain_screen',seq[:120])] if role.startswith('tcr_') and len(seq)>120 else [])+([('first_180_residues_peptide_binding_domain_screen',seq[:180])] if role=='mhc' and len(seq)>180 else []):
            hits=[]
            for target,recordlist in trainseq.items():
                if query==target:
                    hits.append(dict(identity=1.0,query_coverage=1.0,target_coverage=1.0,records=recordlist));continue
                if len(query)<20 or min(len(query),len(target))/max(len(query),len(target))<.7:continue
                al=aligner.align(query,target)[0]
                qb,tb=al.aligned
                if len(qb)==0:continue
                qc=int(qb[-1,1]-qb[0,0])/len(query);tc=int(tb[-1,1]-tb[0,0])/len(target)
                if min(qc,tc)<.7:continue
                matched=sum(sum(a==b for a,b in zip(query[int(qs):int(qe)],target[int(ts):int(te)])) for (qs,qe),(ts,te) in zip(qb,tb))
                a,b=al.coordinates
                length=sum(max(int(a[i+1]-a[i]),int(b[i+1]-b[i])) for i in range(len(a)-1))
                identity=matched/length
                if identity>=.30:hits.append(dict(identity=identity,query_coverage=qc,target_coverage=tc,records=recordlist))
            rows.append(dict(pdb=p.stem,source_path=str(p),source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                chain=ch,role_from_reference_chain_mapping=role,scope=scope,sequence=query,length=len(query),
                exact_train_matches=trainseq.get(query,[]),local_30pct_70pct_both_training_hits=hits))
    print(p.stem,[(c,len(s)) for c,s in chains.items()],flush=True)
out=dict(status='passed' if not any(r['exact_train_matches'] or r['local_30pct_70pct_both_training_hits'] for r in rows) else 'overlap_detected',
    actual_raw_template_files=len(paths),full_chains=sum(r['scope']=='full_chain' for r in rows),domain_screens=sum(r['scope']!='full_chain' for r in rows),
    methods='Read-only first-model amino-acid PDB sequence; canonical MSE->Met. Local alignment match2/mismatch-1/gap-open-3/gap-extension-0.5; integer-span coverage at least0.70 both chains, identity including all gap columns at least0.30. Domain screens additionally test first120 TCR/first180 MHC residues. Roles are reference chain mappings and must be checked in final prepared structure manifests.',
    final_prepared_structures_pending=True,rows=rows)
(ROOT/'raw_template_overlap_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k!='rows'},ensure_ascii=False))
