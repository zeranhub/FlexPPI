from pathlib import Path
import json,csv
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
CORRECTED=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
chains=json.loads((CORRECTED/'sequence_audit/chain_inventory.json').read_text())
itrace=json.loads((EXP/'tcr_case_v40_20261017/local_pair_screen/sequence_inventory.json').read_text())
candidates=['A0053','A0055','A000B','A00PD','1k5n_A','A002W','A003I','A0061','A00PP']
rows=[]
for sample in candidates:
    records=[r for r in chains if r['sample_id']==sample]
    for peptide in [r for r in records if len(r['sequence'])<=25]:
        for pdb,obj in itrace.items():
            pep=obj.get('roles',{}).get('peptide')
            if pep and pep['sequence']==peptide['sequence']:
                mhc=obj['roles'].get('mhc',{})
                matching=[r for r in records if len(r['sequence'])>150]
                rows.append({'sample':sample,'pdb':records[0]['biological_system'],'free_peptide_chain':peptide['chain'],'peptide':peptide['sequence'],
                    'itrace_pdb':pdb,'itrace_peptide_chain':pep['chain'],'candidate_non_peptide_chains':[{k:r[k] for k in ['chain','sequence','length']} for r in matching],
                    'itrace_mhc_sequence':mhc.get('sequence'),
                    'exact_mhc_chain_matches':[r['chain'] for r in matching if r['sequence']==mhc.get('sequence')]})
(ROOT/'cross_source_same_peptide_candidates.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
print(json.dumps([{k:v for k,v in r.items() if k not in ['candidate_non_peptide_chains','itrace_mhc_sequence']} for r in rows],indent=2))
