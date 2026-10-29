"""Numerical verification and descriptive summaries, no hypothesis p-values."""
from pathlib import Path
import json
import re
import hashlib
import csv
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/mechanistic_case_studies_20260913'
OUT=ROOT/'output/flexppi_mechanistic_case_studies_20260913'

def bmrb_mapping(specs):
    """Map NMR-STAR author numbering via deposited sequence, never fixed offsets."""
    text=(BASE/'tcr/sources/bmr27441_3.str').read_text()
    seqs=re.findall(r'_Entity\.Polymer_seq_one_letter_code\s*\n;\s*\n(.*?)\n;',text,re.S)
    seq=''.join(seqs[0].split())
    assert len(seq)==245
    block=text.split('_Atom_chem_shift.ID',1)[1].split('stop_',1)[0]
    cols=['ID'];rows=[]
    for line in block.splitlines():
        stripped=line.strip()
        if stripped.startswith('_Atom_chem_shift.'):
            cols.append(stripped.split('.')[-1])
        elif stripped and stripped[0].isdigit():
            vals=stripped.split()
            assert len(vals)==len(cols)
            rows.append(dict(zip(cols,vals)))
    chem=pd.DataFrame(rows)
    chem=chem[chem.Entity_ID=='1']
    authors={int(r.Auth_seq_ID):int(r.Seq_ID) for r in chem.itertuples() if r.Auth_seq_ID not in ['.','?']}
    st=specs['tcr_bound_wtseq']['residues']['E'];pseq=st['sequence'];ids=st['author_ids']
    al=PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,open_gap_score=-3,extend_gap_score=-.5)
    a=al.align(seq,pseq)[0];mapped={}
    for (b,e),(bb,ee) in zip(*a.aligned):
        for i,j in zip(range(b,e),range(bb,ee)):
            if seq[i]==pseq[j]:mapped[i+1]=ids[j]
    # Selected shared sites explicitly named in He2020, NOT the full positive mask.
    selected={'K119':'intensity','Q226':'intensity','F129':'perturbed','S132':'perturbed','E135':'perturbed','T143':'perturbed','R196':'perturbed'}
    result=[]
    for site,readout in selected.items():
        num=int(site[1:]);sid=authors.get(num)
        if sid is None or seq[sid-1]!=site[0]:
            result.append(dict(literature_site=site,status='not verified in deposited author numbering'))
            continue
        result.append(dict(literature_site=site,status='sequence and author numbering verified',bmrb_sequence_index=sid,
            pdb_chain='E',pdb_author_id=mapped.get(sid),experimental_readout=readout,
            reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC8356565/',
            limitation='Selected textual examples only; no negative labels or enrichment test.'))
    (BASE/'tcr/selected_nmr_site_mapping.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result

def main():
    d=pd.read_csv(BASE/'residue_predictions.csv',dtype={'author_id':str,'before_author':str})
    r=pd.read_csv(BASE/'region_predictions.csv')
    states=json.loads((BASE/'states.json').read_text());specs={s['id']:s for s in states}
    tests=[]
    for s in states:
        assert hashlib.sha256(Path(s['path']).read_bytes()).hexdigest()==s['sha256']
        for f in range(5):
            with np.load(BASE/'prediction_cache'/f"{s['id']}_fold{f}.npz") as z:
                pv=json.loads(str(z['provenance']))
                assert pv['input']==s['sha256']
                assert z['prediction'].shape==(sum(c['n'] for c in s['residues'].values()),2)
    tests.append('102 inputs and510 caches: shapes and SHA256 consistent')
    assert not d.duplicated(['contrast','fold','chain','author_id']).any()
    assert d.groupby('contrast').fold.nunique().eq(5).all()
    assert np.allclose(d.delta_rmsf,d.after_rmsf-d.before_rmsf,atol=1e-6,equal_nan=True)
    assert np.allclose(d.delta_chi_entropy,d.after_chi_entropy-d.before_chi_entropy,atol=1e-6,equal_nan=True)
    tests.append('Unique residue keys; fivefold coverage; after-minus-before arithmetic')
    k=['fold','chain','author_id']
    total=d[d.contrast=='tcr_apo_holo'].set_index(k)
    geo=d[d.contrast=='tcr_geometry'].set_index(k)
    context=d[d.contrast=='tcr_context'].set_index(k)
    common=total.index.intersection(geo.index).intersection(context.index)
    residual=total.loc[common,'delta_rmsf']-geo.loc[common,'delta_rmsf']-context.loc[common,'delta_rmsf']
    assert abs(residual).max()<1e-6
    tests.append('TCR decomposition on exact common residue keys: max residual '+str(abs(residual).max()))
    # Altered residue excluded from mutation comparison by exact identity.
    assert not ((d.case=='bap1')&(d.author_id=='49')).any()
    tests.append('BAP1 mutant residue49 excluded from delta comparisons')
    assert all('60I' in s['residues']['M']['author_ids'] for s in states if s['case']=='thrombin')
    for s in states:
        if s['case']=='thrombin':
            q=s['residues']['M'];assert q['sequence'][q['author_ids'].index('60I')]=='T'
    tests.append('Thrombin author60I reverted correctly, insertion code preserved')
    rr=r[(r.case=='kix')&(r.region=='L12 loop')&~r.contrast.str.contains('full31')].copy()
    rr['partner']=rr.contrast.str.extract(r'm\d+_(\w+)_effect')[0].values
    summaries=[]
    for metric in ['rmsf','chi_entropy']:
        w=rr[rr.metric==metric].pivot(index=['source_model','fold'],columns='partner',values='delta')
        assert len(w)==100
        diff=w.mll-w.myb
        byfold=diff.groupby('fold').mean()
        summaries.append(dict(metric=metric,test='L12 MLL minus cMyb',mean=float(diff.mean()),checkpoint_sd=float(byfold.std()),
            checkpoint_means=byfold.to_dict(),n_conformers=20,n_checkpoints=5,n_negative_pairs=int((diff<0).sum()),
            n_negative_conformer_means=int((diff.groupby('source_model').mean()<0).sum()),
            n_negative_checkpoint_means=int((byfold<0).sum())))
    tests.append('KIX includes all20 conformers ×5checkpoints for each principal contrast')
    numerical=r.groupby(['case','contrast','region','metric']).agg(mean=('delta','mean'),checkpoint_sd=('delta','std'),residues=('n','min')).reset_index()
    numerical.to_csv(OUT/'quantitative_summary.csv',index=False)
    sites=bmrb_mapping(specs)
    report=dict(tests=tests,findings=summaries,nmr_site_mapping=sites,
                statistical_scope='Descriptive only; no independence assumption across folds/conformers, no p-values.',
                exclusions='Hchi region summaries exclude A/G/P; mutated sites excluded; unavailable atoms or labels never imputed.')
    (OUT/'numerical_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=True,indent=2))

if __name__=='__main__':
    main()
