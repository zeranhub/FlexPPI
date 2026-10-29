"""Independent scope/value audit against latest float64 source exports."""
from pathlib import Path
import json,sys,hashlib,importlib
ROOT=Path(__file__).resolve().parent;ATLAS=ROOT.parent;EXP=ATLAS.parent
sys.path.insert(0,str(EXP/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
import audit_cross_case_figures as prior

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csv(p):return pd.read_csv(p,float_precision='round_trip',dtype={'id_a':str,'id_b':str})
def close(a,b,where,tol=1e-13):
    aa=np.asarray(a,float);bb=np.asarray(b,float)
    assert aa.shape==bb.shape,(where,aa.shape,bb.shape)
    assert np.allclose(aa,bb,atol=tol,rtol=0,equal_nan=True),(where,np.nanmax(abs(aa-bb)))
    return float(np.nanmax(abs(aa-bb))) if np.isfinite(aa-bb).any() else 0
def compare_frame(actual,expected,where):
    keys=['role','alignment_column']
    assert list(map(tuple,actual[keys].to_numpy()))==list(map(tuple,expected[keys].to_numpy())),(where,'scope keys')
    n=0;maxd=0
    for metric in ['rmsf','hchi']:
        for suffix in ['full_a','full_b','md_a','md_b','var_a','var_b','full_delta','md_delta','delta_sem']:
            k=metric+'_'+suffix
            maxd=max(maxd,close(actual[k],expected[k],(where,k)));n+=len(actual)
    for k in ['same_aa','paired','gap_a','gap_b','rmsf_valid','hchi_valid']:
        assert np.array_equal(actual[k],expected[k]),(where,k)
    return n,maxd
CAT=json.loads((ATLAS/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
bycase={c['case_id']:c for c in CAT['cases']}
results=list(prior.results)

# Recheck generalized figure source subsets, plotted entropy artists, all MD points,
# and bind every artifact to the current hash of the complete alignment source.
for case in results:
    cid=case['case_id']
    if int(cid[1:])<18:continue
    c=bycase[cid];a,b=c['contrasts'][0];pdir=ROOT/cid/'figure'
    origpath=ROOT/cid/'data'/f'{a}_{b}_all_chain_alignment.csv';f=csv(origpath)
    p=csv(pdir/'panel_b_complete_aligned_peptide.csv');h=csv(pdir/'panel_c_complete_comparable_pMHC_interface.csv')
    n,err=compare_frame(p,f[f.role.eq('peptide')],(cid,'b'))
    nh,erh=compare_frame(h,f[f.role.isin(['mhc','peptide'])&f.interface_union&f.paired&f.hchi_valid&f.same_aa],(cid,'c'))
    audit=json.loads((pdir/f'{cid}_case_figure.audit.json').read_text(encoding='utf8'))
    artists=pd.DataFrame(audit['hchi_all_values'])
    close(artists['first'],h.hchi_full_a,(cid,'Hchi first artists'));close(artists['second'],h.hchi_full_b,(cid,'Hchi second artists'))
    for m in ['rmsf','hchi']:
        actual=csv(pdir/f'panel_d_{m}_all_contrasts.csv')
        expected=[]
        for pa,pb in c['contrasts']:
            raw=csv(ROOT/cid/'data'/f'{pa}_{pb}_all_chain_alignment.csv')
            mask=(raw.role.eq('peptide')&raw.paired&raw.rmsf_valid) if m=='rmsf' else (raw.role.isin(['mhc','peptide'])&raw.interface_union&raw.paired&raw.hchi_valid&raw.same_aa)
            expected.append(raw[mask])
        expected=pd.concat(expected,ignore_index=True)
        nn,dd=compare_frame(actual,expected,(cid,'d',m));n+=nn;err=max(err,dd)
    case.update(source_alignment_SHA256=sha(origpath),source_subset_SHA256={p.name:sha(p) for p in pdir.glob('panel_*.csv')},
                numeric_values_compared=n+nh,max_latest_source_difference=max(err,erh),latest_float64_data_match=True,
                generation_script_SHA256=sha(ROOT/'plot_general_cases.py'),all_gap_and_undefined_chi_rows_retained_in_master=True)

SCOPES={'C05':'interface','C06':'cdr','C07':'cdr','C08':'peptide','C09':'cdr','C10':'interface',
        'C11':'cdr','C12':'cdr','C13':'cdr','C14':'parallel','C15':'cdr','C16':'cdr','C17':'whole'}
def bselect(f,s):
    if s=='cdr':return f[f.region_a.fillna('').str.startswith('CDR')|f.region_b_projected.fillna('').str.startswith('CDR')]
    if s in ['interface','parallel']:return f[f.role.isin(['mhc','peptide'])&f.interface_union]
    if s=='peptide':return f[f.role.eq('peptide')]
    return f[f.paired&f.same_aa]
def cselect(f,s):
    if s=='whole':m=f.paired&f.same_aa
    elif s=='peptide':m=f.role.eq('peptide')&f.paired&f.same_aa
    else:m=f.role.isin(['mhc','peptide'])&f.interface_union&f.paired&f.same_aa
    return f[m&f.hchi_valid]
bioman=ATLAS/'figure_design/biological_series_manifest.json'
if bioman.exists():
    bio=json.loads(bioman.read_text(encoding='utf8'))
    for entry in bio['figures']:
        cid=entry['case_id'];c=bycase[cid];scope=SCOPES[cid];pdir=Path(entry['png']).parent
        audit=json.loads(Path(entry['audit']).read_text(encoding='utf8'))
        assert sha(entry['png'])==audit['PNG_SHA256'],(cid,'PNG hash')
        assert audit['biological_contrasts']==c['contrasts']
        for path,h in audit['source_CSV_SHA256'].items():assert sha(path)==h,(cid,'source CSV hash')
        for path,h in audit['source_audits'].items():assert sha(path)==h,(cid,'NPZ source audit hash')
        counts=[];nvalues=0;maxdiff=0
        for j,(a,b) in enumerate(c['contrasts']):
            raw=csv(ROOT/cid/'data'/f'{a}_{b}_all_chain_alignment.csv');br=bselect(raw,scope);hc=cselect(raw,scope)
            # Paired strip panels contain the predeclared valid-pair scope.
            if scope in ['interface','parallel']:br=br[br.paired&br.rmsf_valid]
            elif scope=='whole':br=br[br.paired&br.same_aa&br.rmsf_valid]
            bp=pdir/f'panel_b_scope_contrast{j+1}.csv'
            hp=pdir/f'panel_c_scope_contrast{j+1}.csv'
            assert bp.exists() and hp.exists(),(cid,'missing contrast',j)
            nv,err=compare_frame(csv(bp),br,(cid,'panel b',j));nvalues+=nv;maxdiff=max(maxdiff,err)
            nv,err=compare_frame(csv(hp),hc,(cid,'panel c',j));nvalues+=nv;maxdiff=max(maxdiff,err)
            counts.append(dict(a=a,b=b,valid_RMSF=int((br.paired&br.rmsf_valid).sum()),valid_same_AA_Hchi=len(hc)))
        assert audit['RMSF_valid_counts_by_contrast']==[x['valid_RMSF'] for x in counts]
        assert audit['Hchi_valid_counts_by_contrast']==[x['valid_same_AA_Hchi'] for x in counts]
        for m in ['rmsf','hchi']:
            actual=csv(pdir/f'panel_d_{m}_all_contrasts.csv');expected=[]
            for a,b in c['contrasts']:
                raw=csv(ROOT/cid/'data'/f'{a}_{b}_all_chain_alignment.csv')
                q=bselect(raw,scope) if m=='rmsf' else cselect(raw,scope)
                if m=='rmsf':q=q[q.paired&q.rmsf_valid]
                expected.append(q)
            nv,err=compare_frame(actual,pd.concat(expected,ignore_index=True),(cid,'panel d',m));nvalues+=nv;maxdiff=max(maxdiff,err)
        assert audit['n_MD_replicas_per_state']==(10 if cid=='C17' else 3)
        results.append(dict(case_id=cid,status='passed',PNG=entry['png'],PNG_SHA256=sha(entry['png']),scope=scope,
            complete_contrast_counts=counts,numeric_values_compared=nvalues,max_latest_source_difference=maxdiff,
            latest_float64_data_match=True,all_defined_contrasts_covered=True,n_MD_replicas_per_state=audit['n_MD_replicas_per_state'],
            source_alignment_SHA256=audit['source_CSV_SHA256'],plotted_source_SHA256={p.name:sha(p) for p in pdir.glob('panel_*.csv')}))

# Non-TCR complete-region figures: separate source protocols, distinct sampling units.
non=ATLAS/'evidence/non_tcr_figures/non_tcr_figure_manifest.json'
if non.exists():
    for entry in json.loads(non.read_text(encoding='utf8'))['figures']:
        cid=entry['caseid'];pdir=Path(entry['png']).parent;audit=json.loads(Path(entry['audit']).read_text(encoding='utf8'))
        for source,h in audit['source_sha256'].items():assert sha(source)==h,(cid,'source hash')
        for ext,h in audit['output_sha256'].items():assert sha(entry[ext])==h,(cid,'output hash')
        points=pd.read_csv(pdir/'all_plotted_points.csv',float_precision='round_trip');err=0;n=0
        if cid=='C31':
            src=pd.read_csv(next(k for k in audit['source_sha256'] if k.endswith('all_residues.csv')),float_precision='round_trip')
            for panel,kind,key in [('c','MD','md_percentile'),('c','model','model_percentile'),('e','MD','md_value'),('e','model','model_value')]:
                pp=points[points.panel.eq(panel)&points.kind.eq(kind)].sort_values(['state','residue'])
                ref=src.sort_values(['state','author_residue_number'])
                assert len(pp)==len(ref)==348
                assert list(zip(pp.state,pp.residue))==list(zip(ref.state,ref.author_residue_number))
                err=max(err,close(pp.y,ref[key],(cid,panel,kind)));n+=len(pp)
            summary=pd.read_csv(next(k for k in audit['source_sha256'] if k.endswith('case_metrics.csv')),float_precision='round_trip').set_index('state')
            for row in points[points.panel.eq('d')].itertuples():
                key='l12_md_mean_percentile' if 'MD' in row.kind else 'l12_model_mean_percentile'
                err=max(err,close([row.y],[summary.loc[row.state,key]],(cid,'d')));n+=1
            scope_summary=dict(states=4,complete_KIX_positions_per_state=87,L12_positions=8,MD_independent_trajectories_per_state=1)
        else:
            src=pd.read_csv(next(k for k in audit['source_sha256'] if k.endswith('all66_residue_summary.csv')),float_precision='round_trip')
            for kind,key in [('MD','frozen_MD_joint_Hchi_nats'),('model','predicted_Hchi_nats')]:
                pp=points[points.panel.eq('b')&points.kind.eq(kind)].sort_values('residue_id');ref=src.sort_values('author_residue_id')
                assert len(pp)==len(ref)==66;err=max(err,close(pp.y,ref[key],(cid,'b',kind)));n+=66
            angle=pd.read_csv(next(k for k in audit['source_sha256'] if k.endswith('defined_angles_long.csv')),float_precision='round_trip')
            pp=points[points.panel.eq('c')].sort_values(['residue_id','chi']);ref=angle.sort_values(['author_residue_id','chi_index'])
            assert len(pp)==len(ref)==151
            err=max(err,close(pp.x,ref.predicted_Hchi_nats_AA_centered,(cid,'c x')))
            err=max(err,close(pp.y,ref.fine_entropy_four_origins_mean_nats_AA_centered,(cid,'c y')));n+=302
            hist=pd.read_csv(next(k for k in audit['source_sha256'] if k.endswith('histograms_501frames.csv')),float_precision='round_trip')
            hist=hist[hist.grid.eq('primary_501frames')&hist.window.eq('pooled')]
            pp=points[points.panel.eq('d')]
            assert len(pp)==len(hist)==432
            for site in ['Q2','Q41']:
                for chi in [1,2,3]:
                    p=pp[(pp.residue_id==int(site[1:]))&(pp.chi==chi)].sort_values('x')
                    h=hist[hist.residue.eq(site)&hist.chi.eq(f'chi{chi}')].sort_values('angle_center_degrees')
                    err=max(err,close(p.x,h.angle_center_degrees,(cid,'d x',site,chi)))
                    err=max(err,close(p.y,h.F_star,(cid,'d y',site,chi)));n+=2*len(p)
            scope_summary=dict(valid_Hchi_positions=66,chi_angle_pairs=151,Fstar_histogram_bins=432,MD_independent_trajectories=1)
        results.append(dict(case_id=cid,status='passed',PNG=entry['png'],PNG_SHA256=sha(entry['png']),numeric_values_compared=n,
                            max_latest_source_difference=err,scope_summary=scope_summary,source_SHA256=audit['source_sha256']))

ids=[r['case_id'] for r in results]
assert len(ids)==len(set(ids))
c14=next((r for r in results if r['case_id']=='C14'),None)
if c14:
    assert [(r['a'],r['b'],r['valid_same_AA_Hchi']) for r in c14['complete_contrast_counts']]==[('3qdg','3qdj',16),('3qdm','3qeq',19)]
report=dict(status='passed' if len(ids)==32 else 'passed_for_completed_figures',review='Independent scope and original-number review, current artifact hashes bound',
    cases_reviewed=len(ids),case_ids=sorted(ids),missing_case_ids=[f'C{i:02}' for i in range(1,33) if f'C{i:02}' not in ids],
    cases=sorted(results,key=lambda r:r['case_id']),
    critical_scope_notes=['C14 DMF5 3QDG→3QDJ has 16 same-AA interface Hchi sites; DMF4 3QDM→3QEQ has 19. These distinct results are not interchangeable.',
        'Changed-aa sites are not pure environment evidence. Unpaired residues and undefined chi remain in complete data; no zero filling.',
        'ITrace MD has3 independent trajectories per state; CatchBond10; KIX and ubiquitin have one available trajectory, no repeat confidence inferred.',
        'All model-primary plots show existing model output; MD comparisons provide corroboration and do not establish causation or affinity.'],
    source_data_manifest_SHA256=sha(ROOT/'data_manifest.json'),general_figure_manifest_SHA256=sha(ROOT/'general_figure_manifest.json'))
finalmap=ATLAS/'ATLAS_FIGURE_MAP.json'
if finalmap.exists():
    bundle=json.loads(finalmap.read_text(encoding='utf8'))
    assert set(bundle['figures'])==set(ids)
    for r in report['cases']:
        bound=bundle['figures'][r['case_id']]
        assert r['PNG_SHA256']==bound['sha256']['png']==sha(bound['png']),(r['case_id'],'Final-map PNG differs')
        for ext in ['pdf','svg','caption']:
            assert sha(bound[ext])==bound['sha256'][ext],(r['case_id'],'Final-map source file changed',ext)
        r['final_map_PNG_source']=bound['png']
        r['final_map_source_SHA256']=bound['sha256']
    report['final_figure_map_SHA256']=sha(finalmap)
    report['all_32_final_map_artifacts_bound']=True
out=ATLAS/'catalog/final_scope_and_numbers_review.json';out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
# Replace the earlier partial handoff report with the consolidated final result.
(ROOT/'cross_case_figure_numeric_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:report[k] for k in ['status','cases_reviewed','missing_case_ids']},ensure_ascii=False,indent=2))
