"""Auditable paired inference analysis; no hypothesis tests on five checkpoints.

Every scientific contrast uses deposited apo/holo coordinates. Artificial
coverage/context perturbations are exported as sensitivity controls only.
Regions are fixed from the literature/geometry before predictions are read.
"""
from pathlib import Path
import json
import re
import hashlib
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'artifacts/tcr_real_states_20260915'
OUT = ROOT / 'output/tcr_real_states_20260915'


def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def read_atoms(path):
    records = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith('ATOM  '):
            continue
        key = (line[21], line[22:27].strip())
        records.setdefault(key, {})[line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return records


def canonical(case, chain):
    return ('E' if chain in 'ACE' else 'F') if case == 'B4.2.3' else chain


def number(author):
    return int(re.match(r'-?\d+', str(author)).group())


def fit_distances(x, y):
    x, y = np.array(x), np.array(y)
    u, _, vt = np.linalg.svd((x-x.mean(0)).T @ (y-y.mean(0)))
    d = np.eye(3); d[-1,-1] = np.linalg.det(u @ vt)
    aligned = (x-x.mean(0)) @ (u @ d @ vt) + y.mean(0)
    return np.linalg.norm(aligned-y, axis=1)


def fixed_regions(case, chain, author, pmhc, cd3):
    n = number(author)
    regions = ['All matched TCR']
    if pmhc <= 5:
        regions.append('pMHC interface')
    elif pmhc <= 12:
        regions.append('pMHC proximal')
    else:
        regions.append('pMHC distal')
    if case == 'B4.2.3':
        if chain == 'F' and 127 <= n <= 138:
            regions.append('H3-associated 127-138')
        if chain == 'F' and n in [133,134]:
            regions.extend(['N133/K134', 'N133' if n == 133 else 'K134'])
        if chain == 'F' and n in [127,130,133,134,138]:
            regions.append('Five literature sites')
        if chain == 'F' and pmhc > 12 and not 127 <= n <= 138:
            regions.append('Other distal beta')
    else:
        if chain == 'A' and n in [104,182]:
            regions.extend(['S104/V182', 'S104' if n == 104 else 'V182'])
        if chain == 'A' and (99 <= n <= 109 or 177 <= n <= 187):
            regions.append('Lock-site windows')
        if cd3 <= 5 and pmhc > 12:
            regions.append('Distal CD3 contact')
        if cd3 > 12 and pmhc > 12 and not (chain == 'A' and (99 <= n <= 109 or 177 <= n <= 187)):
            regions.append('Other distal noncontact')
    return regions


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    states_list = json.loads((BASE/'states.json').read_text())
    specs = {s['id']:s for s in states_list}
    cons = json.loads((BASE/'contrasts.json').read_text())
    prep = json.loads((BASE/'preparation_audit.json').read_text())
    # This definition is exported BEFORE any prediction file is opened.
    save_json(OUT/'region_definitions.json', {
        'selection': 'Literature-defined author-number windows and fixed holo heavy-atom distances; not prediction peak selection.',
        'B4.2.3': {'holo_chain_beta':'F', 'H3-associated 127-138':[127,138], 'individual_functional_sites':[133,134],
            'five_literature_sites':[127,130,133,134,138], 'window_caveat':'Predefined sequence window; not a DSSP claim that all twelve residues form a helix.',
            'source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5440810/'},
        '1G4-CD3': {'holo_chain_alpha':'A', 'lock_sites':[104,182], 'lock_site_windows':[[99,109],[177,187]],
            'window_caveat':'Each functional anchor plus/minus five author-numbered residues, no fabricated mutant.',
            'source':'https://www.nature.com/articles/s41467-025-66939-7'},
        'geometry': {'interface_max_heavy_distance_A':5, 'proximal_max_A':12, 'distal_min_A':12,
            'distance_reference':'Deposited full holo protein context, fixed for every contrast, even apo/apo controls.',
            'contact_count':'Number of other CA atoms within 8 A, all observed input protein chains.'},
        'regional_primary_mask':'B4: all four experimental TCR copies shared 421 positions; 1G4: all four experimental structures shared 366 positions.',
        'hchi_primary_filter':'Exclude A/G/P from regional Hchi aggregation; all predictions retained in residue-level CSV.',
        'uncertainty':'Mean and sample SD across five checkpoints; not biological replication and not a confidence interval.'})
    atoms = {s['id']:read_atoms(s['path']) for s in states_list}
    contact_counts = {}
    for sid, ats in atoms.items():
        keys = list(ats)
        xyz = np.array([ats[k]['CA'] for k in keys])
        counts = cKDTree(xyz).query_ball_point(xyz, 8, return_length=True)-1
        contact_counts[sid] = dict(zip(keys, counts))
    common = {}
    common['B4.2.3'] = {('E' if c=='alpha' else 'F',a) for c,ids in prep['common_ids']['b4'].items() for a in ids}
    common['1G4-CD3'] = {(c,a) for c,ids in prep['common_ids']['one_g4'].items() if c in 'AB' for a in ids}
    fixed, memberships = {}, []
    for case, sid, pmhc_ch, cd3_ch in [('B4.2.3','b4_5ivx_full_full','ABP',''),('1G4-CD3','1g4_9c3e_full_full','H','DEFGXY')]:
        ats = atoms[sid]
        pm = cKDTree(np.array([x for (c,r),a in ats.items() if c in pmhc_ch for x in a.values()]))
        cd = cKDTree(np.array([x for (c,r),a in ats.items() if c in cd3_ch for x in a.values()])) if cd3_ch else None
        for (ch,a), coordinates in ats.items():
            if ch not in specs[sid]['target_chains']:
                continue
            dp = float(pm.query(np.array(list(coordinates.values())))[0].min())
            dc = float(cd.query(np.array(list(coordinates.values())))[0].min()) if cd is not None else np.inf
            regs = fixed_regions(case, ch, a, dp, dc)
            fixed[(case,ch,a)] = (dp, dc, regs)
            memberships.append(dict(case=case,chain=ch,author_id=a,common_all=(ch,a) in common[case],pmhc_distance_A=dp,
                cd3_distance_A=dc if np.isfinite(dc) else None,regions=';'.join(regs)))
    pd.DataFrame(memberships).to_csv(OUT/'residue_region_membership.csv', index=False)
    cached = {}
    for s in states_list:
        for fold in range(5):
            p = BASE/'prediction_cache'/(s['id']+f'_fold{fold}.npz')
            with np.load(p, allow_pickle=False) as z:
                assert z['mask'].all() and z['prediction'].shape == (s['n_residues'],2)
                assert np.isfinite(z['prediction']).all()
                cache = {(str(c),str(a)):(str(aa),v) for c,a,aa,v in zip(z['chain'],z['author'],z['aa'],z['prediction'])}
                assert len(cache) == s['n_residues']
                cached[(s['id'],fold)] = cache
    rows, geometry = [], []
    for co in cons:
        before, after = co['reference_state'], co['comparison_state']
        valid_for_fit = [m for m in co['mapping'] if (canonical(co['case'],m['comparison_chain']),m['comparison_author_id']) in common[co['case']]]
        keys_a = [(m['reference_chain'],m['reference_author_id']) for m in valid_for_fit]
        keys_b = [(m['comparison_chain'],m['comparison_author_id']) for m in valid_for_fit]
        displacement = fit_distances([atoms[before][k]['CA'] for k in keys_a],[atoms[after][k]['CA'] for k in keys_b])
        static_map = dict(zip(keys_b,displacement))
        geometry.append(dict(contrast=co['id'],case=co['case'],kind=co['kind'],n_fit=len(displacement),
            static_CA_RMSD_A=float(np.sqrt(np.mean(displacement**2))),
            note='One whole-TCR rigid fit, common across all experimental states; NOT RMSF.'))
        for m in co['mapping']:
            ka=(m['reference_chain'],m['reference_author_id']);kb=(m['comparison_chain'],m['comparison_author_id'])
            cc=canonical(co['case'],kb[0]); ca=kb[1]
            dp,dc,regs = fixed.get((co['case'],cc,ca),(np.nan,np.nan,[]))
            for fold in range(5):
                aa,va=cached[(before,fold)][ka];ab,vb=cached[(after,fold)][kb]
                assert aa == ab == m['amino_acid'], 'Residue identity mismatch'
                rows.append(dict(case=co['case'],contrast=co['id'],kind=co['kind'],fold=fold,reference_state=before,comparison_state=after,
                    reference_chain=ka[0],reference_author_id=ka[1],comparison_chain=kb[0],chain=cc,author_id=ca,author_number=number(ca),amino_acid=aa,
                    common_all=(cc,ca) in common[co['case']],hchi_eligible=aa not in 'AGP',
                    reference_rmsf=float(va[0]),comparison_rmsf=float(vb[0]),delta_rmsf=float(vb[0]-va[0]),
                    reference_hchi=float(va[1]),comparison_hchi=float(vb[1]),delta_hchi=float(vb[1]-va[1]),
                    pmhc_distance_A=dp,cd3_distance_A=dc if np.isfinite(dc) else np.nan,
                    delta_CA_contacts=int(contact_counts[after][kb]-contact_counts[before][ka]),
                    static_CA_displacement_A=static_map.get(kb,np.nan),regions=';'.join(regs)))
    df=pd.DataFrame(rows)
    df.to_csv(OUT/'residue_predictions.csv', index=False)
    pd.DataFrame(geometry).to_csv(OUT/'static_geometry_controls.csv',index=False)
    region_rows=[]
    d=df[df.common_all].copy()
    d['region']=d.regions.str.split(';');d=d.explode('region')
    for metric in ['rmsf','hchi']:
        dm=d if metric=='rmsf' else d[d.hchi_eligible]
        for (case,co,kind,region,fold),g in dm.groupby(['case','contrast','kind','region','fold']):
            if not region:
                continue
            region_rows.append(dict(case=case,contrast=co,kind=kind,region=region,fold=int(fold),metric=metric,n_residues=len(g),
                mean_reference=float(g['reference_'+metric].mean()),mean_comparison=float(g['comparison_'+metric].mean()),
                mean_delta=float(g['delta_'+metric].mean()),median_delta=float(g['delta_'+metric].median()),
                mean_abs_delta=float(g['delta_'+metric].abs().mean())))
    rd=pd.DataFrame(region_rows)
    rd.to_csv(OUT/'region_fold_statistics.csv',index=False)
    summary=rd.groupby(['case','contrast','kind','region','metric']).agg(mean_delta=('mean_delta','mean'),sd_delta=('mean_delta','std'),
        minimum_fold_delta=('mean_delta','min'),maximum_fold_delta=('mean_delta','max'),n_residues=('n_residues','first'),
        n_folds=('fold','nunique'),mean_reference=('mean_reference','mean'),mean_comparison=('mean_comparison','mean')).reset_index()
    assert (summary.n_folds==5).all()
    summary['same_sign_all_five']=(summary.minimum_fold_delta>0)|(summary.maximum_fold_delta<0)
    summary.to_csv(OUT/'region_summary.csv',index=False)
    # Coverage sensitivity is assessed at exactly the same residue IDs; no paired
    # region-count differences can explain this diagnostic.
    sensitivity=[]
    bio=df[df.common_all & df.kind.isin(['real_apo_holo','detergent_apo_holo_environment_confound'])]
    for cid,full in bio[bio.contrast.str.endswith('_full')].groupby('contrast'):
        other=bio[bio.contrast==cid[:-5]+'_common']
        merged=full.merge(other,on=['chain','author_id','fold'],suffixes=('_full','_common'),validate='one_to_one')
        assert len(merged)==len(full)
        for metric in ['rmsf','hchi']:
            m=merged if metric=='rmsf' else merged[merged.hchi_eligible_full]
            r=m.groupby(['chain','author_id'])[[f'delta_{metric}_full',f'delta_{metric}_common']].mean()
            x,y=r.iloc[:,0].to_numpy(),r.iloc[:,1].to_numpy()
            sensitivity.append(dict(contrast=cid,metric=metric,n_residues=len(r),
                spearman_full_common=float(spearmanr(x,y).statistic),mean_absolute_coverage_effect=float(np.mean(abs(x-y))),
                mean_absolute_biological_delta=float(np.mean(abs(x))),same_sign_fraction=float(np.mean(x*y>0)),
                no_statistical_test=True))
    pd.DataFrame(sensitivity).to_csv(OUT/'coverage_robustness.csv',index=False)
    ca_rows=[]
    for cid,g in bio[bio.contrast.str.endswith('_full')].groupby('contrast'):
        r=g.groupby(['chain','author_id']).agg(delta_rmsf=('delta_rmsf','mean'),delta_hchi=('delta_hchi','mean'),
            delta_CA_contacts=('delta_CA_contacts','first'),hchi_eligible=('hchi_eligible','first'))
        for metric in ['rmsf','hchi']:
            q=r if metric=='rmsf' else r[r.hchi_eligible]
            ca_rows.append(dict(contrast=cid,metric=metric,spearman_delta_vs_CA_contacts=float(spearmanr(q['delta_'+metric],q.delta_CA_contacts).statistic),
                n_residues=len(q),note='Geometry association, not experimental accuracy or proof of improvement over baseline.'))
    pd.DataFrame(ca_rows).to_csv(OUT/'contact_baseline_association.csv',index=False)
    save_json(OUT/'analysis_audit.json',dict(n_states=len(specs),n_checkpoints=5,n_predictions=len(cached),n_contrasts=len(cons),
        n_residue_fold_contrast_rows=len(df),common_residues={k:len(v) for k,v in common.items()},
        no_imputation=True,no_new_training=True,no_mutant_structure_created=True,no_significance_tests=True,
        primary_difference='holo minus independently deposited apo; negative RMSF/Hchi means lower model-predicted flexibility.',
        units={'rmsf':'angstrom','hchi':'nats (joint discretized chi state entropy)'},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print(summary[(summary.contrast.str.endswith('_full')) & summary.region.isin(['N133','K134','H3-associated 127-138','S104','V182','Lock-site windows'])].to_string(index=False),flush=True)


if __name__=='__main__':
    main()
