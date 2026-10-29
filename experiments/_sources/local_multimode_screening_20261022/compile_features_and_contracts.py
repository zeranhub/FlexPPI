"""Compile full unchanged-residue profiles; no outcome-dependent region choice."""
from pathlib import Path
from datetime import datetime,timezone
import csv,json,hashlib,sys
import numpy as np
from Bio.PDB import PDBParser,ShrakeRupley

ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
def readcsv(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def csvout(p,rr):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)

def main():
    pred=ROOT/'predictions/full_pool';state=json.loads((pred/'progress.json').read_text())
    assert state['status']=='passed' and state['completed_this_invocation']+state['resumed_existing']==1211
    refpath=ROOT/'runtime/training_reference/same_AA_training_ECDF.npz'
    ref=np.load(refpath,allow_pickle=False)
    stats=json.loads((ROOT/'runtime/training_reference/summary.json').read_text())['same_AA_stats']
    labels=readcsv(ROOT/'antigen_labels/curated/curated_candidates.csv')
    original=json.loads((ROOT/'antigen_labels/curated/panels.json').read_text())['panels']
    static={r['candidate_id']:r for r in readcsv(ROOT/'structures/static_features.csv')}
    manifest=json.loads((ROOT/'structures/input_manifest.json').read_text());templates={t['template_id']:t for t in manifest['templates']}
    out=ROOT/'features';out.mkdir(exist_ok=True)
    contracts={k:[] for k in ['primary','top_decile','proteinmpnn002','three_dose_AUC']}
    site_rows=[];profile_audit=[]
    # Tien et al. (2013) theoretical maximal whole-residue ASA, Angstrom squared.
    maxasa=dict(zip('ARNDCQEGHILKMFPSTWYV',[129,274,195,193,167,225,223,104,224,197,201,236,224,240,159,155,172,285,263,174]))
    for panel in original:
        pid=panel['panel_id'];tid=pid+'__AF3';t=templates[tid]
        rr=[r for r in labels if r['panel_id']==pid];ids=[r['candidate_id'] for r in rr]
        p=np.load(pred/(pid+'__parent')/'prediction.npz',allow_pickle=False)
        rp=np.load(pred/(pid+'__TCR_removed_parent')/'prediction.npz',allow_pickle=False)
        seq=p['amino_acid'];chain=p['chain'];rid=p['author_residue_id'];valid=p['geometry_valid']
        U=p['environment_mask']&valid&(chain!='B');UH=U&p['hchi_chemically_defined']
        assert U.sum()>0 and UH.sum()>0
        lookup={(c,r):i for i,(c,r) in enumerate(zip(rp['chain'],rp['author_residue_id']))}
        # Same geometrical parent for contacts and SASA, no candidate label used.
        st=PDBParser(QUIET=True).get_structure(pid,t['input_pdb']);model=st[0]
        sr=ShrakeRupley(probe_radius=1.4,n_points=100);sr.compute(model,level='A')
        partners=[r for c in model if c.id!='B' for r in c.get_residues() if r.id[0]==' ']
        for j in np.flatnonzero(chain=='B'):
            aa=str(seq[j]);pos=int(rid[j]);res=model['B'][(' ',pos,' ')]
            sideasa=sum(a.sasa for a in res if a.name not in {'N','CA','C','O','OXT'} and a.element not in {'H','D'})
            contacts=0;xyz=np.array([a.coord for a in res if a.element not in {'H','D'}])
            for r in partners:
                yy=np.array([a.coord for a in r if a.element not in {'H','D'}])
                if np.min(np.sum((xyz[:,None]-yy[None,:])**2,axis=2))<=4.5**2:contacts+=1
            R=float(p['mean'][j,0]);H=float(p['mean'][j,1]);a0=stats[aa+'__rmsf']['RMS_scale']
            percentileR=float(np.searchsorted(ref[aa+'__rmsf'],R,side='right')/len(ref[aa+'__rmsf']))
            percentileH=.5 if aa in {'G','A'} else float(np.searchsorted(ref[aa+'__chi_entropy'],H,side='right')/len(ref[aa+'__chi_entropy']))
            k=lookup['B',str(pos)];assert rp['amino_acid'][k]==aa
            context=[abs(R-float(rp['mean'][k,0]))/a0]
            if aa not in {'G','A'}:context.append(abs(H-float(rp['mean'][k,1]))/stats[aa+'__chi_entropy']['RMS_scale'])
            site_rows.append(dict(panel_id=pid,position_key=f'B:{pos}',position=pos,parent_aa=aa,parent_R=R,parent_H=H,
                parent_R_percentile=percentileR,parent_H_percentile=percentileH,H_physically_defined=aa not in {'G','A'},
                parent_partner_context_response=float(np.mean(context)),static_contact_degree=contacts,
                static_relative_SASA=sideasa/maxasa[aa],sidechain_SASA_A2=sideasa,
                ASA_reference_definition='sidechain area divided by Tien2013 maximal WHOLE-residue ASA; exposure score, not conventional total-residue RSA',
                context_definition='Full complex minus computational TCR-chain removal; same peptide/HLA coordinates; not apo dynamics'))
        rpattern=[];hpattern=[];bindings=[]
        for c in ids:
            folder=pred/c;a=json.loads((folder/'candidate_audit.json').read_text());assert a['status']=='passed'
            for name,digest in a['output_SHA256'].items():assert sha(folder/name)==digest
            z=np.load(folder/'prediction.npz',allow_pickle=False)
            assert np.array_equal(z['chain'],chain) and np.array_equal(z['author_residue_id'],rid)
            assert np.array_equal(z['geometry_valid'],valid) and np.array_equal(z['environment_mask'],p['environment_mask'])
            assert np.all(z['amino_acid'][U]==seq[U]) and np.all(z['same_aa_environment_mask'][U])
            assert np.all(chain[U]!='B') and not np.any(z['changed_vs_parent'][U])
            rpattern.append(z['delta_mean_vs_parent'][U,0]);hpattern.append(z['delta_mean_vs_parent'][UH,1])
            bindings.append(dict(candidate_id=c,prediction_sha256=sha(folder/'prediction.npz'),parent_prediction_sha256=a['parent_prediction_sha256']))
        R=np.asarray(rpattern);H=np.asarray(hpattern)
        Rscale=np.array([stats[str(a)+'__rmsf']['RMS_scale'] for a in seq[U]])
        Hscale=np.array([stats[str(a)+'__chi_entropy']['RMS_scale'] for a in seq[UH]])
        profile_audit.append(dict(panel_id=pid,N=len(ids),R_dimensions=int(U.sum()),H_dimensions=int(UH.sum()),
            R_max_abs_delta=float(np.max(abs(R))),H_max_abs_delta=float(np.max(abs(H))),
            all_inputs_fixed_backbone=True,unchanged_chemical_identity_at_every_feature=True,
            no_spatial_means_or_selected_good_residues=True))
        common=dict(candidate_ids=np.asarray(ids),sequences=np.asarray([r['sequence'] for r in rr]),R_pattern=R,H_pattern=H,
            R_training_scale=Rscale,H_training_scale=Hscale,checkpoint_sha256=p['checkpoint_sha256'],
            unchanged_spectator_U_verified=np.asarray(True),all_library_variable_positions_excluded=np.asarray(True),
            zero_shot_score_names=np.asarray(['HERMES_fixed','ProteinMPNN']),zero_shot_scores_oriented_higher_better=np.asarray(True),
            template_sha256=np.asarray(t['input_sha256']),regions_sha256=np.asarray(sha(ROOT/'structures'/f'{pid}__regions.json')),
            training_reference_sha256=np.asarray(sha(refpath)),R_chain=chain[U],R_author_residue_id=rid[U],H_chain=chain[UH],H_author_residue_id=rid[UH])
        for setting in contracts:
            idx=np.flatnonzero([r['complete_three_dose']=='True' for r in rr]) if setting=='three_dose_AUC' else np.arange(len(rr))
            selected=[rr[j].copy() for j in idx]
            if setting=='three_dose_AUC':
                for r in selected:r['value_oriented_higher_better']=r['log_dose_response_area'];r['endpoint']='three_dose_log_concentration_normalized_response_area'
            fp=out/setting/f'{pid}.npz';fp.parent.mkdir(exist_ok=True)
            fields={key:(val[idx] if key in {'candidate_ids','sequences','R_pattern','H_pattern'} else val) for key,val in common.items()}
            mpnn='ProteinMPNN002' if setting=='proteinmpnn002' else 'ProteinMPNN'
            fields['static_values']=np.asarray([[float(static[ids[j]]['HERMES_fixed_delta']),float(static[ids[j]][mpnn+'_delta'])] for j in idx])
            fields['zero_shot_scores']=np.asarray([[float(static[ids[j]]['HERMES_fixed_raw']),float(static[ids[j]][mpnn+'_raw'])] for j in idx])
            np.savez_compressed(fp,**fields)
            lp=fp.with_suffix('.labels.csv');csvout(lp,selected)
            pp=dict(panel,labels_csv=str(lp),labels_sha256=sha(lp),features_npz=str(fp),features_sha256=sha(fp),
                site_labels_csv=str(lp),site_labels_sha256=sha(lp),site_scores_csv=str(out/'parent_site_scores.csv'))
            if setting in {'top_decile','three_dose_AUC'}:
                pp.update(hit_kind='relative_top_decile',declared_hit_threshold=None,
                    hit_definition='Relative top-decile retrieval within the source-defined finite cohort; not above-parent improvement')
            if setting=='three_dose_AUC':
                pp.update(endpoint='three_dose_log_concentration_normalized_response_area',units='normalized response area on equal log-dose intervals',
                    mode1_eligible=len(idx)==171,mode2_eligible=len(idx)>=60,mode3_eligible=len(idx)>=100,
                    complete_candidate_pool=True,cohort_scope='Only source-complete three-dose curves, defined without predictions; biased measurement-complete cohort sensitivity, not full1197-pool validation',
                    candidate_count=len(idx))
            contracts[setting].append(pp)
        dump(out/f'{pid}__source_bindings.json',dict(profile_audit=profile_audit[-1],prediction_bindings=bindings))
    csvout(out/'parent_site_scores.csv',site_rows)
    for setting,pp in contracts.items():
        for panel in pp:panel['site_scores_sha256']=sha(out/'parent_site_scores.csv')
        dump(out/f'{setting}__contracts.json',dict(setting=setting,panels=pp,
            labels_not_used_in_feature_definition=True,training_reference_sha256=sha(refpath)))
    dump(out/'feature_compilation_audit.json',dict(status='passed',panels=profile_audit,
        compiled_at_utc=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),all1197primary_variants=True,
        weights_not_updated=True,endpoint_values_never_used_for_regions_or_dimension_choice=True))
    print(json.dumps(dict(status='features_compiled',panels=profile_audit)))

if __name__=='__main__':main()
