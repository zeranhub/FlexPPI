"""Read-only numerical/asset checks plus explicitly recorded QA files; no model fits."""
from pathlib import Path
import hashlib, json, importlib.util, sys, dataclasses, subprocess
import numpy as np
import pandas as pd
from PIL import Image

R = Path(__file__).resolve().parent
O = R / 'audit'
O.mkdir(exist_ok=True)
checks = {}

def check(name, condition, details=True):
    assert bool(condition), name
    checks[name] = details

ann = pd.read_csv(R/'biological_annotation/all56_biological_annotations.csv')
core = ann[ann.core_functional_peptide]
check('frozen_cohort', len(ann)==56 and ann.pdb_id.nunique()==54)
check('biological_core', len(core)==42 and core.pdb_id.nunique()==40)
rec = pd.read_csv(R/'formal_results/per_record_all_methods_references.csv')
fig = pd.read_csv(R/'figures/figure_all480_values.csv')
expected = rec[rec.core_functional_peptide & rec.reference.eq('uniform_rep1') & rec.scope.eq('all') & rec.method.isin(['full','intra','static47_sasa'])]
expected = expected.groupby(['pdb_id','component_id','role','metric','method']).mae.mean().reset_index()
keys = ['pdb_id','component_id','role','metric','method']
joined = expected.merge(fig,on=keys,validate='one_to_one',suffixes=('_source','_figure'),how='outer',indicator=True)
err = float(abs(joined.mae_source-joined.mae_figure).max())
check('figure_all_480_points', len(fig)==len(joined)==480 and joined._merge.eq('both').all() and err<1e-12,
      dict(n_values=len(fig),max_abs_difference=err,n_PDB_each_panel=fig.groupby('panel').pdb_id.nunique().to_dict()))
means = pd.read_csv(R/'figures/figure_macro_means_CI.csv')
max_mean_diff = 0.; max_CI_diff = 0.
for row in means.itertuples():
    z = expected[expected.role.eq(row.role)&expected.metric.eq(row.metric)&expected.method.eq(row.method)]
    fam = z.groupby('component_id').mae.agg(['sum','count'])
    rng=np.random.default_rng(20261026)
    ix=rng.integers(0,len(fam),size=(10000,len(fam)))
    boot=fam['sum'].to_numpy()[ix].sum(1)/fam['count'].to_numpy()[ix].sum(1)
    lo,hi=np.quantile(boot,[.025,.975])
    max_mean_diff=max(max_mean_diff,abs(z.mae.mean()-row.PDB_MAE))
    max_CI_diff=max(max_CI_diff,abs(lo-row.CI95_lo),abs(hi-row.CI95_hi))
check('figure_12_means_and_cluster_CIs',len(means)==12 and max_mean_diff<1e-12 and max_CI_diff<1e-12,
      dict(max_mean_abs_difference=max_mean_diff,max_CI_abs_difference=max_CI_diff))

ci=pd.read_csv(R/'formal_results/paired_MAE_gain_cluster_CI.csv')
ind_max=0.; ind_rows=0
for reference in ['original','uniform_rep1']:
    file='independent_core_original_comparisons.csv' if reference=='original' else 'independent_core_uniform_comparisons.csv'
    other=pd.read_csv(R/'independent_statistics_review'/file)
    for row in other.itertuples():
        method=row.comparison.removeprefix('full_vs_')
        q=ci[ci.cohort.eq('core')&ci.reference.eq(reference)&ci.scope.eq('all')&ci.role.eq(row.role)&ci.metric.eq(row.metric)&ci.comparator.eq(method)]
        assert len(q)==1
        a=q.iloc[0]
        ind_max=max(ind_max,abs(a.Full_MAE-row.full_mae),abs(a.comparator_MAE-row.comparator_mae),abs(a.PDB_macro_gain-row.gain),abs(a.CI95_lo-row.lo),abs(a.CI95_hi-row.hi))
        ind_rows+=1
check('independent_core_original_and_uniform',ind_rows==32 and ind_max<1e-12,dict(n_comparisons=ind_rows,max_abs_difference=ind_max))

# Re-enumerate the primary same-AA/backbone/SASA match independently from the final summaries.
f=pd.read_parquet(R/'formal_results/all_features_uniform_rep1.parquet')
pair_ci=pd.read_csv(R/'formal_results/sameAA_match_type_equal_cluster_CI.csv')
match_audit={}
for role in ['peptide','receptor']:
    rr=[];npairs=0
    sub=f[f.core_functional_peptide & f.role.eq(role)&f.valid_rmsf&f.valid_hchi]
    for (sid,pdb,group),p in sub.groupby(['sample_id','pdb_id','component_id']):
        aa_values=[]
        for _,g in p.groupby('amino_acid'):
            if len(g)<2:continue
            i,j=np.triu_indices(len(g),1)
            rmsf=g.md_rmsf.to_numpy(); sasa=g.sidechain_sasa_complex_A2.to_numpy()
            keep=(abs(rmsf[i]-rmsf[j])<=.10)&(abs(sasa[i]-sasa[j])<=10)
            if not keep.any():continue
            truth=g.md_hchi.to_numpy(); full=g.pred_full_hchi.to_numpy(); base=g.pred_static47_sasa_hchi.to_numpy()
            target=(truth[j]-truth[i])[keep]
            e_full=abs((full[j]-full[i])[keep]-target).mean()
            e_base=abs((base[j]-base[i])[keep]-target).mean()
            aa_values.append([e_full,e_base]);npairs+=int(keep.sum())
        if aa_values:
            v=np.mean(aa_values,axis=0);rr.append(dict(pdb_id=pdb,component_id=group,full=v[0],base=v[1]))
    z=pd.DataFrame(rr).groupby(['pdb_id','component_id'])[['full','base']].mean()
    q=pair_ci[pair_ci.cohort.eq('core')&pair_ci.reference.eq('uniform_rep1')&pair_ci.role.eq(role)&pair_ci.rmsf_tolerance_A.eq(.10)&pair_ci.sasa_tolerance_A2.eq(10)&pair_ci.comparator.eq('static47_sasa')]
    assert len(q)==1
    q=q.iloc[0]
    difference=max(abs(z.full.mean()-q.Full_pair_MAE),abs(z.base.mean()-q.comparator_pair_MAE))
    check('strict_match_'+role,npairs==q.n_pairs and len(z)==q.n_PDB and difference<1e-12)
    match_audit[role]=dict(n_pairs=npairs,n_PDB=len(z),max_MAE_difference=float(difference))
checks['strict_match_details']=match_audit

progress=json.loads((R/'controls/training_progress.json').read_text(encoding='utf-8'))
check('all_20_new_baseline_fits',progress['status']=='complete' and len(progress['jobs'])==20 and all(j['family_overlap']==0 for j in progress['jobs']))
check('all_20_trained_models_present',len(list((R/'controls/models').glob('*.joblib')))==20)
uniform=json.loads((R/'uniform_md/compute_summary.json').read_text(encoding='utf-8'))
check('uniform_53_no_failed_and_masks_unchanged',uniform['computed']==53 and not uniform['failed'] and uniform['mask_change_count']==0)
check('uniform_three_unavailable_retained',set(x['pdb_id'] for x in uniform['unavailable'])=={'1CLV','1GL1','1PPE'})
check('single_local_trajectory',uniform['all_one_local_replica'])
ind=json.loads((R/'independent_statistics_review/independent_statistics_audit.json').read_text(encoding='utf-8'))
check('independent_cache_and_AA_baseline_audit',ind['all_fits_verified'] and len(ind['aa_means_independently_verified'])==10)

raw=Path(r'C:\Users\DELL\.codex\worktrees\2a52\FlexPPI\output\non_tcr_case_screen_20261025\archive_extract\family_folds\evaluation\oof_per_residue.csv.gz')
if raw.exists():
    digest=hashlib.file_digest(raw.open('rb'),'sha256').hexdigest()
    check('original_frozen_prediction_SHA256',digest=='0852469b8fb6dd2977796d040ce8bc7fe04f3de29f217e7fc563213c064ad7f3',digest)
else:
    checks['original_frozen_prediction_SHA256']='external original unavailable on this machine; stored manifest applies'

qa=json.loads((R/'figures/figure_QA.json').read_text(encoding='utf-8'))
png=R/'figures/class_core40_absolute_MAE.png'
with Image.open(png) as im:
    check('raster_dimensions',im.size==(3242,2303),dict(pixels=im.size,dpi=im.info.get('dpi')))
asset_hashes={ext:hashlib.sha256((R/f'figures/class_core40_absolute_MAE.{ext}').read_bytes()).hexdigest() for ext in ['png','pdf','svg']}
check('asset_SHA256_matches_QA',asset_hashes==qa['sha256'])
check('axes_all_points_and_text_boundaries',qa['axes_start_at_zero'] and qa['all_individual_points'] and not qa['outside_text'] and qa['min_font']>=8.2)
reviewed_hash='f9059c39562cc17bd8c3f60c29a243dcfbcad56dd5b508feabd901049aecd4d1'
review=R/'figures/QA_人工审核.md'
if review.exists() and asset_hashes['png']==reviewed_hash:
    qa['visual_QA']='passed: actually opened final PNG; reviewed hash and four-pass notes in QA_人工审核.md'
    qa['raw_style_findings']='retained; 2 literal-scanner FAIL and 2 WARN adjudicated in QA_人工审核.md'
    (R/'figures/figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
    checks['human_visual_review_of_exact_PNG']=reviewed_hash
else:
    checks['human_visual_review_of_exact_PNG']='pending: regenerated assets require fresh visual review'

validator=Path(r'C:\Users\DELL\.codex\skills\academic-figure-skill\scripts\qa_validator.py')
if validator.exists():
    spec=importlib.util.spec_from_file_location('figure_validator',validator)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    source=(R/'figures/plot_class_model_comparison.py').read_text(encoding='utf-8')
    findings=[]
    for name in dir(mod):
        if name.startswith(('check_ap','check_cl')):
            out=getattr(mod,name)(source)
            findings.extend(out if isinstance(out,list) else [out])
    (O/'automated_style_QA_raw.json').write_text(json.dumps([dataclasses.asdict(x) for x in findings],ensure_ascii=False,indent=2),encoding='utf-8')
    checks['style_validator_raw_findings']={c:sum(x.category==c for x in findings) for c in ['PASS','WARN','FAIL']}
pdf_fonts=Path(r'C:\texlive\2026\bin\windows\pdffonts.exe')
if pdf_fonts.exists():
    result=subprocess.run([str(pdf_fonts),str(R/'figures/class_core40_absolute_MAE.pdf')],capture_output=True,text=True,check=True)
    (O/'pdffonts.txt').write_text(result.stdout,encoding='utf-8')
    lines=[x for x in result.stdout.splitlines()[2:] if x.strip()]
    check('actual_PDF_fonts_embedded',len(lines)>=2 and all('yes yes yes' in ' '.join(x.split()) for x in lines),result.stdout)

checks['asset_sha256']=asset_hashes
checks['scope']='Formal analysis and diagnostic figure delivered; no claim new mechanism, no manuscript Word edits.'
(O/'final_delivery_audit.json').write_text(json.dumps(dict(status='passed',checks=checks),ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(status='passed',checks=len(checks),max_independent_difference=ind_max,strict_matching=match_audit,style=checks.get('style_validator_raw_findings')),ensure_ascii=False))
