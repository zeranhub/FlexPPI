from pathlib import Path
import pandas as pd, numpy as np, json, hashlib
from scipy.stats import spearmanr
R=Path(__file__).resolve().parent.parent
a=pd.read_csv(R/"all_contrasts_complete_per_position.csv")
s=pd.read_csv(R/"all_contrasts_predefined_regions.csv")
m=pd.read_csv(R/"exact_input_predictions_and_MD_labels.csv")
reported=pd.read_csv(R/"MD_absolute_prediction_metrics.csv")
nmr=pd.read_csv(R/"all_named_NMR_sites_prediction_check.csv")
st=pd.read_csv(R/"all_states_per_residue_predictions.csv")
checks=[]
def selected_rows(contrast,metric,region=None,chain=None,same_type=False):
 x=a[(a.contrast==contrast)&(a.variant=="full")].copy()
 if region is not None:x=x[x.imgt_region_after==region]
 if chain is not None:x=x[x.chain==chain]
 if same_type:x=x[x.same_aa]
 return x[x["valid_"+metric]]
for contrast in ["flex_joint_state_change_full","flex_joint_state_change_fv","flex_crystal_copy_sensitivity"]:
 for metric in ["rmsf","hchi"]:
  for scope in ["H_CDR1","antigen_loop_73_78"]:
   x=a[(a.contrast==contrast)&(a.variant=="full")].copy()
   if scope=="H_CDR1":x=x[(x.chain=="H")&(x.imgt_region_after=="CDR1")]
   else:x=x[(x.chain=="A")&(x.author_number_after.between(73,78))]
   x=x[x["valid_"+metric]]
   v=x["delta_"+metric]
   y=s[(s.contrast==contrast)&(s.variant=="full")&(s.scope==scope)&(s.metric==metric)]
   assert len(y)==1,(contrast,metric,scope)
   err=abs(v.mean()-float(y.iloc[0].mean_delta))
   assert err<1e-12,(contrast,metric,scope,err,len(x),x.imgt_region_after.unique().tolist())
   checks.append(dict(contrast=contrast,scope=scope,metric=metric,n=len(x),mean_delta=float(v.mean()),mean_abs_delta=float(v.abs().mean()),n_increase=int((v>0).sum()),n_decrease=int((v<0).sum()),aggregation_difference=err))
md=[]
for sample in m.sample_id.unique():
 x=m[(m.sample_id==sample)&(m.variant=="full")]
 for metric,pc,tc,vc in [("rmsf","rmsf_pred","rmsf_MD","rmsf_valid"),("hchi","hchi","chi_entropy","chi_entropy_valid")]:
  z=x[x[vc]&x.geometry_valid_pred&x.model_mask].copy()
  if metric=="hchi":z=z[z.hchi_chemical_valid]
  err=float(np.abs(z[pc]-z[tc]).mean());rho=float(spearmanr(z[pc],z[tc]).statistic)
  y=reported[(reported.sample_id==sample)&(reported.variant=="full")&(reported.scope=="all")&(reported.metric==metric)].iloc[0]
  assert len(z)==int(y.n) and abs(err-y.mae)<1e-12 and abs(rho-y.spearman)<1e-12,(sample,metric,len(z),y.n,err,y.mae,rho,y.spearman)
  md.append(dict(sample_id=sample,metric=metric,n=len(z),mae=err,spearman=rho))
ns=[]
for c in nmr.model_contrast.unique():
 x=nmr[(nmr.model_contrast==c)&(nmr.variant=="full")&nmr.exact_reference_identity_eligible]
 ns.append(dict(contrast=c,eligible=len(x),mean_delta_rmsf=float(x.delta_rmsf.mean()),rmsf_decrease=int((x.delta_rmsf<0).sum()),rmsf_increase=int((x.delta_rmsf>0).sum()),scope="Free millisecond NMR named sites; qualitative only, not a numeric MD target"))
summary=dict(status="passed",aggregation_checks=checks,MD_join_independent_recalculation=md,NMR_named_site_summary=ns,limits=["No inference rerun. No changed source CSVs.","Family audit pending; no strict external-family validation claim.","Fixed-backbone mutations are sensitivity inputs, not observed mutant structures or measured dynamic effects.","Published FlexHEL MD concerns free antigen; observed complex rigid/flex comparisons change both antibody and antigen.","Crystal-copy and full/Fv stability are sensitivity checks, not biological-repeat uncertainty.","Hchi differs with side-chain chemistry; compare only same-type chemically valid sites. Crystal rotamer differences are not Hchi ground truth."])
summary["source_sha256"]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [R/"all_contrasts_complete_per_position.csv",R/"all_contrasts_predefined_regions.csv",R/"MD_absolute_prediction_metrics.csv",R/"exact_input_predictions_and_MD_labels.csv",R/"all_named_NMR_sites_prediction_check.csv"]}
(R/"reference/independent_prediction_review.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(summary,ensure_ascii=False,indent=2))
