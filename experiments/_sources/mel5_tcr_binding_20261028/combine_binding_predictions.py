from pathlib import Path
import csv,json,hashlib,math
root=Path(__file__).resolve().parent
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
panel=read(root/'binding_panel_inputs.csv')
ids=[r['case_id'] for r in panel]
ergo=read(root/'ergo/panel_results/all_ergo_scores.csv')
net=read(root/'nettcr/nettcr22_pan_predictions.csv')
method_order=['ERGO-II McPAS','ERGO-II VDJdb','NetTCR2.2 known-overlap-excluded','NetTCR2.2 all20','ImRex downsample','ImRex full']
all_scores={}
for name,tag in zip(method_order[:2],['ERGO-II_mcpas','ERGO-II_vdjdb']):
 rows=[r for r in ergo if r['model']==tag]
 assert len(rows)==12
 all_scores[name]={r['case_id']:float(r['binding_score']) for r in rows}
for name,col in zip(method_order[2:4],['NetTCR22_pan_known_overlap_excluded_4_mean','NetTCR22_pan_all20_mean']):
 assert [r['case_id'] for r in net]==ids
 all_scores[name]={r['case_id']:float(r[col]) for r in net}
for name,suffix in zip(method_order[4:],['downsample','full']):
 rows=read(root/f'imrex/predictions_{suffix}.csv')
 assert [r['case_id'] for r in rows]==ids
 for p,r in zip(panel,rows):assert p['CDR3_beta']==r['cdr3'] and p['peptide']==r['antigen.epitope']
 all_scores[name]={r['case_id']:float(r['prediction_score']) for r in rows}
for d in all_scores.values():assert set(d)==set(ids) and all(math.isfinite(v) and 0<=v<=1 for v in d.values())
wide=[];long=[];sums=[]
for row in panel:
 w=row.copy()
 for method in method_order:
  score=all_scores[method][row['case_id']]
  w[method]=score
  long.append(dict(case_id=row['case_id'],tcr=row['tcr'],peptide=row['peptide'],experimental_binding=int(row['experimental_binding']),experimental_KD_nM=row['experimental_KD_nM'],model_configuration=method,binding_score=score,reference_class_at_0p5=int(score>=.5)))
 wide.append(w)
for method in method_order:
 d=all_scores[method];ela=d['4jff'];recs=[r for r in long if r['model_configuration']==method]
 sums.append(dict(model_configuration=method,ELA_score=ela,G4A_score=d['4jfd'],I7A_score=d['4jfe'],G4A_delta=d['4jfd']-ela,I7A_delta=d['4jfe']-ela,G4A_lower_than_ELA=d['4jfd']<ela,I7A_lower_than_ELA=d['4jfe']<ela,I5A_NB_score=d['mature_I5A'],reference_TP=sum(r['experimental_binding']==1 and r['reference_class_at_0p5']==1 for r in recs),reference_FN=sum(r['experimental_binding']==1 and r['reference_class_at_0p5']==0 for r in recs),reference_TN=sum(r['experimental_binding']==0 and r['reference_class_at_0p5']==0 for r in recs),reference_FP=sum(r['experimental_binding']==0 and r['reference_class_at_0p5']==1 for r in recs)))
for name,rows in [('combined_binding_scores_wide.csv',wide),('combined_binding_scores_long.csv',long),('descriptive_model_contrasts.csv',sums)]:
 with (root/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(root/'combined_results_audit.json').write_text(json.dumps(dict(rows=12,configurations=6,model_families=3,all72_scores_finite_in_0_1=True,panel_input_sha256=hashlib.sha256((root/'binding_panel_inputs.csv').read_bytes()).hexdigest(),score_summary=sums,interpretation='Only within-method score order is compared; 0.5 is a fixed uncalibrated reference, not validated KD threshold. All original core four states experimentally bind. No method is required to output KD; correct mutation score direction is explicitly retained.'),ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(sums,indent=2))
