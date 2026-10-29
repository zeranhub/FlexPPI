from pathlib import Path
import json,hashlib,subprocess,shutil
import numpy as np,pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parent;OLD=R.parent/'revision_v93_dual_case_entropy';BASE=R.parents[2]
O=R/'audit';O.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
plot=pd.read_csv(R/'figure3/source_data/all_plotted_values.csv');source=pd.read_csv(R/'figure3/source_data/defined_angles_long.csv',float_precision='round_trip')
stats=pd.read_csv(R/'figure3/source_data/per_chi_correlations.csv',float_precision='round_trip');stats=stats[stats.predictor.eq('predicted_joint_Hchi')&stats.metric.eq('fine_entropy')]
counts={};err=[];rh=[]
for chi in range(1,5):
    z=source[source.chi.eq('chi'+str(chi))];d=plot[plot.panel.eq('f')&plot.chi.eq(chi)]
    assert len(d)==len(z)*2
    for kind,col in [('predicted_joint_Hchi_AA_centered','predicted_Hchi_nats_AA_centered'),('MD_single_chi_fine_entropy_AA_centered','fine_entropy_four_origins_mean_nats_AA_centered')]:
        merged=d[d.kind.eq(kind)].assign(resid=lambda a:a.identity.astype(int)).merge(z,left_on='resid',right_on='author_residue_id',validate='one_to_one')
        assert len(merged)==len(z);err.extend(abs(merged.value-merged[col]).to_list())
    rho=float(spearmanr(z.predicted_Hchi_nats_AA_centered,z.fine_entropy_four_origins_mean_nats_AA_centered).statistic)
    given=float(stats.loc[stats.chi.eq('chi'+str(chi)),'rho_AA_centered'].item());rh.append(abs(rho-given));counts[chi]=dict(n=len(z),rho=rho)
assert max(err)<1e-12 and max(rh)<1e-12
previous=pd.read_csv(OLD/'figure3/source_data/all_plotted_values.csv')
for panel in ['b','c','e']:
    a=plot[plot.panel.eq(panel)].sort_values(['kind','identity','series']);b=previous[previous.panel.eq(panel)].sort_values(['kind','identity','series'])
    assert np.allclose(a.value,b.value,atol=1e-14,rtol=0)
current=plot[plot.panel.eq('g')].sort_values(['metric','identity','series']);old=previous[previous.panel.eq('f')].sort_values(['metric','identity','series']);assert np.allclose(current.value,old.value,atol=1e-14,rtol=0)
s4=pd.read_csv(R/'supplement4/figS4_all_plotted_points.csv')
selection=pd.read_csv(R/'supplement4/retained_prediction_groups_all498.csv')
assert s4.groupby('panel').size().to_dict()=={'a':6126,'b':498,'c':332,'d':432}
for row in s4[s4.panel.eq('c')].itertuples():
    value=selection.loc[selection.metric_key.eq(row.metric)&selection.sample_id.eq(row.sample_id)&selection.predictor.eq('full_hchi'),row.kind].item();assert abs(value-row.y)<1e-12
for stem,folder in [('fig3_dual_case_entropy','figure3'),('figS4_joint_readouts_angular_occupancy','supplement4')]:
    path=R/folder;fonts=subprocess.run([r'C:\texlive\2026\bin\windows\pdffonts.exe',str(path/(stem+'.pdf'))],check=True,capture_output=True,text=True).stdout
    assert 'ArialMT' in fonts and 'yes yes yes' in fonts
    (O/(folder+'_embedded_fonts.txt')).write_text(fonts,encoding='utf-8')
    if folder=='figure3':
        qa=json.loads((path/'figure_QA.json').read_text(encoding='utf-8'))
        assert qa['min_font']==9 and qa['print_font_at_164mm']>=8 and not qa['outside_text']
        qa['all83_records_in_two_cohort_descriptors']=True;qa.pop('all83_records_in_both_main_diagnostics',None)
        qa.update(visual_review='Root opened final full PNG and actual Word page11; all restored scatter panels visible.',reviewed_png_sha256=sha(path/(stem+'.png')))
        (path/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
    else:
        qa=json.loads((path/'figS4_audit.json').read_text(encoding='utf-8'));qa.update(status='passed',visual_QA='Root opened final PNG and actual Word page29; complete paired selection diagnostics visible.',reviewed_png_sha256=sha(path/(stem+'.png')))
        (path/'figS4_audit.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
    for ext in ['png','pdf','svg']:
        target=BASE/(('figure3_v94_20261026_恢复逐角关联.' if folder=='figure3' else 'figureS4_v94_20261026_构象占据补充.')+ext)
        shutil.copy2(path/(stem+'.'+ext),target)
audit=dict(status='passed',restored_f_all151_points=True,total_main_quantitative_values=len(plot),per_chi_counts_and_correlations=counts,source_value_max_difference=max(err),correlation_max_difference=max(rh),other_main_data_unchanged=True,previous_main_selection_now_S4c_all83=True,no_new_predictions_or_MD=True,actual_PDF_fonts_embedded=True)
(O/'figure_numerical_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(audit,ensure_ascii=False))
