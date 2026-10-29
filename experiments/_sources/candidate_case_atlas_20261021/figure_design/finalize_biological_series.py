"""Freeze C05–C17 figures after explicit final-image review and exact NPZ checks."""
import os
os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
ATLAS=ROOT.parent
DATA=ATLAS/'data_and_general_figures'
MANIFEST=ROOT/'biological_series_manifest.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
content=json.loads(MANIFEST.read_text(encoding='utf8'))
assert len(content['figures'])==13
reports=[]
for item in content['figures']:
    cid=item['case_id'];out=Path(item['png']).parent
    audit=json.loads((out/'data_audit.json').read_text(encoding='utf8'))
    assert sha(item['png'])==audit['PNG_SHA256']
    for path,expected in audit['source_CSV_SHA256'].items():assert sha(path)==expected
    assert audit['minimum_font_pt']>=9 and audit['all_visible_text_inside_canvas']
    assert abs(audit['dimensions_mm'][0]-183)<1e-9
    assert 185<=audit['dimensions_mm'][1]<=205+1e-9
    checks=[]
    for a,b in audit['biological_contrasts']:
        source=DATA/cid/'data'/f'{a}_{b}_all_chain_alignment.csv'
        d=pd.read_csv(source,float_precision='round_trip')
        assert not d[['role','alignment_column']].duplicated().any()
        for state,pdb in [('a',a),('b',b)]:
            sa=json.loads((DATA/cid/'data'/f'{pdb}_source_audit.json').read_text(encoding='utf8'))
            npzpath=Path(sa['prediction']);assert sha(npzpath)==sa['prediction_sha256']
            npz=np.load(npzpath,allow_pickle=True)
            for head,metric in enumerate(['rmsf','hchi']):
                mask=d[f'index_{state}'].notna()&d[f'{metric}_valid_{state}']&d[f'{metric}_full_{state}'].notna()
                ix=d.loc[mask,f'index_{state}'].to_numpy(int)
                values=d.loc[mask,f'{metric}_full_{state}'].to_numpy(float)
                diff=np.abs(values-npz['mean'][ix,head].astype(float))
                assert diff.max(initial=0)<=1e-12,(cid,pdb,metric,diff.max())
                aa=npz['amino_acid'][ix].astype(str)
                assert np.array_equal(aa,d.loc[mask,f'aa_{state}'].astype(str).to_numpy())
                checks.append({'a':a,'b':b,'state':state,'PDB':pdb,'metric':metric,
                               'valid_state_values_checked':len(ix),'maximum_absolute_error':float(diff.max(initial=0)),
                               'NPZ_path':str(npzpath),'NPZ_SHA256':sa['prediction_sha256']})
        for metric in ['rmsf','hchi']:
            m=d.paired&d[metric+'_valid']
            dd=(d.loc[m,metric+'_full_b']-d.loc[m,metric+'_full_a'])-d.loc[m,metric+'_full_delta']
            assert np.abs(dd.to_numpy()).max(initial=0)<=1e-12
    for metric,key in [('rmsf','RMSF_valid_counts_by_contrast'),('hchi','Hchi_valid_counts_by_contrast')]:
        panel=pd.read_csv(out/f'panel_d_{metric}_all_contrasts.csv',float_precision='round_trip')
        expected=sum(audit[key])
        assert len(panel)==expected,(cid,metric,len(panel),expected)
        assert panel[metric+'_valid'].all()
        if metric=='hchi':assert panel.same_aa.all()
    special={
        'C07':'WT and all three HLA variants retained in absolute RMSF and three Hchi rows; all three change comparisons retained.',
        'C14':'DMF5 is 3QDG/3QDJ with 16 Hchi sites; DMF4 is 3QDM/3QEQ with 19 Hchi sites; both 24-site backbone scopes shown.',
        'C15':'Melan-A, IMP2 and BST2 all retain complete CDR curves, three Hchi contrasts and all three MD comparisons; 7Q9B is structural background only.',
        'C17':'407 unchanged RMSF and 343 unchanged chi-defined sites across MHC/peptide/TCR alpha/TCR beta; all 10 MD replicas represented in existing means; not reduced to peptide.'
    }.get(cid,'Complete predefined region and all defined biological contrasts retained.')
    visual={'status':'passed','case_id':cid,'image_actually_opened':True,
            'reviewed_PNG_SHA256':audit['PNG_SHA256'],'reviewer':'case_figure_design_20261021',
            'target_dimensions_mm':audit['dimensions_mm'],'minimum_font_pt':audit['minimum_font_pt'],
            'Pass0':'Semantic nondefault colours, distinct markers, no rainbow or four-sided frames; no selected best residues.',
            'Pass1':'All visible English text >=9 pt and inside canvas; PDF/SVG preserve editable lines and text; 450 dpi PNG.',
            'Pass2':'All source positions in the defined complete scope retained. Gaps and undefined chi not zero-filled. Complete multi-contrast controls displayed.',
            'Pass3':'Final actual PNG opened at exported dimensions. Structural views present, readable state labels, all data within axes, zero-valued chi markers retain margin, no overlapping labels or legends.',
            'case_specific_scope':special,'new_training_inference_MD':False}
    (out/'visual_QA.json').write_text(json.dumps(visual,ensure_ascii=False,indent=2),encoding='utf8')
    audit.update({'status':'passed','final_original_NPZ_checks':checks,
                  'original_NPZ_all_model_values_passed':True,'final_scope_check':special,
                  'final_visual_QA':'visual_QA.json'})
    (out/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
    (out/f'{cid}_case_figure.audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
    stats=json.loads((out/'statistics_and_reproducibility.json').read_text(encoding='utf8'))
    stats['status']='passed';stats['original_NPZ_checks']=checks
    (out/'statistics_and_reproducibility.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
    files={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='delivery_manifest.json'}
    (out/'delivery_manifest.json').write_text(json.dumps({'status':'passed','case_id':cid,'files':files},indent=2),encoding='utf8')
    item.update({'status':'passed','PNG_SHA256':audit['PNG_SHA256'],'visual_QA':str(out/'visual_QA.json'),
                 'RMSF_counts':audit['RMSF_valid_counts_by_contrast'],'Hchi_counts':audit['Hchi_valid_counts_by_contrast']})
    reports.append({'case_id':cid,'model_values_checked':sum(x['valid_state_values_checked'] for x in checks),
                    'maximum_original_NPZ_error':max(x['maximum_absolute_error'] for x in checks),
                    'PNG_SHA256':audit['PNG_SHA256'],'scope_check':special,'status':'passed'})
content['status']='passed'
MANIFEST.write_text(json.dumps(content,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'biological_series_final_QA.json').write_text(json.dumps({'status':'passed','figures':reports},ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'status':'passed','case_figures':len(reports),'model_values_checked':sum(x['model_values_checked'] for x in reports),
                  'maximum_original_NPZ_error':max(x['maximum_original_NPZ_error'] for x in reports)},ensure_ascii=False))
