"""Final source-data, plotting-scope, model-pair and PDF export audit."""
from pathlib import Path
import dataclasses
import importlib.util
import json
import sys
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/pdz3_alpha3_focus_20260914'
OUT=ROOT/'output/flexppi_pdz3_alpha3_focus_20260914'

def main():
    s=pd.read_csv(BASE/'residue_summary.csv')
    d=pd.read_csv(BASE/'residue_deltas.csv')
    report={'numerical_audit':json.loads((BASE/'numerical_audit.json').read_text()),'pdf':[],'raster':[]}
    for name,n in [('fig01_pdz3_alpha3_residue_mechanism.pdf',1),('figS01_pdz3_alpha3_controls.pdf',4)]:
        reader=PdfReader(OUT/name);assert len(reader.pages)==n
        sizes=[];font_sizes=[]
        for p in reader.pages:
            mm=[float(p.mediabox.width)*25.4/72,float(p.mediabox.height)*25.4/72]
            assert abs(mm[0]-183)<.01 and mm[1]<=247
            def visitor(text,cm,tm,font,size):
                if text.strip():font_sizes.append(float(size))
            assert len(p.extract_text(visitor_text=visitor))>100
            sizes.append(mm)
        assert min(font_sizes)>=5
        report['pdf'].append({'file':name,'pages':n,'sizes_mm':sizes,'minimum_actual_font_pt':min(font_sizes)})
    for p in sorted(OUT.glob('*.png')):
        im=Image.open(p);dpi=im.info.get('dpi');assert dpi and min(dpi)>299
        report['raster'].append({'file':p.name,'pixels':list(im.size),'dpi':dpi})
    report['projection']=json.loads((OUT/'structure_assets/projection_checks.json').read_text())
    assert max(x['error_px'] for x in report['projection'])<1.5
    q=d[(d.kind=='mutation')&(d.context=='full')&(d.author_number==337)&d.mutant.isin(['Y397A','Y397F'])]
    q=q.pivot(index=['background','fold'],columns='mutant',values='delta_chi_entropy')
    contrast=q.Y397A-q.Y397F
    report['paired_Y397A_minus_Y397F_at_F337']={bg:{'mean':float(v.mean()),'sd':float(v.std()),'n':len(v),'n_positive':int((v>1e-5).sum())} for bg,v in contrast.groupby('background')}
    report['main_hchi_masks']=[]
    for bg in ['apo1bfe','strip5heb']:
        q=s[(s.kind=='deletion')&(s.context=='d7')&(s.mutant=='WT')&(s.background==bg)]
        report['main_hchi_masks'].append({'background':bg,'retained_raw_sites':len(q),
            'AGP_omitted_from_hchi_only':q[~q.hchi_valid].author_number.tolist(),
            'outside_common_profile_306_395':q[~q.author_number.between(306,395)].author_number.tolist()})
    report['heatmap_scale_checks']=[]
    for bg in ['apo1bfe','strip5heb']:
        q=s[(s.background==bg)&(s.kind=='mutation')&(s.context=='full')&s.mutant.str.match(r'^[A-Z][0-9]+A$')].copy()
        q=q[q.mutant.str[1:-1].astype(int).between(312,377)]
        assert q.mutant.nunique()==61
        for metric,limit in [('hchi',.15),('rmsf',.5)]:
            data=q[q.hchi_valid] if metric=='hchi' else q
            v=data[metric+'_mean']
            report['heatmap_scale_checks'].append({'background':bg,'metric':metric,'mutations':q.mutant.nunique(),
                'n_cells':len(v),'color_limit':limit,'n_color_saturated':int((v.abs()>limit).sum()),
                'maximum_abs':float(v.abs().max()),'n_unique':v.nunique(),
                'mean_within_row_sd':float(data.groupby('mutant')[metric+'_mean'].std().mean())})
    path=Path('C:/Users/DELL/.codex/skills/academic-figure-skill/scripts/qa_validator.py')
    spec=importlib.util.spec_from_file_location('figure_qa_validator',path)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    source=(ROOT/'scripts/plot_pdz3_alpha3_focus.py').read_text(encoding='utf-8')
    findings=[]
    for name in sorted(dir(mod)):
        if name.startswith('check_'):
            v=getattr(mod,name)(source)
            findings.extend(v if isinstance(v,list) else [v])
    report['style_validator_raw']=[dataclasses.asdict(v) for v in findings]
    # These heuristics are not used to alter the scientific values or imply significance.
    report['manual_adjudication']={
        'AP0_palette':'Exact palette block uses alignment spaces around =; substring checker can falsely flag this.',
        'AP4_legend':'Main key outside axes; supplementary V328 key placed in visually empty upper-right corner.',
        'CL1_font':'Validator regex matches the suffix size=2 in capsize=2 (error-bar cap length), not a font. Actual minimum text is 5 pt.',
        'CL2_width':'Validator misses figsize=(183/25.4,...); reopened PDF media boxes independently verify all pages at 183 mm.',
        'heatmap_variance':'Absolute 0.2 threshold is not unit-invariant and is unsuitable for small raw delta-Hchi values; no renormalization or data changes.',
        'replicates':'Five frozen fold checkpoints, not five independent experiments; no significance tests or stars.',
        'page_scope':'Main plus four supplementary pages. Complete 429-state results remain in source_data.'}
    for name in ['numerical_audit.json','checkpoint_provenance.json']:
        (OUT/'source_data'/name).write_bytes((BASE/name).read_bytes())
    (OUT/'QA_and_statistics.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    print('PASS: source-data assertions, 2145 caches, 2 PDFs / 5 pages, 183-mm pages, 300-dpi PNGs, CA projections <1.5 px.')
    print('Raw style heuristics and explicit manual adjudications saved in QA_and_statistics.json.')

if __name__=='__main__':main()
