from pathlib import Path
import json,hashlib,zipfile
base=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010'); root=base/'section33_expansion_20261025/model_first_story_search_20261025'; v=root/'revision_v81_absolute_profiles'; old=root/'revision_v80_redesign'
f=v/'figure3'; m=v/'manuscript'
qa=json.loads((f/'figure_QA.json').read_text(encoding='utf-8'))
qa['visual_review']='passed: root and independent numerical/visual audit'
qa['five_window_ranges_checked']={'contact_groups':[1.4945415632520695,2.0640260378519693],'all_nine_Met':[1.3599,2.9024],'all_intervals_within_axes':True}
qa['pdf_fonts']='Arial and Arial Bold: embedded CID TrueType; confirmed by pdffonts'
qa['source_data_identical_to_v80']={}
for name in ['cam_md_state_contrast_all_common_residues.csv','cam_md_state_contrast_all_residues_five_blocks.csv']:
    data=(f/'source_data'/name).read_bytes(); assert data==(old/'figure3/source_data'/name).read_bytes(); qa['source_data_identical_to_v80'][name]=hashlib.sha256(data).hexdigest()
(f/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
for name in ['fig3_CaM_integrated.png','fig3_CaM_integrated.pdf','fig3_CaM_integrated.svg','fig3_body.py']:
    (f/name).unlink()
(f/'source_data/all_plotted_quantities.csv').unlink()
prep=v/'prepare_v81.py'; s=prep.read_text(encoding='utf-8').replace("(m/'v80_text_plan.json')", "(old/'manuscript/v80_text_plan.json')").replace("(m/'build_v80_manuscript.py')", "(old/'manuscript/build_v80_manuscript.py')"); prep.write_text(s,encoding='utf-8')
for name in ['audit_v80_scientific_sources.py','build_v80_manuscript.py','v80_build_audit.json','v80_final_text.txt','v80_scientific_source_audit.json','v80_source_audit.json','v80_text_plan.json']:
    (m/name).unlink()
audit=json.loads((m/'v81_build_audit.json').read_text(encoding='utf-8')); audit['render_visual_QA_pending']=False; audit['render_visual_QA']='28 pages; changed pages 10,11,12 visually checked; remaining 25 pages pixel-identical to reviewed v80'
(m/'v81_build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
report=json.loads((v/'render/v81_page_comparison.json').read_text(encoding='utf-8')); report['changed_pages_visual_review']='passed: 10,11,12';report['overall_review']='passed';(v/'render/v81_page_comparison.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
archive=base/'paper_v81_20261025_图3双状态绝对值源文件.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(f.rglob('*')):
        if p.is_file(): z.write(p,p.relative_to(v))
with zipfile.ZipFile(archive) as z: assert z.testzip() is None
print(json.dumps({'docx':audit['output'],'source_zip':str(archive),'zip_bytes':archive.stat().st_size,'render_review':'passed','pages':28},ensure_ascii=False))