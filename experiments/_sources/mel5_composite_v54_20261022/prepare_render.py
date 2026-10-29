from pathlib import Path
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'mel5_contact_v53_20261022'
old_doc='paper_v53_20261022_中文修订_MEL5残基接触与双动态联合图_三线表.docx'
new_doc='paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx'
for name in ['export_v53_word.ps1','render_v53.py']:
    s=(OLD/name).read_text(encoding='utf-8')
    s=s.replace(old_doc,new_doc).replace('v53','v54')
    (ROOT/name.replace('v53','v54')).write_text(s,encoding='utf-8')
print('Prepared owned read-only Word export and bundled render_docx.py wrapper.')
